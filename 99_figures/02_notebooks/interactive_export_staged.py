"""Two-stage interactive image export pipeline.

Stage 1 – TissueAligner
    Rotate the raw image and broadly crop to a single tissue piece.
    Saves an aligned .npy array + _stage1_params.json for Stage 2.

Stage 2 – RegionExporter
    Load the aligned tile, draw a refined ROI and optional zoom inset,
    then export publication-quality figures (PDF / TIFF / PNG).

Usage
-----
# ── Stage 1 ──────────────────────────────────────────────────────────────────
aligner = TissueAligner(czi_path="/path/to/slide.czi",
                        output_filename="slide_A")
aligner.show()
# → saves slide_A_stage1.npy and slide_A_stage1_params.json

# ── Stage 2 ──────────────────────────────────────────────────────────────────
exporter = RegionExporter(
    stage1_npy="slide_A_stage1.npy",
    stage1_params="slide_A_stage1_params.json",
    output_filename="figure_1A",
)
exporter.show()

# ── Resume a previous Stage 2 session ────────────────────────────────────────
exporter = RegionExporter(
    stage1_npy="slide_A_stage1.npy",
    stage1_params="slide_A_stage1_params.json",
    params_file="figure_1A_parameters.json",
    output_filename="figure_1A",
)
exporter.show()

# ── Rebuild a missing Stage 1 .npy from the JSON ─────────────────────────────
from interactive_export_staged import rebuild_stage1_npy

rebuild_stage1_npy("slide_A_stage1_params.json")
# If the CZI has moved:
rebuild_stage1_npy("slide_A_stage1_params.json", czi_path="/new/path/slide.czi")
# Save to a custom location instead of overwriting the original:
rebuild_stage1_npy("slide_A_stage1_params.json", output_npy="slide_A_stage1_v2.npy")
"""

import napari
from napari.utils.notifications import show_info
import numpy as np
from scipy import ndimage
import json
import czifile
from pathlib import Path
from datetime import datetime
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib_scalebar.scalebar import ScaleBar
from qtpy.QtWidgets import (
    QPushButton, QVBoxLayout, QWidget, QLabel, QSlider,
    QLineEdit, QHBoxLayout, QFileDialog,
)
from qtpy.QtCore import Qt


# ─────────────────────────────────────────────────────────────────────────────
# Shared utilities
# ─────────────────────────────────────────────────────────────────────────────

def extract_pixel_size_from_czi(czi_path):
    """Return the X pixel size in µm from CZI metadata (defaults to 1.0)."""
    with czifile.CziFile(czi_path) as czi:
        metadata_xml = czi.metadata()
    try:
        root = ET.fromstring(metadata_xml)
        for distance in root.findall(".//Distance"):
            if distance.get('Id') == 'X':
                value = float(distance.find('Value').text)
                pixel_size_um = value * 1e6
                print(f"  Pixel size: {pixel_size_um:.4f} µm")
                return pixel_size_um
    except Exception as e:
        print(f"  Could not extract pixel size ({e}). Using 1.0 µm.")
    return 1.0


def load_czi_image(czi_path):
    """Load a CZI file, squeeze extra dims, reorder to (H, W, C), fix black bg."""
    print(f"Loading CZI: {Path(czi_path).name}")
    with czifile.CziFile(czi_path) as czi:
        arr = np.squeeze(czi.asarray())
    print(f"  Shape: {arr.shape}")

    if arr.ndim == 3 and arr.shape[0] in (3, 4):
        if arr.shape[0] < min(arr.shape[1], arr.shape[2]):
            arr = np.moveaxis(arr, 0, -1)
            print(f"  Transposed to: {arr.shape}")

    if arr.ndim == 3 and arr.shape[-1] == 3:
        black = (arr[:, :, 0] == 0) & (arr[:, :, 1] == 0) & (arr[:, :, 2] == 0)
        n = int(np.sum(black))
        if n > 0:
            arr[black] = 255
            print(f"  Converted {n:,} black background pixels to white")

    return arr


def rotate_image(img, angle):
    """Rotate *img* by *angle* degrees, using cv2 when available.

    Fills empty regions with white (255).  Expands the canvas to fit the
    rotated image (no clipping).
    """
    try:
        import cv2
        h, w = img.shape[:2]
        rad = np.radians(angle)
        ca, sa = np.abs(np.cos(rad)), np.abs(np.sin(rad))
        new_w = int(h * sa + w * ca)
        new_h = int(h * ca + w * sa)
        if max(h, w, new_h, new_w) >= 32767:
            raise OverflowError("Image too large for cv2.warpAffine")
        M = cv2.getRotationMatrix2D((w / 2, h / 2), angle, 1.0)
        M[0, 2] += (new_w - w) / 2
        M[1, 2] += (new_h - h) / 2
        bv = (255,) * img.shape[2] if img.ndim == 3 else 255
        return cv2.warpAffine(img, M, (new_w, new_h),
                              flags=cv2.INTER_LINEAR,
                              borderMode=cv2.BORDER_CONSTANT,
                              borderValue=bv)
    except (ImportError, OverflowError):
        orig_dtype = img.dtype
        cval = np.iinfo(orig_dtype).max if np.issubdtype(orig_dtype, np.integer) else 1.0
        if img.ndim == 3:
            result = np.stack([
                ndimage.rotate(img[:, :, c], angle,
                               reshape=True, order=1, cval=cval)
                for c in range(img.shape[2])
            ], axis=-1)
        else:
            result = ndimage.rotate(img, angle, reshape=True, order=1, cval=cval)
        # scipy always returns float64; cast back to the original dtype
        if np.issubdtype(orig_dtype, np.integer):
            info = np.iinfo(orig_dtype)
            result = np.clip(result, info.min, info.max)
        return result.astype(orig_dtype)


def _make_scalebar(pixel_size_um, fixed_um=None):
    """Return a matplotlib_scalebar ScaleBar artist."""
    kwargs = dict(location='lower right', box_alpha=0.7, color='black',
                  frameon=True, scale_loc='top')
    if fixed_um is not None and fixed_um > 0:
        return ScaleBar(pixel_size_um * 1e-6, 'm',
                        fixed_value=fixed_um, fixed_units='um', **kwargs)
    return ScaleBar(pixel_size_um * 1e-6, 'm', length_fraction=0.25, **kwargs)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 – TissueAligner
# ─────────────────────────────────────────────────────────────────────────────

class TissueAligner:
    """Stage 1: rotate the raw image and broadly crop to one tissue piece.

    Workflow
    --------
    1. Use the rotation slider to level the tissue.
    2. Click **"Preview Rotation"** – the viewer will display a downsampled
       rotated version of the image.
    3. Draw a rectangle on the **"Broad Crop"** layer (green) to encompass the
       tissue.  The rectangle should be drawn *after* previewing the rotation so
       that coordinates are in the rotated image space.
    4. Click **"Save Aligned (Stage 1)"**.

    Outputs (in *output_dir*)
    -------------------------
    ``<stem>_stage1.npy``          – rotated, cropped image array
    ``<stem>_stage1_params.json``  – metadata for Stage 2

    Notes
    -----
    * If rotation is 0 you may skip the preview and draw directly on the
      original image.
    * When rotation ≠ 0 the crop rectangle **must** be drawn on the rotated
      preview, not the original.
    """

    def __init__(self, image_data=None, czi_path=None, pixel_size_um=None,
                 output_dir=None, output_filename=None):
        self.czi_path = czi_path

        if czi_path is not None:
            image_data = load_czi_image(czi_path)
            if pixel_size_um is None:
                pixel_size_um = extract_pixel_size_from_czi(czi_path)

        if image_data is None:
            raise ValueError("Provide either image_data or czi_path.")

        # Normalise to (H, W[, C])
        if image_data.ndim == 3 and image_data.shape[0] in (3, 4):
            if image_data.shape[0] < min(image_data.shape[1], image_data.shape[2]):
                image_data = np.moveaxis(image_data, 0, -1)

        self.original_image = image_data
        self.pixel_size_um = pixel_size_um if pixel_size_um is not None else 1.0
        self.rotation = 0.0
        self.output_dir = output_dir if output_dir is not None else str(Path.cwd())
        self.output_filename = output_filename if output_filename is not None else "tissue_aligned"
        self._preview_mode = False
        self._preview_ds = 1  # downsampling factor used in the current preview

        self.viewer = napari.Viewer()
        self._setup_viewer()

    # ── Napari setup ──────────────────────────────────────────────────────────

    def _setup_viewer(self):
        is_rgb = self.original_image.ndim == 3 and self.original_image.shape[-1] in (3, 4)
        scale = (self.pixel_size_um, self.pixel_size_um)
        self.image_layer = self.viewer.add_image(
            self.original_image, name='Image (Original)',
            rgb=is_rgb, scale=scale)

        self.crop_layer = self.viewer.add_shapes(
            name='Broad Crop',
            face_color='transparent',
            edge_color='lime',
            edge_width=4)

        # Image must be below the shapes layer so rectangles are visible
        self.viewer.layers.move(self.viewer.layers.index(self.image_layer), 0)

        self._build_controls()
        # Activate drawing mode after controls are built so Qt doesn't reset the selection
        self._activate_crop_drawing()

    def _build_controls(self):
        w = QWidget()
        lay = QVBoxLayout()

        lay.addWidget(QLabel("─── Stage 1: Tissue Alignment ───"))

        # Rotation slider
        lay.addWidget(QLabel("Rotation (°)"))
        self.rot_slider = QSlider(Qt.Horizontal)
        self.rot_slider.setMinimum(-180)
        self.rot_slider.setMaximum(180)
        self.rot_slider.setValue(0)
        self.rot_slider.setTickInterval(10)
        self.rot_slider.setTickPosition(QSlider.TicksBelow)
        self.rot_slider.valueChanged.connect(self._on_rot_changed)
        self.rot_label = QLabel("0°")
        lay.addWidget(self.rot_slider)
        lay.addWidget(self.rot_label)

        # Preview / reset buttons
        btn_row = QHBoxLayout()
        prev_btn = QPushButton("Preview Rotation")
        prev_btn.clicked.connect(self._preview_rotation)
        btn_row.addWidget(prev_btn)
        reset_btn = QPushButton("Reset Preview")
        reset_btn.clicked.connect(self._reset_preview)
        btn_row.addWidget(reset_btn)
        btn_w = QWidget(); btn_w.setLayout(btn_row)
        lay.addWidget(btn_w)

        # Drawing / auto-detect buttons
        suggest_btn = QPushButton("Suggest Crop from Tissue")
        suggest_btn.setToolTip(
            "Auto-detect non-white pixels in the current view and draw\n"
            "a suggested bounding rectangle.  Run after 'Preview Rotation'.")
        suggest_btn.clicked.connect(self._suggest_crop)
        lay.addWidget(suggest_btn)

        draw_btn = QPushButton("Draw Crop Rectangle (manual)")
        draw_btn.setToolTip(
            "Click this to activate the Broad Crop layer and switch\n"
            "to rectangle-draw mode before drawing on the image.")
        draw_btn.clicked.connect(self._activate_crop_drawing)
        lay.addWidget(draw_btn)

        lay.addWidget(QLabel("\n─── Output ───"))

        # Output directory
        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel("Output Dir:"))
        self.dir_input = QLineEdit(self.output_dir)
        self.dir_input.textChanged.connect(lambda t: setattr(self, 'output_dir', t))
        dir_row.addWidget(self.dir_input)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse_dir)
        dir_row.addWidget(browse)
        dir_w = QWidget(); dir_w.setLayout(dir_row)
        lay.addWidget(dir_w)

        # Filename stem
        fn_row = QHBoxLayout()
        fn_row.addWidget(QLabel("Filename stem:"))
        self.fn_input = QLineEdit(self.output_filename)
        self.fn_input.textChanged.connect(lambda t: setattr(self, 'output_filename', t))
        fn_row.addWidget(self.fn_input)
        fn_w = QWidget(); fn_w.setLayout(fn_row)
        lay.addWidget(fn_w)

        # Instructions
        lay.addWidget(QLabel(
            "\nInstructions:\n"
            "1. Adjust rotation slider to level tissue.\n"
            "2. Click 'Preview Rotation' to apply.\n"
            "3. Draw rectangle on 'Broad Crop' layer\n"
            "   (green) around the tissue piece.\n"
            "   (Draw AFTER previewing rotation.)\n"
            "4. Click 'Save Aligned' to save output.\n\n"
            "Outputs:\n"
            "  <stem>_stage1.npy\n"
            "  <stem>_stage1_params.json"
        ))

        save_btn = QPushButton("Save Aligned (Stage 1)")
        save_btn.clicked.connect(self._save_aligned)
        lay.addWidget(save_btn)

        self.info_label = QLabel(f"Pixel size: {self.pixel_size_um:.4f} µm")
        lay.addWidget(self.info_label)

        w.setLayout(lay)
        self.viewer.window.add_dock_widget(w, name='Stage 1 Controls', area='right')

    # ── Rotation preview ──────────────────────────────────────────────────────

    def _suggest_crop(self):
        """Auto-detect non-white tissue in the current image and draw a bounding crop."""
        img = self.image_layer.data
        print("Detecting tissue (non-white pixels)…")

        # Threshold: a pixel is "tissue" if any channel is darker than 240
        if img.ndim == 3:
            not_white = np.any(img < 240, axis=2)
        else:
            not_white = img < 240

        rows = np.where(np.any(not_white, axis=1))[0]
        cols = np.where(np.any(not_white, axis=0))[0]

        if len(rows) == 0 or len(cols) == 0:
            show_info(
                "No tissue detected — the visible image is entirely white.\n\n"
                "If rotation is set, click 'Preview Rotation' first, then\n"
                "click 'Suggest Crop from Tissue' again."
            )
            return

        r1, r2 = int(rows[0]), int(rows[-1])
        c1, c2 = int(cols[0]), int(cols[-1])
        scale_y, scale_x = self.image_layer.scale[0], self.image_layer.scale[1]

        # Convert pixel bounds → world coordinates (µm) for the shapes layer
        rect = np.array([
            [r1 * scale_y, c1 * scale_x],
            [r1 * scale_y, c2 * scale_x],
            [r2 * scale_y, c2 * scale_x],
            [r2 * scale_y, c1 * scale_x],
        ])

        # Clear old shapes and add new suggestion
        if len(self.crop_layer.data) > 0:
            self.crop_layer.selected_data = set(range(len(self.crop_layer.data)))
            self.crop_layer.remove_selected()
        self.crop_layer.add_rectangles(rect)

        # Switch to pan/zoom so the user can inspect, then refine if needed
        self.crop_layer.mode = 'pan_zoom'
        self.viewer.layers.selection.active = self.crop_layer

        w_px = int((c2 - c1) * scale_x / self.pixel_size_um)
        h_px = int((r2 - r1) * scale_y / self.pixel_size_um)
        print(f"  Tissue bounds (preview): rows {r1}–{r2}, cols {c1}–{c2}")
        print(f"  Full-res crop size: {w_px}×{h_px} px")
        show_info(
            f"Tissue detected!\n"
            f"Suggested crop: {w_px}×{h_px} full-res pixels.\n\n"
            "Adjust the rectangle if needed, then click\n"
            "'Save Aligned (Stage 1)'."
        )

    def _activate_crop_drawing(self):
        """Select the Broad Crop layer and switch it to rectangle-draw mode."""
        self.viewer.layers.selection.active = self.crop_layer
        self.crop_layer.mode = 'add_rectangle'

    def _on_rot_changed(self, value):
        self.rotation = float(value)
        suffix = " (PREVIEW ACTIVE)" if self._preview_mode else ""
        self.rot_label.setText(f"{value}°{suffix}")

    def _preview_rotation(self):
        if self.rotation == 0:
            show_info("Rotation is 0° — nothing to preview.")
            return

        print(f"Generating rotation preview ({self.rotation}°)…")
        h, w = self.original_image.shape[:2]
        ds = max(1, max(h, w) // 4096)
        if ds > 1:
            src = (self.original_image[::ds, ::ds, :]
                   if self.original_image.ndim == 3
                   else self.original_image[::ds, ::ds])
            print(f"  Downsampled {ds}× for speed ({src.shape[1]}×{src.shape[0]})")
        else:
            src = self.original_image

        rotated = rotate_image(src, self.rotation)
        contrast = self.image_layer.contrast_limits
        self.viewer.layers.remove(self.image_layer)

        is_rgb = rotated.ndim == 3 and rotated.shape[-1] in (3, 4)
        # Scale accounts for downsampling so shapes stay in correct world coords
        preview_scale = (self.pixel_size_um * ds, self.pixel_size_um * ds)
        self.image_layer = self.viewer.add_image(
            rotated, name='Image (Rotated Preview)',
            rgb=is_rgb, scale=preview_scale)
        if not is_rgb:
            self.image_layer.contrast_limits = contrast

        # Push image back behind the shapes layer
        self.viewer.layers.move(self.viewer.layers.index(self.image_layer), 0)

        self._preview_ds = ds
        self._preview_mode = True
        self.rot_label.setText(f"{self.rotation}° (PREVIEW ACTIVE)")

        # Re-select the shapes layer — adding the image layer above it deselects it
        self._activate_crop_drawing()

        show_info(
            "Rotation preview active.\n"
            "Draw your broad crop rectangle, then click 'Save Aligned'.\n"
            "Click 'Reset Preview' to return to the original image."
        )

    def _reset_preview(self):
        if not self._preview_mode:
            show_info("No preview active.")
            return
        contrast = self.image_layer.contrast_limits
        self.viewer.layers.remove(self.image_layer)
        is_rgb = self.original_image.ndim == 3 and self.original_image.shape[-1] in (3, 4)
        self.image_layer = self.viewer.add_image(
            self.original_image, name='Image (Original)',
            rgb=is_rgb, scale=(self.pixel_size_um, self.pixel_size_um))
        if not is_rgb:
            self.image_layer.contrast_limits = contrast

        # Push image back behind the shapes layer
        self.viewer.layers.move(self.viewer.layers.index(self.image_layer), 0)

        self._preview_mode = False
        self._preview_ds = 1
        self.rot_label.setText(f"{self.rotation}°")

        # Re-select the shapes layer after swapping the image layer
        self._activate_crop_drawing()

        show_info("Reset to original image.")

    def _browse_dir(self):
        d = QFileDialog.getExistingDirectory(None, "Select Output Directory", self.output_dir)
        if d:
            self.output_dir = d
            self.dir_input.setText(d)

    # ── Save ──────────────────────────────────────────────────────────────────

    def _get_crop_rect_px(self):
        """Return (x1, y1, x2, y2) in full-res pixels, or None.

        Shapes layer stores world coordinates (µm).  Dividing by pixel_size_um
        gives full-resolution pixel coordinates regardless of whether a
        downsampled preview is displayed.
        """
        if len(self.crop_layer.data) == 0:
            return None
        shape = self.crop_layer.data[-1]  # (4, 2) in layer data space
        scale = np.asarray(self.crop_layer.scale)[-2:]
        translate = np.asarray(self.crop_layer.translate)[-2:]
        shape = shape * scale + translate  # world coords (µm)
        rows = shape[:, 0] / self.pixel_size_um
        cols = shape[:, 1] / self.pixel_size_um
        return (int(np.floor(cols.min())),   # x1
                int(np.floor(rows.min())),   # y1
                int(np.ceil(cols.max())),    # x2
                int(np.ceil(rows.max())))    # y2

    def _save_aligned(self):
        """Rotate full image, crop to broad rect, save NPY + params JSON."""
        if self.rotation != 0 and not self._preview_mode:
            show_info(
                "Rotation is set but no preview is active.\n\n"
                "Please click 'Preview Rotation' first, then draw your\n"
                "crop rectangle on the rotated image before saving."
            )
            return

        crop_rect = self._get_crop_rect_px()
        if crop_rect is None:
            show_info("Please draw a rectangle on the 'Broad Crop' layer first!")
            return

        print(f"Saving aligned image (rotation={self.rotation}°)…")

        # Rotate full-resolution image
        if self.rotation != 0:
            print("  Rotating full-resolution image (may take a moment)…")
            rotated = rotate_image(self.original_image, self.rotation)
        else:
            rotated = self.original_image

        # Clamp crop rect to image bounds
        h, w = rotated.shape[:2]
        x1, y1, x2, y2 = crop_rect
        x1 = max(0, min(x1, w)); x2 = max(0, min(x2, w))
        y1 = max(0, min(y1, h)); y2 = max(0, min(y2, h))

        cropped = rotated[y1:y2, x1:x2]
        print(f"  Cropped shape: {cropped.shape}")

        # Save
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        base = Path(self.output_dir) / self.output_filename
        npy_path = str(base) + "_stage1.npy"
        params_path = str(base) + "_stage1_params.json"

        np.save(npy_path, cropped)
        print(f"  Saved: {npy_path}")

        params = {
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "czi_path": str(self.czi_path) if self.czi_path else None,
                "original_shape": list(self.original_image.shape),
                "pixel_size_um": self.pixel_size_um,
                "stage": 1,
            },
            "parameters": {
                "rotation": float(self.rotation),
                # Crop rect in rotated full-res pixel coords (x1, y1, x2, y2)
                "broad_crop_bbox": [x1, y1, x2, y2],
            },
            "output": {
                "npy_path": npy_path,
                "aligned_shape": list(cropped.shape),
            },
        }
        with open(params_path, 'w') as f:
            json.dump(params, f, indent=2)
        print(f"  Saved: {params_path}")

        show_info(
            f"Stage 1 complete!\n\n"
            f"Saved: {Path(npy_path).name}\n"
            f"       {Path(params_path).name}\n"
            f"Shape: {cropped.shape[1]}×{cropped.shape[0]} px\n\n"
            f"Load in Stage 2:\n"
            f"  RegionExporter(\n"
            f"    stage1_npy='{npy_path}',\n"
            f"    stage1_params='{params_path}'\n"
            f"  ).show()"
        )

    def show(self):
        """Start the Napari viewer (blocking)."""
        print("Starting Stage 1 Napari viewer…")
        try:
            napari.run()
        except RuntimeError as e:
            print(f"Note: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 – headless rebuild
# ─────────────────────────────────────────────────────────────────────────────

def rebuild_stage1_npy(params_file, czi_path=None, output_npy=None):
    """Rebuild a Stage 1 .npy file from the saved _stage1_params.json.

    Useful when the .npy was deleted or moved, or when you want to regenerate
    it on a different machine.

    Parameters
    ----------
    params_file : str or Path
        Path to the ``*_stage1_params.json`` file produced by TissueAligner.
    czi_path : str or Path, optional
        Override the CZI path stored in the JSON.  Required if the original
        CZI has been moved since Stage 1 was run.
    output_npy : str or Path, optional
        Where to write the rebuilt .npy.  Defaults to the ``npy_path`` stored
        inside the JSON (i.e. the original location).

    Returns
    -------
    npy_path : str
        Path to the written .npy file.
    """
    params_file = Path(params_file)
    print(f"[rebuild_stage1_npy] Reading params from {params_file}")

    with open(params_file) as f:
        saved = json.load(f)

    meta   = saved["metadata"]
    p      = saved["parameters"]
    out_   = saved.get("output", {})

    # Resolve CZI path
    czi = czi_path or meta.get("czi_path")
    if czi is None:
        raise ValueError(
            "No czi_path found in params and none provided. "
            "Pass czi_path= explicitly."
        )
    czi = Path(czi)
    if not czi.exists():
        raise FileNotFoundError(f"CZI file not found: {czi}")

    rotation       = float(p.get("rotation", 0.0))
    broad_crop     = p.get("broad_crop_bbox")   # [x1, y1, x2, y2] in rotated px
    pixel_size_um  = float(meta.get("pixel_size_um", 1.0))

    if broad_crop is None:
        raise ValueError("'broad_crop_bbox' not found in params — cannot crop.")

    x1, y1, x2, y2 = broad_crop

    # Resolve output path
    if output_npy is None:
        output_npy = out_.get("npy_path")
    if output_npy is None:
        # Fall back to placing it next to the params file
        output_npy = params_file.with_name(
            params_file.name.replace("_stage1_params.json", "_stage1.npy")
        )
    output_npy = Path(output_npy)

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"  Loading CZI: {czi.name}")
    image = load_czi_image(str(czi))
    print(f"  Original shape: {image.shape}")

    # ── Rotate ────────────────────────────────────────────────────────────────
    if rotation != 0:
        print(f"  Rotating by {rotation}° (full resolution, may take a moment)…")
        image = rotate_image(image, rotation)
        print(f"  Rotated shape: {image.shape}")
    else:
        print("  No rotation.")

    # ── Crop ──────────────────────────────────────────────────────────────────
    h, w = image.shape[:2]
    x1c = max(0, min(x1, w));  x2c = max(0, min(x2, w))
    y1c = max(0, min(y1, h));  y2c = max(0, min(y2, h))
    cropped = image[y1c:y2c, x1c:x2c]
    print(f"  Cropped: ({x1c},{y1c})→({x2c},{y2c})  shape={cropped.shape}")

    # ── Save ──────────────────────────────────────────────────────────────────
    output_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(str(output_npy), cropped)
    print(f"  Saved: {output_npy}")

    # Warn if the rebuilt shape differs from what was saved originally
    orig_shape = out_.get("aligned_shape")
    if orig_shape is not None and list(cropped.shape) != orig_shape:
        print(
            f"  WARNING: rebuilt shape {list(cropped.shape)} differs from "
            f"original {orig_shape}.  This can happen if cv2 vs scipy was used "
            f"for rotation, or if a different CZI was passed."
        )

    return str(output_npy)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 – RegionExporter
# ─────────────────────────────────────────────────────────────────────────────

class RegionExporter:
    """Stage 2: draw a refined ROI and optional zoom inset on the aligned tile.

    Loads the .npy + _stage1_params.json produced by TissueAligner.  No further
    rotation is applied here (already handled in Stage 1).

    Workflow
    --------
    1. Adjust contrast via the Napari layer controls (left panel).
    2. Draw a rectangle on the **"ROI"** layer (cyan) for the main region.
    3. (Optional) Draw a rectangle on the **"Zoom Region"** layer (yellow) for
       a zoom inset.  The zoom rectangle should fall inside the main ROI.
    4. Set output directory and filename.
    5. Click **"Export Image"**.

    Outputs
    -------
    ``<filename>.pdf/.tiff/.png``          – overview figure
    ``<filename>_zoom.pdf/.tiff/.png``     – zoom inset (if drawn)
    ``<filename>_parameters.json``         – reproducibility params

    Reload a previous session
    -------------------------
    Pass ``params_file="<filename>_parameters.json"`` to restore the ROI /
    zoom rectangles and contrast from a previous run.

    Pre-populate rectangles
    -----------------------
    ``roi_rect``  / ``zoom_rect`` accept either:

    * ``(width_um, height_um)``          – rectangle centred on the image
    * ``(x1_um, y1_um, x2_um, y2_um)``  – explicit bounding box

    Both are in µm.  A ``params_file`` will override these if provided.
    """

    def __init__(self, stage1_npy, stage1_params,
                 params_file=None,
                 output_dir=None, output_filename=None,
                 scalebar_um=None, zoom_scalebar_um=None,
                 roi_rect=None, zoom_rect=None):
        self.stage1_npy = stage1_npy
        self.stage1_params_path = stage1_params

        with open(stage1_params, 'r') as f:
            s1 = json.load(f)
        self.pixel_size_um = s1["metadata"]["pixel_size_um"]
        self._stage1_meta = s1

        print(f"Loading Stage 1 image: {Path(stage1_npy).name}")
        self.image = np.load(stage1_npy)
        print(f"  Shape: {self.image.shape}  dtype: {self.image.dtype}  |  pixel size: {self.pixel_size_um:.4f} µm")

        # If saved as float64 (scipy rotation fallback), convert to uint8 so
        # Napari displays it correctly (float images are expected in 0–1 range).
        if self.image.dtype == np.float64 or self.image.dtype == np.float32:
            print(f"  Converting float image (range {self.image.min():.1f}–{self.image.max():.1f}) to uint8")
            self.image = np.clip(self.image, 0, 255).astype(np.uint8)

        # Export settings
        self.export_dpi = 300
        self.figure_width_cm = 15.0
        self.output_dir = output_dir if output_dir is not None else str(Path.cwd())
        self.output_filename = output_filename if output_filename is not None else "publication_export"
        self.scalebar_um = scalebar_um
        self.zoom_scalebar_um = zoom_scalebar_um

        # ROI state
        self.roi_bbox = None
        self.zoom_bbox = None
        self.vmin = None
        self.vmax = None

        # Pre-populate rectangles from roi_rect / zoom_rect (µm).
        # params_file takes precedence if both are given.
        if roi_rect is not None:
            self.roi_bbox = self._parse_rect_um(roi_rect, self.image.shape)
        if zoom_rect is not None:
            self.zoom_bbox = self._parse_rect_um(zoom_rect, self.image.shape)

        if params_file is not None and Path(params_file).exists():
            self._load_parameters(params_file)

        self.viewer = napari.Viewer()
        self._setup_viewer()

    # ── Rect helpers ──────────────────────────────────────────────────────────

    def _parse_rect_um(self, rect, shape):
        """Convert a µm rect spec to a pixel bbox (x1, y1, x2, y2).

        Parameters
        ----------
        rect : (w_um, h_um) or (x1_um, y1_um, x2_um, y2_um)
            Two-element tuple  → width × height, centred on the image.
            Four-element tuple → explicit bounding box.
        shape : (H, W[, C])
        """
        ps = self.pixel_size_um
        img_h, img_w = shape[:2]
        if len(rect) == 2:
            w_px = int(round(rect[0] / ps))
            h_px = int(round(rect[1] / ps))
            cx, cy = img_w // 2, img_h // 2
            x1 = max(0, cx - w_px // 2)
            y1 = max(0, cy - h_px // 2)
            x2 = min(img_w, x1 + w_px)
            y2 = min(img_h, y1 + h_px)
        elif len(rect) == 4:
            x1 = max(0, int(round(rect[0] / ps)))
            y1 = max(0, int(round(rect[1] / ps)))
            x2 = min(img_w, int(round(rect[2] / ps)))
            y2 = min(img_h, int(round(rect[3] / ps)))
        else:
            raise ValueError("roi_rect / zoom_rect must be (w, h) or (x1, y1, x2, y2) in µm")
        return (x1, y1, x2, y2)

    # ── Napari setup ──────────────────────────────────────────────────────────

    def _setup_viewer(self):
        is_rgb = self.image.ndim == 3 and self.image.shape[-1] in (3, 4)
        scale = (self.pixel_size_um, self.pixel_size_um)

        self.image_layer = self.viewer.add_image(
            self.image, name='Aligned Image',
            rgb=is_rgb, scale=scale)

        if self.vmin is not None and self.vmax is not None:
            self.image_layer.contrast_limits = (self.vmin, self.vmax)

        # ROI layer (cyan)
        self.roi_layer = self.viewer.add_shapes(
            name='ROI',
            face_color='transparent',
            edge_color='cyan',
            edge_width=3)

        if self.roi_bbox is not None:
            x1, y1, x2, y2 = self.roi_bbox
            rect = np.array([[y1, x1], [y1, x2], [y2, x2], [y2, x1]]) * self.pixel_size_um
            self.roi_layer.add_rectangles(rect)

        # Zoom layer (yellow)
        self.zoom_layer = self.viewer.add_shapes(
            name='Zoom Region',
            face_color='transparent',
            edge_color='yellow',
            edge_width=2)

        if self.zoom_bbox is not None:
            zx1, zy1, zx2, zy2 = self.zoom_bbox
            zrect = np.array([[zy1, zx1], [zy1, zx2], [zy2, zx2], [zy2, zx1]]) * self.pixel_size_um
            self.zoom_layer.add_rectangles(zrect)

        # Image must be below both shapes layers so rectangles are visible
        self.viewer.layers.move(self.viewer.layers.index(self.image_layer), 0)

        self._build_controls()
        # Activate ROI drawing after controls are built so Qt doesn't reset the selection
        self._activate_roi_drawing()

    def _build_controls(self):
        w = QWidget()
        lay = QVBoxLayout()

        lay.addWidget(QLabel("─── Stage 2: Region Export ───"))

        # Drawing activation buttons
        draw_roi_btn = QPushButton("Draw ROI Rectangle (cyan)")
        draw_roi_btn.setToolTip("Activate the ROI layer in rectangle-draw mode.")
        draw_roi_btn.clicked.connect(self._activate_roi_drawing)
        lay.addWidget(draw_roi_btn)

        draw_zoom_btn = QPushButton("Draw Zoom Rectangle (yellow)")
        draw_zoom_btn.setToolTip("Activate the Zoom Region layer in rectangle-draw mode.")
        draw_zoom_btn.clicked.connect(self._activate_zoom_drawing)
        lay.addWidget(draw_zoom_btn)

        lay.addWidget(QLabel("\n─── Export Settings ───"))

        # DPI
        lay.addWidget(QLabel(f"Export DPI"))
        self.dpi_slider = QSlider(Qt.Horizontal)
        self.dpi_slider.setMinimum(150)
        self.dpi_slider.setMaximum(600)
        self.dpi_slider.setValue(self.export_dpi)
        self.dpi_slider.setTickInterval(50)
        self.dpi_slider.setTickPosition(QSlider.TicksBelow)
        self.dpi_label = QLabel(f"{self.export_dpi} DPI")
        self.dpi_slider.valueChanged.connect(self._on_dpi_changed)
        lay.addWidget(self.dpi_slider)
        lay.addWidget(self.dpi_label)

        # Figure width
        lay.addWidget(QLabel("Figure Width"))
        self.width_slider = QSlider(Qt.Horizontal)
        self.width_slider.setMinimum(50)    # 5.0 cm (stored as tenths)
        self.width_slider.setMaximum(500)   # 50.0 cm
        self.width_slider.setValue(int(self.figure_width_cm * 10))
        self.width_slider.setTickInterval(10)
        self.width_slider.setTickPosition(QSlider.TicksBelow)
        self.width_label = QLabel(f"{self.figure_width_cm:.1f} cm")
        self.width_slider.valueChanged.connect(self._on_width_changed)
        lay.addWidget(self.width_slider)
        lay.addWidget(self.width_label)

        # Scale bar (overview)
        sb_row = QHBoxLayout()
        sb_row.addWidget(QLabel("Scale bar (µm):"))
        self.sb_input = QLineEdit("" if self.scalebar_um is None else str(self.scalebar_um))
        self.sb_input.setPlaceholderText("auto")
        self.sb_input.setToolTip("Fixed scale bar length in µm for the overview image.\nLeave blank to let matplotlib_scalebar choose automatically.")
        self.sb_input.textChanged.connect(self._on_scalebar_changed)
        sb_row.addWidget(self.sb_input)
        sb_w = QWidget(); sb_w.setLayout(sb_row)
        lay.addWidget(sb_w)

        # Scale bar (zoom)
        zsb_row = QHBoxLayout()
        zsb_row.addWidget(QLabel("Zoom scale bar (µm):"))
        self.zsb_input = QLineEdit("" if self.zoom_scalebar_um is None else str(self.zoom_scalebar_um))
        self.zsb_input.setPlaceholderText("auto")
        self.zsb_input.setToolTip("Fixed scale bar length in µm for the zoom image.\nLeave blank to let matplotlib_scalebar choose automatically.")
        self.zsb_input.textChanged.connect(self._on_zoom_scalebar_changed)
        zsb_row.addWidget(self.zsb_input)
        zsb_w = QWidget(); zsb_w.setLayout(zsb_row)
        lay.addWidget(zsb_w)

        lay.addWidget(QLabel("\n─── Output Settings ───"))

        # Output directory
        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel("Output Dir:"))
        self.dir_input = QLineEdit(self.output_dir)
        self.dir_input.textChanged.connect(lambda t: setattr(self, 'output_dir', t))
        dir_row.addWidget(self.dir_input)
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse_dir)
        dir_row.addWidget(browse)
        dir_w = QWidget(); dir_w.setLayout(dir_row)
        lay.addWidget(dir_w)

        # Filename
        fn_row = QHBoxLayout()
        fn_row.addWidget(QLabel("Filename (no ext):"))
        self.fn_input = QLineEdit(self.output_filename)
        self.fn_input.textChanged.connect(lambda t: setattr(self, 'output_filename', t))
        fn_row.addWidget(self.fn_input)
        fn_w = QWidget(); fn_w.setLayout(fn_row)
        lay.addWidget(fn_w)

        # Instructions
        lay.addWidget(QLabel(
            "\nInstructions:\n"
            "1. Adjust contrast via layer controls.\n"
            "2. Click 'Draw ROI Rectangle' then draw\n"
            "   on the image for the main region.\n"
            "3. (Optional) Click 'Draw Zoom Rectangle'\n"
            "   then draw a zoom inset inside the ROI.\n"
            "4. Set output directory and filename.\n"
            "5. Click 'Export Image'.\n\n"
            "Exports: .pdf / .tiff / .png\n"
            "(+ _zoom if zoom region drawn)\n"
            "and _parameters.json"
        ))

        export_btn = QPushButton("Export Image")
        export_btn.clicked.connect(self._export)
        lay.addWidget(export_btn)

        self.info_label = QLabel(f"Pixel size: {self.pixel_size_um:.4f} µm")
        lay.addWidget(self.info_label)

        w.setLayout(lay)
        self.viewer.window.add_dock_widget(w, name='Stage 2 Controls', area='right')

    # ── Drawing activation ────────────────────────────────────────────────────

    def _activate_roi_drawing(self):
        """Select the ROI layer and switch to rectangle-draw mode."""
        self.viewer.layers.selection.active = self.roi_layer
        self.roi_layer.mode = 'add_rectangle'

    def _activate_zoom_drawing(self):
        """Select the Zoom Region layer and switch to rectangle-draw mode."""
        self.viewer.layers.selection.active = self.zoom_layer
        self.zoom_layer.mode = 'add_rectangle'

    # ── Control callbacks ─────────────────────────────────────────────────────

    def _on_dpi_changed(self, value):
        self.export_dpi = value
        self.dpi_label.setText(f"{value} DPI")

    def _on_width_changed(self, value):
        self.figure_width_cm = value / 10.0
        self.width_label.setText(f"{self.figure_width_cm:.1f} cm")

    def _on_scalebar_changed(self, text):
        try:
            self.scalebar_um = float(text) if text.strip() else None
        except ValueError:
            pass  # keep previous value while user is mid-type

    def _on_zoom_scalebar_changed(self, text):
        try:
            self.zoom_scalebar_um = float(text) if text.strip() else None
        except ValueError:
            pass  # keep previous value while user is mid-type

    def _browse_dir(self):
        d = QFileDialog.getExistingDirectory(None, "Select Output Directory", self.output_dir)
        if d:
            self.output_dir = d
            self.dir_input.setText(d)

    # ── Shape extraction ──────────────────────────────────────────────────────

    @staticmethod
    def _layer_world_coords(layer):
        """Return the last shape's vertices in world coordinates (µm).

        Accounts for any layer-level scale/translate that Napari may have
        applied when the user moved a shape.
        """
        shape = layer.data[-1]
        scale = np.asarray(layer.scale)
        translate = np.asarray(layer.translate)
        if scale.ndim == 0:
            scale = np.array([float(scale), float(scale)])
        else:
            scale = scale[-2:]
        if translate.ndim == 0:
            translate = np.array([float(translate), float(translate)])
        else:
            translate = translate[-2:]
        return shape * scale + translate

    def _rect_info(self, layer):
        """Return {'corners', 'bbox'} for the last shape in *layer*, or None."""
        if len(layer.data) == 0:
            return None
        shape = self._layer_world_coords(layer)  # (4, 2) in µm world coords
        corners_px = shape / self.pixel_size_um  # (row, col) in pixels
        return {
            'corners': corners_px,
            'bbox': (
                int(np.floor(corners_px[:, 1].min())),   # x1
                int(np.floor(corners_px[:, 0].min())),   # y1
                int(np.ceil(corners_px[:, 1].max())),    # x2
                int(np.ceil(corners_px[:, 0].max())),    # y2
            ),
        }

    # ── Export ────────────────────────────────────────────────────────────────

    def _export(self):
        roi_info = self._rect_info(self.roi_layer)
        if roi_info is None:
            show_info("Please draw a rectangle on the 'ROI' layer first!")
            return

        zoom_info = self._rect_info(self.zoom_layer)

        # Store bboxes
        self.roi_bbox = roi_info['bbox']
        self.zoom_bbox = zoom_info['bbox'] if zoom_info else None

        vmin, vmax = self.image_layer.contrast_limits

        # Crop ROI from aligned image
        x1, y1, x2, y2 = self.roi_bbox
        h, w = self.image.shape[:2]
        x1 = max(0, min(x1, w)); x2 = max(0, min(x2, w))
        y1 = max(0, min(y1, h)); y2 = max(0, min(y2, h))
        roi = self.image[y1:y2, x1:x2]

        # Zoom rect in ROI-crop coordinates
        zoom_rect = None
        if zoom_info is not None:
            zc = zoom_info['corners']  # (row, col) pixels
            zoom_rect = (
                zc[:, 0].min() - y1,
                zc[:, 1].min() - x1,
                zc[:, 0].max() - y1,
                zc[:, 1].max() - x1,
            )

        print(f"ROI crop: {roi.shape}  contrast: {vmin:.0f}–{vmax:.0f}")

        # Apply contrast stretch to uint8
        roi = np.clip(roi, vmin, vmax)
        roi = ((roi - vmin) / (vmax - vmin) * 255).astype(np.uint8)

        # Prepare output paths
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        base = Path(self.output_dir) / self.output_filename

        self._export_overview(roi, zoom_rect, base)
        if zoom_rect is not None:
            self._export_zoom(roi, zoom_rect, base)

        params_path = str(base) + "_parameters.json"
        self._save_parameters(params_path, vmin, vmax)

        print(f"Exported to: {self.output_dir}")
        show_info(f"Export complete!\nSaved to: {self.output_dir}")

    def _export_overview(self, roi, zoom_rect, base):
        """Save overview image (+ optional zoom-box overlay)."""
        aspect = roi.shape[0] / roi.shape[1]
        fig_w = self.figure_width_cm / 2.54
        fig_h = fig_w * aspect

        fig, ax = plt.subplots(figsize=(fig_w, fig_h),
                               dpi=self.export_dpi, facecolor='white')
        ax.imshow(roi, interpolation='none')
        ax.axis('off')
        ax.set_facecolor('white')

        if zoom_rect is not None:
            r1, c1, r2, c2 = zoom_rect
            ax.add_patch(Rectangle(
                (c1, r1), c2 - c1, r2 - r1,
                linewidth=2, edgecolor='black', facecolor='none'))

        ax.add_artist(_make_scalebar(self.pixel_size_um, self.scalebar_um))
        plt.tight_layout(pad=0)

        for ext in ('pdf', 'tiff', 'png'):
            kw = {'dpi': self.export_dpi, 'facecolor': 'white'}
            if ext == 'pdf':
                kw['format'] = 'pdf'
            fig.savefig(f"{base}.{ext}", **kw)
        plt.close(fig)
        print(f"  Saved: {base}.pdf / .tiff / .png")

    def _export_zoom(self, roi, zoom_rect, base):
        """Save the zoomed detail image."""
        r1, c1, r2, c2 = zoom_rect
        zr1 = max(0, int(np.floor(r1))); zr2 = min(roi.shape[0], int(np.ceil(r2)))
        zc1 = max(0, int(np.floor(c1))); zc2 = min(roi.shape[1], int(np.ceil(c2)))
        zoom_roi = roi[zr1:zr2, zc1:zc2]

        if zoom_roi.size == 0:
            print("  Zoom region is empty – skipping zoom export.")
            show_info("Zoom region is empty. No zoom image exported.")
            return

        aspect = zoom_roi.shape[0] / zoom_roi.shape[1]
        fig_w = self.figure_width_cm / 2.54

        fig, ax = plt.subplots(figsize=(fig_w, fig_w * aspect),
                               dpi=self.export_dpi, facecolor='white')
        ax.imshow(zoom_roi, interpolation='none')
        ax.axis('off')
        ax.set_facecolor('white')
        ax.add_artist(_make_scalebar(self.pixel_size_um, self.zoom_scalebar_um))
        plt.tight_layout(pad=0)

        for ext in ('pdf', 'tiff', 'png'):
            kw = {'dpi': self.export_dpi, 'facecolor': 'white'}
            if ext == 'pdf':
                kw['format'] = 'pdf'
            fig.savefig(f"{base}_zoom.{ext}", **kw)
        plt.close(fig)
        print(f"  Saved: {base}_zoom.pdf / .tiff / .png")

    # ── Parameters ────────────────────────────────────────────────────────────

    def _save_parameters(self, path, vmin, vmax):
        params = {
            "metadata": {
                "timestamp": datetime.now().isoformat(),
                "stage1_npy": str(self.stage1_npy),
                "stage1_params": str(self.stage1_params_path),
                "pixel_size_um": self.pixel_size_um,
                "stage": 2,
            },
            "parameters": {
                "roi_bbox": list(self.roi_bbox) if self.roi_bbox else None,
                "zoom_bbox": list(self.zoom_bbox) if self.zoom_bbox else None,
                "contrast_min": float(vmin),
                "contrast_max": float(vmax),
                "scalebar_um": self.scalebar_um,
                "zoom_scalebar_um": self.zoom_scalebar_um,
                "figure_width_cm": self.figure_width_cm,
            },
        }
        with open(path, 'w') as f:
            json.dump(params, f, indent=2)
        print(f"  Saved params: {path}")

    def _load_parameters(self, params_file):
        print(f"Loading Stage 2 parameters: {params_file}")
        with open(params_file, 'r') as f:
            data = json.load(f)
        p = data["parameters"]
        if p.get("roi_bbox"):
            self.roi_bbox = tuple(p["roi_bbox"])
            print(f"  ROI bbox:  {self.roi_bbox}")
        if p.get("zoom_bbox"):
            self.zoom_bbox = tuple(p["zoom_bbox"])
            print(f"  Zoom bbox: {self.zoom_bbox}")
        if p.get("contrast_min") is not None:
            self.vmin = p["contrast_min"]
        if p.get("contrast_max") is not None:
            self.vmax = p["contrast_max"]
        if p.get("scalebar_um") is not None:
            self.scalebar_um = p["scalebar_um"]
        if p.get("zoom_scalebar_um") is not None:
            self.zoom_scalebar_um = p["zoom_scalebar_um"]
        if p.get("figure_width_cm") is not None:
            self.figure_width_cm = p["figure_width_cm"]

    def show(self):
        """Start the Napari viewer (blocking)."""
        print("Starting Stage 2 Napari viewer…")
        try:
            napari.run()
        except RuntimeError as e:
            print(f"Note: {e}")
