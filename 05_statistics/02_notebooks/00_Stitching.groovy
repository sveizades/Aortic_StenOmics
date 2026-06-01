/**
 * Convert TIFF fields of view to a pyramidal OME-TIFF, unmodified stitching logic from Pete Bankhead.
 *
 * This version:
 *   1) Lets you pick multiple input folders (each containing *.tif/*.tiff).
 *   2) Lets you pick ONE output directory where all stitched .ome.tif files will go.
 *   3) For each folder, uses Pete's logic to parse metadata, build a SparseImageServer & write OME-TIFF.
 *
 * Stitching logic, parseRegion(), parseRegionFromTIFF(), OMEPyramidWriter settings, etc. remain unchanged.
 *
 * @author 
 *   Original: Pete Bankhead
 *   Modified solely for input/output prompting: [Your Name]
 * QuPath version: 0.5.x
 */

import qupath.fx.dialogs.FileChoosers
import qupath.lib.common.GeneralTools
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.images.servers.ImageServers
import qupath.lib.images.servers.SparseImageServer
import qupath.lib.images.writers.ome.OMEPyramidWriter
import qupath.lib.regions.ImageRegion

import javax.imageio.ImageIO
import javax.imageio.plugins.tiff.BaselineTIFFTagSet
import javax.imageio.plugins.tiff.TIFFDirectory
import java.awt.image.BufferedImage

import static qupath.lib.gui.scripting.QPEx.*

// ------------------------------------------------------------------------------------
// Only changed: We now prompt for multiple input folders instead of multiple files
// ------------------------------------------------------------------------------------
List<File> inputFolders = []
boolean keepChoosingFolders = true
while (keepChoosingFolders) {
    def chosenFolder = FileChoosers.promptForDirectory('Select a folder with .tif/.tiff (Cancel to finish)', null)
    if (chosenFolder == null) {
        keepChoosingFolders = false
    } else if (!inputFolders.contains(chosenFolder)) {
        inputFolders << chosenFolder
    }
}

if (!inputFolders) {
    print 'No input folders selected. Stopping.'
    return
}

// ------------------------------------------------------------------------------------
// Only changed: Now we prompt once for the output directory where .ome.tif files go
// ------------------------------------------------------------------------------------
def outputDir = FileChoosers.promptForDirectory('Select output directory for .ome.tif files', null)
if (outputDir == null) {
    print 'No output directory selected. Stopping.'
    return
}

// ------------------------------------------------------------------------------------
// For each input folder, gather .tif files, then run EXACT Pete’s stitching logic
// ------------------------------------------------------------------------------------
for (folder in inputFolders) {

    // This is how we gather the files that Pete’s script would have gotten from promptForMultipleFiles
    List<File> files = folder.listFiles()?.findAll { f ->
        f.isFile() &&
        (f.name.toLowerCase().endsWith('.tif') || f.name.toLowerCase().endsWith('.tiff')) &&
        !f.name.toLowerCase().endsWith('.ome.tif') // skip existing OME-TIFFs
    } ?: []

    // If no .tif, skip
    if (!files) {
        print "No TIFF files found in folder: ${folder}"
        continue
    }

    // This 'baseName' logic is the same style as Pete’s default
    // but naming the output after the folder
    String baseName = folder.name

    // Create a unique .ome.tif file in the output directory
    File fileOutput = new File(outputDir, baseName + '.ome.tif')
    int count = 1
    while (fileOutput.exists()) {
        fileOutput = new File(outputDir, baseName + '-' + count + '.ome.tif')
        count++
    }

    print "Parsing regions from ${files.size()} files in '${folder.name}'..."

    // -------------------
    // PETE’S STITCHING LOGIC, UNCHANGED
    // -------------------
    def builder = new SparseImageServer.Builder()
    files.parallelStream().forEach { f ->
        def region = parseRegion(f)
        if (region == null) {
            print 'WARN: Could not parse region for ' + f
            return
        }
        def serverBuilder = ImageServerProvider
                .getPreferredUriImageSupport(BufferedImage.class, f.toURI().toString())
                .getBuilders().get(0)
        builder.jsonRegion(region, 1.0, serverBuilder)
    }
    print 'Building server...'
    def server = builder.build()
    server = ImageServers.pyramidalize(server)

    long startTime = System.currentTimeMillis()
    String pathOutput = fileOutput.getAbsolutePath()
    new OMEPyramidWriter.Builder(server)
            .downsamples(server.getPreferredDownsamples())
            .tileSize(512)
            .channelsInterleaved()
            .parallelize()
            .losslessCompression()
            .build()
            .writePyramid(pathOutput)
    long endTime = System.currentTimeMillis()
    print('Image written to ' + pathOutput + ' in ' +
          GeneralTools.formatNumber((endTime - startTime) / 1000.0, 1) + ' s')
    server.close()
}

// ------------------------------------------------------------------------------------
// PETE’S parseRegion, parseRegionFromTIFF & supporting methods, UNCHANGED
// ------------------------------------------------------------------------------------
static ImageRegion parseRegion(File file, int z = 0, int t = 0) {
    if (checkTIFF(file)) {
        try {
            return parseRegionFromTIFF(file, z, t)
        } catch (Exception e) {
            print e.getLocalizedMessage()
        }
    }
    return null
}

/**
 * Check for TIFF 'magic number'.
 */
static boolean checkTIFF(File file) {
    file.withInputStream {
        def bytes = it.readNBytes(4)
        short byteOrder = toShort(bytes[0], bytes[1])
        int val
        if (byteOrder == 0x4949) {
            // Little-endian
            val = toShort(bytes[3], bytes[2])
        } else if (byteOrder == 0x4d4d) {
            val = toShort(bytes[2], bytes[3])
        } else
            return false
        return val == 42 || val == 43
    }
}

/**
 * Combine two bytes to create a short, in the given order
 */
static short toShort(byte b1, byte b2) {
    return (b1 << 8) + (b2 << 0)
}

/**
 * Parse an ImageRegion from a TIFF image, using the metadata.
 */
static ImageRegion parseRegionFromTIFF(File file, int z = 0, int t = 0) {
    int x, y, width, height
    file.withInputStream {
        def reader = ImageIO.getImageReadersByFormatName("TIFF").next()
        reader.setInput(ImageIO.createImageInputStream(it))
        def metadata = reader.getImageMetadata(0)
        def tiffDir = TIFFDirectory.createFromMetadata(metadata)

        double xRes = getRational(tiffDir, BaselineTIFFTagSet.TAG_X_RESOLUTION)
        double yRes = getRational(tiffDir, BaselineTIFFTagSet.TAG_Y_RESOLUTION)

        double xPos = getRational(tiffDir, BaselineTIFFTagSet.TAG_X_POSITION)
        double yPos = getRational(tiffDir, BaselineTIFFTagSet.TAG_Y_POSITION)

        width = tiffDir.getTIFFField(BaselineTIFFTagSet.TAG_IMAGE_WIDTH).getAsLong(0) as int
        height = tiffDir.getTIFFField(BaselineTIFFTagSet.TAG_IMAGE_LENGTH).getAsLong(0) as int

        x = Math.round(xRes * xPos) as int
        y = Math.round(yRes * yPos) as int
    }
    return ImageRegion.createInstance(x, y, width, height, z, t)
}

/**
 * Helper for parsing rational from TIFF metadata.
 */
static double getRational(TIFFDirectory tiffDir, int tag) {
    long[] rational = tiffDir.getTIFFField(tag).getAsRational(0)
    return rational[0] / (double)rational[1]
}
