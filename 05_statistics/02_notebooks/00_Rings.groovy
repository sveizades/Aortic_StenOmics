import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.precision.GeometryPrecisionReducer
import org.locationtech.jts.geom.PrecisionModel

import qupath.lib.objects.PathObjects
import qupath.lib.roi.GeometryTools
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.RoiTools

import org.locationtech.jts.geom.Geometry
import org.locationtech.jts.geom.GeometryCollection
import org.locationtech.jts.geom.MultiPolygon
import org.locationtech.jts.geom.Polygon
import org.locationtech.jts.geom.GeometryFactory

// import org.locationtech.jts.geom.util.GeometryFixer
// import org.locationtech.jts.operation.valid.IsValidOp

// import org.locationtech.jts.geom.TopologyException

// PrecisionModel pm = new PrecisionModel(PrecisionModel.FIXED)

// Geometry clean(Geometry g) {
  //   if (g == null || g.isEmpty()) return g
    // g = GeometryPrecisionReducer.reduce(g, pm)
    // g = g.buffer(0)                    // drop NaN‑Z, snap nodes, fix rings
  //   if (!g.isValid()) {
   //      g = GeometryFixer.fix(g)       // JTS ≥ 1.18
 //    }
  //   return g
// }


Geometry flattenToMultiPolygon(Geometry geom) {
    if (geom instanceof Polygon || geom instanceof MultiPolygon) {
        return geom
    }
    if (geom instanceof GeometryCollection) {
        def factory = geom.getFactory()
        def polygons = []
        for (int i = 0; i < geom.getNumGeometries(); i++) {
            def g = geom.getGeometryN(i)
            if (g instanceof Polygon) {
                polygons << g
            }
        }
        if (polygons.isEmpty()) {
            return null
        }
        return factory.createMultiPolygon(polygons.toArray(new Polygon[0]))
    }
    return null
}

// Geometry safeUnion(Geometry a, Geometry b) {
//     a = clean(a)
//     b = clean(b)
// //     if (a == null || a.isEmpty()) return b
//     if (b == null || b.isEmpty()) return a
//     try {
//         return a.union(b)
//     } catch (TopologyException te) {        // last‑resort fallback
//         return GeometryCombiner.combine(a, b).buffer(0)
//     }
// }


//Geometry safeUnion(Geometry geom1, Geometry geom2) {
//   GeometryFactory factory = geom1.getFactory()
//
//    Geometry f1 = flattenToMultiPolygon(geom1)
//    Geometry f2 = flattenToMultiPolygon(geom2)
//
//    if (f1 == null) return f2
//    if (f2 == null) return f1
//
//    return f1.union(f2)
//}


def imageData = getCurrentImageData()
def server = imageData.getServer()
def cal = server.getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) {
    print 'We need the pixel size information here!'
    return
}

double expandMarginMicrons = 100.0
int howManyTimes = 1
double expandPixels = expandMarginMicrons / cal.getAveragedPixelSizeMicrons()
def plane = ImagePlane.getDefaultPlane()
PrecisionModel PM = new PrecisionModel(PrecisionModel.FIXED)

def allAnnotations = getAnnotationObjects()
def nodularTissues = allAnnotations.findAll { it.getPathClass() == getPathClass("Nodular Tissue") }
def tissueMasks = allAnnotations.findAll { it.getPathClass() == getPathClass("Tissue After Subtraction") }

def newAnnotations = []

nodularTissues.each { nodular ->
    def tumorGeom = nodular.getROI().getGeometry()
    def overlappingTissue = tissueMasks.find { it.getROI().getGeometry().intersects(tumorGeom) }

    if (overlappingTissue == null) {
        println "No matching tissue found for nodular region: ${nodular.getName()}"
        return
    }

    def tissueGeom = overlappingTissue.getROI().getGeometry()
    def cleanTumorGeom = tissueGeom.intersection(tumorGeom)
    def tumorROIClean = GeometryTools.geometryToROI(cleanTumorGeom, plane)
    def cleanTumor = PathObjects.createAnnotationObject(tumorROIClean, getPathClass("Nodular Tissue"))
    cleanTumor.setName(nodular.getName())  // <-- Keep original name here

    def annotationsToAdd = [cleanTumor]

    // === Add Inner Rim (100 µm inside the nodule) ===
    //def innerRimGeom = cleanTumorGeom.buffer(-expandPixels)
    //def innerBand = cleanTumorGeom.difference(innerRimGeom)
    //innerBand = GeometryPrecisionReducer.reduce(innerBand, PM)

    //if (!innerBand.isEmpty()) {
    //    def roiInner = GeometryTools.geometryToROI(innerBand, plane)
    //    def innerAnn = PathObjects.createAnnotationObject(roiInner, getPathClass("Inner Rim"))
    //    innerAnn.setName(nodular.getName() + " Inner Rim 100 µm")
    //    innerAnn.setColorRGB(getColorRGB(200, 100, 100))
    //    annotationsToAdd << innerAnn
    //}

    // === Add Outer Rim (100 µm outside the nodule, within tissue) ===
    def outerBuffer = cleanTumorGeom.buffer(expandPixels)
    def outerBand = outerBuffer.difference(cleanTumorGeom).intersection(tissueGeom)
    

    // === Add Inner Rim (100 µm inside the nodule) ===
    def innerRimGeom = outerBuffer.buffer(-expandPixels*1)
    def innernodule = outerBuffer.buffer(-expandPixels*2)
    //def innerRimGeom = outerBand.buffer(-expandPixels*2)
    //def innerBand = innerRimGeom.difference(outerBand).intersection(tissueGeom)
    def innerBand = innerRimGeom.buffer(0).difference(innernodule).intersection(tissueGeom)
    innerBand = GeometryPrecisionReducer.reduce(innerBand, PM)
    
    outerBand = outerBand.buffer(0).difference(innerBand.buffer(0)).difference(innernodule.buffer(0))
    outerBand = GeometryPrecisionReducer.reduce(outerBand, PM)
    
    if (!outerBand.isEmpty()) {
        def roiOuter = GeometryTools.geometryToROI(outerBand, plane)
        def outerAnn = PathObjects.createAnnotationObject(roiOuter, getPathClass("Outer Rim"))
        outerAnn.setName(nodular.getName().replace(' Nodular Tissue', '') + " Outer Rim 100 µm")
        outerAnn.setColor(100, 200, 100)
        annotationsToAdd << outerAnn
    }

    if (!innerBand.isEmpty()) {
        def roiInner = GeometryTools.geometryToROI(innerBand, plane)
        def innerAnn = PathObjects.createAnnotationObject(roiInner, getPathClass("Inner Rim"))
        innerAnn.setName(nodular.getName().replace(' Nodular Tissue', '') + " Inner Rim 100 µm")
        innerAnn.setColor(200, 100, 100)
        annotationsToAdd << innerAnn
    }
    
     if (!innernodule.isEmpty()) {
        def roiInnernod = GeometryTools.geometryToROI(innernodule, plane)
        def innerAnnnod = PathObjects.createAnnotationObject(roiInnernod, getPathClass("Inner Nodule"))
        innerAnnnod.setName(nodular.getName().replace(' Nodular Tissue', '') + " Inner nodule 100 µm")
        innerAnnnod.setColor(96, 96, 96)
        annotationsToAdd << innerAnnnod
    }
    
    
    // === Combined Rim (inner + outer) ===
//    if (innerBand != null && outerBand != null && !innerBand.isEmpty() && !outerBand.isEmpty()) {
//        def combinedRim = innerRimGeom.intersection(tissueGeom)
//if (combinedRim == null) {
 //   println "Union resulted in null geometry."
//    return
//    }
//    
//
//        combinedRim = GeometryPrecisionReducer.reduce(combinedRim, PM)
//        def roiCombined = GeometryTools.geometryToROI(combinedRim, plane)
//        def combinedAnn = PathObjects.createAnnotationObject(roiCombined, getPathClass("Combined Rim"))
//        combinedAnn.setName(nodular.getName() + " Combined Rim 100 µm")
 //       combinedAnn.setColor(150, 150, 255)
 //       annotationsToAdd << combinedAnn
//    }
    // === Continue with expanded margins (like before) ===
    for (i = 0; i < howManyTimes; i++) {
        def currentArea = annotationsToAdd[-1].getROI().getGeometry()
        def areaExpansion = currentArea.buffer(expandPixels)
        areaExpansion = areaExpansion.intersection(tissueGeom)
        areaExpansion = areaExpansion.buffer(0).difference(cleanTumorGeom)

        if (i >= 1) {
            (1..i).each { k ->
                def remove = annotationsToAdd[-k].getROI().getGeometry()
                areaExpansion = areaExpansion.difference(remove)
            }
        }

      // areaExpansion = GeometryPrecisionReducer.reduce(areaExpansion, PM)
      //  def roiExpansion = GeometryTools.geometryToROI(areaExpansion, plane)
        //int nameValue = (i + 1) * expandMarginMicrons
        //def annotationExpansion = PathObjects.createAnnotationObject(roiExpansion, getPathClass(nameValue.toString()))
        //annotationExpansion.setName("Margin " + nameValue + " microns")
        //annotationExpansion.setColorRGB(getColorRGB(50 * (i + 1), 40 * (i + 1), 200 - 30 * (i + 1)))
        //annotationsToAdd << annotationExpansion
    }

    def remainingTissueGeom = tissueGeom.difference(cleanTumorGeom)
    annotationsToAdd.each {
        remainingTissueGeom = remainingTissueGeom.buffer(0).difference(it.getROI().getGeometry().buffer(0))
    }
    def remainingTissueROI = GeometryTools.geometryToROI(remainingTissueGeom, plane)
    def remainingTissue = PathObjects.createAnnotationObject(remainingTissueROI, getPathClass("Other Tissue"))
    remainingTissue.setName(nodular.getName().replace(' Nodular Tissue', '') + " Other Tissue")
    annotationsToAdd << remainingTissue

    newAnnotations.addAll(annotationsToAdd)
    removeObject(nodular, true)
}

addObjects(newAnnotations)
resetSelection()
fireHierarchyUpdate()
println("Done processing all nodular tissues!")
