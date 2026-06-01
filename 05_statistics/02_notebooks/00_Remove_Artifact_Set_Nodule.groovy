// REMOVE ARTIFACT FROM TISSUE OJBECTS

newAnnotations = []

// Find all Tissue annotations
tissueAnnotations = getAnnotationObjects().findAll { it.getPathClass() == getPathClass("Tissue") }

// Get all Artifact* geometries
artifactGeoms = getAnnotationObjects().findAll {
    pc = it.getPathClass()
    return pc != null && pc.getName().startsWith("Artifact")
}.collect { it.getROI().getGeometry() }

// For each tissue annotation, subtract all artifact geometries
tissueAnnotations.each { tissueAnnotation ->
    tissueGeom = tissueAnnotation.getROI().getGeometry()
    
    artifactGeoms.each { artGeom ->
        tissueGeom = tissueGeom.difference(artGeom)
    }
    
    if (!tissueGeom.isEmpty()) {
        newROI = GeometryTools.geometryToROI(tissueGeom, ImagePlane.getDefaultPlane())
        newAnno = PathObjects.createAnnotationObject(newROI, getPathClass("Tissue After Subtraction"))
        newAnno.setName(tissueAnnotation.getName() + " Tissue After Subtraction")
        newAnnotations << newAnno
    }
}

// Add all new annotations to the image
addObjects(newAnnotations)
fireHierarchyUpdate()

// ASSIGN NODULAR AND NON-NODULAR REGIONS OF TISSUE 

nonNodularAnnotations = []
nodularAnnotations = []

// Get all tissue annotations
tissueAnnotations = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("Tissue After Subtraction")}

// Get all nodule geometries
noduleGeoms = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("Nodule")}
                   .collect{ it.getROI().getGeometry() }

for (tissue in tissueAnnotations) {
    tissueGeom = tissue.getROI().getGeometry()
    
    // Subtract nodules → Non-nodular tissue
    nonNodularGeom = tissueGeom
    noduleGeoms.each { noduleGeom ->
        nonNodularGeom = nonNodularGeom.difference(noduleGeom)
    }
    if (!nonNodularGeom.isEmpty()) {
        roi = GeometryTools.geometryToROI(nonNodularGeom, ImagePlane.getDefaultPlane())
        obj = PathObjects.createAnnotationObject(roi, getPathClass("Non-Nodular Tissue"))
        def base = tissue.getName()
                    .replace(' Tissue After Subtraction', '') 
                    .trim()
        obj.setName(base + ' Non‑Nodular Tissue')
        // obj.setName(tissueAnnotation.getName() + " Non-Nodular Tissue")
        nonNodularAnnotations << obj
    }

    // Intersect nodules → Nodular tissue
    nodularGeom = null
    noduleGeoms.each { noduleGeom ->
        if (tissueGeom.intersects(noduleGeom)) {
            intersection = tissueGeom.intersection(noduleGeom)
            nodularGeom = nodularGeom == null ? intersection : nodularGeom.union(intersection)
        }
    }
    if (nodularGeom != null && !nodularGeom.isEmpty()) {
        roi = GeometryTools.geometryToROI(nodularGeom, ImagePlane.getDefaultPlane())
        obj = PathObjects.createAnnotationObject(roi, getPathClass("Nodular Tissue"))
        def base = tissue.getName()
                    .replace(' Tissue After Subtraction', '') 
                    .trim()
        obj.setName(base + ' Nodular Tissue')
        nodularAnnotations << obj
    }
}

// Add new objects
addObjects(nonNodularAnnotations)
addObjects(nodularAnnotations)

fireHierarchyUpdate()
