// **************** NDVI ****************
// The script shows three methods to
// compute NDVI for an AOI.
// **************************************


var geometry = ee.Geometry.Point(77.27898093808096,28.56234251503353)

// Compute NDVI 3 ways.
var landsat = ee.ImageCollection("LANDSAT/LC08/C01/T1")
    .filterDate('2016-01-01', '2017-01-01')
    .filterBounds(geometry)

var composite = ee.Algorithms.Landsat.simpleComposite({
  collection: landsat,
  asFloat: true
})

// Method 1)
var b5 = composite.select("B5")
var b4 = composite.select("B4")
var ndvi_1 = b5.subtract(b4).divide(b5.add(b4))

// Method 2)
var ndvi_2 = composite.normalizedDifference(["B5", "B4"])

// Method 3)
var ndvi_3 = composite.expression("(b5 - b4) / (b5 + b4)", {
    b5: composite.select("B5"),
    b4: composite.select("B4")
})

Map.addLayer(ndvi_1, {min:0, max:1} , "NDVI")
Map.centerObject(geometry, 11) 


// NDVI over a collection
var ndvi = landsat.map(function(image) {
  var result = image.normalizedDifference(["B5", "B4"]).rename("ndvi")
  return image.addBands(result);
})

print(ndvi)

var maxNDVI = ndvi.max();
print(maxNDVI)
Map.addLayer(maxNDVI, {bands: "ndvi", min:0, max:1}, "Max-NDVI")
Map.centerObject(geometry, 11)

// A chart of NDVI over time.
print(ui.Chart.image.doySeries(ndvi.select('ndvi'), geometry, ee.Reducer.mean(), 30))