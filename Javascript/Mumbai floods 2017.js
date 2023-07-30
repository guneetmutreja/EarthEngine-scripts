// **************** 2017 Mumbai floods ****************
// On August 29, severe floods affected Mumbai.
// Mumbai recorded 468 mm of rainfall in twelve hours,
// the highest in a day in August since 1997, according
// to data from the India Meteorological Department.
// ****************************************************

// adapted from Flood Mapping Guided Project
// by Ujaval Gandhi, Spatial Thoughts

var s1 = ee.ImageCollection("COPERNICUS/S1_GRD");
var gsw = ee.Image("JRC/GSW1_4/GlobalSurfaceWater");
var dem = ee.Image("WWF/HydroSHEDS/03VFDEM")
var districts = ee.FeatureCollection("users/guneetmutreja96/districts")

var before_start = '2017-07-20'
var before_end = '2017-08-10'
var after_start = '2017-08-25'
var after_end = '2017-09-5'

//Map.addLayer(districts)
var mumbai_city = districts.filter(ee.Filter.eq('ID_3', 1177));
var geometry = mumbai_city.geometry()
Map.addLayer(geometry, {color: 'grey'}, 'Mumbai city')
Map.centerObject(geometry, 12)

var filtered = s1
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
  .filter(ee.Filter.eq('resolution_meters', 10))
  .filter(ee.Filter.bounds(geometry))
  .select(['VV', 'VH'])
  

var beforeCollection = filtered.filter(ee.Filter.date(before_start, before_end))
var afterCollection = filtered.filter(ee.Filter.date(after_start, after_end))

var before = beforeCollection.mosaic().clip(geometry)
var after = afterCollection.mosaic().clip(geometry)

var addRatioBand = function(image){
  var ratioBand = image.select('VV').divide(image.select('VH')).rename('VV/VH')
  return image.addBands(ratioBand)
}

var beforeRGB = addRatioBand(before)
var afterRGB = addRatioBand(after)

var visParams = {
  min: [-25, -25, 0],
  max: [0, 0, 2]
}

Map.addLayer(beforeRGB, visParams, 'Before')
Map.addLayer(afterRGB, visParams, 'After')

var beforeFiltered = ee.Image(toDB(RefinedLee(toNatural(before))))
var afterFiltered = ee.Image(toDB(RefinedLee(toNatural(after))))

Map.addLayer(beforeFiltered, {min:-25, max:0}, 'BeforeFiltered')
Map.addLayer(afterFiltered, {min:-25, max:0}, 'AfterFiltered')

var smoothingRadius = 50;
var beforeSmoothened = beforeFiltered.focal_median(smoothingRadius, 'circle', 'meters')
var afterSmoothened = afterFiltered.focal_median(smoothingRadius, 'circle', 'meters')

Map.addLayer(beforeSmoothened, {min:-25, max:0}, 'BeforeSmoothened')
Map.addLayer(afterSmoothened, {min:-25, max:0}, 'AfterSmoothened')

var difference = afterSmoothened.divide(beforeSmoothened)
var diffThreshold = 1.25

var flooded = difference.gt(diffThreshold).rename("water").selfMask()
Map.addLayer(flooded, {min:0, max:1, palette:['green']}, 'Initial flood estimate')

var permanentWater = gsw.select('seasonality').gte(5).clip(geometry)
var flooded = flooded.updateMask(permanentWater)

var slopeThresh = 5
var terrain = ee.Algorithms.Terrain(dem);
var slope = terrain.select('slope')
var flooded = flooded.updateMask(slope.lt(slopeThresh))

var connectedPixelsThresh = 8
var connections = flooded.connectedPixelCount(25)
var flooded = flooded.updateMask(connections.gt(connectedPixelsThresh))

Map.addLayer(flooded, {min:0, max:1, palette:['blue']}, 'Flooded areas')

var totalArea = geometry.area().divide(10000)
print(totalArea)

var stats = flooded.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: geometry,
  scale: 10,
  maxPixels: 1e10,
  tileScale: 16
})

var floodedArea = ee.Number(stats.get('water')).divide(10000)
print(floodedArea)


// ********** Speckle filtering function ****************

// Functions to convert from/to dB
function toNatural(img) {
  return ee.Image(10.0).pow(img.select(0).divide(10.0));
}

function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}


// Refined Lee speckle filter as coded in the SNAP 3.0 S1TBX:
//
//   https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/RefinedLee.java
//
// by Guido Lemoine
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}

