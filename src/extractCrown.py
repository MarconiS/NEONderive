
def h5refl2array(refl_filename, epsg):
    refl_filename = full_path
    hdf5_file = h5py.File(refl_filename, 'r')
    file_attrs_string = str(list(hdf5_file.items()))
    file_attrs_string_split = file_attrs_string.split("'")
    sitename = file_attrs_string_split[1]

    #Extract the reflectance & wavelength datasets
    reflArray = hdf5_file[sitename]['Reflectance']
    wavelengths =reflArray['Reflectance_Data'].value

    # Create dictionary containing relevant metadata information
    metadata = {}
    metadata['mapInfo'] = reflArray['Metadata']['Coordinate_System']['Map_Info'].value
    metadata['wavelength'] = reflArray['Metadata']['Spectral_Data']['Wavelength'].value
    metadata['shape'] = wavelengths.shape
    #Extract no data value & scale factor
    metadata['noDataVal'] = float(reflArray['Reflectance_Data'].attrs['Data_Ignore_Value'])
    metadata['scaleFactor'] = float(reflArray['Reflectance_Data'].attrs['Scale_Factor'])
    #metadata['interleave'] = reflData.attrs['Interleave']
    metadata['bad_band_window1'] = np.array([1340, 1445])
    metadata['bad_band_window2'] = np.array([1790, 1955])
    metadata['epsg'] = str(epsg)

    mapInfo_string = str(metadata['mapInfo']);
    mapInfo_split = mapInfo_string.split(",")
    mapInfo_split

    # Extract the resolution & convert to floating decimal number
    metadata['res'] = {}
    metadata['res']['pixelWidth'] = float(mapInfo_split[5])
    metadata['res']['pixelHeight'] = float(mapInfo_split[6])
    # Extract the upper left-hand corner coordinates from mapInfo
    xMin = float(mapInfo_split[3])  # convert from string to floating point number
    yMax = float(mapInfo_split[4])

    # Calculate the xMax and yMin values from the dimensions
    xMax = xMin + (metadata['shape'][1] * metadata['res']['pixelWidth'])  # xMax = left edge + (# of columns * resolution)",
    yMin = yMax - (metadata['shape'][0] * metadata['res']['pixelHeight'])  # yMin = top edge - (# of rows * resolution)",
    metadata['extent'] = (xMin, xMax, yMin, yMax)  # useful format for plotting
    metadata['ext_dict'] = {}
    metadata['ext_dict']['xMin'] = xMin
    metadata['ext_dict']['xMax'] = xMax
    metadata['ext_dict']['yMin'] = yMin
    metadata['ext_dict']['yMax'] = yMax
    hdf5_file.close

    return metadata, wavelengths


# In[2]:

def stack_subset_bands(reflArray, reflArray_metadata, bands, clipIndex):
    subArray_rows = clipIndex['yMax'] - clipIndex['yMin']
    subArray_cols = clipIndex['xMax'] - clipIndex['xMin']

    stackedArray = np.zeros((subArray_rows, subArray_cols, len(bands)), dtype=np.int16)
    band_clean_dict = {}
    band_clean_names = []

    for i in range(len(bands)):
        band_clean_names.append("b" + str(bands[i]) + "_refl_clean")
        band_clean_dict[band_clean_names[i]] = subset_clean_band(reflArray, reflArray_metadata, clipIndex, bands[i])
        stackedArray[..., i] = band_clean_dict[band_clean_names[i]]

    return stackedArray


# In[3]:

def subset_clean_band(reflArray, reflArray_metadata, clipIndex, bandIndex):
    bandCleaned = reflArray[clipIndex['yMin']:clipIndex['yMax'], clipIndex['xMin']:clipIndex['xMax'],
                  bandIndex - 1].astype(np.int16)

    return bandCleaned


# In[4]:

def array2raster(newRaster, reflBandArray, reflArray_metadata, extent, ras_dir, epsg):
    NP2GDAL_CONVERSION = {
        "uint8": 1,
        "int8": 1,
        "uint16": 2,
        "int16": 3,
        "uint32": 4,
        "int32": 5,
        "float32": 6,
        "float64": 7,
        "complex64": 10,
        "complex128": 11,
    }

    pwd = os.getcwd()
    os.chdir(ras_dir)
    cols = reflBandArray.shape[1]
    rows = reflBandArray.shape[0]
    bands = reflBandArray.shape[2]
    pixelWidth = float(reflArray_metadata['res']['pixelWidth'])
    pixelHeight = -float(reflArray_metadata['res']['pixelHeight'])
    originX = extent['xMin']
    originY = extent['yMax']

    driver = gdal.GetDriverByName('GTiff')
    gdaltype = NP2GDAL_CONVERSION[reflBandArray.dtype.name]
    outRaster = driver.Create(newRaster, cols, rows, bands, gdaltype)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # outband = outRaster.GetRasterBand(1)
    # outband.WriteArray(reflBandArray[:,:,x])
    for band in range(bands):
        outRaster.GetRasterBand(band + 1).WriteArray(reflBandArray[:, :, band])

    outRasterSRS = osr.SpatialReference()
    #outRasterSRS.ImportFromEPSG(reflArray_metadata['epsg'])
    #outRasterSRS.ExportToWkt()
    outRasterSRS.ImportFromEPSG(epsg)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outRaster.FlushCache()
    os.chdir(pwd)



# In[5]:

def calc_clip_index(clipExtent, h5Extent, xscale=1, yscale=1):
    h5rows = h5Extent['yMax'] - h5Extent['yMin']
    h5cols = h5Extent['xMax'] - h5Extent['xMin']

    ind_ext = {}
    ind_ext['xMin'] = round((clipExtent['xMin'] - h5Extent['xMin']) / xscale)
    ind_ext['xMax'] = round((clipExtent['xMax'] - h5Extent['xMin']) / xscale)
    ind_ext['yMax'] = round(h5rows - (clipExtent['yMin'] - h5Extent['yMin']) / yscale)
    ind_ext['yMin'] = round(h5rows - (clipExtent['yMax'] - h5Extent['yMin']) / yscale)

    return ind_ext




import numpy as np
import h5py
import gdal, osr
#import matplotlib.pyplot as plt
import sys
import ogr, os

full_path =sys.argv[1]
chm_path = sys.argv[2]
itc_id = sys.argv[3]
itc_xmin = float(sys.argv[4])
itc_xmax = float(sys.argv[5])
itc_ymin = float(sys.argv[6])
itc_ymax = float(sys.argv[7])
epsg = str(sys.argv[8])
wd = sys.argv[9]

print(itc_id, itc_xmin, itc_xmax, itc_ymin, itc_ymax, epsg)

refl_md, refl = h5refl2array(full_path, epsg = epsg)

rgb = np.r_[0:425]
rgb = np.delete(rgb, np.r_[419:425])
rgb = np.delete(rgb, np.r_[283:315])
rgb = np.delete(rgb, np.r_[192:210])
xmin, xmax, ymin, ymax = refl_md['extent']
print(xmin, xmax, ymin, ymax)

clipExtent = {}
clipExtent['xMin'] = itc_xmin
clipExtent['yMin'] = itc_ymin
clipExtent['yMax'] = itc_ymax
clipExtent['xMax'] = itc_xmax
print(clipExtent)

subInd = calc_clip_index(clipExtent, refl_md['ext_dict'])
subInd['xMax'] = int(subInd['xMax'])
subInd['xMin'] = int(subInd['xMin'])
subInd['yMax'] = int(subInd['yMax'])
subInd['yMin'] = int(subInd['yMin'])
print(subInd)

refl = refl[(subInd['yMin']):subInd['yMax'], (subInd['xMin']):subInd['xMax'], :]
refl.shape
print(refl.shape)

chmExtent = {}
chmExtent['xMin'] = int(itc_xmin/1000) * 1000
chmExtent['yMin'] = int(itc_ymin/1000) * 1000
chmExtent['yMax'] = int(itc_ymax/1000 +1) * 1000
chmExtent['xMax'] = int(itc_xmax/1000+1) * 1000
print(chmExtent)

chmInd = calc_clip_index(clipExtent, chmExtent)
chmInd['xMax'] = int(chmInd['xMax'])
chmInd['xMin'] = int(chmInd['xMin'])
chmInd['yMax'] = int(chmInd['yMax'])
chmInd['yMin'] = int(chmInd['yMin'])

print(chmInd)
chm = gdal.Open(chm_path).ReadAsArray()
chm = chm[chmInd['yMin']:chmInd['yMax'], chmInd['xMin']:chmInd['xMax']]
chm = (chm *100).astype(np.int16)

subArray_rows = subInd['yMax'] - subInd['yMin']
subArray_cols = subInd['xMax'] - subInd['xMin']
hcp = np.zeros((subArray_rows, subArray_cols, len(rgb)+1), dtype=np.int16)

band_clean_dict = {}
band_clean_names = []
for i in range(len(rgb)):
    if i == 0:
        band_clean_names.append("b" + 'chm' + "_refl_clean")
        band_clean_dict[band_clean_names[i]] = chm.astype(np.int16)
        hcp[..., i] = band_clean_dict[band_clean_names[i]]
    else:
        band_clean_names.append("b" + str(rgb[i]) + "_refl_clean")
        band_clean_dict[band_clean_names[i]] = refl[:, :, rgb[i]].astype(np.int16)
        hcp[..., i] = band_clean_dict[band_clean_names[i]]

sub_meta = refl_md
ii = str(itc_id) + '.tif'
ras_dir = wd+"/out/AOP/plot/itcTiff/"
array2raster(ii, hcp, sub_meta, clipExtent, ras_dir, epsg)


