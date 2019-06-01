import rasterio
import numpy as np
from PIL import Image
from rasterio.enums import Resampling
import cv2
import glob
from rasterio.mask import mask
import sys

import geopandas as gpd
from fiona.crs import from_epsg
from shapely.geometry import box

def apply_brightness_contrast(input_img, brightness = 0, contrast = 0):

    if brightness != 0:
        if brightness > 0:
            shadow = brightness
            highlight = 255
        else:
            shadow = 0
            highlight = 255 + brightness
        alpha_b = (highlight - shadow)/255
        gamma_b = shadow

        buf = cv2.addWeighted(input_img, alpha_b, input_img, 0, gamma_b)
    else:
        buf = input_img.copy()

    if contrast != 0:
        f = 131*(contrast + 127)/(127*(131-contrast))
        alpha_c = f
        gamma_c = 127*(1-f)

        buf = cv2.addWeighted(buf, alpha_c, buf, 0, gamma_c)

    return buf

EPSG_NEON = {
    "SCBI" : "32617",
    "STEI" : "32616",
    "CHEQ" : "32615",
    "KONZ" : "32614",
    "GRSM" : "32616",
    "ORNL" : "32616",
    "TALL" : "32616",
    "MOAB" : "32612",
    "OSBS" : "32617",
    "MLBS" : "32617",
    "LENO" : "32616",
    "UKFS" : "32615",
    "SERC" : "32618", 
    "DSNY" : "32617",
    "HARV" : "32618"
    
}
    
path = sys.argv[1]
site = sys.argv[2]
list_chm = (glob.glob(path+"CHM/*"+site+"*.tif"))

for ii in range(len(list_chm)):
    try:
        chm = rasterio.open(list_chm[ii])
        crs = rasterio.crs.CRS({"init": "epsg:"+EPSG_NEON.get(site)})    # or whatever CRS you know the image is in    
        data_chm = chm.read(1)
        hps_file = (glob.glob(path+"itcTiff/*"+list_chm[ii][-12:-4]+"*.tif"))
        hps = rasterio.open(hps_file[0])
        hps_b1 = hps.read(90)
        hps_b2 = hps.read(19)
        hps_b3 = hps.read(58)
        if data_chm.shape[1] > data_chm.shape[0]:
            data_chm = data_chm[:,:-1]           
        elif data_chm.shape[1] < data_chm.shape[0]:
            data_chm = data_chm[:-1,:]

        with rasterio.open(path+"temp/band_90.tif") as dataset:
            band_90 = dataset.read(
                #out_shape=(dataset.count, dataset.height * 2, dataset.width * 2),
                out_shape=(dataset.count, data_chm.shape[1], data_chm.shape[0]),
                resampling=Resampling.bilinear
            ) 
        dataset.close()


        with rasterio.open(path+"temp/band_19.tif") as dataset:
            band_19 = dataset.read(
                out_shape=(dataset.count,  data_chm.shape[1], data_chm.shape[0]),
                resampling=Resampling.bilinear
            )
        dataset.close()


        with rasterio.open(path+"temp/band_58.tif") as dataset:
            band_58 = dataset.read(
                out_shape=(dataset.count,  data_chm.shape[1], data_chm.shape[0]),
                resampling=Resampling.bilinear
            )
        dataset.close()

        new_dataset = rasterio.open(
            path+"temp/band_90.tif",
            'w',
            driver='GTiff',
            height=hps_b1.shape[0],
            width=hps_b1.shape[1],
            count=1,
            dtype=hps_b1.dtype,
            crs=chm.crs,
            transform=hps.transform,)
        new_dataset.write(hps_b1, 1)
        new_dataset.close()

        new_dataset = rasterio.open(
            path+"temp/band_19.tif",
            'w',
            driver='GTiff',
            height=hps_b2.shape[0],
            width=hps_b2.shape[1],
            count=1,
            dtype=hps_b2.dtype,
            crs=chm.crs,
            transform=hps.transform,)
        new_dataset.write(hps_b2, 1)
        new_dataset.close()

        new_dataset = rasterio.open(
            path+"temp/band_58.tif",
            'w',
            driver='GTiff',
            height=hps_b3.shape[0],
            width=hps_b3.shape[1],
            count=1,
            dtype=hps_b3.dtype,
            crs=chm.crs,
            transform=hps.transform,)
        new_dataset.write(hps_b3, 1)
        new_dataset.close()

        up_sample90 = ((band_90 - band_90.min())/(band_90.max() - band_90.min()) * 255).astype('uint8')
        up_sample90 = np.reshape(up_sample90, (band_90.shape[2], band_90.shape[1]))

        up_sample19 = ((band_19 - band_19.min())/(band_19.max() - band_19.min()) * 255).astype('uint8')
        up_sample19 = np.reshape(up_sample19, (band_90.shape[2], band_90.shape[1]))

        up_sample58 = ((band_58 - band_58.min())/(band_58.max() - band_58.min()) * 255).astype('uint8')
        up_sample58 = np.reshape(up_sample58, (band_90.shape[2], band_90.shape[1]))
        red = Image.fromarray(up_sample90).convert('L')
        blue =  Image.fromarray(up_sample19).convert('L')
        green = Image.fromarray(up_sample58).convert('L')
        hue = Image.fromarray(data_chm).convert('L')

        out = Image.merge("RGBA", (red, green, blue,hue))
        out.save(path+"/false_col/4channels_"+list_chm[ii][-12:-4]+".tif")
        im_gray = cv2.imread(path+"/false_col/4channels_"+list_chm[ii][-12:-4]+".tif", cv2.IMREAD_GRAYSCALE)
        u=apply_brightness_contrast(im_gray, 30, 30)
        new_dataset = rasterio.open(
            path+"/false_col/grayscale_"+list_chm[ii][-12:-4]+".tif",
            'w',
            driver='GTiff',
            height=u.shape[0],
            width=u.shape[1],
            count=1,
            dtype=u.dtype,
            crs=chm.crs,
            transform=chm.transform,)
        new_dataset.write(u, 1)
        new_dataset.close()
    except:
        print(ii)
