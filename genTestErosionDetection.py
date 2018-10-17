import numpy as np
import random
import math
import copy
import scipy
from math import ceil, floor
import sys
#import tensorflow as tf
import datetime
import os
import ogr
import osr
import gdal
import gdalnumeric
import gisutils as utils
from skimage.draw import polygon
from skimage import io
# from tensorflow.python.framework import ops

from PIL import Image, ImageOps
from os import listdir
from skimage import img_as_float
from scipy import stats

def getVertices(geometry, geotransform):
    ring = geometry.GetGeometryRef(0)
    pX = []
    pY = []
    for i in range(ring.GetPointCount()): 
        lon, lat, z = ring.GetPoint(i)
        p = utils.CoordinateToPixel(geotransform, (lon,lat))
        pX.append(p[0])
        pY.append(p[1])

    return pX, pY

def FillImage(img, geometry, geoTransform, fill=1):
    c, r = getVertices(geometry,geoTransform)
    rows,cols = polygon(r,c)

    img[rows, cols] = fill

def createBoundingBox(geometry, anchor, size, geotransform):
    cols,rows = getVertices(geometry, geotransform)
    cols.sort()
    rows.sort()

    xmin, ymin, xmax, ymax = max(0, cols[0]-anchor[1]), max(0,rows[0]-anchor[0]), \
                            min(cols[-1]-anchor[1], size), min(rows[-1]-anchor[0], size)
    return (xmin,xmax,ymin,ymax)

def cropImg(data, point, outputfile, size, debug=False):
    if len(data.shape) > 2:
        patch = np.zeros((data.shape[0]-1, size, size), dtype=data.dtype) # Ignora a ultima banda (Infravermelho)
        
        #O ponto de referencia eh o centro do crop
        wd = int(floor(size/2))
        lMin = max(point[0] - wd, 0)
        lMax = min(point[0] + wd, data.shape[1])
        cMin = max(point[1] - wd, 0)
        cMax = min(point[1] + wd, data.shape[2])

        for b in range(3):
            patch[b, 0:lMax-lMin, 0:cMax-cMin] = data[b, lMin:lMax, cMin:cMax]
        patch = np.moveaxis(patch, 0, -1)

    else:
        patch = np.zeros((size, size), dtype=np.uint8) # Ignora a ultima banda (Infravermelho)
        
        #O ponto de referencia eh o centro do crop
        wd = int(floor(size/2))
        lMin = max(point[0] - wd, 0)
        lMax = min(point[0] + wd, data.shape[0])
        cMin = max(point[1] - wd, 0)
        cMax = min(point[1] + wd, data.shape[1])
        
        patch[0:lMax-lMin, 0:cMax-cMin] = data[lMin:lMax, cMin:cMax]*255

    if debug:
        print("##############")
        print(data.shape)
        print(point)
        print((lMin, lMax))
        print((cMin, cMax))
        print("##############")


    scipy.misc.imsave(outputfile, patch)

def createData(imgfile, gt, outputfolder, references, geotransform, size):
    print("Opening img: {}".format(imgfile))
    data = gdalnumeric.LoadFile(imgfile)
    print("Done opening img: {}".format(imgfile))

    if not os.path.exists(os.path.join(outputfolder,'JPEGImages')):
        os.mkdir(os.path.join(outputfolder,'JPEGImages'))
    if not os.path.exists(os.path.join(outputfolder,'Masks')):
        os.mkdir(os.path.join(outputfolder,'Masks'))

    imgfilesplit = os.path.splitext(os.path.split(imgfile)[-1])
    for num, point in enumerate(references):
        outputfileimg = os.path.join(outputfolder,'JPEGImages', imgfilesplit[0] + '_crop' + str(num) + '.jpg')
        outputfilemask = os.path.join(outputfolder,'Masks', imgfilesplit[0] + '_crop' + str(num) + '_mask.jpg')
        print("Creating files: \t{0}{1}".format(outputfileimg, outputfilemask))
        cropImg(data, point, outputfileimg, size)
        cropImg(gt, point, outputfilemask, size)
        
def _addPoints(list, p1, p2, step, geotransform):
    for t in np.arange(0., 1.0, step):
        px = (p2[0] - p1[0])*t + p1[0]
        py = (p2[1] - p1[1])*t + p1[1]

        pixel = utils.CoordinateToPixel(geotransform, (px,py)) # Col, Row
        list.append((pixel[1], pixel[0]))

def addPoints(list, obj, geotransform):

    def dist(p1, p2):
        return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

    if obj.GetGeometryName() == 'MULTILINESTRING':
        for l in range(obj.GetGeometryCount()):
            line = obj.GetGeometryRef(l)
            for i in range(line.GetPointCount() - 1):
                p1p = line.GetPoint(i)
                p2p = line.GetPoint(i+1)

                p1 = (p1p[0], p1p[1])
                p2 = (p2p[0], p2p[1])

                ppl = max(dist(p1,p2) / 200, 1.0)

                _addPoints(list, p1, p2, 1.0/float(ppl), geotransform)
    else:
        for i in range(obj.GetPointCount() - 1):
            p1p = obj.GetPoint(i)
            p2p = obj.GetPoint(i+1)

            p1 = (p1p[0], p1p[1])
            p2 = (p2p[0], p2p[1])

            ppl = max(dist(p1,p2) / 200, 1.0)

            _addPoints(list, p1, p2, 1.0/float(ppl), geotransform)
            
def createDetectionData(shapefile, shapelabels, railwayshapefile, imglist, outputfolder, size):
    print("\n\n ------------- Creating points of interest ------------------- \n\n")
    imgs = imglist

    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    rlwds = ogr.Open(railwayshapefile)
    rlwlayer = rlwds.GetLayer(0)
    
    imgs_files = {}


    shpref = layer.GetSpatialRef()
    rlwref = rlwlayer.GetSpatialRef()
    rlwtransform = osr.CoordinateTransformation(rlwref, shpref)

    # Determina os poligonos de erosao com base em um csv
    # erosionFeatures = []
    # with open(shapelabels) as file:
    #     file.readline()
    #     text = file.read().split('\r\n')
    #     for line in text:
    #         lsplit = line.split(',')
    #         if len(line) < 2:
    #             break

    #         fid, label = lsplit[0], lsplit[1]
    #         if label == 'Erosao':
    #             erosionFeatures.append(int(fid))

    # Checa os poligonos que estao a menos de 5,1m do shp da ferrovia
    closetorlw = [False] * len(layer)
    for fid, f in enumerate(layer):
        geometry = f.GetGeometryRef()
        d = float('inf')
        rlwlayer.ResetReading()
        for feature in rlwlayer:
            railway = feature.GetGeometryRef()
            railway.Transform(rlwtransform)
            d = min(d, geometry.Distance(railway))
            
            if geometry.Distance(railway) <= 5.1:
                closetorlw[fid] = True
                break

    for num, file in enumerate(imgs):
        points_img = []
        references = []
        imgs_files[num] = file
        img = gdal.Open(file)
        geot  = img.GetGeoTransform()
        xAxis = img.RasterXSize # Max columns
        yAxis = img.RasterYSize # Max rows
        
        # Para fazer transformacoes de coordenadas nos shapefiles
        imgref = osr.SpatialReference(wkt = img.GetProjectionRef())
        transform = osr.CoordinateTransformation(shpref, imgref)
        rlwtransform = osr.CoordinateTransformation(rlwref, imgref)    

        # Extent da Imagem // Usado para verificar se poligonos estao sobre uma imagem
        ext = utils.GetExtentGeometry(geot, xAxis, yAxis)
        ext.FlattenTo2D()
        
        # Ground Truth // Matriz e Arquivos
        groundTruth = np.zeros((yAxis, xAxis), dtype='bool8')
        gtname = os.path.splitext(os.path.split(file)[-1])
        if not os.path.isdir(os.path.join(outputfolder, "GroundTruths")):
            os.mkdir(os.path.join(outputfolder, "GroundTruths"))
        gtfile = os.path.join(outputfolder, "GroundTruths", "mask_" + gtname[0] + ".png") 

        # Percorre as linhas do shape da ferrovia e adiciona os pontos centrais que irao gerar os crops 
        rlwlayer.ResetReading()
        for feature in rlwlayer:
            railway = feature.GetGeometryRef()
            railway.Transform(rlwtransform)
            if ext.Intersect(railway):
                intersection = ext.Intersection(railway)
                addPoints(references, intersection, geot)

        # Adiciona os poligonos de erosao no groundTruth
        layer.ResetReading()
        for fid, feature in enumerate(layer):
            geometry = feature.GetGeometryRef()
            geometry.Transform(transform)
            if ext.Intersect(geometry):
                intersection = geometry.Intersection(ext)
                FillImage(groundTruth, intersection, geot)

        # Cria os crops tanto da imagem quanto do groundTruth e salva o arquivo de groundTruth
        createData(file, groundTruth, outputfolder, references, geot, size)
        scipy.misc.imsave(gtfile, groundTruth*255)
        # np.save(os.path.join(outputfolder,"referencepoints.npy"), np.asarray(references))

def printParams(listParams):
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for i in range(1, len(sys.argv)):
        print(listParams[i - 1] + '= ' + sys.argv[i])
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

def main():
    list_params = ['shapefilePath', 'shapefileLabels', 'railwayShp', 'outsize', 
                   'output_path(for images, xml)', 'image list'
                   ]
    if len(sys.argv) < len(list_params) + 1:
        sys.exit('Usage: ' + sys.argv[0] + ' ' + ' '.join(list_params))
    printParams(list_params)

    # Input data 
    index = 1
    # shapefile path
    shapefilepath = sys.argv[index]
    index += 1
    # shapefile labels
    shapefilelabels = sys.argv[index]
    index += 1
    # railway shapefile path
    railwayshape = sys.argv[index]
    index += 1
    # output imgsize
    size = int(sys.argv[index])
    index += 1
    # output path
    outputpath = sys.argv[index]
    index += 1
    # image list
    #imglist = sys.argv[index].split(',')
    #index += 1

    folderlist = sys.argv[index].split(',')
    for i,a in enumerate(folderlist):
        imglist = []
        app = os.path.split(os.path.split(a)[0])[-1]
        for f in os.listdir(a):
            #print(os.path.join(a, f))
            _imglist = os.listdir(os.path.join(a,f))
            #print(_imglist)
            for x in _imglist:
                if ".img" in x and ".xml" not in x:
                    imglist.append(os.path.join(a,f,x))
        #print("Imagens : {}".format(imglist))
        if not os.path.exists(os.path.join(outputpath,app,'NewDetectionDataset')):
            os.mkdir(os.path.join(outputpath,app,'NewDetectionDataset'))
        createDetectionData(shapefilepath, shapefilelabels, railwayshape, imglist, os.path.join(outputpath,app,'NewDetectionDataset'), size)


    #createDetectionData(shapefilepath, shapefilelabels, railwayshape, imglist, outputpath, size)
    #createDetectionData(railwayshape, imglist, outputpath, ips, size)

if __name__ == '__main__':
    main()