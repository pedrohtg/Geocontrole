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

def createBoundingBox(geometry, anchor, size, geotransform):
    cols,rows = getVertices(geometry, geotransform)
    cols.sort()
    rows.sort()

    xmin, ymin, xmax, ymax = max(0, cols[0]-anchor[1]), max(0,rows[0]-anchor[0]), \
                            min(cols[-1]-anchor[1], size-1), min(rows[-1]-anchor[0], size-1)
    return (xmin,xmax,ymin,ymax)

def cropImg(data, point, outputfile, size):
    patch = np.zeros((data.shape[0]-1, size, size), dtype=data.dtype) # Ignora a ultima banda (Infravermelho)
    
    #O ponto de referencia eh o centro do crop
    wd = int(floor(size/2))
    lMin = max(point[0] - wd, 0)
    lMax = min(point[0] + wd, data.shape[1])
    cMin = max(point[1] - wd, 0)
    cMax = min(point[1] + wd, data.shape[2])
    print "##############"
    print data.shape
    print point
    print (lMin, lMax)
    print (cMin, cMax)
    print "##############"
    # Faz reflexao se os limites saem da img
    for b in range(data.shape[0]-1):
        for l in range(lMin, lMax):
            for c in range(cMin, cMax):
                patch[b][l - lMin][c - cMin] = data[b][l][c]

    patch = np.moveaxis(patch, 0, -1)
    print patch[:3,:3,:]
    scipy.misc.imsave(outputfile, patch[:,:,:]) 
    return (max(0, lMin), max(0, cMin))

def genXml(outputfilename, objects, imgname, size):
    text = "<annotation>\n"
    text+= "\t<folder>XXX</folder>\n\t<filename>{0}</filename>\n".format(imgname)
    text+= "\t<source>\n\t\t<database>Erosao Ferrovia TransNordestina</database>\n\t\t<annotation>PATREO</annotation>\n\t\t<image>googleearth</image>\n\t</source>\n"
    text+= "\t<owner>\n\t\t<name>patreo</name>\n\t</owner>\n"
    text+= "\t<size>\n\t\t<width>{0}</width>\n\t\t<height>{1}</height>\n\t\t<depth>3</depth>\n\t</size>\n".format(size, size)
    text+= "\n\t<segmented>1</segmented>\n"

    for obj in objects:
        objText = "\t<object>\n\t\t<name>bridge</name>\n\t\t<pose>Front</pose>\n\t\t<truncated>0</truncated>\n\t\t<difficult>0</difficult>\n\t\t"
        objText+= "<bndbox>\n\t\t\t<xmin>{0}</xmin>\n\t\t\t<ymin>{1}</ymin>\n\t\t\t<xmax>{2}</xmax>\n\t\t\t<ymax>{3}</ymax>\n\t\t</bndbox>\n\t</object>".format(obj[0],obj[1],obj[2],obj[3])
        text += objText

    text += "\n</annotation>"

    out = open(outputfilename, 'w')
    out.write(text)
    out.close()

def createData(imgfile, outputfolder, references, geotransform, size):
    print "Opening img: {}".format(imgfile)
    data = gdalnumeric.LoadFile(imgfile)
    print "Done opening img: {}".format(imgfile)

    imgfilesplit = os.path.splitext(os.path.split(imgfile)[-1])
    for num, ref in enumerate(references):
        outputfileimg = os.path.join(outputfolder,'JPEGImages', imgfilesplit[0] + '_crop' + str(num) + '.jpg')
        outputfilexml = os.path.join(outputfolder,'Annotations', imgfilesplit[0] + '_crop' + str(num) + '.xml')
        print "Creating files: \t{0}\n\t{1}".format(outputfileimg,outputfilexml)
        point = ref[0]
        obj = ref[1]
        upleftCorner = cropImg(data, point, outputfileimg, size)
        bb = createBoundingBox(obj, upleftCorner, size, geotransform)        
        genXml(outputfilexml, [bb], imgfilesplit[1]+'_crop'+str(num), size)

def createDetectionData(shapefile, shapelabels, railwayshapefile, imglist, outputfolder, size):
    print "\n\n ------------- Creating points of interest ------------------- \n\n"
    imgs = imglist

    ds = ogr.Open(shapefile)
    rlwds = ogr.Open(railwayshapefile)
    layer = ds.GetLayer(0)
    rwllayer = rlwds.GetLayer(0)
    references = []
    imgs_files = {}

    shpref = layer.GetSpatialRef()
    rlwref = rwllayer.GetSpatialRef()
    rlwtransform = osr.CoordinateTransformation(rlwref, shpref)

    # Determina os poligonos de erosao com base em um csv
    erosionFeatures = []
    with open(shapelabels) as file:
        file.readline()
        text = file.read().split('\r\n')
        for line in text:
            lsplit = line.split(',')
            if len(line) < 2:
                break

            fid, label = lsplit[0], lsplit[1]
            if label == 'Erosao':
                erosionFeatures.append(int(fid))

    #print erosionFeatures
    # Checa os poligonos que estao a menos de 5,1m do shp da ferrovia
    closetorlw = [False] * len(layer)
    for fid, f in enumerate(layer):
        if fid in erosionFeatures:
            geometry = f.GetGeometryRef()
            d = float('inf')
            rwllayer.ResetReading()
            for feature in rwllayer:
                railway = feature.GetGeometryRef()
                railway.Transform(rlwtransform)
                d = min(d, geometry.Distance(railway))
                
                if geometry.Distance(railway) <= 5.1:
                    closetorlw[fid] = True
                    break
            #print "{0}: dist {1}".format(fid, d)
    k = set([])
    for num, file in enumerate(imgs):
        points_img = []
        references = []
        imgs_files[num] = file
        img = gdal.Open(file)
        geot  = img.GetGeoTransform()
        xAxis = img.RasterXSize # Max columns
        yAxis = img.RasterYSize # Max rows
        
        imgref = osr.SpatialReference(wkt = img.GetProjectionRef())
        transform = osr.CoordinateTransformation(shpref, imgref)
        rlwtransform = osr.CoordinateTransformation(rlwref, imgref)    

        ext = utils.GetExtentGeometry(geot, xAxis, yAxis)
        ext.FlattenTo2D()
        
        # Adiciona os centroids dos poligonos de erosao a lista de referencias
        layer.ResetReading()
        for fid, feature in enumerate(layer):
            if fid in erosionFeatures and closetorlw[fid]:
                geometry = feature.GetGeometryRef()
                geometry.Transform(transform)
                if ext.Intersect(geometry):
                    k.add(fid)
    print "#####"
    print size
    print len(k)
    print k
    print "#####"

def printParams(listParams):
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    for i in xrange(1, len(sys.argv)):
        print listParams[i - 1] + '= ' + sys.argv[i]
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

def main():
    list_params = ['shapefilePath', 'shapefileLabels', 'railwayShp', 'outsize', 
                   'output_path(for images, xml)', 'imgfolders',
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

    folderlist = sys.argv[index].split(',')
    for i,a in enumerate(folderlist):
        imglist = []
        print (i,a)
        for f in os.listdir(a):
            #print os.path.join(a, f)
            _imglist = os.listdir(os.path.join(a,f))
            #print _imglist
            for x in _imglist:
                if ".img" in x and ".xml" not in x:
                    imglist.append(os.path.join(a,f,x))
        #print "Imagens : {}".format(imglist)
        createDetectionData(shapefilepath, shapefilelabels, railwayshape, imglist, outputpath, i)
    index += 1

    #createDetectionData(shapefilepath, shapefilelabels, railwayshape, imglist, outputpath, size)

if __name__ == '__main__':
    main()