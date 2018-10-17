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
                            min(cols[-1]-anchor[1], size), min(rows[-1]-anchor[0], size)
    return (xmin,xmax,ymin,ymax)

def cropImg(data, point, outputfile, size):
    patch = np.zeros((data.shape[0]-1, size, size), dtype=data.dtype) # Ignora a ultima banda (Infravermelho)
    
    #O ponto de referencia eh o centro do crop
    wd = int(floor(size/2))
    lMin = max(point[0] - wd, 0)
    lMax = min(point[0] + wd, data.shape[1])
    cMin = max(point[1] - wd, 0)
    cMax = min(point[1] + wd, data.shape[2])

    patch = data[0:3, lMin:lMax, cMin:cMax]

    patch = np.moveaxis(patch, 0, -1)
    scipy.misc.imsave(outputfile, patch[:,:,:])

def createData(imgfile, outputfolder, references, geotransform, size):
    print "Opening img: {}".format(imgfile)
    data = gdalnumeric.LoadFile(imgfile)
    print "Done opening img: {}".format(imgfile)

    imgfilesplit = os.path.splitext(os.path.split(imgfile)[-1])
    for num, point in enumerate(references):
        outputfileimg = os.path.join(outputfolder,'ValidationImages', imgfilesplit[0] + '_crop' + str(num) + '.jpg')
        print "Creating files: \t{0}".format(outputfileimg)
        cropImg(data, point, outputfileimg, size)
        
def _addPoints(list, p1, p2, step, geotransform):
    for t in np.arange(0., 1.0, step):
        px = (p2[0] - p1[0])*t + p1[0]
        py = (p2[1] - p1[1])*t + p1[1]

        pixel = utils.CoordinateToPixel(geotransform, (px,py)) # Col, Row
        list.append((pixel[1], pixel[0]))

def addPoints(list, obj, ppl, geotransform):
    if obj.GetGeometryName() == 'MULTILINESTRING':
        for l in range(obj.GetGeometryCount()):
            line = obj.GetGeometryRef(l)
            for i in range(line.GetPointCount() - 1):
                p1p = line.GetPoint(i)
                p2p = line.GetPoint(i+1)

                p1 = (p1p[0], p1p[1])
                p2 = (p2p[0], p2p[1])

                _addPoints(list, p1, p2, 1.0/float(ppl), geotransform)
    else:
        for i in range(obj.GetPointCount() - 1):
            p1p = obj.GetPoint(i)
            p2p = obj.GetPoint(i+1)

            p1 = (p1p[0], p1p[1])
            p2 = (p2p[0], p2p[1])

            _addPoints(list, p1, p2, 1.0/float(ppl), geotransform)
            
def createDetectionData(railwayshapefile, imglist, outputfolder, imagespersegment, size):
    print "\n\n ------------- Creating points of interest ------------------- \n\n"
    imgs = imglist

    rlwds = ogr.Open(railwayshapefile)
    rlwlayer = rlwds.GetLayer(0)
    references = []
    imgs_files = {}

    rlwref = rlwlayer.GetSpatialRef()

    for num, file in enumerate(imgs):
        points_img = []
        imgs_files[num] = file
        img = gdal.Open(file)
        geot  = img.GetGeoTransform()
        xAxis = img.RasterXSize # Max columns
        yAxis = img.RasterYSize # Max rows
        
        imgref = osr.SpatialReference(wkt = img.GetProjectionRef())
        rlwtransform = osr.CoordinateTransformation(rlwref, imgref)    

        ext = utils.GetExtentGeometry(geot, xAxis, yAxis)
        ext.FlattenTo2D()
        
        rlwlayer.ResetReading()
        for feature in rlwlayer:
            railway = feature.GetGeometryRef()
            railway.Transform(rlwtransform)
            if ext.Intersect(railway):
                intersection = ext.Intersection(railway)
                addPoints(references, intersection, imagespersegment, geot)

        createData(file, outputfolder, references, geot, size)

def printParams(listParams):
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    for i in xrange(1, len(sys.argv)):
        print listParams[i - 1] + '= ' + sys.argv[i]
    print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

def main():
    list_params = ['railwayShp', 'imgs per segment', 'outsize', 'output_path(for images, xml)',
                   'imgList(image files separated by comma)',
                  ]
    if len(sys.argv) < len(list_params) + 1:
        sys.exit('Usage: ' + sys.argv[0] + ' ' + ' '.join(list_params))
    printParams(list_params)

    # Input data 
    index = 1
    # railway shapefile path
    railwayshape = sys.argv[index]
    index += 1
    # imgs per segment
    ips = int(sys.argv[index])
    index += 1
    # output imgsize
    size = int(sys.argv[index])
    index += 1
    # output path
    outputpath = sys.argv[index]
    index += 1
    # image list
    imglist = sys.argv[index].split(',')
    index += 1

    createDetectionData(railwayshape, imglist, outputpath, ips, size)

if __name__ == '__main__':
    main()