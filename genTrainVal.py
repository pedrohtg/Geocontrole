import os
import sys
import hashlib
import shutil
import random
import gdal
import numpy as np
import ogr
import osr
import scipy.misc
from skimage.draw import polygon
import gisutils as utils
import math
import numpy as np
from skimage import io
from skimage.measure import label
from glob import glob


POINTS_DICT = {}


crop_count = 0

#python genTrainVal.py D:\tcu\ImagensFerrovia\Parte1\IMGS\ D:\tcu\MissoesProntas\Shapes\LinhaFerrovia_Playades.shp D:\tcu\MissoesProntas\Shapes\erosao\Problemas_Ferrovia.shp 448 train D:\tcu\gabriel\dataset\Parte1\ D:\tcu\MissoesProntas\Shapes\erosao\Problemas_Ferrovia.csv    

class BatchColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def print_params(list_params):
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    for i in range(1, len(sys.argv)):
        print((list_params[i - 1] + '= ' + sys.argv[i]))
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')


def save_xml(output_path, img_name, width, height, list_rect, offset_x, offset_y):
    f = open(os.path.join(output_path, img_name + ".xml"), "w+")

    f.write("<annotation> \n\
            <folder>dataset_pontes</folder> \n\
            <filename>{0}.jpg</filename> \n\
            <source> \n\
                <database>XXX</database> \n\
                <annotation>patreo</annotation> \n\
                <image>googleearth</image> \n\
            </source> \n\
            <owner> \n\
                <name>patreo</name> \n\
            </owner> \n\
            <size> \n\
                <width>{1}</width> \n\
                <height>{2}</height> \n\
                <depth>3</depth> \n\
            </size> \n\
            <segmented>0</segmented> \n".format(img_name, width, height))

    for it, rect in enumerate(list_rect):
         f.write("	<object> \n\
                <name>bridge</name> \n\
                <pose>front</pose> \n\
                <truncated>0</truncated> \n\
                <difficult>0</difficult> \n\
                <bndbox> \n\
                    <xmin>{0}</xmin> \n\
                    <ymin>{1}</ymin> \n\
                    <xmax>{2}</xmax> \n\
                    <ymax>{3}</ymax> \n\
                </bndbox> \n\
            </object> \n".format(int(rect[0]-offset_x+1), int(rect[1]-offset_y+1), int(rect[2]-offset_x+1), int(rect[3]-offset_y+1))
                )

    f.write("</annotation> \n")


def is_between(a, b, c):
    # y = 1
    # x = 0
    # compare versus epsilon for floating point values, or != 0 if using integers
    crossproduct = (c[1] - a[1]) * (b[0] - a[0]) - (c[0] - a[0]) * (b[1] - a[1])
    if abs(crossproduct) > 0.001:
        return False

    dotproduct = (c[0] - a[0]) * (b[0] - a[0]) + (c[1] - a[1])*(b[1] - a[1])
    if dotproduct < 0:
        return False

    squaredlengthba = (b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1])
    if dotproduct > squaredlengthba:
        return False

    return True


def bounds_comparison(a):
    if a[1][0] > a[1][0] and a[1][1] == a[1][0]:
        return 1
    elif a[1][0] <= a[1][0]:
        return -1
    else:
        return -1


def sort_bounds(bounds):
    bounds = np.squeeze(np.asarray(bounds))

    bounds_list = []
    for i in range(bounds.shape[0]):
        # print bounds[i]
        bounds_list.append(bounds[i])

    bounds_list.sort(key=bounds_comparison)
    bounds_list = np.asarray(bounds_list)
    bounds_list = np.reshape(bounds_list, bounds.shape)
    # print bounds_list

    return bounds_list


def extract_railway_line_limits_from_shapefile(layer_railway, transform_railway_img, ext):
    # reset layer reading for new reading
    layer_railway.ResetReading()

    lims = []
    for feature in layer_railway:
        railway = feature.GetGeometryRef()
        railway.Transform(transform_railway_img)
        if ext.Intersect(railway):
            intersection = ext.Intersection(railway)
            if intersection.GetGeometryName() == 'MULTILINESTRING':
                for l in range(intersection.GetGeometryCount()):
                    line = intersection.GetGeometryRef(l)
                    for i in range(line.GetPointCount() - 1):
                        pnts = []
                        p1p = line.GetPoint(i)
                        p2p = line.GetPoint(i + 1)

                        p1 = (p1p[0], p1p[1])
                        p2 = (p2p[0], p2p[1])

                        pnts.append(p1)
                        pnts.append(p2)
                        lims.append(pnts)
            else:
                for i in range(intersection.GetPointCount() - 1):
                    pnts = []
                    p1p = intersection.GetPoint(i)
                    p2p = intersection.GetPoint(i + 1)

                    p1 = (p1p[0], p1p[1])
                    p2 = (p2p[0], p2p[1])

                    pnts.append(p1)
                    pnts.append(p2)
                    lims.append(pnts)
    return lims


def get_extent(gt, cols, rows):
    """
    Return list of corner coordinates from a geotransform

    :param gt: geotransform -- C{tuple/list}
    :param cols: number of columns in the dataset -- C{int}
    :param rows: number of rows in the dataset -- C{[float,...,float]}
    :return: coordinates of each corner
    """
    ext = []
    xarr = [0, cols]
    yarr = [0, rows]

    for px in xarr:
        for py in yarr:
            x = gt[0] + (px*gt[1])+(py*gt[2])
            y = gt[3] + (px*gt[4])+(py*gt[5])
            ext.append([x, y])
        yarr.reverse()
    return ext


def get_extent_geometry(gt, cols, rows):
    ext = get_extent(gt, cols, rows)

    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(ext[0][0], ext[0][1])
    ring.AddPoint(ext[0][0], ext[1][1])
    ring.AddPoint(ext[2][0], ext[1][1])
    ring.AddPoint(ext[2][0], ext[0][1])
    ring.AddPoint(ext[0][0], ext[0][1])
    raster_geometry = ogr.Geometry(ogr.wkbPolygon)
    raster_geometry.AddGeometry(ring)

    return raster_geometry


def coordinate_to_pixel(gt, c):
    """
    Return the pixel position from coordinates using a geotransform

    :param gt: geotransform -- C{tuple/list}
    :param c: the coordinate -- C{[float, float]}
    :return: position x,y of the pixel
    """
    c1 = c[0] - gt[0] if gt[1] > 0 else gt[0] - c[0]
    c2 = c[1] - gt[3] if gt[5] < 0 else gt[3] - c[1]

    if gt[2] == 0:
        x = c1/gt[1]
    if gt[4] == 0:
        y = c2/gt[5]

    if gt[2]*gt[4] != 0:
        # print "entrei"
        y = ((c2*gt[1])  - (gt[4]*c1)) / (-gt[2]*gt[4])
        x = (c1 - gt[2]) / gt[1]

    # print "gt: "+str(gt)
    # print "c: "+str(c)
    # print "x: "+str(x)
    # print "y: "+str(y)

    return (int(x), int(y))


def pixel_to_coordinate(gt, p):
    x = gt[0]+(p[0]*gt[1])+(p[1]*gt[2])
    y = gt[3]+(p[0]*gt[4])+(p[1]*gt[5])

    return (x, y)


def get_vertices(geometry, geo_transform):
    ring = geometry.GetGeometryRef(0)
    print(geometry.GetGeometryName)
    p_x = []
    p_y = []
    for i in range(ring.GetPointCount()):
        lon, lat, z = ring.GetPoint(i)
        p = coordinate_to_pixel(geo_transform, (lon, lat))
        p_x.append(p[0])
        p_y.append(p[1])

    return p_x, p_y


def check_boundary(ext, bounds):
    for k in list(POINTS_DICT.keys()):
        if POINTS_DICT[k] != tuple(bounds[1]) and POINTS_DICT[k] != tuple(bounds[0]) and \
                is_between(tuple(bounds[1]), tuple(bounds[0]), POINTS_DICT[k]):

            point1 = ogr.Geometry(ogr.wkbPoint)
            point1.AddPoint(k[2], k[3])
            if ext.Intersect(point1):
                return np.array([POINTS_DICT[k], bounds[0]])
                # print '1', POINTS_DICT[k], bounds[i, 0]
            else:
                return np.array([bounds[1], POINTS_DICT[k]])
                # print '2', bounds[i, 1], POINTS_DICT[k]
    return bounds


def fill_ground_truth(img, geometry, geo_transform, fill=255):
    c, r = get_vertices(geometry, geo_transform)
    rows, cols = polygon(r, c)

    img[rows, cols] = fill


def generate_mask_for_img(layer_erosion, transform_erosion_img, ext, geot_img, y_axis, x_axis):
    ground_truth = np.zeros((y_axis, x_axis), dtype=np.uint8)

    # Adiciona os pontos de erosao a pool de pontos
    layer_erosion.ResetReading()
    for fid, feature in enumerate(layer_erosion):
        geometry = feature.GetGeometryRef()
        geometry.Transform(transform_erosion_img)
        if ext.Intersect(geometry):
            # gera lista de pontos do poligono que pertencem a imagem
            # AddPolygon(geometry, geot, points_img, num, xAxis, yAxis)
            intersection = ext.Intersection(geometry)
            # p_x, p_y = get_vertices(intersection, geot_img)
            fill_ground_truth(ground_truth, intersection, geot_img)


def check_intersection(x, y, p_x_final, p_y_final, min_x, min_y, max_x, max_y):
    if x < max_x and p_x_final > min_x and p_y_final > min_y and y < max_y:
        return True
    return False


def print_from_ring(ring):
    for j in range(ring.GetPointCount()):
        lon, lat, z = ring.GetPoint(j)
        print((lon, lat))


def get_intersection_points(ring, geot_img):
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = -float('inf'), -float('inf')

    for j in range(ring.GetPointCount()):
        lon, lat, z = ring.GetPoint(j)
        min_x = min(min_x, lon)
        min_y = min(min_y, lat)
        max_x = max(max_x, lon)
        max_y = max(max_y, lat)


    if not(math.isinf(min_x) or math.isinf(min_y) or math.isinf(max_x) or math.isinf(max_y)):
        p_min_x, p_min_y = coordinate_to_pixel(geot_img, (min_x, min_y))
        p_max_x, p_max_y = coordinate_to_pixel(geot_img, (max_x, max_y))

    p_min_x = p_min_y = float('inf')
    p_max_x = p_max_y = float('inf')
    # print 'inter points', p_min_x, p_max_y, p_max_x, p_min_y
    # imagem e latitude tem ponto de referencia diferentes
    # enquanto a latitude eh inferior esquerdo, a imagem eh superior esquerdo
    # logo, min lat eh na verdade o maior em termos de pixel
    # tem que inverter o y pq o maximo y na latitude eh o minimo na imagem
    return (p_min_x, p_max_y, p_max_x, p_min_y)

def getVertices(geometry, geotransform):
    if geometry.GetGeometryName() == "GEOMETRYCOLLECTION":
        methodList = [method for method in dir(geometry) if callable(getattr(geometry, method))]
        print(methodList)
        for k in range(0, geometry.GetGeometryCount()):
            g = geometry.GetGeometryRef(k)
            if g is not None:
                if g.GeometryName() == "POLYGON":
                    return getVertices(g, geoTransform)
    
    ring = geometry.GetGeometryRef(0)
    if ring is None:
        return 0,0
    pX = []
    pY = []
    for i in range(ring.GetPointCount()):
        lon, lat, z = ring.GetPoint(i)
        p = utils.CoordinateToPixel(geotransform, (lon,lat))
        pX.append(p[0])
        pY.append(p[1])

    return pX, pY

def FillImage(img, geometry, geoTransform, img_name, dataset_path, fill=1, crop_count=0):
    c, r = getVertices(geometry,geoTransform)
    # print "c: "+str(c)
    # print "r: "+str(r)
    print(("Image Fill: "+str(img_name) + " " + str(crop_count)))
    if c and r:
        rows,cols = polygon(r,c)

        img[rows, cols] = fill

        # name = random.randint(0,50000)
        # imges
        # resto,images_path  = images_path.rsplit("/",1)
        #img_name, resto = img_name.rsplit(".",1)
        #scipy.misc.imsave(os.path.join(dataset_path, 'Masks', os.path.basename(img_name+"_crop"+str(crop_count)+"_mask.jpg")),img*255 )
        # scipy.misc.imsave(dataset_path+img+"_crop"+str(crop_count)+".png", img*255)
        

def create_ground_truth(shapefile, img, shapelabels, imgref, geot, img_name, dataset_path, crop_count=0):
    ds = ogr.Open(shapefile)
    layer = ds.GetLayer(0)
    shpref = layer.GetSpatialRef()
    # imgref = osr.SpatialReference(wkt = img.GetProjectionRef())
    transform = osr.CoordinateTransformation(shpref, imgref)
    # geot  = img.GetGeoTransform()
   

    xAxis = img.shape[1] # Max columns
    yAxis = img.shape[0] # Max rows
    groundTruth = np.zeros((yAxis, xAxis), dtype='bool8')

    ext = utils.GetExtentGeometry(geot, xAxis, yAxis)
    ext.FlattenTo2D()


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

    references = []
    # Adiciona os centroids dos poligonos de erosao a lista de referencias
    layer.ResetReading()
    for fid, feature in enumerate(layer):
        #if fid in erosionFeatures:
        geometry = feature.GetGeometryRef()
        geometry.Transform(transform)
        if ext.Intersect(geometry):
            intersection = geometry.Intersection(ext)
            cntr = intersection.Centroid()
            cntrpxl = utils.CoordinateToPixel(geot,(cntr.GetX(), cntr.GetY()))
            r = cntrpxl[1]
            c = cntrpxl[0]
            references.append(((r,c), intersection))
            FillImage(groundTruth, intersection, geot, img_name, dataset_path, crop_count=crop_count)

    img_name, resto = img_name.rsplit(".",1)
    scipy.misc.imsave(os.path.join(dataset_path, 'Masks', os.path.basename(img_name+"_crop"+str(crop_count)+"_mask.jpg")),groundTruth*255 )
        

def create_patches(img, geot_img, x, y, output_window_size, shapefile, shapelabels, img_name, dataset_path, crop_count=0):
    # get info from the geoTransform of the image
    x_origin = geot_img[0]
    y_origin = geot_img[3]
    pixel_width = geot_img[1]
    print(("Image Patch: "+str(img_name)))
    offset_x, offset_y = 0, 0
    # convert central coordinate to pixel
    p_x_c, p_y_c = coordinate_to_pixel(geot_img, (x, y))

    # get initial pixel - left upper
    p_x, p_y = p_x_c - int(output_window_size / 2), p_y_c - int(output_window_size / 2)
    if p_x < 0:
        offset_x = 0 - p_x
        p_x = 0
    if p_y < 0:
        offset_y = 0 - p_y
        p_y = 0

    # get final pixel - right down
    p_x_size, p_y_size = p_x_c + int(output_window_size / 2) + offset_x, p_y_c + int(output_window_size / 2) + offset_y

    # transform pixels back to coordinates
    x_begin, y_begin = pixel_to_coordinate(geot_img, (p_x, p_y))
    x_final, y_final = pixel_to_coordinate(geot_img, (p_x_size, p_y_size))

    # create polygon (or patch) based on the coordinates
    poly = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(x_begin, y_begin)
    ring.AddPoint(x_begin, y_final)
    ring.AddPoint(x_final, y_final)
    ring.AddPoint(x_final, y_begin)
    ring.AddPoint(x_begin, y_begin)
    ring.CloseRings()
    poly.AddGeometry(ring)

    # create patch array
    xoff = int((x_begin - x_origin) / pixel_width)
    yoff = int((y_origin - y_begin) / pixel_width)
    xcount = int(abs(x_final - x_begin) / pixel_width)
    ycount = int(abs(y_final - y_begin) / pixel_width)
    np_image_array = np.moveaxis(img.ReadAsArray(xoff, yoff, xcount, ycount), 0, -1)[:, :, 0:3]


    imgref = osr.SpatialReference(wkt = img.GetProjectionRef())
    geot  = img.GetGeoTransform()
    xoff_coord, yoff_coord = pixel_to_coordinate(geot, (xoff, yoff))
    # geot[0]= xoff_coord
    # geot[3]= yoff_coord
    geot_coord = (xoff_coord, geot[1], geot[2], yoff_coord, geot[4], geot[5])

    create_ground_truth(shapefile,np_image_array, shapelabels, imgref, geot_coord, img_name, dataset_path, crop_count=crop_count)

    return poly, np_image_array



def from_pixel_to_coordinate(img_path, p):
    img = gdal.Open(img_path)
    geot_img = img.GetGeoTransform()
    return pixel_to_coordinate(geot_img, p)


def create_sets(layer_erosion, percentage=0.2):
    total = len(layer_erosion)
    fids_validation = np.asarray(random.sample(list(range(total)), int(total*percentage)))
    fids_train = np.array(set(np.arange(total)).difference(set(fids_validation)))
    return list(fids_train.tolist()), list(fids_validation)


def create_segmentation_object(path):
    for img in sorted(glob(path+"/*")):
        mask  = io.imread(img)
        mask_bn=np.where(mask>127, 255, 0)
        labeled = label(mask_bn)
        resto,img  = img.rsplit("/",1)
        pathmin,resto = path.rsplit("SegmentationClass",1)
        pathmin = pathmin + "SegmentationObject/"
        if not os.path.exists(pathmin):
            os.mkdir(pathmin)
        scipy.misc.imsave(pathmin+img, labeled)

def comments_cod():

    '''
    ring = geometry.GetGeometryRef(0)
    min_x, min_y = float('inf'), float('inf')
    max_x, max_y = -float('inf'), -float('inf')
    for j in range(ring.GetPointCount()):
        lon, lat, z = ring.GetPoint(j)
        min_x = min(min_x, lon)
        min_y = min(min_y, lat)
        max_x = max(max_x, lon)
        max_y = max(max_y, lat)
    print 'end of ring', min_x, min_y, max_x, max_y, check_intersection(x, y, p_x_final, p_y_final, min_x, min_y, max_x, max_y)
    if check_intersection(x, y, p_x_final, p_y_final, min_x, min_y, max_x, max_y):
        raw_input('press')'''


    '''
    python create_voc_style_from_tcu.py IMAGES_PATH_BELOW /datasets/shps/LinhaFerrovia_Playades.shp /datasets/shps/pontes.shp /datasets/bridge/ 448
    '''

    '''
    '/media/tcu/ImagensFerrovia/Parte1/IMGS/1/DIM_PHR1A_PMS_201608151315351_ORT_1938228101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/10/dim_phr1b_pms_201608141324266_ort_1938237101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/3/DIM_PHR1A_PMS_201608151315544_ORT_1938230101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/4/DIM_PHR1A_PMS_201608151316024_ORT_1938231101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/5/DIM_PHR1A_PMS_201608151316100_ORT_1938232101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/6/DIM_PHR1A_PMS_201608151316199_ORT_1938233101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/7/DIM_PHR1A_PMS_201608151316323_ORT_1938234101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/8/DIM_PHR1A_PMS_201608151316396_ORT_1938235101-001.img', '/media/tcu/ImagensFerrovia/Parte1/IMGS/9/DIM_PHR1B_PMS_201608141324168_ORT_1938236101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/1/DIM_PHR1A_PMS_201608151316508_ORT_1938238101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/10/DIM_PHR1A_PMS_201608221313003_ORT_1938247101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/11/DIM_PHR1B_PMS_201608211320281_ORT_1938248101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/2/DIM_PHR1A_PMS_201608171301204_ORT_1938239101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/3/DIM_PHR1A_PMS_201608171301303_ORT_1938240101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/4/DIM_PHR1A_PMS_201608221312129_ORT_1938241101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/5/DIM_PHR1A_PMS_201608221312219_ORT_1938242101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/6/DIM_PHR1A_PMS_201608221312291_ORT_1938243101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/7/DIM_PHR1A_PMS_201608221312368_ORT_1938244101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/8/DIM_PHR1A_PMS_201608221312448_ORT_1938245101-001.img', '/media/tcu/ImagensFerrovia/Parte2/IMGS/9/DIM_PHR1A_PMS_201608221312528_ORT_1938246101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/1/DIM_PHR1A_PMS_201608171301035_ORT_2000410101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/2/DIM_PHR1A_PMS_201608171301111_ORT_2000411101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/3/DIM_PHR1A_PMS_201608171301204_ORT_2000412101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/4/DIM_PHR1A_PMS_201608171301399_ORT_2000413101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/5/DIM_PHR1A_PMS_201608221312009_ORT_2000414101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/6/DIM_PHR1B_PMS_201608161308354_ORT_2000415101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/7/DIM_PHR1B_PMS_201609251301193_ORT_2000416101-001.img', '/media/tcu/ImagensFerrovia/Parte3/IMGS/8/DIM_PHR1B_PMS_201609251301435_ORT_2000417101-001.img', '/media/tcu/ImagensFerrovia/Parte4/IMGS/1/DIM_PHR1B_PMS_201608161307419_ORT_1939468101-001.img', '/media/tcu/ImagensFerrovia/Parte4/IMGS/2/DIM_PHR1B_PMS_201608161307514_ORT_1939469101-001.img', '/media/tcu/ImagensFerrovia/Parte4/IMGS/3/DIM_PHR1B_PMS_201608161308001_ORT_1939470101-001.img', '/media/tcu/ImagensFerrovia/Parte4/IMGS/4/DIM_PHR1B_PMS_201608161308149_ORT_1939471101-001.img', '/media/tcu/ImagensFerrovia/Parte4/IMGS/5/DIM_PHR1B_PMS_201608301300400_ORT_1939472101-001.img', '/media/tcu/ImagensFerrovia/Parte4/IMGS/6/DIM_PHR1B_PMS_201608301301511_ORT_1939473101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/1/DIM_PHR1A_PMS_201608121250273_ORT_1958903101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/10/DIM_PHR1B_PMS_201608181254140_ORT_1958906101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/11/DIM_PHR1B_PMS_201608181254214_ORT_1958907101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/12/DIM_PHR1B_PMS_201608301301066_ORT_1958908101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/13/DIM_PHR1B_PMS_201608301301140_ORT_1958909101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/14/DIM_PHR1B_PMS_201608301301225_ORT_1958910101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/15/DIM_PHR1B_PMS_201608301301435_ORT_1958911101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/16/DIM_PHR1B_PMS_201608301301511_ORT_1958912101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/17/DIM_PHR1B_PMS_201608301301591_ORT_1958913101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/18/DIM_PHR1B_PMS_201609061258119_ORT_1958914101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/2/DIM_PHR1A_PMS_201608171300374_ORT_1958904101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/3/DIM_PHR1A_PMS_201608171300450_ORT_1958889101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/4/DIM_PHR1B_PMS_201608181253113_ORT_1958890101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/5/DIM_PHR1B_PMS_201608181253189_ORT_1958891101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/6/DIM_PHR1B_PMS_201608181253279_ORT_1958892101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/7/DIM_PHR1B_PMS_201608181253353_ORT_1958905101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/8/DIM_PHR1B_PMS_201608181253446_ORT_1958893101-001.img', '/media/tcu/ImagensFerrovia/Parte5/IMGS/9/DIM_PHR1B_PMS_201608181254043_ORT_1958894101-001.img'
    '''

    '''
    python create_voc_style_from_tcu_online.py '/media/tcu/ImagensFerrovia/Parte1/IMGS/1/DIM_PHR1A_PMS_201608151315351_ORT_1938228101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/10/dim_phr1b_pms_201608141324266_ort_1938237101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/3/DIM_PHR1A_PMS_201608151315544_ORT_1938230101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/4/DIM_PHR1A_PMS_201608151316024_ORT_1938231101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/5/DIM_PHR1A_PMS_201608151316100_ORT_1938232101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/6/DIM_PHR1A_PMS_201608151316199_ORT_1938233101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/7/DIM_PHR1A_PMS_201608151316323_ORT_1938234101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/8/DIM_PHR1A_PMS_201608151316396_ORT_1938235101-001.img','/media/tcu/ImagensFerrovia/Parte1/IMGS/9/DIM_PHR1B_PMS_201608141324168_ORT_1938236101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/1/DIM_PHR1A_PMS_201608151316508_ORT_1938238101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/10/DIM_PHR1A_PMS_201608221313003_ORT_1938247101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/11/DIM_PHR1B_PMS_201608211320281_ORT_1938248101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/2/DIM_PHR1A_PMS_201608171301204_ORT_1938239101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/3/DIM_PHR1A_PMS_201608171301303_ORT_1938240101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/4/DIM_PHR1A_PMS_201608221312129_ORT_1938241101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/5/DIM_PHR1A_PMS_201608221312219_ORT_1938242101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/6/DIM_PHR1A_PMS_201608221312291_ORT_1938243101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/7/DIM_PHR1A_PMS_201608221312368_ORT_1938244101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/8/DIM_PHR1A_PMS_201608221312448_ORT_1938245101-001.img','/media/tcu/ImagensFerrovia/Parte2/IMGS/9/DIM_PHR1A_PMS_201608221312528_ORT_1938246101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/1/DIM_PHR1A_PMS_201608171301035_ORT_2000410101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/2/DIM_PHR1A_PMS_201608171301111_ORT_2000411101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/3/DIM_PHR1A_PMS_201608171301204_ORT_2000412101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/4/DIM_PHR1A_PMS_201608171301399_ORT_2000413101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/5/DIM_PHR1A_PMS_201608221312009_ORT_2000414101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/6/DIM_PHR1B_PMS_201608161308354_ORT_2000415101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/7/DIM_PHR1B_PMS_201609251301193_ORT_2000416101-001.img','/media/tcu/ImagensFerrovia/Parte3/IMGS/8/DIM_PHR1B_PMS_201609251301435_ORT_2000417101-001.img','/media/tcu/ImagensFerrovia/Parte4/IMGS/1/DIM_PHR1B_PMS_201608161307419_ORT_1939468101-001.img','/media/tcu/ImagensFerrovia/Parte4/IMGS/2/DIM_PHR1B_PMS_201608161307514_ORT_1939469101-001.img','/media/tcu/ImagensFerrovia/Parte4/IMGS/3/DIM_PHR1B_PMS_201608161308001_ORT_1939470101-001.img','/media/tcu/ImagensFerrovia/Parte4/IMGS/4/DIM_PHR1B_PMS_201608161308149_ORT_1939471101-001.img','/media/tcu/ImagensFerrovia/Parte4/IMGS/5/DIM_PHR1B_PMS_201608301300400_ORT_1939472101-001.img','/media/tcu/ImagensFerrovia/Parte4/IMGS/6/DIM_PHR1B_PMS_201608301301511_ORT_1939473101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/1/DIM_PHR1A_PMS_201608121250273_ORT_1958903101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/10/DIM_PHR1B_PMS_201608181254140_ORT_1958906101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/11/DIM_PHR1B_PMS_201608181254214_ORT_1958907101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/12/DIM_PHR1B_PMS_201608301301066_ORT_1958908101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/13/DIM_PHR1B_PMS_201608301301140_ORT_1958909101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/14/DIM_PHR1B_PMS_201608301301225_ORT_1958910101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/15/DIM_PHR1B_PMS_201608301301435_ORT_1958911101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/16/DIM_PHR1B_PMS_201608301301511_ORT_1958912101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/17/DIM_PHR1B_PMS_201608301301591_ORT_1958913101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/18/DIM_PHR1B_PMS_201609061258119_ORT_1958914101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/2/DIM_PHR1A_PMS_201608171300374_ORT_1958904101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/3/DIM_PHR1A_PMS_201608171300450_ORT_1958889101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/4/DIM_PHR1B_PMS_201608181253113_ORT_1958890101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/5/DIM_PHR1B_PMS_201608181253189_ORT_1958891101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/6/DIM_PHR1B_PMS_201608181253279_ORT_1958892101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/7/DIM_PHR1B_PMS_201608181253353_ORT_1958905101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/8/DIM_PHR1B_PMS_201608181253446_ORT_1958893101-001.img','/media/tcu/ImagensFerrovia/Parte5/IMGS/9/DIM_PHR1B_PMS_201608181254043_ORT_1958894101-001.img' /datasets/shps/LinhaFerrovia_Playades.shp /datasets/shps/pontes.shp 448 train 1
    '''


def main():
    list_params = ['images_path', 'railway_shapefile', 'erosion_shapefile',
                   'output_window_size', 'process (train|test)', 'outputpath', 'shapelabels']
    if len(sys.argv) < len(list_params) + 1:
        sys.exit('Usage: ' + sys.argv[0] + ' ' + ' '.join(list_params))
    print_params(list_params)


    index = 1
    folderlist = sys.argv[index].split(',')
    #images_path = sys.argv[index].split(',')
    #images_path.sort()

    index = index + 1
    railway_shapefile = sys.argv[index]
    index = index + 1
    erosion_shapefile = sys.argv[index]

    index = index + 1
    output_window_size = int(sys.argv[index])
    half_window_h = int(output_window_size / 2)
    half_window_w = int(output_window_size / 2)

    index = index + 1
    process = sys.argv[index]

    index = index + 1
    outputpath = sys.argv[index]
    
    index = index + 1
    shapelabels = sys.argv[index]

    # index = index + 1
    # cache = bool(sys.argv[index]) -- incluir no final do processo

    # m = hashlib.md5()
    # m.update(process + '_' + '_'.join(base))
    # dataset_path = os.path.join('/datasets', m.hexdigest())  # TODO trocar para: /tmp
    dataset_path = outputpath
    # if the path is a folder
    # if os.path.isdir(dataset_path):
    #     shutil.rmtree(dataset_path)
    #     # print("Cheguei aqui if")
    #     # check if this folder has the other subfolders with expected organization (according Pascal VOC dataset)
    #     if os.path.isdir(os.path.join(dataset_path, 'JPEGImages')) and os.path.isdir(os.path.join(dataset_path, 'JPEGImages')):
    #         # if so, then it is ok to go
    #         # return and start processing the main detection/segmentation algorithm
    #         # print("Cheguei aqui if if")
    #         # return True
    #         pass
    #     else:
    #         shutil.rmtree(dataset_path)
    #         # print("Cheguei aqui if else")
    #         # if not, delete the whole folder and subfolders
    #         # this case can only happen if the user delete the data by himself or if the code crashes or is stopped
            
    # else:
    #     # if the path is not a folder, then create the exptected tree of folder and subfolders (according Pascal VOC)
    #     os.makedirs(dataset_path)
    #     os.makedirs(os.path.join(dataset_path, 'JPEGImages'))
    #     os.makedirs(os.path.join(dataset_path, 'Masks'))


    # read bridge shapefile
    shp_erosion = ogr.Open(erosion_shapefile)
    layer_erosion = shp_erosion.GetLayer(0)
    shp_ref_erosion = layer_erosion.GetSpatialRef()

    # read railway shapefile
    shp_railway = ogr.Open(railway_shapefile)
    layer_railway = shp_railway.GetLayer(0)
    shp_ref_railway = layer_railway.GetSpatialRef()

    # if it is train, then split the dataset (based on the fids of the polygons) for training (80%) and validation (20%)
    # if process == 'train':
    #     fids_train, fids_validation = create_sets(layer_erosion)
    #     np.save(os.path.join(dataset_path, 'ImageSets', 'Segmentation', 'fids_train.npy'), fids_train)
    #     np.save(os.path.join(dataset_path, 'ImageSets', 'Segmentation', 'fids_validation.npy'), fids_validation)

    # f_erosion_test = open(os.path.join(dataset_path, 'ImageSets', 'Segmentation', 'train.txt'), 'w')
    # f_test = open(os.path.join(dataset_path, 'ImageSets', 'Segmentation', 'val.txt'), 'w')
    # f_train = open(os.path.join(dataset_path, 'ImageSets', 'Segmentation', 'trainval.txt'), 'w')
    # if process == 'train' and use_google_dataset:
    #     # TODO fazer essa parte automatica: baixar dataset, extrair, e copiar para as pastas certas
    #     # atualmente esse processo eh feito manualmente
    #     for i in range(501):
    #         f_train.write(str(i) + '\n')

    #step = 0.1
    step = 2000
    unique_names = set()
    # print("path: "+str(images_path))
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
        if not os.path.exists(os.path.join(outputpath,'NewDetectionDataset')):
            os.mkdir(os.path.join(outputpath,'NewDetectionDataset'))

        dataset_path = os.path.join(outputpath,'NewDetectionDataset')
        if not os.path.exists(os.path.join(dataset_path, 'JPEGImages')):
            os.makedirs(os.path.join(dataset_path, 'JPEGImages'))
        if not os.path.exists(os.path.join(dataset_path, 'Masks')):
            os.makedirs(os.path.join(dataset_path, 'Masks'))

        base = [os.path.basename(img).split('.')[0] for img in imglist]
        print(dataset_path)
    
        for img_name in imglist:
            if not os.path.exists(img_name.replace(".img", ".ige")):
                continue
            nomecaminho,resto = img_name.rsplit('\\',1)
            # print("caminho: ",str(nomecaminho))
            # print("path: "+str(images_path))
            print((BatchColors.WARNING + "Running image: " + img_name + BatchColors.ENDC))
            # read image, extract geoTranform and extent
            extents = []
            geotransformations = []
            img = gdal.Open(os.path.join(nomecaminho, img_name))
            # img = gdal.Open(img_name)
            geot_img = img.GetGeoTransform()
            geotransformations.append(geot_img)
            x_axis = img.RasterXSize  # Max columns
            y_axis = img.RasterYSize  # Max rows
            extents.append((x_axis, y_axis))

            # translate coordinates of shapefiles into images'
            img_ref = osr.SpatialReference(wkt=img.GetProjectionRef())
            transform_erosion_img = osr.CoordinateTransformation(shp_ref_erosion, img_ref)
            transform_railway_img = osr.CoordinateTransformation(shp_ref_railway, img_ref)

            # convert gdal image file to array to be used further on -- problem of memory constraints with large image
            # img = np.moveaxis(img.ReadAsArray(), 0, -1)

            # extent
            ext = get_extent_geometry(geot_img, x_axis, y_axis)
            ext.FlattenTo2D()

            # extract line limits of railway that intersect the current image
            bounds = []
            bounds.append(extract_railway_line_limits_from_shapefile(layer_railway, transform_railway_img, ext))
            bounds = sort_bounds(bounds)

            # generate black and white mask pointing out regions of interest
            # generate_mask_for_img(layer_erosion, transform_erosion_img, ext, geot_img, y_axis, x_axis)
            
            cropCount = 0
            for i in range(bounds.shape[0]):
                # if point already exist, then it had been processed before -> continue
                if tuple(bounds[i, 1]) + tuple(bounds[i, 0]) in POINTS_DICT:
                    continue
                else:
                    cur_point = check_boundary(ext, bounds[i])

                t = 0
                num_samples = (1 / float(step)) + 1
                p1 = cur_point[1]
                p2 = cur_point[0]

                for j in range(int(num_samples)):
                    x = (p2[0] - p1[0]) * t + p1[0]
                    y = (p2[1] - p1[1]) * t + p1[1]
                    t = t + step

                    x_p, y_p = coordinate_to_pixel(geot_img, (x, y))
                    if (y_p - half_window_h >= 0 and y_p + half_window_h <= y_axis and
                            x_p - half_window_w >= 0 and x_p + half_window_w <= x_axis):
                        # create georeferenced and array patches
                        geo_patch, np_image = create_patches(img, geot_img, x, y, output_window_size, erosion_shapefile, shapelabels, img_name, dataset_path, crop_count=cropCount)
                        # check if the patch has more than 20% of pixels with black color
                        # this was created for a specific case when the raster is larger but most of it is composed of black
                       
                        cropName = os.path.basename(img_name).replace(".img", "") + "_crop" + str(cropCount) + ".jpg"
                        cropName = os.path.join(dataset_path, 'JPEGImages', cropName)
                        cropCount += 1
                        #print(cropName)
                        scipy.misc.imsave(cropName, np_image)

                        if np.bincount(np_image.astype(int).flatten())[0] < output_window_size * output_window_size * 0.2:
                            fids = []

                            # if it is train, then check for bridge features that intersect with current patch
                            # if it is test, nothing to do here!
                            # if process == 'train':
                            #     layer_erosion.ResetReading()
                            #     inter_ps = []
                            #     for fid, feature in enumerate(layer_erosion):
                            #         geometry = feature.GetGeometryRef()
                            #         geometry.Transform(transform_erosion_img)
                            #         if geo_patch.Intersect(geometry):
                            #             intersection = geo_patch.Intersection(geometry)
                            #             inter_p = get_intersection_points(intersection.GetGeometryRef(0), geot_img)
                            #             if abs(inter_p[2] - inter_p[0]) * abs(inter_p[1] - inter_p[3]) > 100:
                            #                 inter_ps.append(inter_p)
                            #                 fids.append(fid)

                            # validation because of intersection of railway lines in the shapefile
                            # this intersection makes a same point to be processed multiple times
                            # this validation implements a naive way to avoid this
                            # if os.path.basename(img_name).replace(".img", "") + "_" + str(x_p) + "_" + str(y_p) not in unique_names:
                            #     scipy.misc.imsave(os.path.join(dataset_path, 'JPEGImages', os.path.basename(img_name).replace(".img", "") +
                            #                                    "_" + str(x_p) + "_" + str(y_p) + ".jpg"), np_image)
                            #     unique_names.add(os.path.basename(img_name).replace(".img", "") + "_" +
                            #                      str(x_p) + "_" + str(y_p))
                            #     # if fids:
                                #     save_xml(os.path.join(dataset_path, 'Annotations'),
                                #              os.path.basename(img_name).replace(".img", "") + "_" + str(x_p) + "_" + str(y_p),
                                #              y_axis, x_axis, inter_ps,
                                #              (x_p - half_window_w), (y_p - half_window_h))
                                # # print "pontons: x:"+ str(x_p) + "y:" + str(y_p)
                                # if process == 'test' or not fids or (process == 'train' and
                                #                                          any(np.isin(fids, fids_validation)) and
                                #                                          not any(np.isin(fids, fids_train))):
                                #     # print 'test', fids, fids_validation
                                #     f_erosion_test.write(os.path.basename(img_name).replace(".img", "") + "_" +
                                #                         str(x_p) + "_" + str(y_p) + ' ' + ('1' if fids else '-1') + '\n')
                                #     f_test.write(os.path.basename(img_name).replace(".img", "") + "_" +
                                #                  str(x_p) + "_" + str(y_p) + '\n')
                                # if process == 'train' and any(np.isin(fids, fids_train)):
                                #     # print 'train', fids, fids_train
                                #     f_train.write(os.path.basename(img_name).replace(".img", "") + "_" +
                                #                   str(x_p) + "_" + str(y_p) + '\n')

                POINTS_DICT[tuple(bounds[i, 1]) + tuple(bounds[i, 0])] = (x, y)


        # create_segmentation_object(os.path.join(dataset_path, 'SegmentationClass'))

        # f_erosion_test.close()
        # f_test.close()
        # f_train.close()


if __name__ == "__main__":
    main()
