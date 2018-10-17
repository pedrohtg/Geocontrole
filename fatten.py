import os
import sys
import gdal
import numpy as np
import ogr
import osr

def extract_railway_line_limits_from_shapefile(layer_railway):
    # reset layer reading for new reading
    layer_railway.ResetReading()

    lims = []
    for feature in layer_railway:
        railway = feature.GetGeometryRef()    
        if railway.GetGeometryName() == 'MULTILINESTRING':
            for l in range(railway.GetGeometryCount()):
                line = railway.GetGeometryRef(l)
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
            for i in range(railway.GetPointCount() - 1):
                pnts = []
                p1p = railway.GetPoint(i)
                p2p = railway.GetPoint(i + 1)

                p1 = (p1p[0], p1p[1])
                p2 = (p2p[0], p2p[1])

                pnts.append(p1)
                pnts.append(p2)
                lims.append(pnts)
    return lims

# assume points contains at least 3 points
def create_polygon(points):
    poly = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    ring.AddPoint(points[0][0], points[0][1])
    for i in range(1, len(points)):
        ring.AddPoint(points[i][0], points[i][1])
    ring.AddPoint(points[0][0], points[0][1])
    ring.CloseRings()
    poly.AddGeometry(ring)

    return poly

def perpendicular_normalized_vector(vector):
    x = vector[1]
    y = -1*vector[0]

    norm = np.sqrt(x**2 + y**2)
    if norm == 0.000:
        norm = 1.0
    return (x/norm, y/norm)

def translate_point(p, direction, alpha=5.1):
    xnew = p[0] + alpha*direction[0]
    ynew = p[1] + alpha*direction[1]

    return (xnew, ynew)

def fatten(rlwshp, outshapename):
	ds = gdal.Open(rlwshp)
	layer = ds.GetLayer(0)
    lims = extract_railway_line_limits_from_shapefile(layer)

    # Save extent to a new Shapefile
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    # Remove output shapefile if it already exists
    if os.path.exists(outshapename):
        outDriver.DeleteDataSource(outshapename)
    
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outshapename)
    outLayer = outDataSource.CreateLayer("railway track polygons", geom_type=ogr.wkbPolygon)
    
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    outLayer.CreateField(idField)

    for id, segment in enumerate(lims):
        vx = segment[1][0] - segment[0][0]
        vy = segment[1][1] - segment[0][1]
        v = (vx, vy)
        perp = perpendicular_normalized_vector(v)
        _perp = (-1*perp[0], -1*perp[1])

        vertices = []
        vertices.append(translate_point(segment[0], perp))
        vertices.append(translate_point(segment[0], _perp))
        vertices.append(translate_point(segment[1], _perp))
        vertices.append(translate_point(segment[1], perp))

        poly = create_polygon(vertices)

        # Create the feature and set values
        featureDefn = outLayer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(poly)
        feature.SetField("id", id)
        outLayer.CreateFeature(feature)
        feature = None

    ds = None
    outDataSource = None

if __name__ == "__main__":
    rlwshp = sys.argv[1]
    outshapename = sys.argv[2]
	fatten(rlwshp, outshapename)