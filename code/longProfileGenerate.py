'''The process is to assemble a river course from the highest upstream source down to its mouth and generate the long profile. 
we merged polylines into one feature then sampled points along the river course based on user given spacing. 
elevations are interpolated for those points and 3D visualization is performed
Coding for Python 3.7.6
Auhtors: Ozan Tuncbilek, Mengyao Gao , Festina Sadiku, Julius Nyonyo
'''
# -*- coding: utf-8 -*-

# Imports #
# used for directory related operations
import os 
# used for shapefile reading and writing
import shapefile
# used for opening raster data
import gdal
# used for calculating bearings and distances
import math 
# used for visualization result
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Shortcut #
pi = math.pi

# Functions #
#%% a) user input control 
def _createSubdir(workspace, subdirList):
    '''create subdirectories if not already existing'''
    for subdir in subdirList:
        if not os.path.isdir(workspace + "/" + subdir):
            os.mkdir(workspace + "/" + subdir)

def _controlExtension(inName, ext):
    '''enforce a user-defined extension to the input data (e.g. '.shp')'''
    if inName.rfind(".") > 0:
        inName[:inName.rfind(".")] + 'ext'
    else:
        inName + 'ext'
    return inName

def _completePath(workspace, subdir, nameList):
    '''join project dir, subdirectory and file names into complete path'''
    for ix in range(len(nameList)):
        nameList[ix] = workspace + "/" + subdir + "/" + nameList[ix]
    return nameList

def _checkExistence(pathList):
    '''check if the specified datasets exist under the entered paths'''
    check = True
    for data in pathList:
        if not os.path.isfile(data):
            check = False
            break
    return check

#%% b) features merging 

def _mergeFeatures(inFile, outFile):
    '''read the shapefile and combine all polylines to one feature
       _mergeSide function is used within this function '''
    sf = shapefile.Reader(inFile,encoding='unicode_escape')
    geotype = sf.shapeType
    if not (geotype == 3):
        sf.shp.close()
        print('!error: segment centres need type POLYLINE')
        return None
    
    w = shapefile.Writer(outFile, shapeType=3)
    w.fields = list(sf.fields)

    # Get lines from shape objects.
    lines = []
    for it in sf.iterShapeRecords():
        sh = it.shape
        lines.append(sh.points)
    # Sort line so the first line will be the left-most and used as the starting line.
    sortedLines = sorted(lines, key=lambda line: line[0][0])

    cntLines = len(sortedLines)    
    # To track which lines were already merged. First one is merged from the beginning.
    merged = [False] * cntLines
    merged[0] = True
    # Start from the left-most one and add all connected features on one side
    sortedLines,merged=_mergeSide(sortedLines,merged)
    # Reverse the points of the already merged feature
    sortedLines[0].reverse()
    # merging all features on the another side
    sortedLines,merged=_mergeSide(sortedLines,merged)
    
    cnt = sum(merged)
    print("Merged lines : ", cnt)
    shapeEl = sf.shapes()[0]
    shapeEl.points = sortedLines[0]
    w.shape(shapeEl)
    w.record(*sf.records()[0])

def _mergeSide(sortedLines,merged):
    '''merge all polylines on one side
       _calcDistance function is used within this function'''
    while True:
        # define a reasonable smallest distance so we won't connect lines which are too far away.
        smallestDist = 0.1 #km
        # choice is used to mark the linked feature; 
        # flag is used to deal with the direction of feature
        params={'choice':-1,'flag':True}
        # Last point of the polyline which is being built.
        lastPnt = sortedLines[0][-1]
        # go through all non-merged lines and find the linked one
        for i in range(0, len(sortedLines)):
            if merged[i]:
                continue
            firstPnt = sortedLines[i][0]
            endPnt=sortedLines[i][-1]
            # the feature follows the direction 
            if _calcDistance(lastPnt, firstPnt) < smallestDist :
                smallestDist =_calcDistance(lastPnt, firstPnt)
                params['choice'] = i
                params['flag']=True
            # the feature follows the reverse direction
            elif _calcDistance(lastPnt, endPnt) < smallestDist :
                smallestDist =_calcDistance(lastPnt, endPnt)
                params['choice'] = i
                params['flag']=False
            
        if params['choice'] == -1:
            print("no more close lines were found")
            break
        
        # Take line with index |choice| to merge into the building polyline.
        print("choice :", params['choice'])
        merged[params['choice']] = True 
        
        # if the feature in reverse direction, reverse all points and append to the left-most one
        if params['flag']==False:
            sortedLines[params['choice']].reverse()
        for point in sortedLines[params['choice']]:
            sortedLines[0].append(point)
    return (sortedLines,merged)

def _calcDistance(startPoint,endPoint):
    '''calculate the distance between two points'''
    dx = endPoint[0]-startPoint[0]
    dy = endPoint[1]-startPoint[1]
    d= math.sqrt(dx**2 + dy**2)
    return d

#%% c) Sampling points based on the given spacing
def _getCoordinates(inFile):
    '''generate utm coordinates list by reading the generated shapefile'''
    r = shapefile.Reader(inFile)
    coords = r.shape(0).points
    utmCoord = []
    for c in range(len(coords)):
        p = r.shape(0).points[c]
        utmCoord.append(p)
    return utmCoord

def _calcBearing(inFile,utmCoord):
    '''calculate the bearing of each segments'''
    orient = []
    p0 = utmCoord[0]
    
    for p in range(1, len(utmCoord)):
        p1 = utmCoord[p]
        if abs(p0[0]-p1[0]) < 0.00001:
            ang = (pi)/2
        else:
            ang = math.atan( (p0[1]-p1[1]) / (p0[0]-p1[0]) )
        if  (p1[0]-p0[0])<0:
            ang=ang+pi
        orient.append(ang)           
        p0 = p1            
    return orient

def _calcSegmentDist(utmCoord):
    '''calculate the distance of each segments'''
    dist = []
    p0 = utmCoord[0]
    #calculating distances of each segment 
    for u in range(1,len(utmCoord)):
        p1 = utmCoord[u]
        dx = p1[0]-p0[0]
        dy = p1[1]-p0[1]
        d= math.sqrt(dx**2 + dy**2)
        p0 = p1
        dist.append(d)
    return dist
    
def _createVertices (inFile,utmCoord):
    '''calculate the vertices based on the given spacing
       _calcBearing function and _calcSegmentDist function are used within this function'''
    # get the user input spacing
    print("Define the distance for vertices below")
    verDist= float(input())

    bAngs = _calcBearing(inFile,utmCoord)
    dist = _calcSegmentDist(utmCoord)
    
    # we don't need the last point and delete it 
    # so that the lenght of utm_coord would be the same with that of bear list
    utmCoord.pop()
    # calculating vertices and starting from the first point
    vertices=[utmCoord[0]]
    # used to store the remining distance of each segment compared with user given spacing
    distance=0
    for ix,p in enumerate(utmCoord):
        # get the angle, length, starting point, end point of each segment
        angle=bAngs[ix]
        length=dist[ix]
        pointX=p[0]
        pointY=p[1]
        # check if there are remining distance of previous segment
        if(distance!=0):
            verDist2=verDist-distance
            # if the length of the segment is larger than the remaining distance getting from previous segment
            # substract the remaining length of spacing from the lenght of the segment
            if(dist[ix]>=verDist2):
                pointX=pointX+verDist2*math.cos(angle)
                pointY=pointY+verDist2*math.sin(angle)
                vertices.append((pointX,pointY))
                length=length-verDist2
            # if the length of the segment is smaller, go to the next segment
            else:
                distance+=length
                continue
        # get the sampling points of each segment based on user given spacing
        while(length>=verDist):
            pointX=pointX+verDist*math.cos(angle)
            pointY=pointY+verDist*math.sin(angle)
            vertices.append((pointX,pointY))
            length=length-verDist
        distance=length
    return vertices

def _generateVertices(vertices,outFile):
    '''generate the sampling points shapefile to check if the result is correct'''
    file = shapefile.Writer(outFile)
    file.field('FIRST_FLD') 
    for v in vertices:
        file.point(v[0],v[1])
        file.record('First', 'Point')
    file.close()
    
# %% d) Interpolation of Z values 
def _cellPos(inGrid, pnts):
    '''transform from geo-coordinates to grid coordinates'''
    ds = gdal.Open(inGrid)         # open DEM
    cols = ds.RasterXSize          # number of columns
    rows = ds.RasterYSize          # number of rows
    trafo = ds.GetGeoTransform()   # transformation parameters
    gpos = []

    # transform geo-coordinates into grid coordinates
    for p in pnts:
        cell = ( (p[0]-trafo[0])/trafo[1], (trafo[3]-p[1])/(-trafo[5]) )
        if cell[0] > 0 and cell[0] < cols and cell[1] > 0 and cell[1] < rows:
            gpos.append(cell)
        else:
            gpos.append( (-1,-1) )
    return gpos

def _fourCells(ppos):
    '''get the four-cell block relevant for a bilinear interpolation'''
    return [(math.floor(ppos[0]), math.floor(ppos[1])),
            (math.ceil(ppos[0]), math.floor(ppos[1])),
            (math.floor(ppos[0]), math.ceil(ppos[1])),
            (math.ceil(ppos[0]), math.ceil(ppos[1])) ]

def _fourDist(ppos):
    '''the x and y distances to the four cell centres'''
    dx0 = ppos[0] - math.floor(ppos[0])
    dy0 = ppos[1] - math.floor(ppos[1])
    return [dx0, 1.-dx0, dy0, 1.-dy0]

def _openGrid(inGrid):
    '''load DEM as a numpy array using gdal'''
    ds = gdal.Open(inGrid)
    cols = ds.RasterXSize          # number of columns
    rows = ds.RasterYSize          # number of rows
    return ds.GetRasterBand(1).ReadAsArray(0, 0, cols, rows)

def _bilinear(data, cells, dist):   
    '''get elevation from grid by bilinear interpolation '''
    # get the four pixel values addressing indexing row first(!)
    vals = [data[cells[0][1],cells[0][0]],
            data[cells[1][1],cells[1][0]],
            data[cells[2][1],cells[2][0]],
            data[cells[3][1],cells[3][0]]]

    # interpolate
    Zrow0 = dist[1]*vals[0] + dist[0]*vals[1] # row-wise
    Zrow1 = dist[1]*vals[2] + dist[0]*vals[3] # row-wise
    return dist[3]*Zrow0 + dist[2]*Zrow1      # column-wise

def _2DTo3D(gprof, inGrid):
    '''interpolate z-values for sampling points using gridded coordinates and DEM
       _openGrid, _fourCells, _fourDist, _bilinear functions are used within this function'''
    data = _openGrid(inGrid)
    elev = []

    # iterate through alll gridded profile points
    for ix, gp in enumerate(gprof):
        c = _fourCells(gp)
        d = _fourDist(gp)
        elev.append( _bilinear(data, c, d) )
        if ix%500 == 0:
            print('--- interpolate z for '+str(ix)+' points')

    return elev

# %% e) Visualization
def _visualization(coordinates,zs):
    '''perform 3D visualization of the long profile of river'''
    xs=[item[0] for item in coordinates]
    ys=[item[1] for item in coordinates]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i,x in enumerate(xs):
        y=ys[i]
        z=zs[i]
        ax.scatter(x, y, z)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('elevation')
    plt.show()

# %% Extra) Calling functions 
    
# start of main
worksp = '..' # adjust to your needs
inFCName = "wurmUTM.shp"        
inGridName = "clip2.tif"
tmpMergeName= "merge.shp"
tmpVerticesName= "vertices.shp"

# file management
_createSubdir(worksp, ["temp", "output"])
inFCName = _controlExtension(inFCName, ".shp")

inFC   = _completePath(worksp, "shape", [inFCName])[0]
inGrid = _completePath(worksp, "grid", [inGridName])[0]
tmpMerge = _completePath(worksp, "temp", [tmpMergeName])[0]
tmpVertices = _completePath(worksp, "output", [tmpVerticesName])[0]

if not _checkExistence([inFC, inGrid]):
    raise ValueError

# perform tasks
_mergeFeatures(inFC, tmpMerge)
utmCoord= _getCoordinates(tmpMerge)
vertices=_createVertices(tmpMerge,utmCoord)
_generateVertices(vertices,tmpVertices)
gvertices=_cellPos(inGrid,vertices)
elevations=_2DTo3D(gvertices,inGrid)
_visualization(gvertices,elevations)



