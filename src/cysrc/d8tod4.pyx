from __future__ import division
import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

def d8tod4(np.ndarray[DTYPE_t,ndim=2] elv, np.ndarray[DTYPE_t,ndim=2] flwdir, dict dirDict):
    """
    Modify elevation based on the flowdirection, in order to force water to flow toward diagonal direction.
    """
    cdef int nLat = elv.shape[0]
    cdef int nLon = elv.shape[1]
    cdef int ilat
    cdef int ilon
    cdef int ilat2 = 0
    cdef int ilon2 = 0
    cdef int ilat3 = 0
    cdef int ilon3 = 0
    cdef int ilat4 = 0
    cdef int ilon4 = 0
    cdef int jlat1 = 0
    cdef int jlon1 = 0
    cdef int jlat2 = 0
    cdef int jlon2 = 0
    cdef int modify
    cdef int nextLat
    cdef int nextLon
    cdef int num1 = 0
    cdef int num2 = 0
    cdef int num3 = 0
    cdef list modLoc = []

    for ilat in range(1, nLat-1):
        for ilon in range(1, nLon-1):
            idir = dirDict[flwdir[ilat, ilon]]
            if idir == "SW" or idir == "SE" or idir == "NW" or idir == "NE":
                if idir == "NE":
                    ilon2 = ilon
                    ilat2 = ilat - 1
                    ilon3 = ilon + 1
                    ilat3 = ilat - 1
                    ilon4 = ilon + 1
                    ilat4 = ilat 
                elif idir == "SE":
                    ilon2 = ilon + 1
                    ilat2 = ilat 
                    ilon3 = ilon + 1
                    ilat3 = ilat + 1
                    ilon4 = ilon
                    ilat4 = ilat + 1
                elif idir == "SW":
                    ilon2 = ilon
                    ilat2 = ilat + 1
                    ilon3 = ilon - 1
                    ilat3 = ilat + 1
                    ilon4 = ilon - 1
                    ilat4 = ilat
                elif idir == "NW":
                    ilon2 = ilon - 1
                    ilat2 = ilat
                    ilon3 = ilon - 1
                    ilat3 = ilat - 1
                    ilon4 = ilon
                    ilat4 = ilat - 1
                if elv[ilat2,ilon2] >= elv[ilat,ilon] and elv[ilat4,ilon4] >= elv[ilat,ilon]:
                    if elv[ilat2,ilon2] >= elv[ilat4,ilon4]:
                        jlon1 = ilon2 #higher
                        jlat1 = ilat2
                        jlon2 = ilon4 #lower
                        jlat2 = ilat4
                    else:
                        jlon1 = ilon4 #higher
                        jlat1 = ilat4
                        jlon2 = ilon2 #lower
                        jlat2 = ilat2

                    modify = 0
            
                    if modify == 0:
                        nextLon, nextLat = nextlonlat(jlon1, jlat1, dirDict[flwdir[jlat1, jlon1]])
                        if (nextLon == ilon and nextLat == ilat) or (nextLon == ilon3 and nextLat == ilat3):
                            modify = 1
                    if modify == 0:
                        nextLon, nextLat = nextlonlat(jlon2, jlat2, dirDict[flwdir[jlat2, jlon2]])
                        if (nextLon == ilon and nextLat == ilat) or (nextLon == ilon3 and nextLat == ilat3):
                            modify = 2

                    if modify == 1:
                        elv[jlat1, jlon1] = elv[ilat, ilon]
                        modLoc.append([[jlat1,jlon1],[ilat,ilon]])
                        num1=num1+1
                    elif modify == 2:
                        elv[jlat2, jlon2] = elv[ilat, ilon]
                        modLoc.append([[jlat2,jlon2],[ilat,ilon]])
                        num2=num2+1
                    else:
                        elv[jlat1, jlon1] = elv[ilat, ilon]
                        modLoc.append([[jlat1,jlon1],[ilat,ilon]])
                        num3=num3+1
                #elif elv[ilat2,ilon2] >= elv[ilat,ilon] and elv[ilat4,ilon4] < elv[ilat,ilon]:
                #    modLoc.append([[ilat4,ilon4],[ilat,ilon]])
                #elif elv[ilat2,ilon2] < elv[ilat,ilon] and elv[ilat4,ilon4] >= elv[ilat,ilon]:
                #    modLoc.append([[ilat2,ilon2],[ilat,ilon]])
                else:
                    if elv[ilat2,ilon2] > elv[ilat4,ilon4]:
                        modLoc.append([[ilat4,ilon4],[ilat,ilon]])
                    else:
                        modLoc.append([[ilat2,ilon2],[ilat,ilon]])
                    
    return elv, modLoc

def nextlonlat(int lon, int lat, str fdir):

    cdef int nextLat
    cdef int nextLon

    if fdir == "N":
        nextLat = lat -1 
        nextLon = lon
    elif fdir == "NE":
        nextLat = lat - 1
        nextLon = lon + 1
    elif fdir == "E":
        nextLat = lat
        nextLon = lon + 1
    elif fdir == "SE":
        nextLat = lat + 1
        nextLon = lon + 1
    elif fdir == "S":
        nextLat = lat + 1
        nextLon = lon
    elif fdir == "SW":
        nextLat = lat + 1
        nextLon = lon - 1
    elif fdir == "W":
        nextLat = lat
        nextLon = lon - 1
    elif fdir == "NW":
        nextLat = lat - 1
        nextLon = lon - 1
    else:
        nextLat = lat
        nextLon = lon

    return nextLon, nextLat

def snapD4(np.ndarray[DTYPE_t,ndim=2] mapData, list modLoc, np.ndarray[DTYPE_t,ndim=2] rivMap):
    """
    Snap givin location (modLoc[i][0]) by assiginig the value of the center pixel (modLoc[i][1]).
    modLoc is the output of d8tod4.d8tod4(). 
    """
    
    cdef int i
    for i in range(len(modLoc)):
        snapPixel = modLoc[i][0]
        centPixel = modLoc[i][1]
        if rivMap[centPixel[0],centPixel[1]]:
            mapData[snapPixel[0],snapPixel[1]] = mapData[centPixel[0],centPixel[1]]
        else:
            print(centPixel)
    
    return mapData

def check(np.ndarray[DTYPE_t,ndim=2] elv, np.ndarray[DTYPE_t,ndim=2] flwdir, dict dirDict):
    
    cdef int ilat
    cdef int ilon
    cdef int num
    cdef int tNum = 0

    for ilat in range(2, len(elv)-2):
        for ilon in range(2, len(elv)-2):
            num = 0
            if elv[ilat,ilon] >= elv[ilat+1,ilon]:
                num = num + 1
            if elv[ilat,ilon] >= elv[ilat-1,ilon]:
                num = num + 1
            if elv[ilat,ilon] >= elv[ilat,ilon+1]:
                num = num + 1
            if elv[ilat,ilon] >= elv[ilat,ilon-1]:
                num = num + 1
            if num == 0:
                print("error at:%d,%d,%s"%(ilat,ilon,dirDict[flwdir[ilat,ilon]]))
                tNum = tNum + 1
    if tNum == 0:
        print("No error found on the modified DEM.")
