import gc
import os

import dask.dataframe as dd
import datetime
import numpy as np
import pandas as pd
import xarray as xr
from numba import jit


class Forcing(object):
    def __init__(self):
        # general settings
        self.undef = -9999.
        self.buf = 5
        self.thsld = 0.05
        self.distance = 0.10
        self.outDir = "./out/"
        if not os.path.exists(self.outDir):
            os.makedirs(self.outDir)

        # advanced settings. development purpose.
        self.upaCache = "./cache/uparea.nc"
        self.cacheType = "nc"

    # Basic IO component.
    def readPoints(self, filename, domain=None, dask=False):
        """
        Read csv file containing points data to map on the model grid cells.
        filename (str): path to the csv file (see below for the format.)
        domain (list,optional): Your domain info [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
        dask (bool,optional): Using dask to perform lazy loading or not.
            Effective for a large csv file (and your domain is only a part of that file).
        Input csv format:
        ex.)
            id,lat,lon,uparea
            0,35.11,-120.24,300.25
        """
        if dask == True:
            df = dd.read_csv(filename).set_index("id")
        else:
            df = pd.read_csv(filename).set_index("id")
        if isinstance(domain, list):
            df = df[df.lat > domain[0]]  # llcrnrlat
            df = df[df.lat < domain[2]]  # urcrnrlat
            df = df[df.lon > domain[1]]  # llcrnrlon
            df = df[df.lon < domain[3]]  # urcrnrlon
        return df.compute() if dask else df

    def readCachedMap(self, filename, kind="nc"):
        """
        Read cached model domain.
        This cache should be created from the lispy.PreProcess.(mf)preprocess()
        """
        if kind == "nc":
            data = (xr.open_dataset(filename)).to_array()
            lats = data.lat.values
            lons = data.lon.values
            vMap = data.values[0,:,:]
        else:
            raise IOError("data type %s is not supported." % kind)
        return vMap, lats, lons

    def readDischarge(self, dschgFile, kind="nc", timeName="time",
                      idName="id", sDate=None, eDate=None):
        """
        Read discharge time series for points. The id of this file should be correspond with pointInfoFile.
        Data Type:
            netCDF:
                values: (time,points)
                dim_1: time
                dim_2: id
        Args: dschgFile(str): a path to the discharge file which contains discharge timesiries over multiple points.
              kind(str): default "nc"; an input data type
              timeName(str): default "time"; a dimension name for a time axis in your data
              idName(str): default "id"; a dimension name for a id (location-id) axis in your data
              sDate(datetime.datetime): default None; starting date to truncate data. 
                                        Use all data length if not specified.
              eDate(datetime.datetime): default None; ending date to truncate data.
                                        Use all data length if not specified.
        """
        if kind == "nc":
            data = (xr.open_dataset(dschgFile)).to_array()
            if all([isinstance(sDate, datetime.datetime),isinstance(eDate, datetime.datetime)]):
                data = data.sel(time=slice(sDate,eDate))
            elif any([isinstance(sDate, datetime.datetime),isinstance(eDate, datetime.datetime)]):
                raise IOError("Both sDate and eDate should be specified.")
            dschg = data.values[0,:,:]
            tmIdx = data[timeName].values
            idIdx = data[idName].values
        else:
            raise IOError("data type %s is not supported." % kind)
        return dschg, tmIdx, idIdx

    def writeBdy(self, ID, discharge, tDiffs, oName):
        """
        Write .bdy file for the model.
        """
        line1 = "inflow_%d\n" % ID
        line2 = "%d seconds\n" % len(tDiffs)
        elm = np.vstack([discharge.reshape(1,-1),tDiffs.reshape(1,-1)])
        df = pd.DataFrame(elm.T)
        df.iloc[:,1] = df.iloc[:,1].values.astype(np.int64)
        with open(oName, "a") as f:
            [f.write(l) for l in [line1, line2]]
            df.to_csv(f, index=False, header=False, sep=" ")

    # main modules
    def locate(self, filename, domain=None, dask=False):
        """
        Locate points on the model grids.
        """
        df = self.readPoints(filename, domain=domain, dask=dask)
        ids = df.index.tolist()
        fLats = df.lat.tolist()
        fLons = df.lon.tolist()
        upareas = df.uparea.tolist()
        del df
        gc.collect()

        upaMap, lats, lons = self.readCachedMap(
            self.upaCache, kind=self.cacheType)

        coords = mapPoints(upareas, fLats, fLons, upaMap, lats, lons,
                           self.thsld, self.buf, self.distance)
        odf = pd.DataFrame(coords)
        odf.index = ids
        odf.columns = ["lat", "lon"]
        return odf.dropna(how="any")

    def makeForcing(self,
                    dschgFile,
                    pointInfoFile,
                    domain=None,
                    sDate=None,
                    eDate=None,
                    dask=False,
                    kind="nc",
                    timeName="time",
                    idName="id",
                    ex="lisfld"):
        """
        Making a input forcing file set (.dci and .bdy) for your domain and forcing data.
        Arguments:
            dschgFile (str): discharge file path. Currently only netCDF format is supported.
                See lispy.Forcing.readDischarge() for the specific format of the netCDF.
            pointInfoFile (str): forcing point information file. Currently only csv format is supported.
                See lispy.Forcing.readPoints() for the specific format of this csv.
            domain (list,optional): Your domain info [llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon]
            dask (bool,optional): Using dask to perform lazy loading or not.
                Effective for a large csv file (and your domain is only a part of that file).
            kind (str, optional): data type of the discharge file. currently only netCDF format is supported.
            timeName (str, optional): the dimension name of time axis in your discharge data.
            idName (str, optional): the dimension name of id axis in your discharge data.
            ex (str, optional): experiment name that will appear on your output files.
        User caution:
            Current version only supports the "P" source identifer and QVAR in a subgrid configuration.
        """
        dschg, tmIdx, idIdx = self.readDischarge(
            dschgFile, kind=kind, timeName=timeName, idName=idName, sDate=sDate, eDate=eDate)

        # convert into the model time step (seconds).
        tmDiffs = [(date - tmIdx[0])/np.timedelta64(1, "s") for date in tmIdx]
        tmDiffs = np.array(tmDiffs).astype(np.int64)
        ### This is a numpy.timedelta object, which is a wrapper for the aactual datetime.timedela object.

        df = self.locate(pointInfoFile, domain=domain, dask=dask)
        # This part can be flexible in a later version.
        bIdentifer = ["P" for i in range(0, len(df))]
        bType = ["QVAR" for i in range(0, len(df))]
        inflwName = ["inflow_%d" % ID for ID in df.index]
        df["identifer"] = bIdentifer
        df["bType"] = bType
        df["inflwName"] = inflwName
        df = df[["identifer", "lon", "lat", "bType", "inflwName"]]
        oName = os.path.join(self.outDir, "setting_%s.bci" % ex)
        df.to_csv(oName, index=False, header=False, sep=" ")
        print(oName)

        # write it into bdy file.
        oName = os.path.join(self.outDir, "inflow_%s.bdy" % ex)
        with open(oName,"w") as f:
            sDate = np.datetime_as_string(tmIdx[0])
            eDate = np.datetime_as_string(tmIdx[-1])
            line1 = "inflow [m3/s] from %s to %s\n" % (sDate, eDate)
            f.write(line1)
        idIdx = idIdx.tolist()
        for ID in df.index:
            i = idIdx.index(ID)
            self.writeBdy(ID, dschg[:,i], tmDiffs, oName)
        print(oName)


# Numba or Cython later.
@jit
def locatePoint(uparea, fLat, fLon, upaMap, lats, lons, thsld, buf, distance):
    """
    Locate a point on a map based on upstream area information.
        uparea: upstream area of a point
        fLat: latitude of a point
        fLon: longitude of a point
        upaMap: upstream area 2d map array of base model grids
        lats: latitude series of base model grids
        lons: longitude series of base model grids
    """
    latIdx = np.argmin(((lats - fLat)**2))
    lonIdx = np.argmin(((lons - fLon)**2))
    distErr = (((lats[latIdx] - fLat))**2 + ((lons[lonIdx] - fLon))**2)**0.5
    if distErr > distance:
        return np.nan, np.nan

    err = 1e+20
    lat = np.nan
    lon = np.nan
    for iLatIdx in range(latIdx - buf, latIdx + buf):
        for iLonIdx in range(lonIdx - buf, lonIdx + buf):
            upa = upaMap[iLatIdx, iLonIdx]
            cErr = ((upa - uparea)**2)**0.5
            if cErr < err:
                err = cErr
                lat = lats[iLatIdx]
                lon = lons[iLonIdx]
    if err > uparea * thsld:
        lat = np.nan
        lon = np.nan

    return lat, lon


# Numba or Cython later.
@jit
def mapPoints(upareas, fLats, fLons, upaMap, lats, lons, thsld, buf, distance):
    """
    Iterate self.locatePoint to create latlon series representing each point on a model grids.
        upareas (list): list of uparea information of points to locate.
        fLats (list): list of latitude information of points to locate.
        fLons (list): list of longitude information of points to locate.
        upaMap: upstream area 2d map array of base model grids
        lats : latitude series of base model grids
        lons: longitude series of base model grids
    """
    mCoords = [
        locatePoint(upareas[i], fLats[i], fLons[i], upaMap, lats, lons, thsld,
                    buf, distance) for i in range(0, len(upareas))
    ]

    return mCoords

if __name__ == "__main__":
    test = Forcing()
    dschgFile = "./deploy/discharge_missouri.nc"
    pointInfoFile = "./deploy/pointInfo_missouri.csv"
    sDate = datetime.datetime(1984,1,1)
    eDate = datetime.datetime(1985,1,1)
    test.makeForcing(dschgFile,pointInfoFile,sDate=sDate,eDate=eDate)