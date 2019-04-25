import gc
import os
import sys

import dask as ds
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr


class PreProcess(object):
    """
    Description: Preprocessing the hydrography data into the LISFLOOD-FP friendly (most of them are ARC ascii format) format.
    Dependency: rasterio > 1.0.0
                xarray > 0.11.3
    Methods: preprocess(*args): serial computation with single input data. This assumes that your model domain is included in a single file.
                                No parallelization is implemented.
             mfpreprocess(*arg): parallel computation with multiple input data. This assumes that your model domain is located over multiple files.
                                 Note that you can use this with single input data but with parallel computation.
                                 If your input file is very large enough to contain your whole domain, this would be a option you may think.
    Required data: Hydrography (upstream area, surface elevation (hydrologically corrected), and river width) information.
                   Now only a geotiff format is supported, but it should be very easy to implement netCDF format.
    User caution: Please make sure that you assure the RAM storage enough to contain 1.5 or 2 x (your data size for one variable (e.g., elevation) for your WHOLE DOMAIN)
                      even if you are using mfpreprocess.
                  The mfpreprocess stores mapped object (lazy loaded) in a format of dask array and xarray as possible as we can,
                      but it will compile your data when it dumps to the txt format.
                  This means that you will need a RAM at least the size of your single variable data.
                  Note that you can also reduce the total size of the output domain by specifying domain
                      ([llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon]) optional argument,
                      by cutting off the unneeded edge of your input tiles.
    """

    def __init__(self):
        self.undef = -9999.

    # Basic IO components
    def readGeoTiff(self,
                    fileName,
                    latName="y",
                    lonName="x",
                    bandNum=0,
                    domain=None):
        """Read geotiff formatted file with rasteio API in xarray"""
        data = xr.open_rasterio(fileName).rename({
            latName: "lat",
            lonName: "lon"
        })
        if isinstance(domain, list):
            data = self.__domainSlice(data, domain)
        lats = data["lat"].values
        lons = data["lon"].values
        cellsize = data.res[0]
        values = data.values[bandNum]
        return values, lats, lons, cellsize

    def mfreadGeoTiff(self,
                      files,
                      nCols,
                      nRows,
                      latName="y",
                      lonName="x",
                      bandNum=0,
                      domain=None,
                      memLat=1,
                      memLon=1):
        """
        Read multiple geotiff files (tiles) and organize them as a dask array.
        Note: the input file list should be aligned in the C-order (cols>rows).
        """
        # lazy load all data as tiles
        idx = 0
        mapTiles = []
        for row in range(0, nRows):
            rowTiles = []
            for col in range(0, nCols):
                print(files[idx])
                data = xr.open_rasterio(files[idx]).rename({
                    latName: "lat",
                    lonName: "lon"
                })
                idx = idx + 1
                rowTiles.append(data)
            rowTiles = xr.concat(rowTiles, dim="lon")
            mapTiles.append(rowTiles)
        nLat_single = len(data["lat"])
        nLon_single = len(data["lon"])
        # get whole domain
        mapTiles = xr.concat(mapTiles, dim="lat")

        # slicing it is specified
        if isinstance(domain, list):
            mapTiles = self.__domainSlice(mapTiles, domain)
        lats = mapTiles["lat"].values
        lons = mapTiles["lon"].values
        cellsize = mapTiles.res[0]
        values = ds.array.from_array(
            mapTiles.values[bandNum],
            chunks=[nLat_single * memLat, nLon_single * memLon])
        return values, lats, lons, cellsize

    def __domainSlice(self, dataset, domain):

        llcrnrlat = domain[0]
        llcrnrlon = domain[1]
        urcrnrlat = domain[2]
        urcrnrlon = domain[3]
        data = dataset.sel(lat=slice(urcrnrlat, llcrnrlat))
        if len(dataset["lat"]) == 0:
            sys.stderr.write("Runtime Warning: Got a latitudinal shape of 0. \
                 Check if the order of the latitudinal axis is NS. \
                 If not, the xarray will return a empty array.")
        data = dataset.sel(lon=slince(llcrnrlon, urlrnrlon))
        return dataset

    def makeHeader(self, lats, lons, cellsize, order="NS"):
        """Make a header for the lisflood-fp model"""
        nrows = len(lats)
        ncols = len(lons)
        if order == "NS":
            yllcorner = lats[-1]
        elif order == "SN":
            yllcorner = lats[0]
        else:
            raise IOError("order %s is not defined." % order)
        xllcorner = lons[0]
        header = "ncols %d\nnrows %d\nxllconer %f\nyllcorner %f\ncellsize %f\nNODATA_value %d\n" % (
            ncols, nrows, xllcorner, yllcorner, cellsize, self.undef)
        return header

    def dump(self, array2D, header, fileName):
        """dump data and header into the file (filename)"""
        df = pd.DataFrame(array2D)
        with open(fileName, "w") as f:
            f.write(header)
        with open(fileName, "a") as f:
            df.to_csv(f, header=False, index=False)

    def daskDump(self, darray2D, header, fileName):
        """
        dump dask array and header into the file (filename).
        Note that you need at least more storage in RAM than the actual your 2D domain size,
        as here dask finally compile your data.
        """
        df = (dd.from_dask_array(darray2D)).compute()
        with open(fileName, "w") as f:
            f.write(header)
        with open(fileName, "a") as f:
            df.to_csv(f, header=False, index=False)
        del df
        gc.colloct()

    # Specific operation toward the hydrography
    def defineRivers(self, upareaMap, thsld):
        """Define boolean river network map"""
        rivMap = upareaMap > thsld  # boolean Array
        return rivMap

    def maskNoRivers(self, array2d, rivMap):
        """Mask out the raster pixel where no river is defined."""
        array2d[rivMap == False] = self.undef
        return array2d

    def lazyMaskNoRivers(self, darray2d, rivMap):
        """
        Mask out the raster pixel where no river is defined.
        Dask parallel computing version.
        """
        array2d = ds.array.where(rivMap, x=darray2d, y=self.undef)
        return array2d

    def cacheAsNc(self, array2d, lats, lons, name):
        data = xr.DataArray(
            array2d, coords={
                "lat": lats,
                "lon": lons
            }, dims=["lat","lon"], name=name)
        if not os.path.exists("./cache/"):
            os.makedirs("./cache/")
        data = data.to_netcdf("./cache/%s.nc" % name)

    # Main modules to make input topography data
    def preprocess(self,
                   upaPath,
                   elvPath,
                   wthPath,
                   thsld,
                   domain=None,
                   prefix="lisfld",
                   latName="y",
                   lonName="x",
                   bandNum=0):
        """
        Main module to process hydrography data into the lisflood format:
        Assuming that each file will contain your whole modeled domain.
        If it is not a case (your domain is within multiple files), use mfpreprocess().
        """
        # make .dem.ascii
        elv, lats, lons, cellsize = self.readGeoTiff(
            elvPath,
            latName=latName,
            lonName=lonName,
            bandNum=bandNum,
            doamin=domain)
        header = self.makeHeader(lats, lons, cellsize)
        oName = "%s.dem.ascii" % prefix
        self.dump(elv, header, oName)
        # make river network
        uparea, *_ = self.readGeoTiff(
            upaPath,
            latName=latName,
            lonName=lonName,
            bandNum=bandNum,
            domain=domain)
        self.cacheAsNc(uparea, lats, lons, "uparea")
        rivMap = self.defineRivers(uparea, thsld)
        # make .width.asc
        width, *_ = self.readGeoTiff(
            wthPath,
            latName=latName,
            lonName=lonName,
            bandNum=bandNum,
            domain=domain)
        width = self.maskNoRivers(width, rivMap)
        oName = "%s.width.asc" % prefix
        self.dump(width, header, oName)
        # make .bank.asc
        banke = self.maskNoRivers(elv, header, thsld)
        oNmae = "%s.bank.asc" % prefix
        self.dump(banke, header, oName)

    def mfpreprocess(self,
                     upaPath,
                     elvPath,
                     wthPath,
                     thsld,
                     nCols,
                     nRows,
                     domain=None,
                     prefix="lisfld",
                     latName="y",
                     lonName="x",
                     bandNum=0,
                     memLat=1,
                     memLon=1):
        """
        Main module to process hydrography data into the lisflood format:
        Same as preprocess but from multiple files. Designed to reduce required memory as small as possible.
        """
        # make .dem.ascii
        elv, lats, lons, cellsize = self.mfreadGeoTiff(
            elvPath,
            nCols,
            nRows,
            latName=latName,
            lonName=lonName,
            bandNum=bandNum,
            domain=domain,
            memLat=memLat,
            memLon=memLon)
        header = self.makeHeader(lats, lons, cellsize)
        print("output file informaton:\n%s"%header)
        oName = "%s.dem.ascii" % prefix
        print("%s" % oName)
        self.daskDump(elv, header, oName)

        #make river network
        uparea, *_ = self.mfreadGeoTiff(
            upaPath,
            nCols,
            nRows,
            latName=latName,
            lonName=lonName,
            bandNum=bandNum,
            domain=domain,
            memLat=memLat,
            memLon=memLon)
        self.cacheAsNc(uparea, lats, lons, "uparea")
        rivMap = self.defineRivers(uparea, thsld)
        del uparea
        gc.collect()

        #make .width.asc
        width, *_ = self.mfreadGeoTiff(
            wthPath,
            nCols,
            nRows,
            latName=latName,
            lonName=lonName,
            bandNum=bandNum,
            domain=domain,
            memLat=memLat,
            memLon=memLon)
        width = self.lazyMaskNoRivers(width, rivMap)
        oName = "%s.width.asc" % prefix
        self.daskDump(width, header, oName)
        print("%s" % oName)

        #make .bank.asc
        banke = self.lazyMaskNoRivers(elv, rivMap)
        oName = "%s.bank.asc" % prefix
        self.daskDump(banke, header, oName)
        print("%s" % oName)


if __name__ == "__main__":
    upaPaths = [
        "./deploy/MERIT_HYDRO/upa/upa_n30w120/n35w100_upa.tif",
        "./deploy/MERIT_HYDRO/upa/upa_n30w120/n35w095_upa.tif"
    ]
    elvPaths = [
        "./deploy/MERIT_HYDRO/elv/elv_n30w120/n35w100_elv.tif",
        "./deploy/MERIT_HYDRO/elv/elv_n30w120/n35w095_elv.tif"
    ]
    wthPaths = [
        "./deploy/MERIT_HYDRO/wth/wth_n30w120/n35w100_wth.tif",
        "./deploy/MERIT_HYDRO/wth/wth_n30w120/n35w095_wth.tif"
    ]
    nCols = 2
    nRows = 1
    thsld = 24.04
    test = PreProcess()
    test.mfpreprocess(upaPaths, elvPaths, wthPaths, thsld, nCols, nRows)
