# -*- coding: utf-8 -*-a
"""
Dependency:
    geoviews
    holoviews

Class:
    Visualize()

Todo:
    * discharge hydrograph
    * visualize multiple output over time
"""
import pandas as pd
import numpy as np
import xarray as xr
import geoviews as gv
import holoviews as hv
from holoviews.operation import datashader
from bokeh.models import WMTSTileSource
from bokeh.io import save
hv.extension("bokeh")

class Visualize(object):
    """Wrapper to visualize output from lisflood-fp"""

    def __init__(self):
        
        self.mapTiles = self.__defineMap()
    
    # basic IO modules
    def __defineMap(self):
        """define a backgound map tile source"""
        from bokeh.models import WMTSTileSource
        url = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{Z}/{Y}/{X}.jpg'
        wmts = WMTSTileSource(url=url)
        mapTiles = gv.WMTS(wmts)
        return mapTiles

    def readData(filename, headerNum=6):
        """read output file from LISFLOOD-FP.
        
        Args:
            filename (str): a file path to the output.
            headerNum (int): a number of header lines.
        
        Returns:
            dataframe containing results (pandas.dataframe)
        """
        df = pd.read_csv(filename, skiprows=headerNum, header=None, sep="\s+")
        return df

    def readCache(self, filename, kind="nc"):
        """
        Read cached model domain.
        This cache should be created from the lispy.PreProcess.(mf)preprocess()

        Args:
            filename (str): a file path of cached file.  
            kind (nc): data format type

        Returns:
            lats (np.ndarray): latitudinal series of your domain.
            lons (np.ndarray): longitudinal series of your domain.
        """
        if kind == "nc":
            data = (xr.open_dataset(filename)).to_array()
            lats = data.lat.values
            lons = data.lon.values
        else:
            raise IOError("data type %s is not supported." % kind)
        return lats, lons

    def constDataArray(self, df, lats, lons, name, undef=-9999):
        """
        constract DataArray from pandas.dataframe. Mask undefined values.
        
        Args:
            df (pandas.dataframe): data frame of your 2d output.
            lats (numpy.ndarray): latitudinal series of your map domain.
            lons (numpy.ndarray): longitudinal series of your map domain.
            name (str): name of you output (e.g., width, elevation, etc.).
            undef (int or float): undefined value in your output.
        
        Returns:
            DataArray of your output (xarray.DataArray)
        """
        data = df.values
        data[data == -9999] = np.nan
        dArray = xr.DataArray(data, coords={"lat":lats, "lon":lons}, dims=["lat","lon"]).rename(name)
        return dArray

    # visualization modules
    def plotMap(self, dArray, name, width=500, height=250, cmap="gist_earth_r", alpha=0,5):
        """
        plot the DataArray onto the map.
        
        Args:
            dArray (xarray.DataArray): DataArray of your output.
            name (str): name of your data
            width (int): width of the output image
            height (int): height of the output image
            cmap (str): colormap
            alpha (float): alpha value
        
        Returns:
            bokeh image object
        """
        dataset = gv.DataSet(dArray)
        img = dataset.to(gv.Image, ["lon","lat"], name)
        img_shaded = datashader.regrid(img)
        img_out = img.opts(width=width, height=height, alpha=alpha, colorbar=True, cmap=cmap, tools=["hover"]) * self.mapTiles
        return img_out

    # higher ranked API for easy use
    def show(self, filename, name, cacheFile, undef=-9999):
        """Higher API for an instant visualization of 2D map output.
        
        Args:
            filename (str): a file path to your output
            name (str): a name of your output (e.g. width, elevation, etc.)
            cacheFile (str): a file path to your cached netcdf data (should be in cache/ directory)
            undef (int or float): undefined value in your output data

        Returns:
            image object
        """
        df = self.readData(filename)
        lats, lons = self.readCache(cacheFile)
        dArray = self.constDataArray(df, lats, lons, name, undef=undef)
        img = self.plotMap(dArray, name)
        return ing

    def saveHtml(self, img, outName):
        """save bokeh object in a outName html file."""
        save(img, outName)

