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

    def __init__(self):
        
        self.mapTiles = self.__defineMap()
    
    # basic IO modules
    def __defineMap(self):
        from bokeh.models import WMTSTileSource
        url = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{Z}/{Y}/{X}.jpg'
        wmts = WMTSTileSource(url=url)
        mapTiles = gv.WMTS(wmts)
        return mapTiles

    def readData(filename, headerNum=6):
        df = pd.read_csv(filename, skiprows=headerNum, header=None, sep="\s+")
        return df

    def readCache(self, filename, kind="nc"):
        """
        Read cached model domain.
        This cache should be created from the lispy.PreProcess.(mf)preprocess()

        Args:
            filename (str): a file path of cached file.  
            kind (nc): data format type
        """
        if kind == "nc":
            data = (xr.open_dataset(filename)).to_array()
            lats = data.lat.values
            lons = data.lon.values
        else:
            raise IOError("data type %s is not supported." % kind)
        return lats, lons

    def constDataArray(self, df, lats, lons, name, undef=-9999):
        data = df.values
        data[data == -9999] = np.nan
        dArray = xr.DataArray(data, coords={"lat":lats, "lon":lons}, dims=["lat","lon"]).rename(name)
        return dArray

    # visualization modules
    def plotMap(self, dArray, name, width=500, height=250, cmap="gist_earth_r", alpha=0,5):
        dataset = gv.DataSet(dArray)
        img = dataset.to(gv.Image, ["lon","lat"], name)
        img_shaded = datashader.regrid(img)
        img_out = img.opts(width=width, height=height, alpha=alpha, colorbar=True, cmap=cmap, tools=["hover"]) * self.mapTiles
        return img_out

    # higher ranked API for easy use
    def show(self, filename, name, cacheFile, undef=-9999):
        df = self.readData(filename)
        lats, lons = self.readCache(cacheFile)
        dArray = self.constDataArray(df, lats, lons, name, undef=undef)
        img = self.plotMap(dArray, name)
        return show

    def saveHtml(self, img, outName):
        save(img, outName)

