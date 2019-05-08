# LISFLOOD-FPy (lifpy)  
lifpy is a Python-based wrapper to prepare data, handle the model, and visualize results from the hydraulic model, LISFLOOD-FP.  
lifpy is developed especially for the use of large scale modeling and benefitted from lazy loading and parallelization function of xarray and dask.  
## Setup the environment  
The dependency of lifpy is:  
- xarray  
- rasterio
- dask  
- numpy
- cython  
- pandas  
- geoviews  
- holoviews  
- matplotlib  
- bokeh  
  
If you are using Anaconda,  

```conda create -n lifpy python=3.6 xarray rasterio geoviews holoviews```  

should create a necessary environment.  
Instead, you can also use enviromnent.yml in this repository to create a environment:  

```conda create -n environment.yml```  

## Prepare for the simulation  
Please also see the ```demo.ipynb``` for the specific example to use this library.  
  
First and foremost, you need to prepare the input files for the simulation. The current version of lifpy needs the following datasets in specific format:
- __Hydrography data (GeoTiff)__
  * Flow direction  
  * River width  
  * Upstream area size  
  * Surface elevation  
- __River discharge data (netCDF4)__
  * shoud contain time and id (river identification number or string) coordinates.  
  
Once your data is ready:
```python
import lifpy as lfp  

# define required arguments  
elevPath = "list of path to the surface elevation files"  
upaPath = "list of path to the surface upstream area size files"  
wthPath = "list of path to the surface river width files"  
thsld = "the minimum size of upstrea area. The pixel whose upsream area size is above this number will be considered as rivers"  
nCols = "number of files for a longitudinal axis"  
nRaws = "number of files for a latitudinal axis"  

lfp.PreProcess.mfpreprocess(upaPath,elvPath,wthPath,thsld,nCols,nRaws)  
```
will create a topography data required for your lisflood-fp simulation.  
  
Next, the input forcing dataset is required. The current version of lifpy needs the following datasets in a specific format:  
- __Discharge data (netCDF)__  
  * discharge time series at each point (2d array)
  * time index (dimension 1)
  * river reach id (dimension 2)
- __Discharge point information (csv)__  
This is a file to define the coordinates (lat/lon) and upstream area of each river id in your river discharge data.  
  * The sample format is:  
  
    |id|lat|lon|uparea|
    |---|---|---|---|  
    |0|35.11|-120.24|300.25|
  
After you prepared the data, you can make forcing data with:
```python
dschgFile = "path to your discharge data"
pointInfoFile = "path to your point info data"

lifpy.Forcing.makeForcing(dschgFile, pointInfoFile)
```
These are the files you need to run the LISFLOOD-FP (subgrid).  
  
# Visualizing results  
Current lifpy only supports the instant visualization of a snapshot (at specific time) of output file (e.g. res-001.txt). Please see the Todo in the document to see the future implementation.  
  
lifpy has a higher API for instant visualization to check your simulation:
```python
fileName = "your results (.txt) from LISFLOOD-FP"
name = "name of the result (e.g. width, elevation, etc.)"
cacheFile = "path to the cached netCDF file that lifpy.PreProcess.mfpreprocess generates."
# in default it is in cache/uparea.nc  
img = lifpy.Visualize.show(fileName, name, cacheFile)
```
This function uses a datashader as a pipeline, and those grids in the figure is dynamically regirdded which enables smooth loading.  
This is useful when you see a big picture of your simulation (which has more than a millions of data points in a large scale hi-res modeling) which is otherwise very expensive to visualize.  

This will output a Bokeh interactive plot, and you can also save as a html link with:  
```python
import bokeh.io
outName = "outputName.html"
bokeh.io.save(img, outName)
```
or just using a simple wrapper of lifpy:
```python
lifpy.Visualize.save(img, outName)
```
You can save as a png or svg plot directly with bokeh library, but usually it is useful to open the html link and adjust the plot whatever you like, and save the plot as png from that interective plot.  
