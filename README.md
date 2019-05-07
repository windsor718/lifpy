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

```conda create -n lifpy python=3.6 xarray rasterio geiviews holoviews```  

should create a necessary environment.  
Instead, you can also use enviromnent.yml in this repository to create a environment:  

```conda create -n environment.yml```  

## Prepare for the simulation  
First and foremost, you need to prepare the input files for the simulation. The current version of lifpy needs the following datasets in specific format:
- __Hydrography data (GeoTiff)__
  * Flow direction  
  * River width  
  * Upstream area size  
  * Surface elevation  
- __River discharge data (netCDF4)__
  * shoud contain time and id (river identification number or string) coordinates.  
Once your data is ready:
```import lifpy.PreProcess as lfp```  
```elevPath = "list of path to the surface elevation files"```  
```upaPath = "list of path to the surface upstream area size files"```
```wthPath = "list of path to the surface river width files"```  
