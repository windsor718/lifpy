# LISFLOOD-FPy (lispy)  
This is a Python toolkit to prepare data for, handle, and visualize the results from the hydraulic model, LISFLOOD-FP.  
## 1. Preparing your data for the LISFLOOD-FP.  
You need to prepare geographical (DEM, river width, etc.) and hydrological (input hydrograph) boundary condition to run LISFLOOD-FP.  
lispy enables to reduce the pain to prepare those inputs.  
### 1.1 Making geographical boundary condition.  
lispy.PreProcess() instance has toolkits to create the geographical forcing to run LISFLOOD-FP. The operation is done by lazy lodeing and parallelization using the Dask and Xarray, resulting in fast and memory efficient processing even with large data.  
#### What you need.  
* HYDROGRAPHY DATA:  
  * DEM (Hydrologically corrected), river width, and upstream area in a GeoTiff format.  
* RAM storage enough to store at least a single map of your whole domain.  
  * Eventhough lispy try to parallelize your domain, the module needs to dump your data into ARC ascii format. This is the time when the module compile and load whole domain onto your RAM. Although there is a way to parallelie this process, curremt version only supports dumping your domain by storing whole domain. (You need to secure more memory storage to run the model anyway!)  
#### Instruction:
lispy.PreProcess.preprocess(upareaFile,elevationFile,widthFile,threshold,args*)  
Create a geographical input files from single tile of hydrography data.  
Positional arguments:  
+ upareaFile(str): file path to the uparea file.  
+ elevationFile(str): file path to the elevation file.  
+ widthFile(str): file path to the width file.  
+ threshold(int or float): upstream area threshold to extract rivers. The raster pixels having upstream area above this threshold will be a river.  
  
Optional arguments:  
+ domain(list): Default None. your model domain [llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon]. The output files will be sliced to fit this domain. Note that nearest value within this domain will be selected and no interpoltion will be happened to secure the hydrography infromation.  
+ latName(str): Default y. The dimension name in the latitudinal axis.  
+ lonName(str): Default x. The dimension name in the longitudinal axis.  
+ bandNum(int): Default 0. The band number that stores your data.  
  
Note that no parallelization is implemented in this function. If your domain is small enough compared to an original hydrography file (usually devided into tiles), using this function may be efficient.  
  
lispy.PreProcess.mfpreprocess(upareaFiles,elevationFiles,widthFiles,threshold,nCols,nRows,args*)  
Same as the preprocess(), but implemented with parallelization. Suitable for the large domain. The mfpreprocess() can also read multiple input sources and concatenate them based on nCols and nRows information. The input file list should be aligned with C order.  
Positional arguments:
+ upareaFiles(list): list pf file paths to the uparea file.
+ elevationFiles(list): list of file paths to the elevation file.
+ widthFiles(list): list of file paths to the width file.
+ threshold(int or float): upstream area threshold to extract rivers. The raster pixels having upstream area above this threshold will be a river.
+ nCols(int): number of tiles (files) in columns (longitudinal axis).  
+ nRows(int): number of tiles (files) in rows (latitudinal axis).  
Optional arguments:
+ domain(list): Default None. your model domain [llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon]. The output files will be sliced to fit this domain. Note that nearest value within this domain will be selected and no interpoltion will be happened to secure the hydrography infromation.
+ latName(str): Default y. The dimension name in the latitudinal axis.
+ lonName(str): Default x. The dimension name in the longitudinal axis.
+ bandNum(int): Default 0. The band number that stores your data.
