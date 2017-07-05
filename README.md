**Disclaimer: I started this projected as a personal exercise and I am still experimenting with the code. I would recommend thoroughly examining this code before adopting it for your personal use.** 


# Description

This is a Python implementation of the Price Weller Pinkel (PWP) ocean mixed layer model. This code is based on the MATLAB version of the PWP model, originally written by [Peter Lazarevich and Scott Stoermer](http://www.po.gso.uri.edu/rafos/research/pwp/) (U. Rhode Island) and later modified by Byron Kilbourne (University of Washington) and Sarah Dewey (University of Washington).

For a detailed description of the theory behind the model, see the original [Price et al. (1986)](http://onlinelibrary.wiley.com/doi/10.1029/JC091iC07p08411/full) paper. A much shorter review of the algorithm is provided in the [HYCOM documentation for the PWP](https://hycom.org/attachments/067_pwp.pdf); a google search may yield produce better sources.

The code presented here is functionally similar to its MATLAB equivalent (see *matlab_files/PWP_Byron.m*), but I have made significant changes to the code organization and flow. One major difference is that this code is split into two files: **PWP.py** and **PWP_helper.py**. 

*PWP.py* contains the core numerical algorithms for the PWP model and is mostly a line-by-line translation of the original MATLAB code. 

*PWP_helper.py* contains helper functions to facilitate model initialization, output analysis and other miscellaneous tasks. Many of these functions were introduced in this implementation.

Another major modification is the addition of a basic thermodynamic sea-ice model. This sea-ice model is coupled to the PWP is able to grow and melt ice over a full annual cycle. These algorithms are stored in *PWP_ice.py*.

Examples for setting up model runs are provided in *PWP_demos.py*.  

# Required modules/libraries
To run this code, you'll need Python 3 and the following libraries:

+ Numpy
+ Scipy
+ Matplotlib
+ xarray
+ seawater

The first three modules are available with the popular python distributions such as [Anaconda](https://www.continuum.io/downloads) and [Canopy](https://store.enthought.com/downloads/#default). You can get the other two modules via the `pip install` command from the unix command line:

```
pip install xarray
pip install seawater
```

Besides the python libraries listed here, this repository should have everything you need to do a model run with the provided datasets.


# PWP ocean model
[See jupyter notebook for details and examples](https://github.com/earlew/pwp_python/blob/with_sea_ice_v2/PWP_model.ipynb).



# Thermodynamic Sea-Ice Model
[See jupyter notebook for details and examples](https://github.com/earlew/pwp_python/blob/with_sea_ice_v2/PWP_ice_model.ipynb)



## Future work
+ Implement a way to incorporate ice-velocity in determining the ocean-ice heat flux.
+ Write a helper function to facilitate comparison between different runs.
+ Give user more control over what plots get created.
+ Finish documenting this model.


