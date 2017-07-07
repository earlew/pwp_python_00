**Disclaimer: I started this projected as a personal exercise and I am still experimenting with the code. I would recommend thoroughly examining this code before adopting it for your personal use.** 


# Description

This is a Python implementation of the Price-Weller-Pinkel (PWP) 1-D ocean mixing model.

The code in this repository has its roots in the original Fortran implementation (c. 1989-1999) written by Jim Price (WHOI), which is available [here](http://www.whoi.edu/science/po/people/jprice/website/projects_upperocean.html). Price adapted the Fortran code into the first MATLAB version of the model (2001). A translation into MATLAB 5.X by Peter Lazarevich and Scott Stoermer (University of Rhode Island) is available [here](http://www.po.gso.uri.edu/rafos/research/pwp/). Subsequent modifications were added by Byron Kilbourne (University of Washington) and Sarah Dewey (UW). This implementation is the work of Earle Wilson (UW), with minor contributions from Ethan Campbell (UW).

For the theory behind the model, see the original description in [Price et al. (1986)](http://onlinelibrary.wiley.com/doi/10.1029/JC091iC07p08411/full). A short review of the algorithm is provided in the [HYCOM documentation for the PWP model](https://hycom.org/attachments/067_pwp.pdf); a Google search may yield better sources.

The code presented here is functionally similar to its MATLAB equivalent (see *matlab_files/PWP_Byron.m*), but I have made significant changes to the code organization and flow. One major difference is that this code is split into two files: **PWP.py** and **PWP_helper.py**. 

*PWP.py* contains the core numerical algorithms for the PWP model and is mostly a line-by-line translation of the original MATLAB code. 

*PWP_helper.py* contains helper functions to facilitate model initialization, output analysis, and other miscellaneous tasks. Many of these functions were introduced in this implementation.

Another major modification is the addition of a basic thermodynamic sea-ice model. This coupled sea-ice model allows ice to grow and melt over a full annual cycle. These algorithms are stored in *PWP_ice.py*.

Examples for setting up model runs are provided in *PWP_demos.py*.  

# Required modules/libraries
To run this code, you'll need Python 3 and the following libraries:

+ Numpy
+ Scipy
+ Matplotlib
+ [xarray](http://xarray.pydata.org/en/stable/)
+ [GSW-Python](https://teos-10.github.io/GSW-Python/)
+ [Python-Seawater](http://pythonhosted.org/seawater/eos80.html)

The first three modules are available with popular Python distributions such as [Anaconda](https://www.continuum.io/downloads) and [Canopy](https://store.enthought.com/downloads/#default). You can get the other three modules via the `pip install` command from the UNIX command line:

```
pip install xarray
pip install gsw
pip install seawater
```

Besides the Python libraries listed here, this repository should have everything you need to do a model run with the provided datasets.


# PWP ocean model
[See Jupyter notebook for details and examples](https://github.com/earlew/pwp_python/blob/with_sea_ice_v2/PWP_model.ipynb).



# Thermodynamic sea-ice model
[See Jupyter notebook for details and examples](https://github.com/earlew/pwp_python/blob/with_sea_ice_v2/PWP_ice_model.ipynb).



## Future work
+ Implement a way to incorporate ice-velocity in determining the ocean-ice heat flux.
+ Write a helper function to facilitate comparison between different runs.
+ Give user more control over what plots get created.
+ Finish documenting this model.


