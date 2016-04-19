### PWP

Python implementation of the Price Weller Pinkel (PWP) mixed layer model. This is based on the MATLAB implementation of the PWP model written by Byron Kilbourne (University of Washington). I did this re-write as a personal exercise, so I would be cautious about using this code to do any serious work. 

The code presented here is functionally similar to */matlab_version/PWP_Byron.m*. However, I made significant changes to the code structure and organization. One big difference is that this code is split into two files. Them main file, *PWP.py*, contains the core numerical algorithms for the PWP model. These numerical algorithms are essentially direct translations of their MATLAB equivalents. The second file, *PWP_helper.py*, contains helper functions to facilitate model initialization, output analysis and other miscellaneous tasks. 

To run the code, you can type `%run PWP.py` from the iPython command line. This calls the PWP.run() function, which is the main function for the script. Alternatively, you can import PWP.py as a module then run the model directly:

```
import PWP
PWP.run()
```

This approach allows you to modify the model settings and parameters.

```
PWP.run(met_data='somewhere_else.nc', overwrite=False, diagnostics=False)
```

To get a feel for how this code/model is organized, the `PWP.run()` would be a good place to start. This function has detailed doc-file, so it should be self-explanatory. 

This repository also contains sample surface forcing and initial profile data files. Therefore, if you copy this repository to your local directory, you should be able to run the model immediately. 

The code was written using Python version 2.7.6.
