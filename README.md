# NMRtoSTL
Converts 2D NMR data to an STL file for 3D printing (Currently only works with TopSpin totxt outputs)

Inspired by https://doi.org/10.1021/acs.jchemed.0c01130 and https://doi.org/10.1002/cmr.a.21470

# Use

Select a 2D Spectum you wish to convert and open it in TopSpin. 
Make sure you are fully zoomed out and use the totxt command. Select the location and save.
Then run `python NMRtoSTL.py <yourOutputTextFile>` from your command line and it will output an STL file with the same name.

# Features to add
- [ ] X and Y limits
- [ ] read from generic files

# Known issues
Makes a mess of spectra that have been zoomed in. 
