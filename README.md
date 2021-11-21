# NMRtoSTL
Converts 2D NMR data to an STL file for 3D printing.

Inspired by https://doi.org/10.1021/acs.jchemed.0c01130 and https://doi.org/10.1002/cmr.a.21470

# Use
1. Find the folder where you NMR data is saved. Currently supported formats are:
  - Bruker 1D processed data (1r) [minimum 3 spectra]
  - Bruker 2D processed data (2rr)
  - NMRpipe 2D processed data (.ft2)

2. Copy the file path. In the case of Bruker data, the file path of the parent folder containing the 'fid' file should be used.
3. Run `python NMRtoSTL.py <yourFilePath>` from your command line and it will output an STL file with the same name.

# Features to add
- [x] X and Y limits
- [ ] read from generic files (in progress)
- [ ] reduce number of mesh points
- [ ] add arguments for limits, size, stacking, smoothing and threshold to command line input

# Known issues
- Read from Agilent/Varian data still work in progress.
- Fusion 360 reports errors: 'Mesh is not closed' and 'Mesh does not have a positive volume' with the .stl file. The current workaround is to run mesh repair with the 'wrap' option in Fusion.
