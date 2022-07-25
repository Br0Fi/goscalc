## Scripts

These scripts have the purpose of automating the generation of input and collection of output for `wavegen` and `goscalc`, permitting to generate the GOS data in bulk for many edges.
Two of these scripts currently depend on the Panta Rhei software developed by CEOS, which provides an EELS Atlas with a list of edges and corresponding energy. 
Since this software package is not yet publicly available, it should be replaced in the future with hyperspy, which contains similar data.

#### goscalc_confgen.py

This is used to generate the configuration files for wavegen and goscalc for a give element or edge.

Dependencies:
numpy
matplotlib
mendeleev - Provides electronic configuration data for the different elements
panta_rhei - Provides EELS Atlas
json - support for the json format used by goscalc


#### compute_all_gos.py

This script is used to run call wavegen and goscalc to compute the GOS for all edges included in the EELS atlas.

Dependencies:
mendeleev - used to associate atomic numbers to symbols
panta_rhei - Provides EELS atlas
goscalc_confgen - generates the configurations for the binary executables
subprocess - used to launch the binary executables


#### assemble_gos5_file.py

This script assembles all computed GOS into a gos5 archive. For more information on the gos5 file format please see https://gitlab.com/gguzzina/gos5

Dependencies:
numpy
h5py - support for hdf5 format used by the gos5 specification
json - used to load the metadata from the metadata.json file
