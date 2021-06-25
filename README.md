# GOScalc

## Introduction

This is a program to calculate the total generalized oscillator strength for a certain (sub-)orbital.
For specific information refer to [1] (German).  
The additional program Wavegen is based on [2] and was provided by Prof. Krüger, who modified it together with M.Frigge [3].

The required class hankel_trafo.cc was excluded in the Gitlab publication because of the restrictive CPC license as it is an adapted version of
NumSBT (aanz_v2), which was written in Fortran90 and published by Peter Koval and J.D. Talman in Computer Physics Communications in 2010.
This was published under the CPC licence which forbids republication even in adapted form unless permission is given (which I didn't request).
The original version is NumSBT (aanz_v2), which was written in Fortran90 and published by Peter Koval and J.D. Talman in Computer Physics Communications in 2010. (CPC license)  
	 v2: [https://data.mendeley.com/datasets/y294ttxyw4/1](https://data.mendeley.com/datasets/y294ttxyw4/1)  
	 v3: [https://data.mendeley.com/datasets/m3fc83rytv/1]( https://data.mendeley.com/datasets/m3fc83rytv/1)  
	 ([corresponding article](https://www.sciencedirect.com/science/article/pii/S0010465508003329 ); here the link to the program is broken)
On request I will gladly provide hankel_trafo.cc or you can use any other program to do the hankel/spherical bessel transformation or incorporate the fortran program above. goscalc will not function without it.

## Installation and compiling (Linux systems):
Libraries used: Armadillo (requiring OpenBLAS or standard BLAS and LAPACK), FFTW, Boost, WignerSymbols. WignerSymbols can be found [here](https://github.com/joeydumont/wignerSymbols)<!--- TODO: is superlu needed?-->  
Make sure you have cmake, make and a fortran and c++ compiler installed (gfortran, gcc).

Create directory build (if it doesn't exist already) and cd into it:
```bash
	mkdir build/
	cd build/
```
Compile with cmake and make:
```bash
    cmake ..
     make
```
Compile wavegen_mod with:
```bash
	gfortran -std=legacy -o wavegen_mod wavegen_mod.f
```

## Usage

+ create the wavegen.dat file:
    + The config file for wavegen should be named wavegen.dat.
    +	It tells the program which exchange-correlation (XC) functional to use,
        the atomic number (Z) of the simulated atom and its electron configuration.
        + Options for the XC functional are LDA (Local-Density Approximation) and GGA (Generalized Gradient Approximation). LDA is recommended.
	+ The electron configuration is given in rows, where each row contains the principal quantum number (n),
	   the azimuthal quantum number (l) and the occupation number (ON) separated by spin direction (ONup and ONdown).
    + Conceptually, the wavegen.dat then looks like this:


            XC-functional
            Z
            n l ONup ONdown
            n l ONup ONdown
            ...............
            n l Onup Ondown

    + With the quantum numbers n,l and the number of electrons with spin up or down for that subshell.
	   Empty shells should not be listed.

+ Example files are provided for copper, silicon, carbon and lead.
	+ These need to be renamed to wavegen.dat for usage.
	+ Notice: The spin is not maximized correctly as according to Hund's second rule.
        + 	This is because goscalc doesn't include spin effects.
        	     So instead spin is equally distributed; with the excess electron for odd numbers of electrons put as spin down because later the spinup file is used
            	(reduces effect of imbalance).
	Electron configurations can be looked up [here](https://sciencenotes.org/list-of-electron-configurations-of-elements/).
+ execute wavegen in the same directory
```bash
    ./wavegen
```
+ place the waveup.dat file in the same directory as the config.json and the goscalc executable
+ fill in config.json with the desired parameters and the output directory name
+ execute goscalc
```bash
    ./goscalc
```
+ alternatively you can pass the path to the config file as a command line argument
```bash
    ./goscalc /path/to/config.json
```

### Output:
	The output directory includes a copy of the config file, the command line log,
    the k values in k.dat, the corresponding generalized oscillator strengths in gos.dat for the
	energy losses, which result from the desired free energies saved in free_energies.dat.

### Additional Directories:
+ element_configs contains some example configuration files to use with wavegen and goscalc for different elements

### Restrictions:
+ GOS for ions can not be calculated, because their atomic potential does not fall off to zero within the mesh given by wavegen (or at all, technically), but contwace.c requires this.
+ the number of mesh points is hard coded into wavegen. It can be changed by changing mmax. It should probably be a power of 2, if not just for efficiency reasons.

## Known Issues:
+ The calculation of the atomic wave functions by wavegen leads to a relatively (relative to the value of the wave function at this r) significant (and unphysical) jump in the potential.
	As this only happens at higher r (>5 Angstrom), where the potential is small anyways, it is assumed, that this doesn't lead to a significant error
+ When compared to the GOS tables used by Gatan's EELS Analysis (2.3.2) there is a significant qualitative difference in the GOSs, while the overall shape is very similar.
	The difference appears to be stronger for higher l, though this hasn't been tested rigorously. It is unclear where this stems from.

## Bibliography
[1] Leonhard Segger, Berechnung generalisierter Oszillatorenstärken für die Quantifizierung von EEL-Spektren, Bachelorarbeit, WWU-Münster 2019 <!---TODO: Add link to AG Website as soon as thesis is up-->
[2] Hamann, Phys.Rev. B 40 (1989), 2980  
[3] Frigge, Kohl, Krüger, Microscopy Conference 2011 (Kiel), IM5.P174, <!--- TODO: Citation from proceedings journal -->
	Calculation of relativistic differential cross-sections for use in microanalysis. Abstract available at [https://www.uni-muenster.de/imperia/md/content/physik_pi/kohl/mc2011/im5_p175.pdf](https://www.uni-muenster.de/imperia/md/content/physik_pi/kohl/mc2011/im5_p175.pdf)
