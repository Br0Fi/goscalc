GOScalc
=======

This is a program to calculate the total generalized oscillator strength for a certain (sub-)orbital.
	For specific information refer to [1] (German).
The additional program Wavegen is based on [2] and was provided by Prof. Krüger, who modified it.

The required class hankel_trafo.cc was excluded in the Gitlab publication because of the restrictive cpc license as it is an adapted version of
NumSBT (aanz_v2), which was written in Fortran90 and published by Peter Koval and J.D. Talman in Computer Physics Communications in 2010.
This was published under the CPC licence which forbids republication even in adapted form unless permission is given (which I didn't request).
The original version is NumSBT (aanz_v2), which was written in Fortran90 and published by Peter Koval and J.D. Talman in Computer Physics Communications in 2010. (CPC license)
	 v2: https://data.mendeley.com/datasets/y294ttxyw4/1
	 v3: https://data.mendeley.com/datasets/m3fc83rytv/1
	 (corresponding article: https://www.sciencedirect.com/science/article/pii/S0010465508003329 ; here the link to the program is broken)
On request I can provide hankel_trafo.cc or you can use any other program to do the hankel/spherical bessel transformation or incorporate the fortran program above.

##Installation and compiling (Linux systems):
+ Libraries used: Armadillo, FFTW, Boost, WignerSymbols
-TODO: meinen Installations-Zettel finden und abschreiben, wenn das noch nötig ist

+ create directory build and cd into it
```bash 
	mkdir build/
	cd build/
```
+ compile with cmake and make
```bash
	cmake && make
```
+ compile wavegen_mod with
```bash
	gfortran -std=legacy -o wavegen_mod wavegen_mod.f
```

## Usage:
+ create the wavegen.dat file:
	The config file for wavegen should be named wavegen.dat.
	It tells the program which exchange-correlation (XC) functional to use,
	the atomic number (Z) of the simulated atom and its electron configuration.
	The electron configuration is described in rows, where each row contains the principal quantum number (n),
	the azimuthal quantum number (l) and the occupation number (ON) separated by spin direction (ONup and ONdown).
	Conceptually, the config file for wavegen looks like this:
	
	XC-functional
	Z
	n l ONup ONdown
	n l ONup ONdown
	...............
	n l Onup Ondown
	
	With the quantum numbers n,l and the number of electrons with spin up or down for that subshell.
	Empty shells should not be listed.
	

	Example files are provided for copper, silicon, carbon and lead.
	These need to be renamed to wavegen.dat for usage.
	Notice: The spin is not maximized correctly as according to Hund's second rule.
          This is because goscalc doesn't include spin effects.
          So instead spin is equally distributed; with the excess electron for
          odd numbers put as spin down because later the spinup file is used
          (reduces effect of imbalance).
    Electron configurations can be looked up here: https://sciencenotes.org/list-of-electron-configurations-of-elements/
+ execute wavegen in the same directory
+ place the waveup.dat file in the same directory as the config.json and the goscalc executable
+ fill in config.json with the desired parameters and the output directory name
+ execute goscalc.
+ alternatively you can pass the path to the config file as a command line argument

## Output:
	The output directory includes a copy of the config file, the command line log,
	the k values in k.dat, the corresponding generalized oscillator strengths in gos.dat for the 
	energy losses, which result from the desired free energies saved in free_energies.dat.
	
## Restrictions:
+ GOS for ions can not be calculated, because their atomic potential does not fall of to zero within the mesh given by wavegen
+ the number of mesh points is hard coded into wavegen. It can be changed by changing mmax. It should probably be a power of 2, if not just for efficiency reasons.

##additional directories:
+ autogen is supposed to generate the electron configuration (wip, not currently usable)
+ reference_code includes code by Frigge, Majert and Talman for the spherical Hankel transformation

##Known Issues:
+ the calculation of the atomic wave functions by wavegen leads to a relatively (relative to the value of the wave function at this r) significant (and unphysical) jump in the potential.
	As this only happens at higher r (>5 Angstrom), where the potential is small anyways, it is assumed, that this doesn't lead to a significant error
+ When compared to the GOS tables used by Gatan's 

# Bibliography
[1] Leonhard Segger, Berechnung generalisierter Oszillatorenstärken für die Quantifizierung von EEL-Spektren, Bachelorarbeit, WWU-Münster 2019
[2] Hamann, Phys.Rev. B 40 (1989) 2980
