#include "prepare_wavefunctions.hpp"

#include"goscalc_constants.hpp"
#include"contwave.c"

/**
 *@file prepare_wavefunctions.cpp
 *@brief Functions to prepare the results of wavegen (source)
 */




void prepare_bound_wave(Goscalc_shared_data& gsd){
	gsd.prepd_bound_wave = gsd.cd.bound_wavefunction;

	//convert atomic units to Angstrom and convert u from wavegen to R = u/r;
	//convert R from 1/sqrt(au^3) to 1/sqrt(Ang^3)
	const double x_conversion = cst::au_to_ang;
	//const double y_conversion = 1./sqrt(pow(cst::au_to_ang,3));
	const double y_conversion = 1./sqrt(cst::au_to_ang);

	Value_pairs& bound_wave = gsd.prepd_bound_wave;
	for (unsigned int i=0; i<bound_wave.size(); i++){
		bound_wave(i,0) = bound_wave(i,0)*x_conversion;
		bound_wave(i,1) = (bound_wave(i,1)/bound_wave(i,0))*y_conversion;
    //[R(r)] = A^{-3/2} = au^{-1/2} * A^-1 * [au_to_A]^{-1/2}
    //  = u(r)/r * [au_to_A]^{-1/2} with [au_to_A] = 1A/1au
	}

}



void prepare_free_waves(Goscalc_shared_data& gsd){
	unsigned int mmax = static_cast<unsigned>(gsd.cd.potential.size());
	const unsigned int emax = gsd.cd.energy_free_steps; //no. of energy steps
	const double einc = gsd.cd.energy_free_increase; //energy increase per step
	const unsigned int max_lfree = gsd.cd.max_considered_lfree;

	//create dummyArray for use with contwave;
  // 3 dimensions: r, free energy, free l
	std::vector<std::vector<std::vector<double>>> dummy_array;
	dummy_array.resize(mmax);
	for (unsigned int i=0; i<mmax; i++){
		dummy_array[i].resize(emax+1);
		for (unsigned int e=0; e<=emax; e++){
			dummy_array[i][e].resize(max_lfree+1);
		}
	}


	//use contwave from \cite Frigge2011
	std::vector<double> pot_0 = gsd.cd.potential.column_as_vector(0);
	std::vector<double> pot_1 = gsd.cd.potential.column_as_vector(1);

	const int max_lfree_int = static_cast<int>(max_lfree);
	const int emax_int = static_cast<int>(emax); //number of steps (not actual energy value)
	const int mmax_int = static_cast<int>(mmax);
	contwave(gsd.cd.amesh,emax_int,max_lfree_int,einc,gsd.cd.energy_free_start,mmax_int,pot_0,pot_1,dummy_array);
  //stores calculated wave functions in dummy_array



  //create 2dim array: energy, free l containing value_pairs objects (giving 3rd dimension: r)
	gsd.free_waves.resize(gsd.cd.energy_free_steps+1);
	for (auto& free_waves_same_energy : gsd.free_waves){
		free_waves_same_energy.resize(max_lfree+1);
		for (Value_pairs& free_wave : free_waves_same_energy){
			free_wave = Value_pairs(gsd.cd.potential.size());
		}
	}



	//parse data to output (for one energy value) and convert atomic units (au)
	//to Angstrom and 1/(sqrt(au^3)*sqrt(hartree)) to 1/(sqrt(Angstrom^3)*sqrt(eV))
	const double x_conversion = cst::au_to_ang;
	const double y_conversion = 1./(sqrt(pow(cst::au_to_ang,3))*sqrt(cst::hart_to_eV));
	for (unsigned int i_e = 0; i_e<gsd.cd.energy_free_steps+1; i_e++){
		for (unsigned int i_lf = 0; i_lf<=max_lfree; i_lf++){
			for (unsigned int i=0; i<gsd.cd.potential.size(); i++){
				gsd.free_waves[i_e][i_lf](i,0) = gsd.cd.potential(i,0)*x_conversion; //use r from potential for lattice
				gsd.free_waves[i_e][i_lf](i,1) = dummy_array[i][i_e][i_lf]*y_conversion;
			}
		}
	}

}
