#ifndef goscalc_shared_data_HPP
#define goscalc_shared_data_HPP

#include<string>
#include<vector>

#include<armadillo>

#include"value_pairs.hpp"

/**
 *@file goscalc_shared_data.hpp
 *@brief Struct to share input data and results throughout goscalc
 */


///Data read from the config file and from the DFT file
struct Config_data{
	std::string config_filename;					///<Name of config file
	std::string dft_filename;						///<Name of file with DFT results from wavegen
	std::string output_dir_name;					///<Name of the output directory


	unsigned int n_bound;                           ///<Principal quantum number of the bound atomic electron
	unsigned int l_bound;                           ///<Azimuthal quantum number of the bound atomic electron
	unsigned int max_considered_lfree;				///<Highest azimuthal quantum number of the free electron considered in the calculation
	double energy_free_start;                       ///<Lowest energy for the free atomic electron for which GOS is calculated in eV
	unsigned int energy_free_steps;                 ///<Amount of energy steps for the free atomic electron for which GOS is calculated
	double energy_free_increase;                    ///<Increase of the energy of the free atomic electron per step in eV

	double max_kvalue;	          			    ///<Maximum k-value to be computed in Angstrom

	double amesh;									///<amesh parameter used in the calculation of the wave function of the free atomic electron (see wavegen_mod)
	unsigned int atomic_number;						///<Amount of protons in the nucleus (from DFT file)
	double bound_energy;             				///<Energy of the bound atomic electron (from DFT file)
	Value_pairs bound_wavefunction;	                ///<xy-value of the wave function of the bound atomic electron (from DFT file)
	Value_pairs potential;					        ///<xy-values of the atomic potential (from DFT file)
};








///Struct for data sharing within the program
struct Goscalc_shared_data{
	///Allow construction of shared storage only with a config file
	Goscalc_shared_data(const Config_data& config) : cd(config){};

	///Delete copy-assignement operator (move equivalent implicitly deleted)
	Goscalc_shared_data & operator = (const Goscalc_shared_data&) = delete;

	///Delete copy operator (move equivalent implicitly deleted)
	Goscalc_shared_data(const Goscalc_shared_data&) = delete;

	Value_pairs prepd_bound_wave;
	std::vector<std::vector<Value_pairs>> free_waves;

	///Data read from config file
	const Config_data cd;


	///All transition potentials calculated within a single program run
	std::vector<double> results_energy_losses;
	std::vector<Value_pairs> results_gos;
};



#endif
