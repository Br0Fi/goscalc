#include<string>

#include<fftw3.h>

#include"goscalc_shared_data.hpp"
#include"logger.hpp"
#include"goscalc_io.hpp"
#include"goscalc_constants.hpp"
#include"prepare_wavefunctions.hpp"
#include"calc_gos.hpp"
#include"hankel_trafo.hpp"



/**
 *@file main.cpp
 *@brief Main file
 */



///The main
int main(int argc, char* argv[]){
	std::string config_filename;
	if (argc==1){
		config_filename = "config.json";
	} else if (argc==2) {
		config_filename = std::string(argv[1]);
	} else {
		throw goscalc_io_exception("More than one config-file passed as command line argument!");
	}

	logg::out()<<"This is GOScalc version "<<VERSION<<".\n";
	auto start_time =  std::chrono::system_clock::now();


	//save on planing time for the Fourier Transforms by importing wisdom
	if (fftw_import_wisdom_from_filename("fftw.wis")==1){logg::out()<<"reading fftw-wisdom...\n";}


	logg::out()<<"reading config data from "<<config_filename<<"...\n";
	Goscalc_shared_data gsd(read_input(config_filename));

	prepare_bound_wave(gsd);
	prepare_free_waves(gsd);

	gsd.results_energy_losses = calc_energy_losses(gsd);

	const unsigned int no_energy_losses = static_cast<unsigned int>(gsd.results_energy_losses.size());
	gsd.results_gos.resize(no_energy_losses);
	const std::string l_string = cst::subshell_abbrv.at(gsd.cd.l_bound);
	const std::string element = cst::elemental_abbrv.at(gsd.cd.atomic_number);
	const std::string shell_name = element + "-" +std::to_string(gsd.cd.n_bound)+l_string;
	const double electrons_in_shell = (2*gsd.cd.l_bound+1)*2;
	logg::out()<<"Assuming fully filled "<<shell_name<<" (sub)shell with "<<electrons_in_shell<<" electrons.\n";
	logg::out()<<"The energy of the bound atomic electron is "<<gsd.cd.bound_energy<<"eV.\n";

	unsigned int max_lambda = gsd.cd.max_considered_lfree + gsd.cd.l_bound;
  Hankel_trafo trafo_plan = Hankel_trafo(gsd.prepd_bound_wave.column_as_vector(0), max_lambda, gsd.cd.max_kvalue);
	auto ini_time = std::chrono::system_clock::now();
	auto ini_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(ini_time-start_time);
	logg::out()<<"Initialization time was "<<ini_elapsed.count()<<" milliseconds.\n";
	for (unsigned int i_e = 0; i_e<no_energy_losses; i_e++){
		const double e_free = gsd.cd.energy_free_start+i_e*gsd.cd.energy_free_increase;
		logg::out()<<"working on free energy "<<e_free<<"eV...\n";
		gsd.results_gos[i_e] =  calculate_gos(i_e, trafo_plan, gsd)*electrons_in_shell;
	}


	auto end_time = std::chrono::system_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time-start_time);
	logg::out()<<"Execution time was "<<elapsed.count()<<" milliseconds.\n";
	logg::out()<<"saving results...\n";
	save_results(gsd);

	//export gathered fftw wisdom
	if (fftw_export_wisdom_to_filename("fftw.wis")==1){logg::out()<<"exporting fftw-wisdom...\n";}






}
