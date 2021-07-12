#include"goscalc_io.hpp"

#include<fstream>
#include<sstream>
#include<math.h>

#include<armadillo>
#include<boost/property_tree/ptree.hpp>
#include<boost/property_tree/json_parser.hpp>
// currently (2021-07-12) causes compiler warning, seems to be issue with boost's json_parser:
// see https://github.com/boostorg/property_tree/issues/51
#include<boost/format.hpp>
#include<boost/math/special_functions/round.hpp>
#include<boost/algorithm/string.hpp>

#include"path_switcher.hpp"
#include"goscalc_constants.hpp"
#include"logger.hpp"

/**
 *@file goscalc_io.cpp
 *@brief Functions for input and output (source)
 */

namespace{
	Config_data read_config_file(std::string config_filename){
		namespace pt = boost::property_tree;

		pt::ptree input;
		pt::read_json(config_filename.c_str(), input);

		Config_data cd;
		cd.config_filename = config_filename;
		cd.dft_filename = input.get<std::string>("dft_filename");
		cd.output_dir_name = input.get<std::string>("output_dir_name");

		cd.n_bound = input.get<unsigned int>("n_bound");
		cd.l_bound = input.get<unsigned int>("l_bound");
		if (cd.n_bound<=cd.l_bound){
			throw goscalc_io_exception("n_bound must be larger than l_bound.");
		}
		cd.max_considered_lfree = input.get<unsigned int>("max_considered_lfree");
		cd.energy_free_start = input.get<double>("energy_free_start");
		cd.energy_free_steps = input.get<unsigned int>("energy_free_steps");
		cd.energy_free_increase = input.get<double>("energy_free_increase");
		cd.max_kvalue = input.get<double>("max_kvalue_Ang");

		return cd;
	}


	unsigned int read_atomic_number(std::string dft_filename){
		std::ifstream file(dft_filename.c_str());
		if (!file.is_open()){
			throw goscalc_io_exception("Wavefunction input file not found.");
		}

		//read atomic number string
		std::string dummy_string;
		std::getline(file, dummy_string);
		file.close();

		//pull atomic number from string
		std::stringstream dummy_stream;
		dummy_stream.str(dummy_string);
		dummy_stream>>dummy_string; //skip text
		unsigned int atomic_number;
		dummy_stream>>atomic_number;
		return atomic_number;
	}


	std::vector<double> read_bound_energies(std::string dft_filename){
		std::ifstream file(dft_filename.c_str());
		if (!file.is_open()){
			throw goscalc_io_exception("Wavefunction input file not found.");
		}

		//read all values of the bound energy
		std::string dummy_string;
		std::getline(file, dummy_string);//skip atomic number
		std::getline(file, dummy_string);
		boost::trim(dummy_string);
		file.close();

		//search through stringstream to find correct energy
		std::stringstream dummy_stream;
		dummy_stream.str(dummy_string);
		dummy_stream>>dummy_string; //skip "eigenvalues_eV"
		std::vector<double> bound_energies;
		double dummy_double;
		while (!dummy_stream.eof()){
			dummy_stream >> dummy_double;
			bound_energies.push_back(dummy_double);
		}
		return bound_energies;
	}


	unsigned int index_nl_from_n_and_l(unsigned int n, unsigned int l){
		unsigned int index_nl = 0;
		for (unsigned int i=1; i<n; i++){
			index_nl += i;
		}

		return index_nl+l;
	}


	void read_wavegen_file(Config_data& cd){
		const std::string dft_filename = cd.dft_filename;
		const unsigned index_nl = index_nl_from_n_and_l(cd.n_bound,cd.l_bound);

		cd.atomic_number = read_atomic_number(dft_filename);
		const std::vector<double> bound_energies = read_bound_energies(dft_filename);
		cd.bound_energy = bound_energies[index_nl];
		cd.potential = Value_pairs(dft_filename,1);
		if (index_nl>=bound_energies.size()){
			std::string err_msg = "n and l not included in file " + dft_filename + ".";
			throw goscalc_io_exception(err_msg);
		}
		const unsigned col_no_bound_wave = index_nl+1+1; //+1 for x-values and +1 to skip potential
		cd.bound_wavefunction = Value_pairs(dft_filename,col_no_bound_wave);
	}



	void set_amesh(Config_data& cd){
    //hard-coded version, dependend on implementation in wavegen.f:
		//const double at_num = static_cast<double>(cd.atomic_number);
		//const double pot_size = static_cast<double>(cd.potential.size());
		//cd.amesh = pow(45.*160.*at_num, 1./pot_size);
    //better, dynamic version:
    const double rmax = cd.bound_wavefunction(static_cast<unsigned int>(cd.potential.size())-1, 0);
    const double rmin = cd.bound_wavefunction(0,0);
    const double pot_size = static_cast<double>(cd.potential.size());
    cd.amesh = std::pow((rmax/rmin), 1/(pot_size-1));
    // = exp(dr) with dr = ln(rmax/rmin)/(mmax-1)
	}


}


Config_data read_input(std::string config_filename){
	Config_data cd = read_config_file(config_filename);
	read_wavegen_file(cd);
	set_amesh(cd);

	return cd;
}

namespace{

	///creates the output directory \p output_dir_name in \p work_path
	/*!An output_directory with the name \p output_dir_name is created
   	   in the work_path \p work_path. If the directory already exists, a one, two, etc. is
   	   appended to the directory name instead of a zero. If more than a hundred directories
   	   exist, an exception is thrown.
   	   @param[in] output_dir_name name of the created output directory (without appended zero)
   	   @param[in] work_path path from which the program is run*/
	boost::filesystem::path create_output_directory(std::string output_dir_name,
													   boost::filesystem::path work_path){
		namespace bfs = boost::filesystem;

		//output directory path if this is the first calculation with the output_dir_name
		unsigned int iteration=0;
		std::string no_output_dir_name = output_dir_name;
		bfs::path out_path= work_path / bfs::path(no_output_dir_name);

		//determine how many simulation result directories with the given output_dir_name are in the
		//working directory and set out_path accordingly
		while (exists(out_path)){
			iteration++;
			no_output_dir_name = output_dir_name+std::to_string(iteration);
			out_path = work_path / bfs::path(no_output_dir_name);
			if (iteration>100){
				std::string error_msg = "More than a hundred result directories of the same name!\n"
										"Clean up your working directory! -> Aborting";
				throw goscalc_io_exception(error_msg);
			}
		}

		//create output directory
		if ( !bfs::create_directory(out_path) ){
			throw goscalc_io_exception("Could not create output directory!->Aborting");
		}

		return out_path;
	}

	void save_k_values(const arma::colvec& k){
		std::ofstream output_file("k.dat");
		output_file<<"#k in inverse Angstrom\n";
		for (unsigned int row = 0; row < k.size(); row++){
			output_file<<k[row]<<"\n";
		}
		output_file.close();
	}

	void save_free_energies(const std::vector<double>& e_losses,double bound_energy){
		std::ofstream output_file("free_energies.dat");
		output_file<<"#free_energies in eV\n";
		for (unsigned int row = 0; row < e_losses.size(); row++){
			output_file<<e_losses[row]+bound_energy<<"\t";
		}
		output_file.close();
	}


	void save_gos(const Goscalc_shared_data& gsd){
		std::ofstream output_file("gos.dat");
		output_file<<"#GOS in 1/eV for e_loss =\n#";
		for (auto energy_loss : gsd.results_energy_losses){
			output_file<<energy_loss<<"eV\t\t";
		}
		output_file<<"\n";
		for (unsigned int row = 0; row < gsd.results_gos[0].size(); row++){
			for (const auto& gos : gsd.results_gos){
				output_file<<gos(row,1)<<"\t\t";
			}
			output_file<<"\n";
		}

		output_file.close();

	}

	void write_logfile(const std::string& filename){
		std::ofstream out_file(filename);
		out_file<<logg::out().get_log();
		out_file.close();
	}

}


void save_results(const Goscalc_shared_data& gsd){
	namespace bfs = boost::filesystem;

	bfs::path work_path = bfs::current_path();
	bfs::path out_path = create_output_directory(gsd.cd.output_dir_name, work_path);

	//switch current path to output directory
	//(switched back upon destruction of out_switcher)
	Path_switcher out_switcher(out_path);

	save_k_values(gsd.results_gos[0].column_as_colvec(0));
	save_free_energies(gsd.results_energy_losses,gsd.cd.bound_energy);
	save_gos(gsd);



	bfs::path config_file_wpath( work_path / bfs::path(gsd.cd.config_filename) );
	bfs::path config_file_opath( out_path / bfs::path(gsd.cd.config_filename) );
	bfs::copy_file(config_file_wpath,config_file_opath.filename(),bfs::copy_option::fail_if_exists);

	write_logfile("goscalc.log");

}
