#ifndef multex_io_HPP
#define multex_io_HPP


#include"goscalc_shared_data.hpp"


/**
 *@file transpot_io.hpp
 *@brief Functions for input and output (header)
 */


///Input/output-exceptions in goscalc
class goscalc_io_exception : public std::ios_base::failure{
 public:
	///Construct with error message \p error_msg
	explicit goscalc_io_exception(const std::string& error_msg)
			 : std::ios_base::failure(error_msg) {}
};




///Read config file \p config_filename
Config_data read_input(std::string config_filename);

///Save results of program run (which are stored in \p tsd)
void save_results(const Goscalc_shared_data& gsd);



#endif
