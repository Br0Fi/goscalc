/* Load data from wavegen and write calculated GOS to files
    Copyright (C) 2019-2021 The GOScalc devolopers

    This file is part of GOScalc.

    GOScalc is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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
