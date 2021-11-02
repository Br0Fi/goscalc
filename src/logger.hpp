/* define a log stream for GOScalc to be communicated with the user via the console (and saved)
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

#ifndef logger_HPP
#define logger_HPP

#include<iostream>
#include<sstream>


/**
 *@file logger.hpp
 *@brief Class for simultaneous logging and console output
 */


///Includes logger and functions for logging
namespace logg{
	///Simultaneous logging and user output
	class logger{
	 public:
		///Write to terminal and save internally
		template<typename T> logger& operator<<(T&& rhs){
			std::cout << rhs;
			log_stream << rhs;
			return *this;
		}

		///Save \p input internally only
		template<typename T> logger& log_only(T&& input){
			log_stream << input;
			return *this;
		}

		///Get internal storage as string
		std::string get_log() const{
			return log_stream.str();
		}

		///Provides access to logger
		friend logger& out();

	 private:
		logger() = default;
		std::stringstream log_stream;

	};

	///Function to access logger object (singleton)
	inline logger& out(){
		static logger l;
		return l;
	}

}


#endif
