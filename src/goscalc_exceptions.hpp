/* Define exceptions for GOScalc
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

#ifndef goscalc_exceptions_HPP
#define goscalc_exceptions_HPP

#include<exception>


/**
 *@file goscalc_exceptions.hpp
 *@brief Exceptions for goscalc
 */

///Exceptions in goscalc
class goscalc_exception : public std::exception{
public:
	///Construct with error message \p error_msg
	explicit goscalc_exception(const std::string& error_msg)
			 : error_message(error_msg) {}	
	///Return error message
	virtual const char* what() const throw(){
		return error_message.c_str();	
	}

private:
	const std::string error_message;

};

///Input/output-exceptions in goscalc
class goscalc_io_exception : public goscalc_exception{
 public:
	///Construct with error message \p error_msg
	explicit goscalc_io_exception(const std::string& error_msg)
			 : goscalc_exception(error_msg) {}
};






#endif
