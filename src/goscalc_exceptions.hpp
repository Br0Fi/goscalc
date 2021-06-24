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
