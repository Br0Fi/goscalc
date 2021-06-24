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
