/* Data structure to contain a function and its corresponding mesh values (~ the x-axis)
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

#ifndef value_pairs_HPP
#define value_pairs_HPP

#include<vector>
#include<string>
#include<exception>

#include<armadillo>
#include<boost/operators.hpp>

/**
 *@file value_pairs.hpp
 *@brief Class to save and operate on two columns of values, i.e. x and f(x) (header)
 */
 
 
 ///Exceptions in Value_pairs
 class value_pairs_exception : public std::exception{
 public:
 	///Construct with error message \p error_msg
 	explicit value_pairs_exception(const std::string& error_msg)
 			 : error_message(error_msg) {}	
 	///Return error message
 	virtual const char* what() const throw(){
 		return error_message.c_str();	
 	}
 
 private:
 	const std::string error_message;
 
 };
 
 
 ///Input/output-exceptions in Value_pairs
 class value_pairs_io_exception : public value_pairs_exception{
  public:
 	///Construct with error message \p error_msg
 	explicit value_pairs_io_exception(const std::string& error_msg)
 			 : value_pairs_exception(error_msg) {}
 };
 


///Holds two columns of values in a \f$x\f$ - \f$f(x)\f$-fashion
class Value_pairs :
		boost::arithmetic<Value_pairs,
		boost::arithmetic<Value_pairs, double,
		boost::arithmetic<Value_pairs, arma::colvec
		>>>{
 public:
	Value_pairs() = default;														///<Default constructor of empty Value_pairs object with zero entries
	explicit Value_pairs(size_t entries);											///<Construct empty Value_pairs object with \p entries entries

	
	///Read two columns from file \p filename
	/* Read two columns form file \p filename. The file is supposed to contain
	 * columns of values. This constructor reads the first column as the x-values and
	 * the column \p col_no as the y-values. If col_no = 0, the y-values will
	 * be equal to the x-values. If col_no = 1, the second column in the file
	 * will be used as y-values and so on... If the file is not found or \p col_no
	 * is out of range of the amount of columns in the file, an exception is thrown.
	 */
	Value_pairs(std::string filename, unsigned int col_no);
	double& operator() (unsigned int index, unsigned int x_or_y);					///<Get x- (\p x_or_y=0) or y-value (\p x_or_y=1) at row \p index
	const double& operator() (unsigned int index, unsigned int x_or_y) const;		///<Const-cast of operator()
	Value_pairs& operator+= (const Value_pairs& rhs);								///<+= second column of \p rhs to second column of this
	Value_pairs& operator-= (const Value_pairs& rhs);								///<-= second column of \p rhs to second column of this
	Value_pairs& operator*= (const Value_pairs& rhs);								///<*= second column of \p rhs to second column of this
	Value_pairs& operator/= (const Value_pairs& rhs);								///</= second column of \p rhs to second column of this
	Value_pairs& operator+= (const arma::colvec& rhs);								///<+= \p rhs -column to second column of this
	Value_pairs& operator-= (const arma::colvec& rhs);								///<-= \p rhs -column to second column of this
	Value_pairs& operator*= (const arma::colvec& rhs);								///<*= \p rhs -column to second column of this
	Value_pairs& operator/= (const arma::colvec& rhs);								///</= \p rhs -column to second column of this
	Value_pairs& operator+= (const double& rhs);									///<+= \p rhs to second column of this
	Value_pairs& operator-= (const double& rhs);									///<-= \p rhs to second column of this
	Value_pairs& operator*= (const double& rhs);									///<*= \p rhs to second column of this
	Value_pairs& operator/= (const double& rhs);									///</= \p rhs to second column of this
	
	std::size_t size() const;														///<Get size of x-column and y-column
	std::vector<double> column_as_vector(unsigned int x_or_y) const;				///<Return x- (\p x_or_y=0) or y-column (\p x_or_y=1) as a std::vector
	const arma::colvec& column_as_colvec(unsigned int x_or_y) const;						///<Return x- (\p x_or_y=0) or y-column (\p x_or_y=1) as an arma::colvec
	void print(const std::string& filename) const;									///<Print to file

 private: 
	arma::colvec x_values;
	arma::colvec y_values;
};




#endif
