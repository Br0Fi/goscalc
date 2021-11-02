/* perform the Hankel transformation
MIT License
Copyright (c) 2009-2021 Peter Koval, J.D. Talman, Leonhard Segger

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef hankel_trafo_HPP
#define hankel_trafo_HPP

#include<exception>
#include<complex>
#include<fftw3.h>

#include"value_pairs.hpp"

/**
 *@file hankel_trafo.hpp
 *@brief Class that performs Hankel transform (header)
 */

 ///Exceptions in hankel_trafo
 class hankel_trafo_exception : public std::exception{
 public:
 	///Construct with error message \p error_msg
 	explicit hankel_trafo_exception(const std::string& error_msg)
 			 : error_message(error_msg) {}
 	///Return error message
 	virtual const char* what() const throw(){
 		return error_message.c_str();
 	}

 private:
 	const std::string error_message;

 };



///Perform Hankel transforms
/*!
 * Perform Hankel transforms from a given dataset. Intermediate results of each
 * Hankel transform are saved to speed up subsequent Hankel transforms.
 * The algorithm works by writing the Hankel transform as a sum over Fourier
 * transforms. For more details see \cite koval2010 and \cite talman 2009
 * (there instead of "Hankel transform", the performed operation is called
 * "spherical Bessel transform").
 */
class Hankel_trafo{
 public:


	///Default copy constructor
	Hankel_trafo(const Hankel_trafo& rhs) = default;
  //logarithmic mesh rr, l_max maximum l value to be transformed with, k_max maximum k-value for output mesh
	Hankel_trafo(const std::vector<double>& rr, const unsigned int l_max, const double k_max);

  //destruktor
  ~Hankel_trafo();

	///Perform \p l -th Hankel transform with function ff
	Value_pairs perform_hankel_trafo(const std::vector<double>& ff, unsigned int l);

 private:

  //the maximum l (or lambda) value that can be computed with this object
  const unsigned int l_max;

  //std::vector<double> ff_fixed; //f(r)
  std::vector<double> rr_fixed;    // r mesh (logarithmic)
  const unsigned int samples;   //number of mesh points
  const unsigned int samples2;  //2*samples
  const double rhomin_fixed, rhomax_fixed; //rho = ln(r)

  std::vector<double> premult;
  std::vector<double> smallr;
  std::vector<double> postdiv;

  std::vector<double> kk; //momentum space mesh (logarithmic)
  std::vector<double> sbt_rr3; //extended r mesh
  std::vector<double> sbt_kk3; //extended k mesh
  //2-dim complex arrays:
  std::vector<std::vector<std::complex<double>>> mult_table1, mult_table2;
  double sbt_rhomin, sbt_kappamin, sbt_rmin, sbt_kmin;

  //following should be threadprivate for multithreading:
  fftw_plan plan12, plan_r2c, plan_c2r;
  std::complex<double> *temp1, *temp2;
  std::complex<double> *r2c_out, *c2r_in;
  double *r2c_in, *c2r_out;
  // !$OMP THREADPRIVATE(plan12,plan_r2c,plan_c2r,temp1,temp2,r2c_in,r2c_out,c2r_in,c2r_out)


};



#endif
