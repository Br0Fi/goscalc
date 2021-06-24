#include"value_pairs.hpp"

#include<fstream>
#include<sstream>


/**
 *@file value_pairs.cc
 *@brief Class to save and operate on two columns of values, i.e. x and f(x) (source)
 */

namespace{
	std::pair<arma::colvec,arma::colvec> read_xy_data (std::string filename, unsigned int col_no){
		std::pair<std::vector<double>,std::vector<double>> xy_data;

		std::ifstream file(filename.c_str());
		if (!file.is_open()){
			throw value_pairs_io_exception("Wavefunction input file not found.");
		}

		//determine amount of y-columns in file
		std::string dummy_string;
		std::getline(file, dummy_string);   //skip atomic number
		std::getline(file, dummy_string);	//skip energy values
		std::getline(file, dummy_string);   //skip u_nl-stuff
		std::getline(file, dummy_string);
		std::stringstream dummy_stream;
		dummy_stream.str(dummy_string);
		unsigned int columns = 0;
		double dummy_double = 0;
		while(dummy_stream){
			dummy_stream>>dummy_double;
			columns++;
		}
		columns -= 1; //substract one to get actual amount of columns

		if ( ! (col_no<columns) ){
			std::cout<<"col_no ="<<col_no<<"\n";
			std::cout<<"columns = "<<columns<<"\n";
			throw value_pairs_io_exception("Desired y-column not in input file.");
		}

		//reset filestream to beginning of file
		file.clear();
		file.seekg(0, std::ios::beg);
		std::getline(file, dummy_string);   //skip atomic number
		std::getline(file, dummy_string);	//skip energy values
		std::getline(file, dummy_string);   //skip u_nl-stuff


		//read xy_data
		double x_value=0, y_value=0;
		while(true){
			file>>x_value;

			if (col_no == 0){y_value = x_value;}
			for (unsigned int i=1; i<columns; i++){
				if (i == col_no){
					file>>y_value;
				} else {
					file>>dummy_double;
				}
			}

			if (file.eof()){break;}
			xy_data.first.push_back(x_value);
			xy_data.second.push_back(y_value);

		}

		file.close();

		return xy_data;
	}

}

Value_pairs::Value_pairs(size_t entries){
	x_values.resize(entries);
	x_values.fill(0);
	y_values.resize(entries);
	y_values.fill(0);
}



Value_pairs::Value_pairs(std::string filename, unsigned int col_no){
	std::pair<arma::colvec,arma::colvec> pairs = read_xy_data(filename, col_no);
	x_values = pairs.first;
	y_values = pairs.second;
}


double& Value_pairs::operator() (unsigned int index, unsigned int x_or_y){
	if (index>=size()){
		throw value_pairs_exception{"Invalid coordinate for Value_pairs class."};
	}
	if ( x_or_y == 0){
		return x_values[index];
	} else if (x_or_y == 1) {
		return y_values[index];
	} else {
		throw value_pairs_exception("Invalid coordinate for Value_pairs class.");
	}
}


const double& Value_pairs::operator() (unsigned int index, unsigned int x_or_y) const{
	return const_cast<Value_pairs&>(*this)(index,x_or_y);
}

Value_pairs& Value_pairs::operator+= (const Value_pairs& rhs){
	if (!arma::approx_equal(x_values,rhs.x_values,"absdiff",0.002)){
		throw value_pairs_exception("Addition of Value_pairs with different x-values.");
	}
	y_values += rhs.y_values;
	return *this;
}

Value_pairs& Value_pairs::operator-= (const Value_pairs& rhs){
	if (!arma::approx_equal(x_values,rhs.x_values,"absdiff",0.002)){
		throw value_pairs_exception("Substraction of Value_pairs with different x-values.");
	}
	y_values -= rhs.y_values;
	return *this;
}

Value_pairs& Value_pairs::operator*= (const Value_pairs& rhs){
	if (!arma::approx_equal(x_values,rhs.x_values,"absdiff",0.002)){
		throw value_pairs_exception("Multiplication of Value_pairs with different x-values.");
	}
	y_values %= rhs.y_values;
	return *this;
}

Value_pairs& Value_pairs::operator/= (const Value_pairs& rhs){
	if (!arma::approx_equal(x_values,rhs.x_values,"absdiff",0.002)){
		throw value_pairs_exception("Division of Value_pairs with different x-values.");
	}
	y_values /= rhs.y_values;
	return *this;
}

Value_pairs& Value_pairs::operator+= (const arma::colvec& rhs){
	y_values += rhs;
	return *this;
}

Value_pairs& Value_pairs::operator-= (const arma::colvec& rhs){
	y_values -= rhs;
	return *this;
}

Value_pairs& Value_pairs::operator*= (const arma::colvec& rhs){
	y_values %= rhs;
	return *this;
}

Value_pairs& Value_pairs::operator/= (const arma::colvec& rhs){
	y_values /= rhs;
	return *this;
}

Value_pairs& Value_pairs::operator+= (const double& rhs){
	y_values += rhs;
	return *this;
}

Value_pairs& Value_pairs::operator-= (const double& rhs){
	y_values -= rhs;
	return *this;
}

Value_pairs& Value_pairs::operator*= (const double& rhs){
	y_values *= rhs;
	return *this;
}

Value_pairs& Value_pairs::operator/= (const double& rhs){
	y_values /= rhs;
	return *this;
}




std::size_t Value_pairs::size() const {
	return x_values.size();
}


std::vector<double> Value_pairs::column_as_vector(unsigned int x_or_y) const{
	std::vector<double> result;
	result.resize(size());
	if (x_or_y == 0){
		for (unsigned int i=0; i<size(); i++){
			result[i]=x_values[i];
		}
	} else if (x_or_y == 1){
		for (unsigned int i=0; i<size(); i++){
			result[i] = y_values[i];
		}
	} else {
		throw value_pairs_exception("Invalid coordinate for Value_pairs class.");
	}

	return result;
}

const arma::colvec& Value_pairs::column_as_colvec(unsigned int x_or_y) const{
	if (x_or_y == 0){
		return x_values;
	} else if (x_or_y == 1){
		return y_values;
	} else {
		throw value_pairs_exception("Invalid coordinate for Value_pairs class.");
	}
}



void Value_pairs::print(const std::string &filename) const{
	std::ofstream file(filename.c_str());
	for (unsigned int i=0; i<size(); i++){
		file<<x_values[i]<<" "<<y_values[i]<<"\n";
	}
	file.close();
}
