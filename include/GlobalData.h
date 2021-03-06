/*
 * GlobalData.h
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "quark.h"

class GlobalData {

private:
	//! A pointer on the class itself
	static GlobalData* instance_;
	//! globally accessible data
	int Lx, Ly, Lz, Lt;
	int number_of_eigen_vec;
	int dim_row;
	int start_config, end_config, delta_config;
	int verbose;
	std::string path_eigenvectors;
	std::string path_perambulators_c;
	std::string path_perambulators_s;
	std::string path_perambulators_u;
	std::string path_output;
	std::vector<quark> quarks;
	void quark_input_data_handling (const std::vector<std::string> quark_configs);


public:
	static GlobalData* Instance ();

	void read_parameters(int ac, char* av[]);

	inline int get_Lx () {
		return Lx;
	}
	inline int get_Ly () {
		return Ly;
	}
	inline int get_Lz () {
		return Lz;
	}
	inline int get_Lt () {
		return Lt;
	}
	inline int get_dim_row () {
		return dim_row;
	}
	inline int get_start_config () {
		return start_config;
	}
	inline int get_end_config () {
		return end_config;
	}
	inline int get_delta_config () {
		return delta_config;
	}
	inline int get_number_of_eigen_vec() {
		return number_of_eigen_vec;
	}
	inline int get_verbose() {
		return verbose;
	}
	inline std::string get_path_eigenvectors() {
		return path_eigenvectors;
	}
	inline std::string get_path_perambulators_c() {
		return path_perambulators_c;
	}
	inline std::string get_path_perambulators_s() {
		return path_perambulators_s;
	}
	inline std::string get_path_perambulators_u() {
		return path_perambulators_u;
	}
	inline std::string get_path_output() {
		return path_output;
	}
	inline std::vector<quark> get_quarks() {
		return quarks;
	}

	//! All con/de-structors are protected to assure that only one instance exists
	//! at once. DO NOT CHANGE!!
protected:
	GlobalData () {
	}
	GlobalData (const GlobalData& other) {
	}
	virtual ~GlobalData () {
	}

};

#endif /* GLOBALDATA_H_ */
