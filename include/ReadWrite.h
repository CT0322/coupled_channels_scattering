/*
 * ReadWrite.h
 *
 *  Created on: Apr 13, 2014
 *      Author: werner
 */

#ifndef _READ_WRITE_H_
#define _READ_WRITE_H_

// TODO: check if they are all necessary. Doesn't matter much though
#include <iostream>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "GlobalData.h"
#include "propagator_io.h"
#include "quark.h"

/***************************Input from files**********************************/

class ReadWrite {

public:
	ReadWrite ();
	virtual ~ReadWrite (); 
	void read_eigenvectors_from_file (const int config_i, const int);
	void read_perambulators_from_file (const int config_i);
	void read_rnd_vectors_from_file (const int config_i);

	Eigen::MatrixXcd* perambulator_c;
        Eigen::MatrixXcd* perambulator_s;
        Eigen::MatrixXcd* perambulator_u;
	Eigen::VectorXcd* rnd_vec_c;
	Eigen::VectorXcd* rnd_vec_s;
	Eigen::VectorXcd* rnd_vec_u;
	Eigen::MatrixXcd V;
};

#endif // _READ_WRITE_H__
