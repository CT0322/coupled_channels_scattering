/*
 * Contrac_twopoint.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef CONTRAC_TWOPOINT_H_
#define CONTRAC_TWOPOINT_H_

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




/******* two point contractions  ************/
	Eigen::VectorXcd  Contractions_TwoPoint_SiSiSoSo(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_TwoPoint_SoSoSiSi(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_TwoPoint_SoSiSiSo(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_TwoPoint_SiSoSoSi(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	std::complex<double>  Contractions_TwoPoint_SiSiSoSo(int, int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	std::complex<double>  Contractions_TwoPoint_SoSoSiSi(int, int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	std::complex<double>  Contractions_TwoPoint_SoSiSiSo(int, int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);
	std::complex<double>  Contractions_TwoPoint_SiSoSoSi(int, int, Eigen::MatrixXcd* , Eigen::MatrixXcd*);

#endif /* CONTRAC_TWOPOINT_H_ */
