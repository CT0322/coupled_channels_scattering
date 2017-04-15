/*
 * Contrac_SiSiSoSo_SoSoSiSi.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef CONTRAC_SISISISI_SOSOSOSO_H_
#define CONTRAC_SISISISI_SOSOSOSO_H_

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


	Eigen::VectorXcd  Contractions_FourPoint_SiSiSiSi_SoSoSoSo(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
#endif /* CONTRAC_SISISISI_SOSOSOSO_H_ */
