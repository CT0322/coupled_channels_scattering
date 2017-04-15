/*
 * Contrac_SoSoSiSi_SiSiSoSo.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef CONTRAC_SISISOSO_SOSI_H_
#define CONTRAC_SISISOSO_SOSI_H_

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


	Eigen::VectorXcd  Contractions_ThreePoint_SiSiSoSo_SoSi(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SiSiSoSo_SoSi_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SiSiSoSo_SoSi_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);

	Eigen::VectorXcd  Contractions_ThreePoint_SiSiSoSo_SoSi_dis(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SiSiSoSo_SoSi_dis_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SiSiSoSo_SoSi_dis_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
#endif /* CONTRAC_SISISOSO_SOSI_H_ */
