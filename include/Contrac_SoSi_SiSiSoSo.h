/*
 * Contrac_SoSoSiSi_SiSiSoSo.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef CONTRAC_SOSI_SISISOSO_H_
#define CONTRAC_SOSI_SISISOSO_H_

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


	Eigen::VectorXcd  Contractions_ThreePoint_SoSi_SiSiSoSo(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSi_SiSiSoSo_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSi_SiSiSoSo_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);

	Eigen::VectorXcd  Contractions_ThreePoint_SoSi_SiSiSoSo_dis(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSi_SiSiSoSo_dis_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSi_SiSiSoSo_dis_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
#endif /* CONTRAC_SOSI_SISISOSO_H_ */
