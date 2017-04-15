/*
 * Contrac_SoSoSiSi_SiSiSoSo.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef CONTRAC_SOSOSISI_SISO_H_
#define CONTRAC_SOSOSISI_SISO_H_

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


	Eigen::VectorXcd  Contractions_ThreePoint_SoSoSiSi_SiSo(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSoSiSi_SiSo_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSoSiSi_SiSo_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);

	Eigen::VectorXcd  Contractions_ThreePoint_SoSoSiSi_SiSo_dis(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSoSiSi_SiSo_dis_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_ThreePoint_SoSoSiSi_SiSo_dis_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*);
#endif /* CONTRAC_SOSOSISI_SISO_H_ */
