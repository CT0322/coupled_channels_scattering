/*
 * Contrac_SiSiSoSo_SoSoSiSi.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef CONTRAC_SISISOSO_SOSOSISI_H_
#define CONTRAC_SISISOSO_SOSOSISI_H_

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


	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_3(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_4(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);

	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_1(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_2(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_3(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_3(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_4(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
	Eigen::VectorXcd  Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_4(int, Eigen::MatrixXcd* , Eigen::MatrixXcd*, Eigen::MatrixXcd*, Eigen::MatrixXcd*);
#endif /* CONTRAC_SISISOSO_SOSOSISI_H_ */
