/*
 * BasicOperator.h
 *
 *  Created on: Mar 26, 2013
 *      Author: knippsch
 */

#ifndef MESONFUNCTIONS_H_
#define MESONFUNCTIONS_H_

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
#include "ReadWrite.h"
#include "propagator_io.h"
#include "quark.h"

// struct for Look-up table in create_gamma and get_operator. To read as
// "in column i the row[i]-element is non-zero and its value is value[i]"
// As Gamma matrices are 4x4 matrices, row and value are 4-vectors

struct gammastruct{
   int index[4];
   std::complex<double> value[4];
};


	void MesonFunc_SinkbarSink(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_diag(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_diag_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_diag_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_col(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_col_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_col_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_row(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_row_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_row_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_diag(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_diag_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_diag_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_col(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_col_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_col_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_row(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_row_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_row_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
        void MesonFunc_SourcebarSource(int t, Eigen::MatrixXcd&, Eigen::VectorXcd*, Eigen::VectorXcd*, Eigen::MatrixXcd*, const int);
        void MesonFunc_SourceSourcebar(int t, Eigen::MatrixXcd&, Eigen::VectorXcd*, Eigen::VectorXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSource(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::VectorXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SourceSink(int t, Eigen::MatrixXcd& , Eigen::VectorXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSourcebar(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::VectorXcd*,  Eigen::MatrixXcd*, const int);
	void MesonFunc_SourcebarSinkbar(int t, Eigen::MatrixXcd& , Eigen::VectorXcd*, Eigen::MatrixXcd*,  Eigen::MatrixXcd*, const int);



	void MesonFunc_SinkbarSink(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  const int);
	void MesonFunc_SinkbarSink_diag(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_diag_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_diag_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_col(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_col_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_col_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_row(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_row_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSink_row_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*,  const int);
	void MesonFunc_SinkSinkbar_diag(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_diag_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_diag_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_col(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_col_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_col_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_row(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_row_shiftup(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkSinkbar_row_shiftdown(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::MatrixXcd*, const int);
        void MesonFunc_SourcebarSource(int t, Eigen::MatrixXcd&, Eigen::VectorXcd*, Eigen::VectorXcd*, const int);
        void MesonFunc_SourceSourcebar(int t, Eigen::MatrixXcd&, Eigen::VectorXcd*, Eigen::VectorXcd*, const int);
	void MesonFunc_SinkSource(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::VectorXcd*, const int);
	void MesonFunc_SourceSink(int t, Eigen::MatrixXcd& , Eigen::VectorXcd*, Eigen::MatrixXcd*, const int);
	void MesonFunc_SinkbarSourcebar(int t, Eigen::MatrixXcd& , Eigen::MatrixXcd*, Eigen::VectorXcd*, const int);
	void MesonFunc_SourcebarSinkbar(int t, Eigen::MatrixXcd& , Eigen::VectorXcd*, Eigen::MatrixXcd*, const int);

        void build_VdaggerV(Eigen::MatrixXcd&, Eigen::MatrixXcd*, const int, const int, const int);

#endif /* MESONFUNCTIONS_H_ */
