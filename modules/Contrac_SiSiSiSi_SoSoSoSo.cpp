/*
 * Contractions.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: Liu
 *   
 */

#include "Contrac_SiSiSiSi_SoSoSoSo.h"


// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

//no time shift

Eigen::VectorXcd Contractions_FourPoint_SiSiSiSi_SoSoSoSo(int t_sink, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*quarks[0].number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_source_dil;

    for(int t_source=0; t_source<Lt; t_source++){

    if(dil_T == "TB"){
      t_source_dil = t_source*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source_dil = t_source%number_of_dilution_T;
    }
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
    }

		corr0 = ((*MesonFuncSinkSink0).block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source].adjoint()).trace();
		corr1 = ((*MesonFuncSinkSink1).block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source].adjoint()).trace();
		corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
    }
	return corr;
}


Eigen::VectorXcd Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(int t_sink, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*quarks[0].number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_source_dil;

    for(int t_source=0; t_source<Lt; t_source++){
    if(dil_T == "TB"){
      t_source_dil = t_source*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source_dil = t_source%number_of_dilution_T;
    }
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
    }

      corr[(t_sink - t_source +Lt)%Lt] += ((*MesonFuncSinkSink0).block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source].adjoint()*(*MesonFuncSinkSink1).block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source].adjoint()).trace();
   }
	return corr;
}

