/*
 * Contractions.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: Liu
 *   
 */

#include "Contrac_SoSoSiSi_SiSiSoSo.h"


// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

//no time shift


Eigen::VectorXcd Contractions_FourPoint_SoSoSiSi_SiSiSoSo_dis(int t_source, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	int t_source_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
//	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*quarks[0].number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0;
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

   for(int t_sink=0; t_sink<Lt; ++t_sink){
    if(dil_T == "TB"){
      t_sink_dil = t_sink*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
    }
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSourceSource0[t_sink]*MesonFuncSinkSink1[t_sink].block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source].adjoint()*(MesonFuncSinkSink0[t_source].block(t_sink_dil*dim_block, 0, dim_block, dim_block)).adjoint()).trace();
  }
	return corr;
}

//time shift  t t+1 t0 t0+1 

Eigen::VectorXcd Contractions_FourPoint_SoSoSiSi_SiSiSoSo_dis_1(int t_source, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
//	int t_source_dil;
	int t_source1_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_source1 = (t_source+1)%Lt;
	int t_sink1, t_sink_dil=0;
//	int t_sink1_dil=0;
    if(dil_T == "TB"){
//      t_source_dil = t_source*number_of_dilution_T/Lt;
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
//         t_source_dil = t_source%number_of_dilution_T;
         t_source1_dil = t_source1%number_of_dilution_T;
    }
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
    }

   for(int t_sink=0; t_sink<Lt; ++t_sink){
	 t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
      t_sink_dil = t_sink*number_of_dilution_T/Lt;
//      t_sink1_dil = t_sink1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
//	 t_sink1_dil = t_sink1%number_of_dilution_T;
    }
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSourceSource0[t_sink]*MesonFuncSinkSink1[t_sink1].block(0, t_source1_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source1].adjoint()*(MesonFuncSinkSink0[t_source].block(t_sink_dil*dim_block, 0, dim_block, dim_block)).adjoint()).trace();
  }
	return corr;
}

//time shift t t+1 t0+1 t0

Eigen::VectorXcd Contractions_FourPoint_SoSoSiSi_SiSiSoSo_dis_2(int t_source, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	int t_source_dil;
//	int t_source1_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_source1 = (t_source+1)%Lt;
	int t_sink1, t_sink_dil=0;
//	int t_sink1_dil=0;
    if(dil_T == "TB"){
      t_source_dil = t_source*number_of_dilution_T/Lt;
//      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source_dil = t_source%number_of_dilution_T;
//         t_source1_dil = t_source1%number_of_dilution_T;
    }
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
    }

   for(int t_sink=0; t_sink<Lt; ++t_sink){
	 t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
      t_sink_dil = t_sink*number_of_dilution_T/Lt;
//      t_sink1_dil = t_sink1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
//	 t_sink1_dil = t_sink1%number_of_dilution_T;
    }
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSourceSource0[t_sink]*MesonFuncSinkSink1[t_sink1].block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source].adjoint()*(MesonFuncSinkSink0[t_source1].block(t_sink_dil*dim_block, 0, dim_block, dim_block)).adjoint()).trace();
  }
	return corr;
}

//time shift t+1 t t0+1 t0

Eigen::VectorXcd Contractions_FourPoint_SoSoSiSi_SiSiSoSo_dis_3(int t_source, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	int t_source_dil;
//	int t_source1_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_source1 = (t_source+1)%Lt;
	int t_sink1;
//	int t_sink_dil=0;
	int t_sink1_dil=0;
    if(dil_T == "TB"){
      t_source_dil = t_source*number_of_dilution_T/Lt;
//      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source_dil = t_source%number_of_dilution_T;
//         t_source1_dil = t_source1%number_of_dilution_T;
    }
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
    }

   for(int t_sink=0; t_sink<Lt; ++t_sink){
	 t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
 //     t_sink_dil = t_sink*number_of_dilution_T/Lt;
      t_sink1_dil = t_sink1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
//	 t_sink_dil = t_sink%number_of_dilution_T;
	 t_sink1_dil = t_sink1%number_of_dilution_T;
    }
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSourceSource0[t_sink1]*MesonFuncSinkSink1[t_sink].block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source].adjoint()*(MesonFuncSinkSink0[t_source1].block(t_sink1_dil*dim_block, 0, dim_block, dim_block)).adjoint()).trace();
  }
	return corr;
}

//time shift t+1 t t0 t0+1

Eigen::VectorXcd Contractions_FourPoint_SoSoSiSi_SiSiSoSo_dis_4(int t_source, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
//	int t_source_dil;
	int t_source1_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_source1 = (t_source+1)%Lt;
	int t_sink1;
//	int t_sink_dil=0;
	int t_sink1_dil=0;
    if(dil_T == "TB"){
//      t_source_dil = t_source*number_of_dilution_T/Lt;
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
//         t_source_dil = t_source%number_of_dilution_T;
         t_source1_dil = t_source1%number_of_dilution_T;
    }
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
    }

   for(int t_sink=0; t_sink<Lt; ++t_sink){
	 t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
 //     t_sink_dil = t_sink*number_of_dilution_T/Lt;
      t_sink1_dil = t_sink1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
//	 t_sink_dil = t_sink%number_of_dilution_T;
	 t_sink1_dil = t_sink1%number_of_dilution_T;
    }
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSourceSource0[t_sink1]*MesonFuncSinkSink1[t_sink].block(0, t_source1_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_source1].adjoint()*(MesonFuncSinkSink0[t_source].block(t_sink1_dil*dim_block, 0, dim_block, dim_block)).adjoint()).trace();
  }
	return corr;
}
