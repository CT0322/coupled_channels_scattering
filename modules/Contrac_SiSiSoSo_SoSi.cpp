/*
 * Contractions.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: Liu
 *   
 */

#include "Contrac_SiSiSoSo_SoSi.h"


// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

//no time shift

Eigen::VectorXcd Contractions_ThreePoint_SiSiSoSo_SoSi(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource, Eigen::MatrixXcd* MesonFuncSourceSink){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source_dil;
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
       corr[(t_sink - t_source +Lt)%Lt] = (MesonFuncSinkSink[t_sink].block(t_source_dil*4*number_of_dilution_E, 0, 4*number_of_dilution_E, 4*number_of_dilution_E)*MesonFuncSourceSource[t_sink]*(MesonFuncSourceSink[t_source].block(0, t_sink_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)).adjoint()).trace();
  }
	return corr;
}

Eigen::VectorXcd Contractions_ThreePoint_SiSiSoSo_SoSi_dis(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource, Eigen::MatrixXcd* MesonFuncSourceSink){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	int t_source_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
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
	corr0 = (MesonFuncSinkSink[t_sink].block(0, t_sink_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)*MesonFuncSourceSource[t_sink]).trace();
	corr1 = ((MesonFuncSourceSink[t_source].block(0, t_source_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)).adjoint()).trace();
      corr((t_sink - t_source + Lt)%Lt) = corr0*corr1;
  }
	return corr;
}

//time shift  t t+1 t0 

Eigen::VectorXcd Contractions_ThreePoint_SiSiSoSo_SoSi_1(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource, Eigen::MatrixXcd* MesonFuncSourceSink){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source_dil;
	int t_sink1;
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
    t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
      t_sink1_dil = t_sink1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink1_dil = t_sink1%number_of_dilution_T;
    }
       corr[(t_sink - t_source +Lt)%Lt] = (MesonFuncSinkSink[t_sink].block(t_source_dil*4*number_of_dilution_E, 0, 4*number_of_dilution_E, 4*number_of_dilution_E)*MesonFuncSourceSource[t_sink1]*(MesonFuncSourceSink[t_source].block(0, t_sink1_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)).adjoint()).trace();
  }
	return corr;
}

Eigen::VectorXcd Contractions_ThreePoint_SiSiSoSo_SoSi_dis_1(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource, Eigen::MatrixXcd* MesonFuncSourceSink){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	int t_source_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_sink1;
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
	t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
      t_sink1_dil = t_sink1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink1_dil = t_sink1%number_of_dilution_T;
    }
	corr0 = (MesonFuncSinkSink[t_sink].block(0, t_sink1_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)*MesonFuncSourceSource[t_sink1]).trace();
	corr1 = ((MesonFuncSourceSink[t_source].block(0, t_source_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)).adjoint()).trace();
      corr((t_sink - t_source + Lt)%Lt) = corr0*corr1;
  }
	return corr;
}


//time shift t+1 t t0

Eigen::VectorXcd Contractions_ThreePoint_SiSiSoSo_SoSi_2(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource, Eigen::MatrixXcd* MesonFuncSourceSink){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source_dil;
	int t_sink1;
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
    t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
      t_sink_dil = t_sink*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
    }
       corr[(t_sink - t_source +Lt)%Lt] = (MesonFuncSinkSink[t_sink1].block(t_source_dil*4*number_of_dilution_E, 0, 4*number_of_dilution_E, 4*number_of_dilution_E)*MesonFuncSourceSource[t_sink]*(MesonFuncSourceSink[t_source].block(0, t_sink_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)).adjoint()).trace();
  }
	return corr;
}

Eigen::VectorXcd Contractions_ThreePoint_SiSiSoSo_SoSi_dis_2(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource, Eigen::MatrixXcd* MesonFuncSourceSink){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	int t_source_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_sink1;
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
	t_sink1 = (t_sink+1)%Lt;
    if(dil_T == "TB"){
      t_sink_dil = t_sink*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
    }
	corr0 = (MesonFuncSinkSink[t_sink1].block(0, t_sink_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)*MesonFuncSourceSource[t_sink]).trace();
	corr1 = ((MesonFuncSourceSink[t_source].block(0, t_source_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)).adjoint()).trace();
      corr((t_sink - t_source + Lt)%Lt) = corr0*corr1;
  }
	return corr;
}


