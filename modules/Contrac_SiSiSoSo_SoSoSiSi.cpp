/*
 * Contractions.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: Liu
 *   
 */

#include "Contrac_SiSiSoSo_SoSoSiSi.h"


// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();

//no time shift

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
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
      corr0 = (MesonFuncSinkSink0[t_sink].block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source].adjoint()).trace();
      corr1 = (MesonFuncSourceSource1[t_sink]*(MesonFuncSinkSink1[t_source].block(0, t_sink_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	int t_source_dil;
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
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
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSinkSink0[t_sink].block(t_source_dil*dim_block, 0, dim_block, dim_block)*MesonFuncSourceSource1[t_sink]*(MesonFuncSinkSink1[t_source].block(0, t_sink_dil*dim_block, dim_block, dim_block)).adjoint()*MesonFuncSourceSource0[t_source].adjoint()).trace();
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
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
      corr0 = (MesonFuncSinkSink0[t_sink].block(0, t_sink_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_sink]).trace();
      corr1 = (MesonFuncSourceSource0[t_source].adjoint()*(MesonFuncSinkSink1[t_source].block(0, t_source_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

//time shift  t t+1 t0 t0+1 
Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_1(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source_dil;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

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
      corr0 = (MesonFuncSinkSink0[t_sink].block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source].adjoint()).trace();
      corr1 = (MesonFuncSourceSource1[t_sink1]*(MesonFuncSinkSink1[t_source1].block(0, t_sink1_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_1(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source_dil=0;
	int t_sink1, t_source1;
	t_source1 = (t_source+1)%Lt;
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
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSinkSink0[t_sink].block(t_source_dil*dim_block, 0, dim_block, dim_block)*MesonFuncSourceSource1[t_sink1]*(MesonFuncSinkSink1[t_source1].block(0, t_sink1_dil*dim_block, dim_block, dim_block)).adjoint()*MesonFuncSourceSource0[t_source].adjoint()).trace();
  }
	return corr;
}


Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_1(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source_dil;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

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
      corr0 = (MesonFuncSinkSink0[t_sink].block(0, t_sink1_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_sink1]).trace();
      corr1 = (MesonFuncSourceSource0[t_source].adjoint()*(MesonFuncSinkSink1[t_source1].block(0, t_source_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}


//time shift t t+1 t0+1 t0
Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_2(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source1_dil=0;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

    if(dil_T == "TB"){
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source1_dil = t_source1%number_of_dilution_T;
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
      corr0 = (MesonFuncSinkSink0[t_sink].block(0, t_source1_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source1].adjoint()).trace();
      corr1 = (MesonFuncSourceSource1[t_sink1]*(MesonFuncSinkSink1[t_source].block(0, t_sink1_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_2(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source1_dil=0;
	int t_sink1, t_source1;
	t_source1 = (t_source+1)%Lt;
    if(dil_T == "TB"){
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source1_dil = t_source1%number_of_dilution_T;
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
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSinkSink0[t_sink].block(t_source1_dil*dim_block, 0, dim_block, dim_block)*MesonFuncSourceSource1[t_sink1]*(MesonFuncSinkSink1[t_source].block(0, t_sink1_dil*dim_block, dim_block, dim_block)).adjoint()*MesonFuncSourceSource0[t_source1].adjoint()).trace();
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_2(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink1_dil=0, t_source1_dil=0;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

    if(dil_T == "TB"){
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
         t_source1_dil = t_source1%number_of_dilution_T;
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
      corr0 = (MesonFuncSinkSink0[t_sink].block(0, t_sink1_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_sink1]).trace();
      corr1 = (MesonFuncSourceSource0[t_source1].adjoint()*(MesonFuncSinkSink1[t_source].block(0, t_source1_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}


//time shift t+1 t t0+1 t0
Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_3(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source1_dil;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

    if(dil_T == "TB"){
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
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
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
    }
      corr0 = (MesonFuncSinkSink0[t_sink1].block(0, t_source1_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source1].adjoint()).trace();
      corr1 = (MesonFuncSourceSource1[t_sink]*(MesonFuncSinkSink1[t_source].block(0, t_sink_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_3(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source1_dil=0;
	int t_sink1, t_source1;
	t_source1 = (t_source+1)%Lt;
    if(dil_T == "TB"){
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
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
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
    }
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSinkSink0[t_sink1].block(t_source1_dil*dim_block, 0, dim_block, dim_block)*MesonFuncSourceSource1[t_sink]*(MesonFuncSinkSink1[t_source].block(0, t_sink_dil*dim_block, dim_block, dim_block)).adjoint()*MesonFuncSourceSource0[t_source1].adjoint()).trace();
  }
	return corr;
}


Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_3(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source1_dil;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

    if(dil_T == "TB"){
      t_source1_dil = t_source1*number_of_dilution_T/Lt;
    }
    else if(dil_T == "TI"){
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
    }
    else if(dil_T == "TI"){
	 t_sink_dil = t_sink%number_of_dilution_T;
    }
      corr0 = (MesonFuncSinkSink0[t_sink1].block(0, t_sink_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_sink]).trace();
      corr1 = (MesonFuncSourceSource0[t_source1].adjoint()*(MesonFuncSinkSink1[t_source].block(0, t_source1_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}


//time shift t+1 t t0 t0+1
Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_4(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source_dil;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

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
      corr0 = (MesonFuncSinkSink0[t_sink1].block(0, t_source_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource0[t_source].adjoint()).trace();
      corr1 = (MesonFuncSourceSource1[t_sink]*(MesonFuncSinkSink1[t_source1].block(0, t_sink_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_dis_4(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source_dil=0;
	int t_sink1, t_source1;
	t_source1 = (t_source+1)%Lt;
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
      corr((t_sink - t_source + Lt)%Lt) = (MesonFuncSinkSink0[t_sink1].block(t_source_dil*dim_block, 0, dim_block, dim_block)*MesonFuncSourceSource1[t_sink]*(MesonFuncSinkSink1[t_source1].block(0, t_sink_dil*dim_block, dim_block, dim_block)).adjoint()*MesonFuncSourceSource0[t_source].adjoint()).trace();
  }
	return corr;
}

Eigen::VectorXcd Contractions_FourPoint_SiSiSoSo_SoSoSiSi_disdis_4(int t_source, Eigen::MatrixXcd* MesonFuncSinkSink0, Eigen::MatrixXcd* MesonFuncSourceSource1, Eigen::MatrixXcd* MesonFuncSourceSource0, Eigen::MatrixXcd* MesonFuncSinkSink1){
	const int Lt = global_data->get_Lt();
	Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
	std::complex<double> corr0, corr1;
	corr0 = std::complex<double>(0.0, 0.0);
	corr1 = std::complex<double>(0.0, 0.0);
	const std::vector<quark> quarks = global_data->get_quarks();
	const int number_of_dilution_E = quarks[0].number_of_dilution_E;
	const int number_of_dilution_T = quarks[0].number_of_dilution_T;
	const int dim_block = quarks[0].number_of_dilution_D*number_of_dilution_E;
	const std::string dil_T = quarks[0].dilution_T;
	int t_sink_dil=0, t_source_dil;
	int t_sink1, t_source1;
	t_source1 = (t_source + 1)%Lt;

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
      corr0 = (MesonFuncSinkSink0[t_sink1].block(0, t_sink_dil*dim_block, dim_block, dim_block)*MesonFuncSourceSource1[t_sink]).trace();
      corr1 = (MesonFuncSourceSource0[t_source].adjoint()*(MesonFuncSinkSink1[t_source1].block(0, t_source_dil*dim_block, dim_block, dim_block)).adjoint()).trace();
      corr[(t_sink - t_source +Lt)%Lt] += corr0*corr1;
  }
	return corr;
}

