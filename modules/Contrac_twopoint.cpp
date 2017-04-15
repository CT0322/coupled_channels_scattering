/*
 * Contractions.cpp
 *
 *  Created on: Jun 13, 2014
 *      Author: Liu
 *   
 */

#include "Contrac_twopoint.h"


// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();


/******************************************************/
/******************************************************/
/***** CONTRACTION TWO POINT FUNCTIONS CONNECTED ******/
/******************************************************/
/******************************************************/

Eigen::VectorXcd Contractions_TwoPoint_SiSiSoSo(int t_sink, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource){
    const int Lt = global_data->get_Lt();
    Eigen::VectorXcd corr = Eigen::VectorXcd::Zero(Lt);
    int t_source_dil =0;
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_dilution_E = quarks[0].number_of_dilution_E;
    const int number_of_dilution_T = quarks[0].number_of_dilution_T;
    const std::string  dil_T = quarks[0].dilution_T;
   for(int t_source=0; t_source<Lt; ++t_source){
   if(dil_T == "TB")
     t_source_dil = t_source*number_of_dilution_T/Lt;
   else{
      if(dil_T == "TI")
	t_source_dil = t_source%number_of_dilution_T;
      else{
	std::cout<<"unknown time dilution scheme."<<std::endl;
	exit(0);
      }
   }
      corr((t_sink - t_source + Lt)%Lt) = ((*MesonFuncSinkSink).block(0, t_source_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)*(MesonFuncSourceSource[t_source]).adjoint()).trace();

  }
   return corr;
}


/********* Two point function contraction***************/
std::complex<double> Contractions_TwoPoint_SiSiSoSo(int t_sink, int t_source, Eigen::MatrixXcd* MesonFuncSinkSink, Eigen::MatrixXcd* MesonFuncSourceSource){
    const int Lt = global_data->get_Lt();
    std::complex<double> corr(0.0, 0.0);
    int t_source_dil =0;
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_dilution_E = quarks[0].number_of_dilution_E;
    const int number_of_dilution_T = quarks[0].number_of_dilution_T;
    const std::string  dil_T = quarks[0].dilution_T;
    if(dil_T == "TB")
      t_source_dil = t_source*number_of_dilution_T/Lt;
    else if(dil_T == "TI")
      t_source_dil = t_source%number_of_dilution_T;
    else{
         std::cout<<"unknown time dilution scheme."<<std::endl;
         exit(0);
       }

      corr = ((*MesonFuncSinkSink).block(0, t_source_dil*4*number_of_dilution_E, 4*number_of_dilution_E, 4*number_of_dilution_E)*(MesonFuncSourceSource[t_source]).adjoint()).trace();

   return corr;
}




