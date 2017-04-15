/*
 * MesonFunctions.cpp
 *
 *  Created on: Mar 25, 2014
 *      Author: Liuming Liu
 * 
 */

#include "MesonFunctions.h"

namespace{

static const std::complex<double> I(0.0, 1.0);

static void create_gamma(gammastruct& GammaMatrix, const int i){
 switch(i) {

//Identity
  case 0:      
      GammaMatrix.index[0] = 0;
      GammaMatrix.index[1] = 5;
      GammaMatrix.index[2] = 10;
      GammaMatrix.index[3] = 15;
      GammaMatrix.value[0] = 1.0;
      GammaMatrix.value[1] = 1.0;
      GammaMatrix.value[2] = 1.0;
      GammaMatrix.value[3] = 1.0;
      break;

//gamma1
  case 1:
      GammaMatrix.index[0] = 3;
      GammaMatrix.index[1] = 6;
      GammaMatrix.index[2] = 9;
      GammaMatrix.index[3] = 12;
      GammaMatrix.value[0] = -I;
      GammaMatrix.value[1] = -I;
      GammaMatrix.value[2] = I;
      GammaMatrix.value[3] = I;
      break;
  
//gamma2
  case 2:
      GammaMatrix.index[0] = 3;
      GammaMatrix.index[1] = 6;
      GammaMatrix.index[2] = 9;
      GammaMatrix.index[3] = 12;
      GammaMatrix.value[0] = -1.0;
      GammaMatrix.value[1] = 1.0;
      GammaMatrix.value[2] = 1.0;
      GammaMatrix.value[3] = -1.0;
      break;

//gamma3
  case 3:
      GammaMatrix.index[0] = 2;
      GammaMatrix.index[1] = 7;
      GammaMatrix.index[2] = 8;
      GammaMatrix.index[3] = 13;
      GammaMatrix.value[0] = -I;
      GammaMatrix.value[1] = I;
      GammaMatrix.value[2] = I;
      GammaMatrix.value[3] = -I;
      break;

//gamma4
  case 4:
      GammaMatrix.index[0] = 2;
      GammaMatrix.index[1] = 7;
      GammaMatrix.index[2] = 8;
      GammaMatrix.index[3] = 13;
      GammaMatrix.value[0] = 1.0;
      GammaMatrix.value[1] = 1.0;
      GammaMatrix.value[2] = 1.0;
      GammaMatrix.value[3] = 1.0;
      break;

//gamma5
  case 5:
      GammaMatrix.index[0] = 0;
      GammaMatrix.index[1] = 5;
      GammaMatrix.index[2] = 10;
      GammaMatrix.index[3] = 15;
      GammaMatrix.value[0] = 1.0;
      GammaMatrix.value[1] = 1.0;
      GammaMatrix.value[2] = -1.0;
      GammaMatrix.value[3] = -1.0;
      break;


//gamma1*gamma4*gamma5 (-gamma2*gamma3)
  case 6:
      GammaMatrix.index[0] = 1;
      GammaMatrix.index[1] = 4;
      GammaMatrix.index[2] = 11;
      GammaMatrix.index[3] = 14;
      GammaMatrix.value[0] = -I;
      GammaMatrix.value[1] = -I;
      GammaMatrix.value[2] = -I;
      GammaMatrix.value[3] = -I;
      break;

//gamma2*gamma4*gamma5 (-gamma3*gamma1)
  case 7:
      GammaMatrix.index[0] = 1;
      GammaMatrix.index[1] = 4;
      GammaMatrix.index[2] = 11;
      GammaMatrix.index[3] = 14;
      GammaMatrix.value[0] = -1.0;
      GammaMatrix.value[1] = 1.0;
      GammaMatrix.value[2] = -1.0;
      GammaMatrix.value[3] = 1.0;
      break;

//gamma3*gamma4*gamma5 (-gamma1*gamma2)
  case 8:
      GammaMatrix.index[0] = 0;
      GammaMatrix.index[1] = 5;
      GammaMatrix.index[2] = 10;
      GammaMatrix.index[3] = 15;
      GammaMatrix.value[0] = -I;
      GammaMatrix.value[1] = I;
      GammaMatrix.value[2] = -I;
      GammaMatrix.value[3] = I;
      break;

//gamma1*gamma4
  case 9:
      GammaMatrix.index[0] = 1;
      GammaMatrix.index[1] = 4;
      GammaMatrix.index[2] = 11;
      GammaMatrix.index[3] = 14;
      GammaMatrix.value[0] = -I;
      GammaMatrix.value[1] = -I;
      GammaMatrix.value[2] = I;
      GammaMatrix.value[3] = I;
      break;

//gamma2*gamma4
  case 10:
      GammaMatrix.index[0] = 1;
      GammaMatrix.index[1] = 4;
      GammaMatrix.index[2] = 11;
      GammaMatrix.index[3] = 14;
      GammaMatrix.value[0] = -1.0;
      GammaMatrix.value[1] = 1.0;
      GammaMatrix.value[2] = 1.0;
      GammaMatrix.value[3] = -1.0;
      break;

//gamma3*gamma4
  case 11:
      GammaMatrix.index[0] = 0;
      GammaMatrix.index[1] = 5;
      GammaMatrix.index[2] = 10;
      GammaMatrix.index[3] = 15;
      GammaMatrix.value[0] = -I;
      GammaMatrix.value[1] = I;
      GammaMatrix.value[2] = I;
      GammaMatrix.value[3] = -I;
      break;

//gamma1*gamma5
  case 12:
      GammaMatrix.index[0] = 3;
      GammaMatrix.index[1] = 6;
      GammaMatrix.index[2] = 9;
      GammaMatrix.index[3] = 12;
      GammaMatrix.value[0] = I;
      GammaMatrix.value[1] = I;
      GammaMatrix.value[2] = I;
      GammaMatrix.value[3] = I;
      break;

//gamma2*gamma5
  case 13:
      GammaMatrix.index[0] = 3;
      GammaMatrix.index[1] = 6;
      GammaMatrix.index[2] = 9;
      GammaMatrix.index[3] = 12;
      GammaMatrix.value[0] = 1.0;
      GammaMatrix.value[1] = -1.0;
      GammaMatrix.value[2] = 1.0;
      GammaMatrix.value[3] = -1.0;
      break;

//gamma3*gamma5
  case 14:
      GammaMatrix.index[0] = 2;
      GammaMatrix.index[1] = 7;
      GammaMatrix.index[2] = 8;
      GammaMatrix.index[3] = 13;
      GammaMatrix.value[0] = I;
      GammaMatrix.value[1] = -I;
      GammaMatrix.value[2] = I;
      GammaMatrix.value[3] = -I;
      break;

//gamma4*gamma5
  case 15:
      GammaMatrix.index[0] = 2;
      GammaMatrix.index[1] = 7;
      GammaMatrix.index[2] = 8;
      GammaMatrix.index[3] = 13;
      GammaMatrix.value[0] = -1.0;
      GammaMatrix.value[1] = -1.0;
      GammaMatrix.value[2] = 1.0;
      GammaMatrix.value[3] = 1.0;
      break;

  default:
      printf("gamma index must be 0~16.");
      exit(0);

  }
 }
}


// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();
/******************************************************************************/
void build_VdaggerV(Eigen::MatrixXcd& VdaggerV, Eigen::MatrixXcd* V, const int Px, const int Py, const int Pz){
   
    clock_t time = clock();
//    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::complex<double> I(0.0, 1.0);

   
    Eigen::VectorXcd momentum = Eigen::VectorXcd::Zero(Lx*Ly*Lz);
    Eigen::MatrixXcd temp = Eigen::MatrixXcd::Zero(Lx*Ly*Lz*3, number_of_eigen_vec);
/*** create momentum ****/
  const double pxunit = 2. * M_PI / (double) Lx;
  const double pyunit = 2. * M_PI / (double) Ly;
  const double pzunit = 2. * M_PI / (double) Lz;
    
   for(int x=0; x<Lx; ++x)
    for(int y=0; y<Ly; ++y)
     for(int z=0; z<Lz; ++z)
       momentum(x*Ly*Lz+y*Lz+z) = exp(-I*(Px*pxunit*x+Py*pyunit*y+Pz*pzunit*z));
/**********************/

    for(int i=0; i<Lx*Ly*Lz; ++i){
      temp.block(i*3, 0, 3, number_of_eigen_vec) = momentum(i)*(*V).block(i*3, 0, 3, number_of_eigen_vec);
    }
  VdaggerV = (*V).adjoint()*temp;

  time = clock() - time;
  printf("\t\t build Vdagger*V success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}



/************************************************/
/******* build sinkbar sink Matrix  ****************/
/************************************************/

void MesonFunc_SinkbarSink(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, number_of_inversion);
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);

    MesonFunc += gammasink.value[i]*temp1.adjoint()*temp2;
  }

//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_diag(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_diag mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkbarSink_diag_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, Lt*dim_block);
  
  const std::string dil_T = quarks[0].dilution_T;
  Eigen::VectorXi dilt1 = Eigen::VectorXi::Zero(Lt);
  Eigen::VectorXi dilt = Eigen::VectorXi::Zero(Lt);
  int ti1;
  if(dil_T == "TB"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti*num_dil_T/Lt;
    ti1 = (ti+1)%Lt;
    dilt1(ti) = ti1*num_dil_T/Lt;    
  }  
  }

  else if(dil_T == "TI"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti%num_dil_T;
    ti1 = (ti+1)%Lt;
    dilt1(ti) = ti1%num_dil_T;    
  }  
  }

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }

for(int ti=0; ti<Lt; ++ti){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, ti*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_diag_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_diag_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, Lt*dim_block);
  
  const std::string dil_T = quarks[0].dilution_T;
  Eigen::VectorXi dilt1 = Eigen::VectorXi::Zero(Lt);
  Eigen::VectorXi dilt = Eigen::VectorXi::Zero(Lt);
  int ti1;
  if(dil_T == "TB"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti*num_dil_T/Lt;
    ti1 = (ti-1+Lt)%Lt;
    dilt1(ti) = ti1*num_dil_T/Lt;    
  }  
  }

  else if(dil_T == "TI"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti%num_dil_T;
    ti1 = (ti-1+Lt)%Lt;
    dilt1(ti) = ti1%num_dil_T;    
  }  
  }

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }

for(int ti=0; ti<Lt; ++ti){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, ti*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_diag_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_col(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

  int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_col mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_col_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);
  int t1=(t+1)%Lt;
  int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_col_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_col_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);
  int t1=(t-1+Lt)%Lt;
  int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_col_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkbarSink_row(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

  int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_row mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_row_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
  int t1=(t+1)%Lt;
  int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_row_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_row_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
  int t1=(t-1+Lt)%Lt;
  int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_row_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

/********************************************/
/********** SinkSinkbar  ********************/
/********************************************/

void MesonFunc_SinkSinkbar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, number_of_inversion);
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);

    MesonFunc += gammasink.value[i]*temp1.adjoint()*temp2;
  }

//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkSinkbar_diag(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_diag mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkSinkbar_diag_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, Lt*dim_block);

  const std::string dil_T = quarks[0].dilution_T;
  Eigen::VectorXi dilt1 = Eigen::VectorXi::Zero(Lt);
  Eigen::VectorXi dilt = Eigen::VectorXi::Zero(Lt);
  int ti1;
  if(dil_T == "TB"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti*num_dil_T/Lt;
    ti1 = (ti+1)%Lt;
    dilt1(ti) = ti1*num_dil_T/Lt;    
  }  
  }

  else if(dil_T == "TI"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti%num_dil_T;
    ti1 = (ti+1)%Lt;
    dilt1(ti) = ti1%num_dil_T;    
  }  
  }

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int ti=0; ti<Lt; ++ti){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, ti*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_diag_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkSinkbar_diag_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, Lt*dim_block);

  const std::string dil_T = quarks[0].dilution_T;
  Eigen::VectorXi dilt1 = Eigen::VectorXi::Zero(Lt);
  Eigen::VectorXi dilt = Eigen::VectorXi::Zero(Lt);
  int ti1;
  if(dil_T == "TB"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti*num_dil_T/Lt;
    ti1 = (ti-1+Lt)%Lt;
    dilt1(ti) = ti1*num_dil_T/Lt;    
  }  
  }

  else if(dil_T == "TI"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti%num_dil_T;
    ti1 = (ti-1+Lt)%Lt;
    dilt1(ti) = ti1%num_dil_T;    
  }  
  }

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int ti=0; ti<Lt; ++ti){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, ti*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
} // ti loop ends here
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_diag_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

void MesonFunc_SinkSinkbar_col(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_col mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_col_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

int t1=(t+1)%Lt;
int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_col_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_col_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

int t1=(t-1+Lt)%Lt;
int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_col_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_row(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_row mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_row_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

int t1=(t+1)%Lt;
int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_row_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_row_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

int t1=(t-1+Lt)%Lt;
int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*VdaggerV)*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_row_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


/**********************************************/
/******* build source source Matrix ***********/
/**********************************************/

void MesonFunc_SourcebarSource(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec1, Eigen::VectorXcd* rnd_vec2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
// clock_t time = clock();
 gammastruct gammasource;
 create_gamma(gammasource, gamma_index);
// const int Lt = global_data->get_Lt();
 const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
 const std::vector<quark> quarks = global_data->get_quarks();
 const int number_of_dilution_E = quarks[0].number_of_dilution_E;
 const std::string dil_E = quarks[0].dilution_E;

 Eigen::MatrixXcd* temp1 = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec1)
 Eigen::MatrixXcd* temp2 = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec2)
 for(int gamma_i=0; gamma_i<4; ++gamma_i){
   temp1[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
   temp2[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
 }
 Eigen::MatrixXcd temp3 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E); //( gamma5 * P^b * rnd_vec1 )
 Eigen::MatrixXcd temp4 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E); // ( VdaggerV * P^b *rnd_vec2 )

 
 if(dil_E == "EB"){
   for(int gamma_i=0; gamma_i<4; ++gamma_i){
//     for(int t=0; t<Lt; ++t)
      for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
         int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
         temp1[gamma_i](ev_i, dil_ev_i) = (*rnd_vec1)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i); 
         temp2[gamma_i](ev_i, dil_ev_i) = (*rnd_vec2)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
      } 
   }
 }
 
 else if(dil_E == "EI"){
   for(int gamma_i=0; gamma_i<4; ++gamma_i){
//     for(int t=0; t<Lt; ++t)
      for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
         int dil_ev_i = ev_i%number_of_dilution_E;
         temp1[gamma_i](ev_i, dil_ev_i) = (*rnd_vec1)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i); 
         temp2[gamma_i](ev_i, dil_ev_i) = (*rnd_vec2)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
      } 
   }
 }

 else{
    std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
    exit(0);
 }

    MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, 4*number_of_dilution_E);
    for(int i=0; i<4; ++i){
       int gamma_row = gammasource.index[i]/4;
       int gamma_col = gammasource.index[i]%4;
       if(gamma_row == 2 || gamma_row == 3)
          temp3 = -1.0*temp1[gamma_row];
       else
          temp3 = temp1[gamma_row];
      temp4 = (*VdaggerV)*temp2[gamma_col];
     MesonFunc.block(gamma_row*number_of_dilution_E, gamma_col*number_of_dilution_E, number_of_dilution_E, number_of_dilution_E) = gammasource.value[i]*temp3.adjoint()*temp4;
   }

delete[] temp1;
delete[] temp2; 
//  time = clock() - time; 
//  printf("\t\tbuild sourcebarsource mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}

void MesonFunc_SourceSourcebar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec1, Eigen::VectorXcd* rnd_vec2, Eigen::MatrixXcd* VdaggerV, const int gamma_index){
// clock_t time = clock();
 gammastruct gammasource;
 create_gamma(gammasource, gamma_index);
// const int Lt = global_data->get_Lt();
 const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
 const std::vector<quark> quarks = global_data->get_quarks();
 const int number_of_dilution_E = quarks[0].number_of_dilution_E;
 const std::string dil_E = quarks[0].dilution_E;

 Eigen::MatrixXcd* temp1 = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec1)
 Eigen::MatrixXcd* temp2 = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec2)
 for(int gamma_i=0; gamma_i<4; ++gamma_i){
   temp1[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
   temp2[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
 }
 Eigen::MatrixXcd temp3 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E); //( gamma5 * P^b * rnd_vec1 )
 Eigen::MatrixXcd temp4 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E); // ( VdaggerV * P^b *rnd_vec2 )

 
 if(dil_E == "EB"){
   for(int gamma_i=0; gamma_i<4; ++gamma_i){
//     for(int t=0; t<Lt; ++t)
      for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
         int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
         temp1[gamma_i](ev_i, dil_ev_i) = (*rnd_vec1)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i); 
         temp2[gamma_i](ev_i, dil_ev_i) = (*rnd_vec2)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
      } 
   }
 }
 
 else if(dil_E == "EI"){
   for(int gamma_i=0; gamma_i<4; ++gamma_i){
//     for(int t=0; t<Lt; ++t)
      for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
         int dil_ev_i = ev_i%number_of_dilution_E;
         temp1[gamma_i](ev_i, dil_ev_i) = (*rnd_vec1)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i); 
         temp2[gamma_i](ev_i, dil_ev_i) = (*rnd_vec2)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
      } 
   }
 }

 else{
    std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
    exit(0);
 }

    MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, 4*number_of_dilution_E);
    for(int i=0; i<4; ++i){
       int gamma_row = gammasource.index[i]/4;
       int gamma_col = gammasource.index[i]%4;
       if(gamma_col == 2 || gamma_col == 3)
          temp4 = -1.0*(*VdaggerV)*temp2[gamma_col];
       else
          temp4 = (*VdaggerV)*temp2[gamma_col];
      temp3 = temp1[gamma_row];
     MesonFunc.block(gamma_row*number_of_dilution_E, gamma_col*number_of_dilution_E, number_of_dilution_E, number_of_dilution_E) = gammasource.value[i]*temp3.adjoint()*temp4;
   }

delete[] temp1;
delete[] temp2; 
//  time = clock() - time; 
//  printf("\t\tbuild sourcebarsource mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}

/************* build SourceSink Matrix *********************/
void MesonFunc_SourceSink(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec, Eigen::MatrixXcd* peram, Eigen::MatrixXcd* VdaggerV, const int gamma_index){

       // clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

        Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }
         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

//        for(int t=0; t<Lt; t++){
           MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, number_of_inversion);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              temp_peram = (*VdaggerV)*(*peram).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
              MesonFunc.block(gamma_row*number_of_dilution_E, 0, number_of_dilution_E, number_of_inversion) = gammaMatrix.value[i]*(temp[gamma_row]).adjoint()*temp_peram;
           }
//        }

//  time = clock() - time;
//  printf("\t\tbuild SourceSink mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/************* build SinkSource Matrix *********************/
void MesonFunc_SinkSource(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram, Eigen::VectorXcd* rnd_vec, Eigen::MatrixXcd* VdaggerV, const int gamma_index){

  //      clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

//      Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }

         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

           MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, 4*number_of_dilution_E);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              MesonFunc.block(0, gamma_col*number_of_dilution_E, number_of_inversion, number_of_dilution_E) = gammaMatrix.value[i]*((*peram).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion)).adjoint()*(*VdaggerV)*temp[gamma_col];
           }

//  time = clock() - time;
//  printf("\t\tbuild Source Sink mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/************* build SourcebarSinkbar Matrix *********************/
void MesonFunc_SourcebarSinkbar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec, Eigen::MatrixXcd* peram, Eigen::MatrixXcd* VdaggerV, const int gamma_index){

 //       clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

        Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }

         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

        temp[2] *= -1.0;
        temp[3] *= -1.0;

           MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, number_of_inversion);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              temp_peram = (*VdaggerV)*(*peram).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
                if(gamma_col==2 || gamma_col==3)
                temp_peram *= -1.0;
              MesonFunc.block(gamma_row*number_of_dilution_E, 0, number_of_dilution_E, number_of_inversion) = gammaMatrix.value[i]*(temp[gamma_row]).adjoint()*temp_peram;
           }

//  time = clock() - time;
//  printf("\t\tbuild Sourcebar Sinkbar mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/************* build Sinkbar Sourcebar Matrix *********************/
void MesonFunc_SinkbarSourcebar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram, Eigen::VectorXcd* rnd_vec, Eigen::MatrixXcd* VdaggerV, const int gamma_index){

 //       clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

        Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }

         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

        temp[2] *= -1.0;
        temp[3] *= -1.0;

           MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, 4*number_of_dilution_E);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              temp_peram = (*peram).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
                if(gamma_row==2 || gamma_row==3)
                temp_peram *= -1.0;
              MesonFunc.block(0, gamma_col*number_of_dilution_E, number_of_inversion, number_of_dilution_E) = gammaMatrix.value[i]*temp_peram.adjoint()*(*VdaggerV)*temp[gamma_col];
           }

//  time = clock() - time;
//  printf("\t\tbuild Sinkbar Sourcebar mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/*************************************************************/
/*************************************************************/
/*** Meson functions for zero momentum   *********************/
/*************************************************************/
/*************************************************************/

/*******build sourcebar source for p=0 ************/
void MesonFunc_SourcebarSource(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec1, Eigen::VectorXcd* rnd_vec2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasource;
  create_gamma(gammasource, gamma_index);
  const std::vector<quark> quarks = global_data->get_quarks(); 
  const int number_of_dilution_E = quarks[0].number_of_dilution_E; 
  const std::string dil_E = quarks[0].dilution_E; 
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  Eigen::MatrixXcd rnd_vec_bar = (*rnd_vec1).block(t*4*number_of_eigen_vec, 0, 4*number_of_eigen_vec, 1);
//  for(int t=0; t<Lt; t++)
    rnd_vec_bar.block(2*number_of_eigen_vec, 0, 2*number_of_eigen_vec, 1) *= -1.0;
  int ev_i;
 
  if(dil_E == "EI"){
       MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, 4*number_of_dilution_E);
    for(int i=0; i<4; ++i){      
       int gamma_row = gammasource.index[i]/4;
       int gamma_col = gammasource.index[i]%4;
          for(int ev=0; ev<number_of_eigen_vec; ++ev){
             ev_i = ev%number_of_dilution_E;
             MesonFunc(gamma_row*number_of_dilution_E  + ev_i, gamma_col*number_of_dilution_E + ev_i) += gammasource.value[i]*std::conj(rnd_vec_bar(gamma_row*number_of_eigen_vec + ev))*(*rnd_vec2)(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec + ev);
          }
       }     
  }
  
  else if(dil_E == "EB"){
       MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, 4*number_of_dilution_E);
    for(int i=0; i<4; ++i){            
       int gamma_row = gammasource.index[i]/4;
       int gamma_col = gammasource.index[i]%4;
          for(int ev=0; ev<number_of_eigen_vec; ++ev){             
             ev_i = ev*number_of_dilution_E/number_of_eigen_vec; 
             MesonFunc(gamma_row*number_of_dilution_E  + ev_i, gamma_col*number_of_dilution_E + ev_i) += gammasource.value[i]*std::conj(rnd_vec_bar(gamma_row*number_of_eigen_vec + ev))*(*rnd_vec2)(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec + ev);
          }
       }     
  }

//  time = clock() - time;
//  printf("\t\tbuild SourcebarSource mesonfunction at zero momentum success -%.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}

void MesonFunc_SourceSourcebar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec1, Eigen::VectorXcd* rnd_vec2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasource;
  create_gamma(gammasource, gamma_index);
  const std::vector<quark> quarks = global_data->get_quarks(); 
  const int number_of_dilution_E = quarks[0].number_of_dilution_E; 
  const std::string dil_E = quarks[0].dilution_E; 
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  Eigen::MatrixXcd rnd_vec_bar = (*rnd_vec2).block(t*4*number_of_eigen_vec, 0, 4*number_of_eigen_vec, 1);
//  for(int t=0; t<Lt; t++)
    rnd_vec_bar.block(2*number_of_eigen_vec, 0, 2*number_of_eigen_vec, 1) *= -1.0;
  int ev_i;
 
  if(dil_E == "EI"){
       MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, 4*number_of_dilution_E);
    for(int i=0; i<4; ++i){      
       int gamma_row = gammasource.index[i]/4;
       int gamma_col = gammasource.index[i]%4;
          for(int ev=0; ev<number_of_eigen_vec; ++ev){
             ev_i = ev%number_of_dilution_E;
             MesonFunc(gamma_row*number_of_dilution_E  + ev_i, gamma_col*number_of_dilution_E + ev_i) += gammasource.value[i]*std::conj((*rnd_vec1)(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec + ev))*rnd_vec_bar(gamma_col*number_of_eigen_vec + ev);
          }
       }     
  }
  
  else if(dil_E == "EB"){
       MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, 4*number_of_dilution_E);
    for(int i=0; i<4; ++i){            
       int gamma_row = gammasource.index[i]/4;
       int gamma_col = gammasource.index[i]%4;
          for(int ev=0; ev<number_of_eigen_vec; ++ev){             
             ev_i = ev*number_of_dilution_E/number_of_eigen_vec; 
             MesonFunc(gamma_row*number_of_dilution_E  + ev_i, gamma_col*number_of_dilution_E + ev_i) += gammasource.value[i]*std::conj((*rnd_vec1)(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec + ev))*rnd_vec_bar(gamma_col*number_of_eigen_vec + ev);
          }
       }     
  }

//  time = clock() - time;
//  printf("\t\tbuild SourcebarSource mesonfunction at zero momentum success -%.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}


/*********build sinkbar sink for p=0 ***************/
void MesonFunc_SinkbarSink(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
// clock_t time = clock();
 gammastruct gammasink;
 create_gamma(gammasink, gamma_index);
 const std::vector<quark> quarks = global_data->get_quarks();
 const int number_of_dilution_E = quarks[0].number_of_dilution_E;
 const int number_of_dilution_T = quarks[0].number_of_dilution_T;
// const int Lt = global_data->get_Lt();
 const int number_of_inversion = 4*number_of_dilution_E*number_of_dilution_T;
 const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
 Eigen::MatrixXcd peram_bar = (*peram1).block(t*4*number_of_eigen_vec, 0, 4*number_of_eigen_vec, number_of_inversion);

 peram_bar.block(2*number_of_eigen_vec, 0, 2*number_of_eigen_vec, number_of_inversion) *= -1.0;

  MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, number_of_inversion);
 for(int i=0; i<4; ++i){
   int gamma_row = gammasink.index[i]/4;
   int gamma_col = gammasink.index[i]%4;
    MesonFunc += gammasink.value[i]*peram_bar.block(gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion).adjoint()*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
 }
// time = clock()-time;
//  printf("\t\tbuild SinkbarSink mesonfunction at zero momentum success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}

void MesonFunc_SinkbarSink_diag(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_diag mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_col(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2,  const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

  int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_col mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_col_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2,  const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);
  int t1=(t+1)%Lt;
  int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_col_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_col_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);
  int t1=(t-1+Lt)%Lt;
  int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_col_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkbarSink_row(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2,  const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

  int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_row mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_row_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
  int t1=(t+1)%Lt;
  int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_row_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkbarSink_row_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int num_dil_T = quarks[0].number_of_dilution_T;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = gamma5 * peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
  int t1=(t-1+Lt)%Lt;
  int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }
  
for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_row == 2 || gamma_row == 3) // (temp = gamma_5 * peram1)
      temp1 = -1.0*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinkbarsink_row_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkSinkbar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
// clock_t time = clock();
 gammastruct gammasink;
 create_gamma(gammasink, gamma_index);
 const std::vector<quark> quarks = global_data->get_quarks();
 const int number_of_dilution_E = quarks[0].number_of_dilution_E;
 const int number_of_dilution_T = quarks[0].number_of_dilution_T;
// const int Lt = global_data->get_Lt();
 const int number_of_inversion = 4*number_of_dilution_E*number_of_dilution_T;
 const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
 Eigen::MatrixXcd peram_bar = (*peram2).block(t*4*number_of_eigen_vec, 0, 4*number_of_eigen_vec, number_of_inversion);

 peram_bar.block(2*number_of_eigen_vec, 0, 2*number_of_eigen_vec, number_of_inversion) *= -1.0;

  MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, number_of_inversion);
 for(int i=0; i<4; ++i){
   int gamma_row = gammasink.index[i]/4;
   int gamma_col = gammasink.index[i]%4;
    MesonFunc += gammasink.value[i]*(*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion).adjoint()*peram_bar.block(gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
 }
// time = clock()-time;
//  printf("\t\tbuild SinkbarSink mesonfunction at zero momentum success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
}

void MesonFunc_SinkSinkbar_diag(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
//  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_diag mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkSinkbar_diag_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, Lt*dim_block);

  const std::string dil_T = quarks[0].dilution_T;
  Eigen::VectorXi dilt1 = Eigen::VectorXi::Zero(Lt);
  Eigen::VectorXi dilt = Eigen::VectorXi::Zero(Lt);
  int ti1;
  if(dil_T == "TB"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti*num_dil_T/Lt;
    ti1 = (ti+1)%Lt;
    dilt1(ti) = ti1*num_dil_T/Lt;    
  }  
  }

  else if(dil_T == "TI"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti%num_dil_T;
    ti1 = (ti+1)%Lt;
    dilt1(ti) = ti1%num_dil_T;    
  }  
  }

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int ti=0; ti<Lt; ++ti){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, ti*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_diag_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

void MesonFunc_SinkSinkbar_diag_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
//  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, Lt*dim_block);

  const std::string dil_T = quarks[0].dilution_T;
  Eigen::VectorXi dilt1 = Eigen::VectorXi::Zero(Lt);
  Eigen::VectorXi dilt = Eigen::VectorXi::Zero(Lt);
  int ti1;
  if(dil_T == "TB"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti*num_dil_T/Lt;
    ti1 = (ti-1+Lt)%Lt;
    dilt1(ti) = ti1*num_dil_T/Lt;    
  }  
  }

  else if(dil_T == "TI"){
  for(int ti=0; ti<Lt; ++ti){
    dilt(ti) = ti%num_dil_T;
    ti1 = (ti-1+Lt)%Lt;
    dilt1(ti) = ti1%num_dil_T;    
  }  
  }

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int ti=0; ti<Lt; ++ti){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt1(ti)*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt(ti)*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, ti*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
} // ti loop ends here
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_diag_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

void MesonFunc_SinkSinkbar_col(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_col mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_col_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

int t1=(t+1)%Lt;
int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_col_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_col_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, dim_block);

int t1=(t-1+Lt)%Lt;
int dil_t_col;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_col = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_col = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dil_t_col*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(dilt*dim_block, 0, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_col_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_row(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_row mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_row_shiftup(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
//  clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

int t1=(t+1)%Lt;
int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_row_shiftup mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }


void MesonFunc_SinkSinkbar_row_shiftdown(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram1, Eigen::MatrixXcd* peram2, const int gamma_index){
 // clock_t time = clock();
  gammastruct gammasink;
  create_gamma(gammasink, gamma_index);
  const int Lt = global_data->get_Lt();
  const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
  const std::vector<quark> quarks = global_data->get_quarks();
  const int number_of_inversion = quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D;
  const int dim_block = quarks[0].number_of_dilution_E*quarks[0].number_of_dilution_D; 
  const int num_dil_T = quarks[0].number_of_dilution_T;
  Eigen::MatrixXcd temp1 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp1 = peram1)
  Eigen::MatrixXcd temp2 = Eigen::MatrixXcd::Zero(number_of_eigen_vec, dim_block); //(temp2 = VdaggerV * gamma5 * peram2)
    MesonFunc = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);

int t1=(t-1+Lt)%Lt;
int dil_t_row;
  const std::string dil_T = quarks[0].dilution_T;
  if(dil_T == "TB"){
    dil_t_row = t1*num_dil_T/Lt;
  }

  else if(dil_T == "TI"){
    dil_t_row = t1%num_dil_T;    
  }  

  else{
  std::cout<<"unknown time dilution scheme."<<std::endl;
  exit(0);
  }


for(int dilt=0; dilt<num_dil_T; ++dilt){
  for(int i=0; i<4; ++i){
    int gamma_row = gammasink.index[i]/4;
    int gamma_col = gammasink.index[i]%4;
    if(gamma_col == 2 || gamma_col == 3) // (temp2 = gamma_5 * VdaggerV * peram2)
      temp2 = -1.0*(*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
    else 
      temp2 = (*peram2).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, dilt*dim_block, number_of_eigen_vec, dim_block);
      temp1 = (*peram1).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, dil_t_row*dim_block, number_of_eigen_vec, dim_block);

    MesonFunc.block(0, dilt*dim_block, dim_block, dim_block) += gammasink.value[i]*temp1.adjoint()*temp2;
  }
}
//  time = clock() - time;
//  printf("\t\tbuild sinksinkbar_row_shiftdown mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
 }

/************* build SourceSink Matrix *********************/
void MesonFunc_SourceSink(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec, Eigen::MatrixXcd* peram, const int gamma_index){

 //       clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

        Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }
         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

//        for(int t=0; t<Lt; t++){
           MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, number_of_inversion);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              temp_peram = (*peram).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
              MesonFunc.block(gamma_row*number_of_dilution_E, 0, number_of_dilution_E, number_of_inversion) = gammaMatrix.value[i]*(temp[gamma_row]).adjoint()*temp_peram;
           }
//        }

//  time = clock() - time;
//  printf("\t\tbuild SourceSink mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/************* build SinkSource Matrix *********************/
void MesonFunc_SinkSource(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram, Eigen::VectorXcd* rnd_vec, const int gamma_index){

  //      clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

//      Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }

         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

           MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, 4*number_of_dilution_E);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              MesonFunc.block(0, gamma_col*number_of_dilution_E, number_of_inversion, number_of_dilution_E) = gammaMatrix.value[i]*((*peram).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion)).adjoint()*temp[gamma_col];
           }

//  time = clock() - time;
//  printf("\t\tbuild Source Sink mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/************* build SourcebarSinkbar Matrix *********************/
void MesonFunc_SourcebarSinkbar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::VectorXcd* rnd_vec, Eigen::MatrixXcd* peram, const int gamma_index){

 //       clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

        Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }

         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

        temp[2] *= -1.0;
        temp[3] *= -1.0;

           MesonFunc = Eigen::MatrixXcd::Zero(4*number_of_dilution_E, number_of_inversion);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              temp_peram = (*peram).block(t*4*number_of_eigen_vec + gamma_col*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
                if(gamma_col==2 || gamma_col==3)
                temp_peram *= -1.0;
              MesonFunc.block(gamma_row*number_of_dilution_E, 0, number_of_dilution_E, number_of_inversion) = gammaMatrix.value[i]*(temp[gamma_row]).adjoint()*temp_peram;
           }

//  time = clock() - time;
//  printf("\t\tbuild Sourcebar Sinkbar mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

/************* build Sinkbar Sourcebar Matrix *********************/
void MesonFunc_SinkbarSourcebar(int t, Eigen::MatrixXcd& MesonFunc, Eigen::MatrixXcd* peram, Eigen::VectorXcd* rnd_vec, const int gamma_index){

  //      clock_t time = clock();
        gammastruct gammaMatrix;
        create_gamma(gammaMatrix, gamma_index);
        const std::vector<quark> quarks = global_data->get_quarks();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_dilution_D = quarks[0].number_of_dilution_D;
        const int number_of_inversion = number_of_dilution_E*number_of_dilution_T*number_of_dilution_D;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const std::string dil_E = quarks[0].dilution_E;
        const std::string dil_T = quarks[0].dilution_T;
        Eigen::MatrixXcd* temp = new Eigen::MatrixXcd[4]; // (P^b*rnd_vec)
        for(int gamma_i=0; gamma_i<4; ++gamma_i){
           temp[gamma_i] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_dilution_E);
         }

        Eigen::MatrixXcd temp_peram = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_inversion);

        if(dil_E == "EB"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i*number_of_dilution_E/number_of_eigen_vec;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
              }
        }

         else if(dil_E == "EI"){
           for(int gamma_i=0; gamma_i<4; ++gamma_i)
              for(int ev_i=0; ev_i<number_of_eigen_vec; ++ev_i){
                 int dil_ev_i = ev_i%number_of_dilution_E;
                 temp[gamma_i](ev_i, dil_ev_i) = (*rnd_vec)(t*4*number_of_eigen_vec + gamma_i*number_of_eigen_vec + ev_i);
           }
         }

        else{
          std::cout<<"unknown eigen vector dilution scheme."<<std::endl;
          exit(0);
        }

        temp[2] *= -1.0;
        temp[3] *= -1.0;

           MesonFunc = Eigen::MatrixXcd::Zero(number_of_inversion, 4*number_of_dilution_E);
           for(int i=0; i<4; i++){
              int gamma_row = gammaMatrix.index[i]/4;
              int gamma_col = gammaMatrix.index[i]%4;
              temp_peram = (*peram).block(t*4*number_of_eigen_vec + gamma_row*number_of_eigen_vec, 0, number_of_eigen_vec, number_of_inversion);
                if(gamma_row==2 || gamma_row==3)
                temp_peram *= -1.0;
              MesonFunc.block(0, gamma_col*number_of_dilution_E, number_of_inversion, number_of_dilution_E) = gammaMatrix.value[i]*temp_peram.adjoint()*temp[gamma_col];
           }

//  time = clock() - time;
//  printf("\t\tbuild Sinkbar Sourcebar mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);

}

