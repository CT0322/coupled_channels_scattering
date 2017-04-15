//============================================================================
// Name        : LapHs.cpp
// Author      : BK
// Version     :
// Copyright   : Copies are prohibited so far
// Description : stochastic LapH code
//============================================================================

#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>
#include <iomanip>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "GlobalData.h"
#include "MesonFunctions.h"
#include "ReadWrite.h"
#include "Contrac_SiSiSiSi_SoSoSoSo.h"
#include "Contrac_twopoint.h"
#include "Momentum_lookup.h"

int main (int ac, char* av[]) {
        
        clock_t time_start = clock();
        clock_t time = clock();
        Eigen::initParallel();
        Eigen::setNbThreads(4);

/*********** global variables ******************/
        GlobalData* global_data = GlobalData::Instance();
        global_data->read_parameters(ac, av);

	const int Lt = global_data->get_Lt();
        const int Lx = global_data->get_Lx();
   //     const int Ly = global_data->get_Ly();
   //     const int Lz = global_data->get_Lz();
        const int end_config = global_data->get_end_config();
        const int delta_config = global_data->get_delta_config();
        const int start_config = global_data->get_start_config();
        const std::vector<quark> quarks = global_data->get_quarks();
        const int n_rndvec_c = quarks[0].number_of_rnd_vec;
        const int n_rndvec_u = quarks[1].number_of_rnd_vec;
        const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
        const int number_of_dilution_E = quarks[0].number_of_dilution_E;
        const int number_of_dilution_T = quarks[0].number_of_dilution_T;
        const int number_of_inversion = quarks[0].number_of_dilution_D*quarks[0].number_of_dilution_T*quarks[0].number_of_dilution_E;
	const int dim_block = 4*number_of_dilution_E;
	const std::string path_output = global_data->get_path_output();
	const int n_rnd = n_rndvec_u*n_rndvec_c;
/*************** end global variables ***************/

        char outfile[400];

/**************** momentum ********************/

	int norm=0;

std::cout<<"start momentum difinition"<<std::endl;
	const int number_of_mom=5;
	const int TP0_dim = 4;
	const int np[5] = {1, 6, 12, 8, 6};
	int p0[1][3] = {{0, 0, 0}}; 
        int p1[6][3] = {{0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0}};
	int p2[12][3] = {{0, 1, 1}, {0, -1, -1}, {1, 1, 0}, {-1, -1, 0}, {1, 0, 1}, {-1, 0, -1}, {0, 1, -1}, {0, -1, 1}, {1, -1, 0}, {-1, 1, 0}, {1, 0, -1}, {-1, 0, 1}};
        int p3[8][3] = {{1, 1, 1}, {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1}, {1, -1, -1}};
        int p4[6][3] = {{0, 0, 2}, {0, 0, -2}, {0, 2, 0}, {0, -2, 0}, {2, 0, 0}, {-2, 0, 0}};
        Eigen::VectorXi** p = new Eigen::VectorXi*[5];
        for(int i=0; i<5; ++i){
                p[i] = new Eigen::VectorXi[np[i]];
                for(int j=0; j<np[i]; ++j){
                        p[i][j] = Eigen::VectorXi::Zero(3);
                }
        }
        for(int j=0; j<1; ++j)
        for(int k=0; k<3; ++k)
                p[0][j][k] = p0[j][k];

        for(int j=0; j<6; ++j)
        for(int k=0; k<3; ++k)
                p[1][j][k] = p1[j][k];

        for(int j=0; j<12; ++j)
        for(int k=0; k<3; ++k)
                p[2][j][k] = p2[j][k];

        for(int j=0; j<8; ++j)
        for(int k=0; k<3; ++k)
                p[3][j][k] = p3[j][k];

        for(int j=0; j<6; ++j)
        for(int k=0; k<3; ++k)
                p[4][j][k] = p4[j][k];
//non-zero total momentum indices
        int TP1_dim = 8;
        int n_TP1[TP1_dim];
        Eigen::Vector2i**** TP1 = new Eigen::Vector2i***[np[1]/2];
        for(int i=0; i<np[1]/2; ++i){
                TP1[i] = new Eigen::Vector2i**[TP1_dim];
                Momentum_lookup(p[1][2*i], 0, 1, TP1[i][0], n_TP1[0]);
                Momentum_lookup(p[1][2*i], 1, 0, TP1[i][1], n_TP1[1]);
                Momentum_lookup(p[1][2*i], 1, 2, TP1[i][2], n_TP1[2]);
                Momentum_lookup(p[1][2*i], 2, 1, TP1[i][3], n_TP1[3]);
                Momentum_lookup(p[1][2*i], 2, 3, TP1[i][4], n_TP1[4]);
                Momentum_lookup(p[1][2*i], 3, 2, TP1[i][5], n_TP1[5]);
                Momentum_lookup(p[1][2*i], 1, 4, TP1[i][6], n_TP1[6]);
                Momentum_lookup(p[1][2*i], 4, 1, TP1[i][7], n_TP1[7]);
        }

        int TP2_dim = 8;
        int n_TP2[TP2_dim];
        Eigen::Vector2i**** TP2 = new Eigen::Vector2i***[np[2]/2];
        for(int i=0; i<np[2]/2; ++i){
                TP2[i] = new Eigen::Vector2i**[TP2_dim];
                Momentum_lookup(p[2][2*i], 1, 1, TP2[i][0], n_TP2[0]);
                Momentum_lookup(p[2][2*i], 0, 2, TP2[i][1], n_TP2[1]);
                Momentum_lookup(p[2][2*i], 2, 0, TP2[i][2], n_TP2[2]);
                Momentum_lookup(p[2][2*i], 2, 2, TP2[i][3], n_TP2[3]);
                Momentum_lookup(p[2][2*i], 1, 3, TP2[i][4], n_TP2[4]);
                Momentum_lookup(p[2][2*i], 3, 1, TP2[i][5], n_TP2[5]);
                Momentum_lookup(p[2][2*i], 2, 4, TP2[i][6], n_TP2[6]);
                Momentum_lookup(p[2][2*i], 4, 2, TP2[i][7], n_TP2[7]);
        }

        int TP3_dim = 6;
        int n_TP3[TP3_dim];
        Eigen::Vector2i**** TP3 = new Eigen::Vector2i***[np[3]/2];
        for(int i=0; i<np[3]/2; ++i){
                TP3[i] = new Eigen::Vector2i**[TP3_dim];
                Momentum_lookup(p[3][2*i], 1, 2, TP3[i][0], n_TP3[0]);
                Momentum_lookup(p[3][2*i], 2, 1, TP3[i][1], n_TP3[1]);
                Momentum_lookup(p[3][2*i], 0, 3, TP3[i][2], n_TP3[2]);
                Momentum_lookup(p[3][2*i], 3, 0, TP3[i][3], n_TP3[3]);
                Momentum_lookup(p[3][2*i], 4, 3, TP3[i][4], n_TP3[4]);
                Momentum_lookup(p[3][2*i], 3, 4, TP3[i][5], n_TP3[5]);
        }


std::cout<<"endl momentum definition"<<std::endl;

        Eigen::MatrixXcd*** VdaggerV = new Eigen::MatrixXcd**[number_of_mom];
	for(int inp=0; inp<number_of_mom; ++inp){
		VdaggerV[inp] = new Eigen::MatrixXcd*[np[inp]];
        	for(int i=0; i<np[inp]; ++i){
                	VdaggerV[inp][i] = new Eigen::MatrixXcd[Lt];
			for(int t=0; t<Lt; ++t)
				VdaggerV[inp][i][t] = Eigen::MatrixXcd::Zero(number_of_eigen_vec, number_of_eigen_vec);
        	}
	}

std::cout<<"end VdaggerV memory allocate"<<std::endl;
/******* meson functions and VdaggerV memory allocate ***********/ 
	Eigen::MatrixXcd**** SinkbarSink_dbarcp = new Eigen::MatrixXcd***[n_rnd];
	Eigen::MatrixXcd***** SourcebarSource_dbarcp = new Eigen::MatrixXcd****[n_rnd];
	Eigen::MatrixXcd**** SinkbarSink_cmbaru = new Eigen::MatrixXcd***[n_rnd];
	Eigen::MatrixXcd***** SourcebarSource_cmbaru = new Eigen::MatrixXcd****[n_rnd];

        Eigen::MatrixXcd**** SinkbarSink_cmbarcp = new Eigen::MatrixXcd***[n_rnd];
        Eigen::MatrixXcd***** SourcebarSource_cmbarcp = new Eigen::MatrixXcd****[n_rnd];

        Eigen::MatrixXcd*** SinkbarSink_dbaru = new Eigen::MatrixXcd**[n_rnd];
        Eigen::MatrixXcd**** SourcebarSource_dbaru = new Eigen::MatrixXcd***[n_rnd];

        for(int rnd=0; rnd<n_rnd; ++rnd){
		SinkbarSink_dbarcp[rnd] = new Eigen::MatrixXcd**[number_of_mom];
		SourcebarSource_dbarcp[rnd] = new Eigen::MatrixXcd***[number_of_mom];
		SinkbarSink_cmbaru[rnd] = new Eigen::MatrixXcd**[number_of_mom];
		SourcebarSource_cmbaru[rnd] = new Eigen::MatrixXcd***[number_of_mom];

                SinkbarSink_cmbarcp[rnd] = new Eigen::MatrixXcd**[number_of_mom];
                SourcebarSource_cmbarcp[rnd] = new Eigen::MatrixXcd***[number_of_mom];

                SinkbarSink_dbaru[rnd] = new Eigen::MatrixXcd*[number_of_mom];
                SourcebarSource_dbaru[rnd] = new Eigen::MatrixXcd**[number_of_mom];

		for(int inp=0; inp<number_of_mom; ++inp){
		SinkbarSink_dbarcp[rnd][inp] = new Eigen::MatrixXcd*[np[inp]];
		SourcebarSource_dbarcp[rnd][inp] = new Eigen::MatrixXcd**[np[inp]];
		SinkbarSink_cmbaru[rnd][inp] = new Eigen::MatrixXcd*[np[inp]];
		SourcebarSource_cmbaru[rnd][inp] = new Eigen::MatrixXcd**[np[inp]];

                SinkbarSink_cmbarcp[rnd][inp] = new Eigen::MatrixXcd*[np[inp]];
                SourcebarSource_cmbarcp[rnd][inp] = new Eigen::MatrixXcd**[np[inp]];

                SinkbarSink_dbaru[rnd][inp] = new Eigen::MatrixXcd[np[inp]];
                SourcebarSource_dbaru[rnd][inp] = new Eigen::MatrixXcd*[np[inp]];

			for(int i=0; i<np[inp]; ++i){
				SinkbarSink_dbarcp[rnd][inp][i] = new Eigen::MatrixXcd[4];
				SourcebarSource_dbarcp[rnd][inp][i] = new Eigen::MatrixXcd*[4];
				SinkbarSink_cmbaru[rnd][inp][i] = new Eigen::MatrixXcd[4];
				SourcebarSource_cmbaru[rnd][inp][i] = new Eigen::MatrixXcd*[4];
       
                                SinkbarSink_cmbarcp[rnd][inp][i] = new Eigen::MatrixXcd[4];
                                SourcebarSource_cmbarcp[rnd][inp][i] = new Eigen::MatrixXcd*[4];

                                SinkbarSink_dbaru[rnd][inp][i] = new Eigen::MatrixXcd[4];
                                SourcebarSource_dbaru[rnd][inp][i] = new Eigen::MatrixXcd[4];

				for(int k=0; k<4; ++k){
					SinkbarSink_dbarcp[rnd][inp][i][k] = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
					SourcebarSource_dbarcp[rnd][inp][i][k] = new Eigen::MatrixXcd[Lt];
					SinkbarSink_cmbaru[rnd][inp][i][k] = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
					SourcebarSource_cmbaru[rnd][inp][i][k] = new Eigen::MatrixXcd[Lt];

					for(int t=0; t<Lt; ++t){
						SourcebarSource_dbarcp[rnd][inp][i][k][t] = Eigen::MatrixXcd::Zero(dim_block, dim_block);
						SourcebarSource_cmbaru[rnd][inp][i][k][t] = Eigen::MatrixXcd::Zero(dim_block, dim_block);
					}
				}

                               for(int k=0; k<4; ++k){
                                        SinkbarSink_cmbarcp[rnd][inp][i][k] = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
                                        SourcebarSource_cmbarcp[rnd][inp][i][k] = new Eigen::MatrixXcd[Lt];

                                        for(int t=0; t<Lt; ++t)
                                                SourcebarSource_cmbarcp[rnd][inp][i][k][t] = Eigen::MatrixXcd::Zero(dim_block, dim_block);
                                }
 
                              for(int k=0; k<4; ++k){
                                        SinkbarSink_dbaru[rnd][inp][i][k] = Eigen::MatrixXcd::Zero(dim_block, number_of_inversion);
                                        SourcebarSource_dbaru[rnd][inp][i][k] = new Eigen::MatrixXcd[Lt];

                                        for(int t=0; t<Lt; ++t)
                                                SourcebarSource_dbaru[rnd][inp][i][k][t] = Eigen::MatrixXcd::Zero(dim_block, dim_block);
                                }


			}
		}	
	} // rnd loop ends here

std::cout<<"end mesonfunctions memory allocate"<<std::endl;


//allocate memory for correlation function TP0	
	Eigen::VectorXcd***** corr_DDstar_DDstar = new Eigen::VectorXcd****[TP0_dim];
	for(int inp1=0; inp1<TP0_dim; ++inp1){
		corr_DDstar_DDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
		for(int inp2 =0; inp2<TP0_dim; ++inp2){
			corr_DDstar_DDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
			for(int i=0; i<np[inp1]*np[inp2]; ++i){
				corr_DDstar_DDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
				for(int k=0; k<3; ++k){
					corr_DDstar_DDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[4];
					for(int n=0; n<4; ++n)
						corr_DDstar_DDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
				}
			}
		}
	}

	Eigen::VectorXcd***** corr_DstarDstar_DstarDstar = new Eigen::VectorXcd****[TP0_dim];
	for(int inp1=0; inp1<TP0_dim; ++inp1){
		corr_DstarDstar_DstarDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
		for(int inp2 =0; inp2<TP0_dim; ++inp2){
			corr_DstarDstar_DstarDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
			for(int i=0; i<np[inp1]*np[inp2]; ++i){
				corr_DstarDstar_DstarDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
				for(int k=0; k<3; ++k){
					corr_DstarDstar_DstarDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[1];
					for(int n=0; n<1; ++n)
						corr_DstarDstar_DstarDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
				}
			}
		}
	}


	Eigen::VectorXcd***** corr_DDstar_DstarDstar = new Eigen::VectorXcd****[TP0_dim];
	for(int inp1=0; inp1<TP0_dim; ++inp1){
		corr_DDstar_DstarDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
		for(int inp2 =0; inp2<TP0_dim; ++inp2){
			corr_DDstar_DstarDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
			for(int i=0; i<np[inp1]*np[inp2]; ++i){
				corr_DDstar_DstarDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
				for(int k=0; k<3; ++k){
					corr_DDstar_DstarDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
					for(int n=0; n<2; ++n)
						corr_DDstar_DstarDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
				}
			}
		}
	}

	Eigen::VectorXcd***** corr_DstarDstar_DDstar = new Eigen::VectorXcd****[TP0_dim];
	for(int inp1=0; inp1<TP0_dim; ++inp1){
		corr_DstarDstar_DDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
		for(int inp2 =0; inp2<TP0_dim; ++inp2){
			corr_DstarDstar_DDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
			for(int i=0; i<np[inp1]*np[inp2]; ++i){
				corr_DstarDstar_DDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
				for(int k=0; k<3; ++k){
					corr_DstarDstar_DDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
					for(int n=0; n<2; ++n)
						corr_DstarDstar_DDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
				}
			}
		}
	}
    
        Eigen::VectorXcd***** corr_JPsiPi_JPsiPi = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_JPsiPi_JPsiPi[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_JPsiPi_JPsiPi[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_JPsiPi_JPsiPi[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_JPsiPi_JPsiPi[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_JPsiPi_JPsiPi[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_EtacRho_EtacRho = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_EtacRho_EtacRho[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_EtacRho_EtacRho[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_EtacRho_EtacRho[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_EtacRho_EtacRho[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_EtacRho_EtacRho[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }



        Eigen::VectorXcd***** corr_DDstar_JPsiPi = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_DDstar_JPsiPi[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_DDstar_JPsiPi[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_DDstar_JPsiPi[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_DDstar_JPsiPi[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<2; ++n)
                                                corr_DDstar_JPsiPi[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_DDstar_EtacRho = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_DDstar_EtacRho[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){ 
                        corr_DDstar_EtacRho[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_DDstar_EtacRho[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_DDstar_EtacRho[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<2; ++n)
                                                corr_DDstar_EtacRho[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }       
                        }
                }
        }

        Eigen::VectorXcd***** corr_DstarDstar_JPsiPi = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_DstarDstar_JPsiPi[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){ 
                        corr_DstarDstar_JPsiPi[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_DstarDstar_JPsiPi[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_DstarDstar_JPsiPi[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_DstarDstar_JPsiPi[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }       
                        }
                }
        }

        Eigen::VectorXcd***** corr_DstarDstar_EtacRho = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_DstarDstar_EtacRho[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_DstarDstar_EtacRho[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_DstarDstar_EtacRho[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_DstarDstar_EtacRho[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_DstarDstar_EtacRho[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_JPsiPi_DDstar = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_JPsiPi_DDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_JPsiPi_DDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_JPsiPi_DDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_JPsiPi_DDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<2; ++n)
                                                corr_JPsiPi_DDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_JPsiPi_DstarDstar = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_JPsiPi_DstarDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_JPsiPi_DstarDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_JPsiPi_DstarDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_JPsiPi_DstarDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_JPsiPi_DstarDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_JPsiPi_EtacRho = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_JPsiPi_EtacRho[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_JPsiPi_EtacRho[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_JPsiPi_EtacRho[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_JPsiPi_EtacRho[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_JPsiPi_EtacRho[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_EtacRho_DDstar = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_EtacRho_DDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_EtacRho_DDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_EtacRho_DDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_EtacRho_DDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<2; ++n)
                                                corr_EtacRho_DDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }

        Eigen::VectorXcd***** corr_EtacRho_DstarDstar = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_EtacRho_DstarDstar[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_EtacRho_DstarDstar[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_EtacRho_DstarDstar[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_EtacRho_DstarDstar[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_EtacRho_DstarDstar[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }  

        Eigen::VectorXcd***** corr_EtacRho_JPsiPi = new Eigen::VectorXcd****[TP0_dim];
        for(int inp1=0; inp1<TP0_dim; ++inp1){
                corr_EtacRho_JPsiPi[inp1] = new Eigen::VectorXcd***[TP0_dim];
                for(int inp2 =0; inp2<TP0_dim; ++inp2){
                        corr_EtacRho_JPsiPi[inp1][inp2] = new Eigen::VectorXcd**[np[inp1]*np[inp2]];
                        for(int i=0; i<np[inp1]*np[inp2]; ++i){
                                corr_EtacRho_JPsiPi[inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                for(int k=0; k<3; ++k){
                                        corr_EtacRho_JPsiPi[inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                        for(int n=0; n<1; ++n)
                                                corr_EtacRho_JPsiPi[inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                }
                        }
                }
        }








//allocate memory for correlation function TP1	
	Eigen::VectorXcd****** corr_DDstar_DDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
	for(int itp=0; itp<np[1]/2; ++itp){
		corr_DDstar_DDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim]; 
		for(int inp1=0; inp1<TP1_dim; ++inp1){
			corr_DDstar_DDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
			for(int inp2 =0; inp2<TP1_dim; ++inp2){
				corr_DDstar_DDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
				for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
					corr_DDstar_DDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DDstar_DDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
						for(int n=0; n<4; ++n)
							corr_DDstar_DDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DstarDstar_DstarDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
	for(int itp=0; itp<np[1]/2; ++itp){
		corr_DstarDstar_DstarDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim]; 
		for(int inp1=0; inp1<TP1_dim; ++inp1){
			corr_DstarDstar_DstarDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
			for(int inp2 =0; inp2<TP1_dim; ++inp2){
				corr_DstarDstar_DstarDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
				for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
					corr_DstarDstar_DstarDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DstarDstar_DstarDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
						for(int n=0; n<1; ++n)
							corr_DstarDstar_DstarDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DDstar_DstarDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
	for(int itp=0; itp<np[1]/2; ++itp){
		corr_DDstar_DstarDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim]; 
		for(int inp1=0; inp1<TP1_dim; ++inp1){
			corr_DDstar_DstarDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
			for(int inp2 =0; inp2<TP1_dim; ++inp2){
				corr_DDstar_DstarDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
				for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
					corr_DDstar_DstarDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DDstar_DstarDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
						for(int n=0; n<2; ++n)
							corr_DDstar_DstarDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DstarDstar_DDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
	for(int itp=0; itp<np[1]/2; ++itp){
		corr_DstarDstar_DDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim]; 
		for(int inp1=0; inp1<TP1_dim; ++inp1){
			corr_DstarDstar_DDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
			for(int inp2 =0; inp2<TP1_dim; ++inp2){
				corr_DstarDstar_DDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
				for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
					corr_DstarDstar_DDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DstarDstar_DDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
						for(int n=0; n<2; ++n)
							corr_DstarDstar_DDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

       
        Eigen::VectorXcd****** corr_JPsiPi_JPsiPi_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_JPsiPi_JPsiPi_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_JPsiPi_JPsiPi_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_JPsiPi_JPsiPi_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_JPsiPi_JPsiPi_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_JPsiPi_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_JPsiPi_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_EtacRho_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_EtacRho_EtacRho_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_EtacRho_EtacRho_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_EtacRho_EtacRho_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_EtacRho_EtacRho_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_EtacRho_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_EtacRho_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }
   

        Eigen::VectorXcd****** corr_DDstar_JPsiPi_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_DDstar_JPsiPi_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_DDstar_JPsiPi_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_DDstar_JPsiPi_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_DDstar_JPsiPi_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DDstar_JPsiPi_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<2; ++n)
                                                        corr_DDstar_JPsiPi_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_DDstar_EtacRho_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_DDstar_EtacRho_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_DDstar_EtacRho_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_DDstar_EtacRho_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_DDstar_EtacRho_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DDstar_EtacRho_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<2; ++n)
                                                        corr_DDstar_EtacRho_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }               
                        }                       
                }                                       
        }

 
        Eigen::VectorXcd****** corr_DstarDstar_JPsiPi_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_DstarDstar_JPsiPi_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_DstarDstar_JPsiPi_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_DstarDstar_JPsiPi_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_DstarDstar_JPsiPi_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DstarDstar_JPsiPi_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_DstarDstar_JPsiPi_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_DstarDstar_EtacRho_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_DstarDstar_EtacRho_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_DstarDstar_EtacRho_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_DstarDstar_EtacRho_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_DstarDstar_EtacRho_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DstarDstar_EtacRho_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_DstarDstar_EtacRho_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_JPsiPi_DDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_JPsiPi_DDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_JPsiPi_DDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_JPsiPi_DDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_JPsiPi_DDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_DDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<2; ++n)
                                                        corr_JPsiPi_DDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_DDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_EtacRho_DDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_EtacRho_DDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_EtacRho_DDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_EtacRho_DDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_DDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<2; ++n)
                                                        corr_EtacRho_DDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }
       

        Eigen::VectorXcd****** corr_JPsiPi_DstarDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_JPsiPi_DstarDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_JPsiPi_DstarDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){ 
                                corr_JPsiPi_DstarDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){ 
                                        corr_JPsiPi_DstarDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_DstarDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_DstarDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }               
                                }       
                        }       
                }       
        }       
        
        
        Eigen::VectorXcd****** corr_EtacRho_DstarDstar_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_EtacRho_DstarDstar_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_EtacRho_DstarDstar_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_EtacRho_DstarDstar_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_EtacRho_DstarDstar_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_DstarDstar_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_DstarDstar_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }               
                                }       
                        }
                }
        }
        

        Eigen::VectorXcd****** corr_JPsiPi_EtacRho_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_JPsiPi_EtacRho_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_JPsiPi_EtacRho_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_JPsiPi_EtacRho_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_JPsiPi_EtacRho_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_EtacRho_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_EtacRho_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_JPsiPi_TP1 = new Eigen::VectorXcd*****[np[1]/2];
        for(int itp=0; itp<np[1]/2; ++itp){
                corr_EtacRho_JPsiPi_TP1[itp] = new Eigen::VectorXcd****[TP1_dim];
                for(int inp1=0; inp1<TP1_dim; ++inp1){
                        corr_EtacRho_JPsiPi_TP1[itp][inp1] = new Eigen::VectorXcd***[TP1_dim];
                        for(int inp2 =0; inp2<TP1_dim; ++inp2){
                                corr_EtacRho_JPsiPi_TP1[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP1[inp1]*n_TP1[inp2]];
                                for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i){
                                        corr_EtacRho_JPsiPi_TP1[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_JPsiPi_TP1[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_JPsiPi_TP1[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }
 
	
//allocate memory for correlation function TP2	
	Eigen::VectorXcd****** corr_DDstar_DDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
	for(int itp=0; itp<np[2]/2; ++itp){
		corr_DDstar_DDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim]; 
		for(int inp1=0; inp1<TP2_dim; ++inp1){
			corr_DDstar_DDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
			for(int inp2 =0; inp2<TP2_dim; ++inp2){
				corr_DDstar_DDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
				for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
					corr_DDstar_DDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DDstar_DDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
						for(int n=0; n<4; ++n)
							corr_DDstar_DDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DstarDstar_DstarDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
	for(int itp=0; itp<np[2]/2; ++itp){
		corr_DstarDstar_DstarDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim]; 
		for(int inp1=0; inp1<TP2_dim; ++inp1){
			corr_DstarDstar_DstarDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
			for(int inp2 =0; inp2<TP2_dim; ++inp2){
				corr_DstarDstar_DstarDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
				for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
					corr_DstarDstar_DstarDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DstarDstar_DstarDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
						for(int n=0; n<1; ++n)
							corr_DstarDstar_DstarDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DDstar_DstarDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
	for(int itp=0; itp<np[2]/2; ++itp){
		corr_DDstar_DstarDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim]; 
		for(int inp1=0; inp1<TP2_dim; ++inp1){
			corr_DDstar_DstarDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
			for(int inp2 =0; inp2<TP2_dim; ++inp2){
				corr_DDstar_DstarDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
				for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
					corr_DDstar_DstarDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DDstar_DstarDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
						for(int n=0; n<2; ++n)
							corr_DDstar_DstarDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DstarDstar_DDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
	for(int itp=0; itp<np[2]/2; ++itp){
		corr_DstarDstar_DDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim]; 
		for(int inp1=0; inp1<TP2_dim; ++inp1){
			corr_DstarDstar_DDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
			for(int inp2 =0; inp2<TP2_dim; ++inp2){
				corr_DstarDstar_DDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
				for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
					corr_DstarDstar_DDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DstarDstar_DDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[2];
						for(int n=0; n<2; ++n)
							corr_DstarDstar_DDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

        Eigen::VectorXcd****** corr_JPsiPi_JPsiPi_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_JPsiPi_JPsiPi_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_JPsiPi_JPsiPi_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_JPsiPi_JPsiPi_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_JPsiPi_JPsiPi_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_JPsiPi_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_JPsiPi_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

        Eigen::VectorXcd****** corr_EtacRho_EtacRho_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_EtacRho_EtacRho_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_EtacRho_EtacRho_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_EtacRho_EtacRho_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_EtacRho_EtacRho_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_EtacRho_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_EtacRho_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

 
     
        Eigen::VectorXcd****** corr_DDstar_JPsiPi_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_DDstar_JPsiPi_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_DDstar_JPsiPi_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_DDstar_JPsiPi_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_DDstar_JPsiPi_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DDstar_JPsiPi_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_DDstar_JPsiPi_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }
  

        Eigen::VectorXcd****** corr_JPsiPi_DDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_JPsiPi_DDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_JPsiPi_DDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_JPsiPi_DDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_JPsiPi_DDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_DDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_JPsiPi_DDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }
        


        Eigen::VectorXcd****** corr_DstarDstar_JPsiPi_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_DstarDstar_JPsiPi_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_DstarDstar_JPsiPi_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_DstarDstar_JPsiPi_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_DstarDstar_JPsiPi_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DstarDstar_JPsiPi_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
                                                for(int n=0; n<1; ++n)
                                                        corr_DstarDstar_JPsiPi_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

        Eigen::VectorXcd****** corr_JPsiPi_DstarDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_JPsiPi_DstarDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_JPsiPi_DstarDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_JPsiPi_DstarDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_JPsiPi_DstarDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_DstarDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_DstarDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

        Eigen::VectorXcd****** corr_DDstar_EtacRho_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_DDstar_EtacRho_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_DDstar_EtacRho_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_DDstar_EtacRho_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_DDstar_EtacRho_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DDstar_EtacRho_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_DDstar_EtacRho_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_DDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_EtacRho_DDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_EtacRho_DDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_EtacRho_DDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_EtacRho_DDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_DDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_EtacRho_DDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }       
                                }               
                        }                               
                }
        }


       Eigen::VectorXcd****** corr_DstarDstar_EtacRho_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_DstarDstar_EtacRho_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_DstarDstar_EtacRho_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_DstarDstar_EtacRho_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_DstarDstar_EtacRho_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DstarDstar_EtacRho_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_DstarDstar_EtacRho_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_DstarDstar_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_EtacRho_DstarDstar_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_EtacRho_DstarDstar_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_EtacRho_DstarDstar_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_EtacRho_DstarDstar_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_DstarDstar_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_DstarDstar_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }



        Eigen::VectorXcd****** corr_EtacRho_JPsiPi_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_EtacRho_JPsiPi_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_EtacRho_JPsiPi_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_EtacRho_JPsiPi_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_EtacRho_JPsiPi_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_JPsiPi_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_JPsiPi_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

        Eigen::VectorXcd****** corr_JPsiPi_EtacRho_TP2 = new Eigen::VectorXcd*****[np[2]/2];
        for(int itp=0; itp<np[2]/2; ++itp){
                corr_JPsiPi_EtacRho_TP2[itp] = new Eigen::VectorXcd****[TP2_dim];
                for(int inp1=0; inp1<TP2_dim; ++inp1){
                        corr_JPsiPi_EtacRho_TP2[itp][inp1] = new Eigen::VectorXcd***[TP2_dim];
                        for(int inp2 =0; inp2<TP2_dim; ++inp2){
                                corr_JPsiPi_EtacRho_TP2[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP2[inp1]*n_TP2[inp2]];
                                for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i){
                                        corr_JPsiPi_EtacRho_TP2[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_EtacRho_TP2[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[1];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_EtacRho_TP2[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }               


//allocate memory for correlation function TP3	
	Eigen::VectorXcd****** corr_DDstar_DDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
	for(int itp=0; itp<np[3]/2; ++itp){
		corr_DDstar_DDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim]; 
		for(int inp1=0; inp1<TP3_dim; ++inp1){
			corr_DDstar_DDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
			for(int inp2 =0; inp2<TP3_dim; ++inp2){
				corr_DDstar_DDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
				for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
					corr_DDstar_DDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DDstar_DDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
						for(int n=0; n<4; ++n)
							corr_DDstar_DDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DstarDstar_DstarDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
	for(int itp=0; itp<np[3]/2; ++itp){
		corr_DstarDstar_DstarDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim]; 
		for(int inp1=0; inp1<TP3_dim; ++inp1){
			corr_DstarDstar_DstarDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
			for(int inp2 =0; inp2<TP3_dim; ++inp2){
				corr_DstarDstar_DstarDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
				for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
					corr_DstarDstar_DstarDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DstarDstar_DstarDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
						for(int n=0; n<4; ++n)
							corr_DstarDstar_DstarDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DDstar_DstarDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
	for(int itp=0; itp<np[3]/2; ++itp){
		corr_DDstar_DstarDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim]; 
		for(int inp1=0; inp1<TP3_dim; ++inp1){
			corr_DDstar_DstarDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
			for(int inp2 =0; inp2<TP3_dim; ++inp2){
				corr_DDstar_DstarDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
				for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
					corr_DDstar_DstarDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DDstar_DstarDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
						for(int n=0; n<4; ++n)
							corr_DDstar_DstarDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}

	Eigen::VectorXcd****** corr_DstarDstar_DDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
	for(int itp=0; itp<np[3]/2; ++itp){
		corr_DstarDstar_DDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim]; 
		for(int inp1=0; inp1<TP3_dim; ++inp1){
			corr_DstarDstar_DDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
			for(int inp2 =0; inp2<TP3_dim; ++inp2){
				corr_DstarDstar_DDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
				for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
					corr_DstarDstar_DDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
					for(int k=0; k<3; ++k){
						corr_DstarDstar_DDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
						for(int n=0; n<4; ++n)
							corr_DstarDstar_DDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
					}
				}	
			}
		}
	}  

        Eigen::VectorXcd****** corr_DDstar_JPsiPi_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_DDstar_JPsiPi_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_DDstar_JPsiPi_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_DDstar_JPsiPi_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_DDstar_JPsiPi_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DDstar_JPsiPi_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_DDstar_JPsiPi_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

        Eigen::VectorXcd****** corr_DstarDstar_JPsiPi_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_DstarDstar_JPsiPi_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_DstarDstar_JPsiPi_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_DstarDstar_JPsiPi_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_DstarDstar_JPsiPi_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DstarDstar_JPsiPi_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_DstarDstar_JPsiPi_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }
 

        Eigen::VectorXcd****** corr_DDstar_EtacRho_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_DDstar_EtacRho_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_DDstar_EtacRho_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_DDstar_EtacRho_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_DDstar_EtacRho_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DDstar_EtacRho_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_DDstar_EtacRho_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }               
                        
        Eigen::VectorXcd****** corr_DstarDstar_EtacRho_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_DstarDstar_EtacRho_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_DstarDstar_EtacRho_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_DstarDstar_EtacRho_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_DstarDstar_EtacRho_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_DstarDstar_EtacRho_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_DstarDstar_EtacRho_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_EtacRho_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_EtacRho_EtacRho_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_EtacRho_EtacRho_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_EtacRho_EtacRho_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_EtacRho_EtacRho_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_EtacRho_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_EtacRho_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_JPsiPi_JPsiPi_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_JPsiPi_JPsiPi_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_JPsiPi_JPsiPi_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_JPsiPi_JPsiPi_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_JPsiPi_JPsiPi_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_JPsiPi_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_JPsiPi_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }



        Eigen::VectorXcd****** corr_JPsiPi_EtacRho_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_JPsiPi_EtacRho_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_JPsiPi_EtacRho_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_JPsiPi_EtacRho_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_JPsiPi_EtacRho_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_EtacRho_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_EtacRho_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_JPsiPi_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_EtacRho_JPsiPi_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_EtacRho_JPsiPi_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_EtacRho_JPsiPi_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_EtacRho_JPsiPi_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_JPsiPi_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_JPsiPi_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_JPsiPi_DDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_JPsiPi_DDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_JPsiPi_DDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_JPsiPi_DDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_JPsiPi_DDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_DDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_JPsiPi_DDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }


        Eigen::VectorXcd****** corr_EtacRho_DDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_EtacRho_DDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_EtacRho_DDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_EtacRho_DDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_EtacRho_DDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_DDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<2; ++n)
                                                        corr_EtacRho_DDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }
                        }
                }
        }

        Eigen::VectorXcd****** corr_JPsiPi_DstarDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_JPsiPi_DstarDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_JPsiPi_DstarDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_JPsiPi_DstarDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_JPsiPi_DstarDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_JPsiPi_DstarDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_JPsiPi_DstarDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }       
                        }               
                }                               
        }                                       
                                                        
                                        
        Eigen::VectorXcd****** corr_EtacRho_DstarDstar_TP3 = new Eigen::VectorXcd*****[np[3]/2];
        for(int itp=0; itp<np[3]/2; ++itp){
                corr_EtacRho_DstarDstar_TP3[itp] = new Eigen::VectorXcd****[TP3_dim];
                for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_EtacRho_DstarDstar_TP3[itp][inp1] = new Eigen::VectorXcd***[TP3_dim];
                        for(int inp2 =0; inp2<TP3_dim; ++inp2){
                                corr_EtacRho_DstarDstar_TP3[itp][inp1][inp2] = new Eigen::VectorXcd**[n_TP3[inp1]*n_TP3[inp2]];
                                for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i){
                                        corr_EtacRho_DstarDstar_TP3[itp][inp1][inp2][i] = new Eigen::VectorXcd*[3];
                                        for(int k=0; k<3; ++k){
                                                corr_EtacRho_DstarDstar_TP3[itp][inp1][inp2][i][k] = new Eigen::VectorXcd[4];
                                                for(int n=0; n<1; ++n)
                                                        corr_EtacRho_DstarDstar_TP3[itp][inp1][inp2][i][k][n] = Eigen::VectorXcd::Zero(Lt);
                                        }
                                }       
                        }               
                }                               
        }             
i




	Eigen::VectorXcd** corr_DDstar_DDstar_sum = new Eigen::VectorXcd*[TP0_dim];
	Eigen::VectorXcd** corr_DstarDstar_DstarDstar_sum = new Eigen::VectorXcd*[TP0_dim];
	Eigen::VectorXcd** corr_DDstar_DstarDstar_sum = new Eigen::VectorXcd*[TP0_dim];
	Eigen::VectorXcd** corr_DstarDstar_DDstar_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_JPsiPi_JPsiPi_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_DDstar_JPsiPi_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_DstarDstar_JPsiPi_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_JPsiPi_DDstar_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_JPsiPi_DstarDstar_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_EtacRho_EtacRho_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_DDstar_EtacRho_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_DstarDstar_EtacRho_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_EtacRho_DDstar_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_EtacRho_DstarDstar_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_JPsiPi_EtacRho_sum = new Eigen::VectorXcd*[TP0_dim];
        Eigen::VectorXcd** corr_EtacRho_JPsiPi_sum = new Eigen::VectorXcd*[TP0_dim];
	for(int i=0; i<TP0_dim; ++i){
		corr_DDstar_DDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
		corr_DstarDstar_DstarDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
		corr_DDstar_DstarDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
		corr_DstarDstar_DDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_JPsiPi_JPsiPi_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_DDstar_JPsiPi_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_DstarDstar_JPsiPi_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_JPsiPi_DDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_JPsiPi_DstarDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_EtacRho_EtacRho_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_DDstar_EtacRho_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_DstarDstar_EtacRho_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_EtacRho_DDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_EtacRho_DstarDstar_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_JPsiPi_EtacRho_sum[i] = new Eigen::VectorXcd[TP0_dim];
                corr_EtacRho_JPsiPi_sum[i] = new Eigen::VectorXcd[TP0_dim];
		for(int j=0; j<TP0_dim; ++j){
			corr_DDstar_DDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
			corr_DstarDstar_DstarDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
			corr_DDstar_DstarDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
			corr_DstarDstar_DDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_JPsiPi_JPsiPi_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_DDstar_JPsiPi_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_DstarDstar_JPsiPi_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_JPsiPi_DDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_JPsiPi_DstarDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_EtacRho_EtacRho_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_DDstar_EtacRho_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_DstarDstar_EtacRho_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_EtacRho_DDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_EtacRho_DstarDstar_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_JPsiPi_EtacRho_sum[i][j] = Eigen::VectorXcd::Zero(Lt);
                        corr_EtacRho_JPsiPi_sum[i][j] = Eigen::VectorXcd::Zero(Lt);

		}
	}

	Eigen::VectorXcd*** corr_DDstar_DDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
	Eigen::VectorXcd*** corr_DstarDstar_DstarDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
	Eigen::VectorXcd*** corr_DDstar_DstarDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
	Eigen::VectorXcd*** corr_DstarDstar_DDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_JPsiPi_JPsiPi_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_DDstar_JPsiPi_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_DstarDstar_JPsiPi_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_JPsiPi_DDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_JPsiPi_DstarDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_EtacRho_EtacRho_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_DDstar_EtacRho_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_DstarDstar_EtacRho_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_EtacRho_DDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_EtacRho_DstarDstar_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_JPsiPi_EtacRho_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
        Eigen::VectorXcd*** corr_EtacRho_JPsiPi_TP1_sum = new Eigen::VectorXcd**[np[1]/2];
	for(int n=0; n<np[1]/2; ++n){
		corr_DDstar_DDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
		corr_DstarDstar_DstarDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
		corr_DDstar_DstarDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
		corr_DstarDstar_DDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_JPsiPi_JPsiPi_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_DDstar_JPsiPi_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_DstarDstar_JPsiPi_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_JPsiPi_DDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_JPsiPi_DstarDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_EtacRho_EtacRho_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_DDstar_EtacRho_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_DstarDstar_EtacRho_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_EtacRho_DDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_EtacRho_DstarDstar_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_JPsiPi_EtacRho_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];
                corr_EtacRho_JPsiPi_TP1_sum[n] = new Eigen::VectorXcd*[TP1_dim];

		for(int inp1=0; inp1<TP1_dim; ++inp1){
			corr_DDstar_DDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
			corr_DstarDstar_DstarDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
			corr_DDstar_DstarDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
			corr_DstarDstar_DDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_JPsiPi_JPsiPi_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_DDstar_JPsiPi_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_DstarDstar_JPsiPi_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_JPsiPi_DDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_JPsiPi_DstarDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_EtacRho_EtacRho_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_DDstar_EtacRho_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_DstarDstar_EtacRho_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_EtacRho_DDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_EtacRho_DstarDstar_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_JPsiPi_EtacRho_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
                        corr_EtacRho_JPsiPi_TP1_sum[n][inp1] = new Eigen::VectorXcd[TP1_dim];
			for(int inp2=0; inp2<TP1_dim; ++inp2){
				corr_DDstar_DDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DstarDstar_DstarDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DDstar_DstarDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DstarDstar_DDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_JPsiPi_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_JPsiPi_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_JPsiPi_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DstarDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_EtacRho_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_EtacRho_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_EtacRho_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DstarDstar_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_EtacRho_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_JPsiPi_TP1_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
			}
		}
	}
	Eigen::VectorXcd** corr_DDstar_DDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
	Eigen::VectorXcd** corr_DstarDstar_DstarDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
	Eigen::VectorXcd** corr_DDstar_DstarDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
	Eigen::VectorXcd** corr_DstarDstar_DDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim]; 
        Eigen::VectorXcd** corr_JPsiPi_JPsiPi_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_DDstar_JPsiPi_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_DstarDstar_JPsiPi_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_JPsiPi_DDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_JPsiPi_DstarDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_EtacRho_EtacRho_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_DDstar_EtacRho_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_DstarDstar_EtacRho_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_EtacRho_DDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_EtacRho_DstarDstar_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_JPsiPi_EtacRho_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
        Eigen::VectorXcd** corr_EtacRho_JPsiPi_TP1_sumsum = new Eigen::VectorXcd*[TP1_dim];
	for(int inp1=0; inp1<TP1_dim; ++inp1){
		corr_DDstar_DDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
		corr_DstarDstar_DstarDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
		corr_DDstar_DstarDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
		corr_DstarDstar_DDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_JPsiPi_JPsiPi_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_DstarDstar_JPsiPi_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_DDstar_JPsiPi_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_JPsiPi_DDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_JPsiPi_DstarDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_EtacRho_EtacRho_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_DDstar_EtacRho_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_DstarDstar_EtacRho_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_EtacRho_DDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_EtacRho_DstarDstar_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_JPsiPi_EtacRho_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
                corr_EtacRho_JPsiPi_TP1_sumsum[inp1] = new Eigen::VectorXcd[TP1_dim];
			for(int inp2=0; inp2<TP1_dim; ++inp2){
				corr_DDstar_DDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DstarDstar_DstarDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DDstar_DstarDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DstarDstar_DDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_JPsiPi_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_JPsiPi_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_JPsiPi_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DstarDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_EtacRho_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_EtacRho_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_EtacRho_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DstarDstar_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_EtacRho_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_JPsiPi_TP1_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);       
			}
	}

	Eigen::VectorXcd*** corr_DDstar_DDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
	Eigen::VectorXcd*** corr_DstarDstar_DstarDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
	Eigen::VectorXcd*** corr_DDstar_DstarDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
	Eigen::VectorXcd*** corr_DstarDstar_DDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_JPsiPi_JPsiPi_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_DstarDstar_JPsiPi_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_DDstar_JPsiPi_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_JPsiPi_DDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_JPsiPi_DstarDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_EtacRho_EtacRho_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_DDstar_EtacRho_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_DstarDstar_EtacRho_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_EtacRho_DDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_EtacRho_DstarDstar_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_JPsiPi_EtacRho_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
        Eigen::VectorXcd*** corr_EtacRho_JPsiPi_TP2_sum = new Eigen::VectorXcd**[np[2]/2];
	for(int n=0; n<np[2]/2; ++n){
		corr_DDstar_DDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
		corr_DstarDstar_DstarDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
		corr_DDstar_DstarDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
		corr_DstarDstar_DDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_JPsiPi_JPsiPi_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_DstarDstar_JPsiPi_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_DDstar_JPsiPi_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_JPsiPi_DDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_JPsiPi_DstarDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_EtacRho_EtacRho_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_DDstar_EtacRho_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_DstarDstar_EtacRho_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_EtacRho_DDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_EtacRho_DstarDstar_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_JPsiPi_EtacRho_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
                corr_EtacRho_JPsiPi_TP2_sum[n] = new Eigen::VectorXcd*[TP2_dim];
		for(int inp1=0; inp1<TP2_dim; ++inp1){
			corr_DDstar_DDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
			corr_DstarDstar_DstarDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
			corr_DDstar_DstarDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
			corr_DstarDstar_DDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_JPsiPi_JPsiPi_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_DstarDstar_JPsiPi_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_DDstar_JPsiPi_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_JPsiPi_DDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_JPsiPi_DstarDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_EtacRho_EtacRho_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_DDstar_EtacRho_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_DstarDstar_EtacRho_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_EtacRho_DDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_EtacRho_DstarDstar_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_JPsiPi_EtacRho_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
                        corr_EtacRho_JPsiPi_TP2_sum[n][inp1] = new Eigen::VectorXcd[TP2_dim];
			for(int inp2=0; inp2<TP2_dim; ++inp2){
				corr_DDstar_DDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DstarDstar_DstarDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DDstar_DstarDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
				corr_DstarDstar_DDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_JPsiPi_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_JPsiPi_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_JPsiPi_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DstarDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_EtacRho_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_EtacRho_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_EtacRho_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DstarDstar_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_EtacRho_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_JPsiPi_TP2_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
			}
		}
	}
      
        Eigen::VectorXcd** corr_DDstar_DDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DstarDstar_DstarDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DDstar_DstarDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DstarDstar_DDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_JPsiPi_JPsiPi_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DDstar_JPsiPi_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DstarDstar_JPsiPi_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_JPsiPi_DDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_JPsiPi_DstarDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_EtacRho_EtacRho_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DDstar_EtacRho_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_DstarDstar_EtacRho_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_EtacRho_DDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_EtacRho_DstarDstar_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_JPsiPi_EtacRho_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];
        Eigen::VectorXcd** corr_EtacRho_JPsiPi_TP2_sumsum = new Eigen::VectorXcd*[TP2_dim];

	for(int inp1=0; inp1<TP2_dim; ++inp1){
                corr_DDstar_DDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DstarDstar_DstarDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DDstar_DstarDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DstarDstar_DDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_JPsiPi_JPsiPi_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DstarDstar_JPsiPi_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DDstar_JPsiPi_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_JPsiPi_DDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_JPsiPi_DstarDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_EtacRho_EtacRho_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DDstar_EtacRho_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_DstarDstar_EtacRho_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_EtacRho_DDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_EtacRho_DstarDstar_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_JPsiPi_EtacRho_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
                corr_EtacRho_JPsiPi_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP2_dim];
			for(int inp2=0; inp2<TP2_dim; ++inp2){        
                                corr_DDstar_DDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_DstarDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_DstarDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_DDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_JPsiPi_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_JPsiPi_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_JPsiPi_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DstarDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_EtacRho_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_EtacRho_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_EtacRho_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DstarDstar_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_EtacRho_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_JPsiPi_TP2_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);

			}
	}

        Eigen::VectorXcd*** corr_DDstar_DDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DstarDstar_DstarDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DDstar_DstarDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DstarDstar_DDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_JPsiPi_JPsiPi_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DstarDstar_JPsiPi_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DDstar_JPsiPi_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_JPsiPi_DDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_JPsiPi_DstarDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_EtacRho_EtacRho_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DDstar_EtacRho_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_DstarDstar_EtacRho_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_EtacRho_DDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_EtacRho_DstarDstar_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_JPsiPi_EtacRho_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
        Eigen::VectorXcd*** corr_EtacRho_JPsiPi_TP3_sum = new Eigen::VectorXcd**[np[3]/2];
	for(int n=0; n<np[3]/2; ++n){
                corr_DDstar_DDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DstarDstar_DstarDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DDstar_DstarDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DstarDstar_DDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_JPsiPi_JPsiPi_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DstarDstar_JPsiPi_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DDstar_JPsiPi_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_JPsiPi_DDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_JPsiPi_DstarDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_EtacRho_EtacRho_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DDstar_EtacRho_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_DstarDstar_EtacRho_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_EtacRho_DDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_EtacRho_DstarDstar_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_JPsiPi_EtacRho_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
                corr_EtacRho_JPsiPi_TP3_sum[n] = new Eigen::VectorXcd*[TP3_dim];
		for(int inp1=0; inp1<TP3_dim; ++inp1){
                        corr_DDstar_DDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DstarDstar_DstarDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DDstar_DstarDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DstarDstar_DDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_JPsiPi_JPsiPi_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DstarDstar_JPsiPi_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DDstar_JPsiPi_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_JPsiPi_DDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_JPsiPi_DstarDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_EtacRho_EtacRho_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DDstar_EtacRho_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_DstarDstar_EtacRho_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_EtacRho_DDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_EtacRho_DstarDstar_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_JPsiPi_EtacRho_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
                        corr_EtacRho_JPsiPi_TP3_sum[n][inp1] = new Eigen::VectorXcd[TP3_dim];
			for(int inp2=0; inp2<TP3_dim; ++inp2){
                                corr_DDstar_DDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_DstarDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_DstarDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_DDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_JPsiPi_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_JPsiPi_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_JPsiPi_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DstarDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_EtacRho_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_EtacRho_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_EtacRho_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DstarDstar_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_EtacRho_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_JPsiPi_TP3_sum[n][inp1][inp2] = Eigen::VectorXcd::Zero(Lt);

			}
		}
	}
        Eigen::VectorXcd** corr_DDstar_DDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DstarDstar_DstarDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DDstar_DstarDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DstarDstar_DDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_JPsiPi_JPsiPi_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DDstar_JPsiPi_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DstarDstar_JPsiPi_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_JPsiPi_DDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_JPsiPi_DstarDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_EtacRho_EtacRho_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DDstar_EtacRho_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_DstarDstar_EtacRho_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_EtacRho_DDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_EtacRho_DstarDstar_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_JPsiPi_EtacRho_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
        Eigen::VectorXcd** corr_EtacRho_JPsiPi_TP3_sumsum = new Eigen::VectorXcd*[TP3_dim];
	for(int inp1=0; inp1<TP3_dim; ++inp1){
                corr_DDstar_DDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DstarDstar_DstarDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DDstar_DstarDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DstarDstar_DDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_JPsiPi_JPsiPi_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DstarDstar_JPsiPi_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DDstar_JPsiPi_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_JPsiPi_DDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_JPsiPi_DstarDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_EtacRho_EtacRho_TP2_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DDstar_EtacRho_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_DstarDstar_EtacRho_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_EtacRho_DDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_EtacRho_DstarDstar_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_JPsiPi_EtacRho_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
                corr_EtacRho_JPsiPi_TP3_sumsum[inp1] = new Eigen::VectorXcd[TP3_dim];
			for(int inp2=0; inp2<TP3_dim; ++inp2){
                                corr_DDstar_DDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_DstarDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_DstarDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_DDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_JPsiPi_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_JPsiPi_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_JPsiPi_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_DstarDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_EtacRho_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DDstar_EtacRho_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_DstarDstar_EtacRho_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_DstarDstar_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_JPsiPi_EtacRho_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
                                corr_EtacRho_JPsiPi_TP3_sumsum[inp1][inp2] = Eigen::VectorXcd::Zero(Lt);
			}
	}

	Eigen::VectorXcd** corr_dbarcp = new Eigen::VectorXcd*[number_of_mom];	
	Eigen::VectorXcd** corr_cmbaru = new Eigen::VectorXcd*[number_of_mom];
        Eigen::VectorXcd** corr_cmbarcp = new Eigen::VectorXcd*[number_of_mom];
        Eigen::VectorXcd** corr_dbaru = new Eigen::VectorXcd*[number_of_mom];	
	for(int i=0; i<number_of_mom; ++i){
		corr_dbarcp[i] = new Eigen::VectorXcd[4];	
		corr_cmbaru[i] = new Eigen::VectorXcd[4];
                corr_cmbarcp[i] = new Eigen::VectorXcd[4];
                corr_dbaru[i] = new Eigen::VectorXcd[4];	
		for(int gi=0; gi<4; ++gi){
			corr_dbarcp[i][gi] = Eigen::VectorXcd::Zero(Lt);	
			corr_cmbaru[i][gi] = Eigen::VectorXcd::Zero(Lt);
                        corr_cmbarcp[i][gi] = Eigen::VectorXcd::Zero(Lt);
                        corr_dbaru[i][gi] = Eigen::VectorXcd::Zero(Lt);
		}
	}



	int n_mom_components = 33;
//define intermediate 2pt functions
	Eigen::VectorXcd**** corr_dbarcp_gammaii = new Eigen::VectorXcd***[n_rnd];	
	Eigen::VectorXcd**** corr_dbarcp_gammai5 = new Eigen::VectorXcd***[n_rnd];	
	Eigen::VectorXcd**** corr_dbarcp_gamma5i = new Eigen::VectorXcd***[n_rnd];	
	Eigen::VectorXcd*** corr_dbarcp_gamma55 = new Eigen::VectorXcd**[n_rnd];
	Eigen::VectorXcd**** corr_cmbaru_gammaii = new Eigen::VectorXcd***[n_rnd];	
	Eigen::VectorXcd**** corr_cmbaru_gammai5 = new Eigen::VectorXcd***[n_rnd];	
	Eigen::VectorXcd**** corr_cmbaru_gamma5i = new Eigen::VectorXcd***[n_rnd];	
	Eigen::VectorXcd*** corr_cmbaru_gamma55 = new Eigen::VectorXcd**[n_rnd];

        Eigen::VectorXcd**** corr_cmbarcp_gammaii = new Eigen::VectorXcd***[n_rnd];
        Eigen::VectorXcd**** corr_cmbarcp_gammai5 = new Eigen::VectorXcd***[n_rnd];
        Eigen::VectorXcd**** corr_cmbarcp_gamma5i = new Eigen::VectorXcd***[n_rnd];
        Eigen::VectorXcd*** corr_cmbarcp_gamma55 = new Eigen::VectorXcd**[n_rnd];
        Eigen::VectorXcd**** corr_dbaru_gammaii = new Eigen::VectorXcd***[n_rnd];
        Eigen::VectorXcd**** corr_dbaru_gammai5 = new Eigen::VectorXcd***[n_rnd];
        Eigen::VectorXcd**** corr_dbaru_gamma5i = new Eigen::VectorXcd***[n_rnd];
        Eigen::VectorXcd*** corr_dbaru_gamma55 = new Eigen::VectorXcd**[n_rnd];

	for(int inr=0; inr<n_rnd; ++inr){
		corr_dbarcp_gammaii[inr] = new Eigen::VectorXcd**[9];
		corr_dbarcp_gammai5[inr] = new Eigen::VectorXcd**[3];
		corr_dbarcp_gamma5i[inr] = new Eigen::VectorXcd**[3];
		corr_dbarcp_gamma55[inr] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
		corr_cmbaru_gammaii[inr] = new Eigen::VectorXcd**[9];
		corr_cmbaru_gammai5[inr] = new Eigen::VectorXcd**[3];
		corr_cmbaru_gamma5i[inr] = new Eigen::VectorXcd**[3];
		corr_cmbaru_gamma55[inr] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];

                corr_cmbarcp_gammaii[inr] = new Eigen::VectorXcd**[9];
                corr_cmbarcp_gammai5[inr] = new Eigen::VectorXcd**[3];
                corr_cmbarcp_gamma5i[inr] = new Eigen::VectorXcd**[3];
                corr_cmbarcp_gamma55[inr] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
                corr_dbaru_gammaii[inr] = new Eigen::VectorXcd**[9];
                corr_dbaru_gammai5[inr] = new Eigen::VectorXcd**[3];
                corr_dbaru_gamma5i[inr] = new Eigen::VectorXcd**[3];
                corr_dbaru_gamma55[inr] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
			
			for(int gi=0; gi<3; ++gi){
				corr_dbarcp_gammai5[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
				corr_dbarcp_gamma5i[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
				corr_cmbaru_gammai5[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
				corr_cmbaru_gamma5i[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];

                                corr_cmbarcp_gammai5[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
                                corr_cmbarcp_gamma5i[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
                                corr_dbaru_gammai5[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
                                corr_dbaru_gamma5i[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];

				
					for(int inp1=0; inp1<number_of_mom; ++inp1)
					for(int inp2=0; inp2<number_of_mom; ++inp2){
						corr_dbarcp_gammai5[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
						corr_dbarcp_gamma5i[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
						corr_cmbaru_gammai5[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
						corr_cmbaru_gamma5i[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
			
                                                corr_cmbarcp_gammai5[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
                                                corr_cmbarcp_gamma5i[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
                                                corr_dbaru_gammai5[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
                                                corr_dbaru_gamma5i[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
	
						for(int i=0; i<np[inp1]*np[inp2]; ++i){
							corr_dbarcp_gammai5[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
							corr_dbarcp_gamma5i[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
							corr_cmbaru_gammai5[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
							corr_cmbaru_gamma5i[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);

                                                        corr_cmbarcp_gammai5[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
                                                        corr_cmbarcp_gamma5i[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
                                                        corr_dbaru_gammai5[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
                                                        corr_dbaru_gamma5i[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
							}
					}
			}

			for(int gi=0; gi<9; ++gi){
				corr_dbarcp_gammaii[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
				corr_cmbaru_gammaii[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
			
                                corr_cmbarcp_gammaii[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];
                                corr_dbaru_gammaii[inr][gi] = new Eigen::VectorXcd*[number_of_mom*number_of_mom];	
					for(int inp1=0; inp1<number_of_mom; ++inp1)
					for(int inp2=0; inp2<number_of_mom; ++inp2){
						corr_dbarcp_gammaii[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
						corr_cmbaru_gammaii[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
				
                                                corr_cmbarcp_gammaii[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
                                                corr_dbaru_gammaii[inr][gi][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
						for(int i=0; i<np[inp1]*np[inp2]; ++i){
							corr_dbarcp_gammaii[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
							corr_cmbaru_gammaii[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);

                                                        corr_cmbarcp_gammaii[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
                                                        corr_dbaru_gammaii[inr][gi][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
							}
					}
			}

		for(int inp1=0; inp1<number_of_mom; ++inp1)
		for(int inp2=0; inp2<number_of_mom; ++inp2){
			corr_dbarcp_gamma55[inr][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
			corr_cmbaru_gamma55[inr][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];

                        corr_cmbarcp_gamma55[inr][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
                        corr_dbaru_gamma55[inr][inp1*number_of_mom+inp2] = new Eigen::VectorXcd[np[inp1]*np[inp2]];
				for(int i=0; i<np[inp1]*np[inp2]; ++i){
					corr_dbarcp_gamma55[inr][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
					corr_cmbaru_gamma55[inr][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);

                                        corr_cmbarcp_gamma55[inr][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
                                        corr_dbaru_gamma55[inr][inp1*number_of_mom+inp2][i] = Eigen::VectorXcd::Zero(Lt);
				}
		}	
	}

	
std::cout<<"end correlator memory allocate"<<std::endl;

	int gammaindex[4] = {12,13,14,5};
	int i_m, j_m;
	int norm_2pt_1=0, norm_2pt_2=0, norm_4pt=0;
	time = clock() - time;
   printf("\t\t pre setting done - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
	time = clock();	

	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1)
		norm_2pt_1++;

	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1){
		if(i0==i1) continue;
		norm_2pt_2++;
	}
		
	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1)
	for(int i2=0; i2<n_rndvec_c; ++i2)
	for(int i3=0; i3<n_rndvec_u; ++i3){
		if((i0==i3) || (i1==i2))  continue;
		norm_4pt++;
	}

 int vec_i[3][2] = {{1,2},{2,0},{0,1}};
 const std::complex<double> ii(0.0, 1.0);
/******************start config loop ********************/
  for(int config_i = start_config; config_i <= end_config; config_i += delta_config){

      std::cout << "\nprozessing configuration: " << config_i << "\n\n";
  
  ReadWrite* rewr = new ReadWrite;

      rewr->read_perambulators_from_file(config_i);
      rewr->read_rnd_vectors_from_file(config_i);
//      rewr->read_eigenvectors_from_file(config_i);
	
	time = clock() - time;
   printf("\t\t read in perambulators and random vectors success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
	time = clock();	

//calculate source mesonfunc


	for(int t=0; t<Lt; t++){
	rewr->read_eigenvectors_from_file(config_i, t);
	for(int inp=1; inp<number_of_mom; ++inp)
        for(int i=0; i<np[inp]/2; ++i){
                build_VdaggerV(VdaggerV[inp][i*2][t], &(rewr->V), p[inp][i*2][0], p[inp][i*2][1], p[inp][i*2][2]);
		VdaggerV[inp][i*2+1][t]=VdaggerV[inp][i*2][t].adjoint();
	}
	}

	time = clock() - time;
   printf("\t\t build VdaggerV success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
	time = clock();
	

	for(int t=0; t<Lt; ++t){
         for(int i0=0; i0<n_rndvec_u; ++i0)
            for(int i1=0; i1<n_rndvec_c; ++i1){
		for(int gi=0; gi<4; ++gi){
			MesonFunc_SourcebarSource(t, SourcebarSource_dbarcp[i0*n_rndvec_c+i1][0][0][gi][t], &(rewr->rnd_vec_u[i0]), &(rewr->rnd_vec_c[i1]), gammaindex[gi]);
			MesonFunc_SourcebarSource(t, SourcebarSource_cmbaru[i1*n_rndvec_u+i0][0][0][gi][t], &(rewr->rnd_vec_c[i1]), &(rewr->rnd_vec_u[i0]), gammaindex[gi]);

			for(int inp=1; inp<number_of_mom; ++inp)
			for(int i=0; i<np[inp]; ++i){
				MesonFunc_SourcebarSource(t, SourcebarSource_dbarcp[i0*n_rndvec_c+i1][inp][i][gi][t], &(rewr->rnd_vec_u[i0]), &(rewr->rnd_vec_c[i1]), &(VdaggerV[inp][i][t]), gammaindex[gi]);
				MesonFunc_SourcebarSource(t, SourcebarSource_cmbaru[i1*n_rndvec_u+i0][inp][i][gi][t], &(rewr->rnd_vec_c[i1]), &(rewr->rnd_vec_u[i0]), &(VdaggerV[inp][i][t]), gammaindex[gi]);
			}
		}
	}
	 
         for(int i0=0; i0<n_rndvec_c; ++i0)
            for(int i1=0; i1<n_rndvec_c; ++i1){
                for(int gi=0; gi<4; ++gi){
                        MesonFunc_SourcebarSource(t, SourcebarSource_cmbarcp[i0*n_rndvec_c+i1][0][0][gi][t], &(rewr->rnd_vec_c[i0]), &(rewr->rnd_vec_c[i1]), gammaindex[gi]);

                        for(int inp=1; inp<number_of_mom; ++inp)
                        for(int i=0; i<np[inp]; ++i){
                                MesonFunc_SourcebarSource(t, SourcebarSource_cmbarcp[i0*n_rndvec_c+i1][inp][i][gi][t], &(rewr->rnd_vec_c[i0]), &(rewr->rnd_vec_c[i1]), &(VdaggerV[inp][i][t]), gammaindex[gi]);
                        }
                }
        }  

         for(int i0=0; i0<n_rndvec_u; ++i0)
            for(int i1=0; i1<n_rndvec_u; ++i1){
                for(int gi=0; gi<4; ++gi){
                        MesonFunc_SourcebarSource(t, SourcebarSource_dbaru[i0*n_rndvec_u+i1][0][0][gi][t], &(rewr->rnd_vec_u[i0]), &(rewr->rnd_vec_u[i1]), gammaindex[gi]);
                        
                        for(int inp=1; inp<number_of_mom; ++inp)
                        for(int i=0; i<np[inp]; ++i){
                                MesonFunc_SourcebarSource(t, SourcebarSource_dbaru[i0*n_rndvec_u+i1][inp][i][gi][t], &(rewr->rnd_vec_u[i0]), &(rewr->rnd_vec_u[i1]), &(VdaggerV[inp][i][t]), gammaindex[gi]);               
                        }       
                }
        }     
	
}
	time = clock() - time;
   printf("\t\t sourcebarsource mesonfunction success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);
	time = clock();



	for(int t_sink=0; t_sink<Lt; ++t_sink){
//	for(int t_sink=0; t_sink<1; ++t_sink){

// start calculate Sinkbar Sink mesonfunction  

 
        for(int i0=0; i0<n_rndvec_c; ++i0)        
               for(int i1=0; i1<n_rndvec_c; ++i1)                                                                                                                      for(int gi=0; gi<4; ++gi){
                           MesonFunc_SinkbarSink_diag(t_sink, SinkbarSink_cmbarcp[i0*n_rndvec_c+i1][0][0][gi], &(rewr->perambulator_c[i0]), &(rewr->perambulator_c[i1]), gammaindex[gi]);
                           for(int inp=1; inp<number_of_mom; ++inp)
                                for(int i=0; i<np[inp]; ++i){
                                     MesonFunc_SinkbarSink_diag(t_sink, SinkbarSink_cmbaru[i0*n_rndvec_c+i1][inp][i][gi], &(rewr->perambulator_c[i0]), &(rewr->perambulator_c[i1]), &(VdaggerV[inp][i][t_sink]), gammaindex[gi]);
                                        }
                }
 	
	for(int i0=0; i0<n_rndvec_u; ++i0)
               for(int i1=0; i1<n_rndvec_u; ++i1)
                       for(int gi=0; gi<4; ++gi){
                           MesonFunc_SinkbarSink_diag(t_sink, SinkbarSink_dbaru[i0*n_rndvec_u+i1][0][0][gi], &(rewr->perambulator_u[i0]), &(rewr->perambulator_u[i1]), gammaindex[gi]);
                           for(int inp=1; inp<number_of_mom; ++inp)
                                for(int i=0; i<np[inp]; ++i){
                                     MesonFunc_SinkbarSink_diag(t_sink, SinkbarSink_dbaru[i0*n_rndvec_u+i1][inp][i][gi], &(rewr->perambulator_u[i0]), &(rewr->perambulator_u[i1]), &(VdaggerV[inp][i][t_sink]), gammaindex[gi]);
                                        }
                }
 
  
	for(int i0=0; i0<n_rndvec_u; ++i0)        
		 for(int i1=0; i1<n_rndvec_c; ++i1) 
                       for(int gi=0; gi<4; ++gi){
                            MesonFunc_SinkbarSink_diag(t_sink, SinkbarSink_cmbaru[i1*n_rndvec_u+i0][0][0][gi], &(rewr->perambulator_c[i1]), &(rewr->perambulator_u[i0]), gammaindex[gi]);
                           for(int inp=1; inp<number_of_mom; ++inp)
                                for(int i=0; i<np[inp]; ++i){
                                     MesonFunc_SinkbarSink_diag(t_sink, SinkbarSink_cmbaru[i1*n_rndvec_u+i0][inp][i][gi], &(rewr->perambulator_c[i1]), &(rewr->perambulator_u[i0]), &(VdaggerV[inp][i][t_sink]), gammaindex[gi]);
                                        }
                }
        
 
	 for(int i0=0; i0<n_rndvec_u; ++i0)
            for(int i1=0; i1<n_rndvec_c; ++i1){
                for(int inp=0; inp<number_of_mom; ++inp){
                for(int i=0; i<np[inp]; ++i){
                        if(inp==0)
                                i_m = i;
                        else
                                i_m = (i%2==0)?(i+1):(i-1);
                        for(int gi=0; gi<4; ++gi){
                                for(int t_dil=0; t_dil<number_of_dilution_T; ++t_dil)
                                SinkbarSink_dbarcp[i0*n_rndvec_c+i1][inp][i][gi].block(0, t_dil*dim_block, dim_block, dim_block) = (SinkbarSink_cmbaru[i1*n_rndvec_u+i0][inp][i_m][gi].block(0, t_dil*dim_block, dim_block, dim_block)).adjoint();
                        }
                }
        }
        }



	time = clock() - time;
   printf("\t\t calculate sink mesonfunc at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();


//start correlation function calculation

//two-point function intermediate

	for(int i0=0; i0<n_rndvec_u; ++i0)
        for(int i1=0; i1<n_rndvec_c; ++i1)
                for(int inp1=0; inp1<number_of_mom; ++inp1)
                for(int inp2=0; inp2<number_of_mom; ++inp2){
                       	for(int i=0; i<np[inp1]; ++i){
                       	for(int j=0; j<np[inp2]; ++j){
			        corr_cmbaru_gamma55[i1*n_rndvec_u+i0][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbaru[i1*n_rndvec_u+i0][inp1][i][3]), SourcebarSource_cmbaru[i1*n_rndvec_u+i0][inp2][j][3]);
				for(int gi1=0; gi1<3; ++gi1){
					 corr_cmbaru_gammai5[i1*n_rndvec_u+i0][gi1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbaru[i1*n_rndvec_u+i0][inp1][i][gi1]), SourcebarSource_cmbaru[i1*n_rndvec_u+i0][inp2][j][3]);
                               		 corr_cmbaru_gamma5i[i1*n_rndvec_u+i0][gi1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbaru[i1*n_rndvec_u+i0][inp1][i][3]), SourcebarSource_cmbaru[i1*n_rndvec_u+i0][inp2][j][gi1]);
                    	      		 for(int gi2=0; gi2<3; ++gi2){
                                     		 corr_cmbaru_gammaii[i1*n_rndvec_u+i0][gi1*3+gi2][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbaru[i1*n_rndvec_u+i0][inp1][i][gi1]), SourcebarSource_cmbaru[i1*n_rndvec_u+i0][inp2][j][gi2]);
                                       		 }
                                        }
                        }
	         	}
        }

        for(int i0=0; i0<n_rndvec_c; ++i0)
        for(int i1=0; i1<n_rndvec_c; ++i1)
                for(int inp1=0; inp1<number_of_mom; ++inp1)
                for(int inp2=0; inp2<number_of_mom; ++inp2){
                        for(int i=0; i<np[inp1]; ++i){
                        for(int j=0; j<np[inp2]; ++j){
                                corr_cmbarcp_gamma55[i0*n_rndvec_c+i1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbarcp[i0*n_rndvec_c+i1][inp1][i][3]), SourcebarSource_cmbarcp[i0*n_rndvec_c+i1][inp2][j][3]);
                                for(int gi1=0; gi1<3; ++gi1){
                                         corr_cmbarcp_gammai5[i0*n_rndvec_c+i1][gi1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbarcp[i0*n_rndvec_c+i1][inp1][i][gi1]), SourcebarSource_cmbarcp[i0*n_rndvec_c+i1][inp2][j][3]);
                                         corr_cmbarcp_gamma5i[i0*n_rndvec_c+i1][gi1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbarcp[i0*n_rndvec_c+i1][inp1][i][3]), SourcebarSource_cmbarcp[i0*n_rndvec_c+i1][inp2][j][gi1]);                                     
				         for(int gi2=0; gi2<3; ++gi2){
                                                 corr_cmbarcp_gammaii[i0*n_rndvec_c+i1][gi1*3+gi2][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_cmbarcp[i0*n_rndvec_c+i1][inp1][i][gi1]), SourcebarSource_cmbarcp[i0*n_rndvec_c+i1][inp2][j][gi2]);                     
                                                       }    
                                    }
                        }                
                        }
      }

        for(int i0=0; i0<n_rndvec_u; ++i0)
        for(int i1=0; i1<n_rndvec_u; ++i1)
                for(int inp1=0; inp1<number_of_mom; ++inp1)
                for(int inp2=0; inp2<number_of_mom; ++inp2){
                        for(int i=0; i<np[inp1]; ++i)
                        for(int j=0; j<np[inp2]; ++j){
                                corr_dbaru_gamma55[i0*n_rndvec_u+i1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_dbaru[i0*n_rndvec_u+i1][inp1][i][3]), SourcebarSource_dbaru[i0*n_rndvec_u+i1][inp2][j][3]);
                                for(int gi1=0; gi1<3; ++gi1){
                                         corr_dbaru_gammai5[i0*n_rndvec_u+i1][gi1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_dbaru[i0*n_rndvec_u+i1][inp1][i][gi1]), SourcebarSource_dbaru[i0*n_rndvec_u+i1][inp2][j][3]);
                                         corr_dbaru_gamma5i[i0*n_rndvec_u+i1][gi1][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_dbaru[i0*n_rndvec_u+i1][inp1][i][3]), SourcebarSource_dbaru[i0*n_rndvec_u+i1][inp2][j][gi1]);
                                         for(int gi2=0; gi2<3; ++gi2){
                                                 corr_dbaru_gammaii[i0*n_rndvec_u+i1][gi1*3+gi2][inp1*number_of_mom+inp2][i*np[inp2]+j] = Contractions_TwoPoint_SiSiSoSo(t_sink, &(SinkbarSink_dbaru[i0*n_rndvec_u+i1][inp1][i][gi1]), SourcebarSource_dbaru[i0*n_rndvec_u+i1][inp2][j][gi2]);
                                                       }
                                    }
                        }
      }
	
     


	time = clock() - time;
   printf("\t\t 2pt intermediate correlation function at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();

	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1)
	for(int gi=0; gi<3; ++gi){
		for(int inp1=0; inp1<number_of_mom; ++inp1){
			for(int inp2=0; inp2<number_of_mom; ++inp2){
				for(int i=0; i<np[inp1]; ++i){
					if(np[inp1] == 1)
                              			i_m = i;
                        		else
                              			i_m = (i%2==0)?(i+1):(i-1);
					for(int j=0; j<np[inp2]; ++j){
						if(np[inp2] == 1)
							j_m=j;
						else
							j_m = (j%2==0)?(j+1):(j-1);
					
						corr_dbarcp_gammai5[i0*n_rndvec_c+i1][gi][inp1*number_of_mom+inp2][i*np[inp2]+j] = corr_cmbaru_gammai5[i1*n_rndvec_u+i0][gi][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m].conjugate(); 
						corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][gi][inp1*number_of_mom+inp2][i*np[inp2]+j] = corr_cmbaru_gamma5i[i1*n_rndvec_u+i0][gi][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m].conjugate(); 
					}
				}
			}
		}
	}	

	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1)
	for(int gi=0; gi<9; ++gi){
		for(int inp1=0; inp1<number_of_mom; ++inp1){
			for(int inp2=0; inp2<number_of_mom; ++inp2){
				for(int i=0; i<np[inp1]; ++i){
					if(np[inp1] == 1)
                              			i_m = i;
                        		else
                              			i_m = (i%2==0)?(i+1):(i-1);
					for(int j=0; j<np[inp2]; ++j){
						if(np[inp2] == 1)
							j_m=j;
						else
							j_m = (j%2==0)?(j+1):(j-1);
					
						corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi][inp1*number_of_mom+inp2][i*np[inp2]+j] = corr_cmbaru_gammaii[i1*n_rndvec_u+i0][gi][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m].conjugate(); 
					}
				}
			}
		}
	}	

	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1)
		for(int inp1=0; inp1<number_of_mom; ++inp1){
			for(int inp2=0; inp2<number_of_mom; ++inp2){
				for(int i=0; i<np[inp1]; ++i){
					if(np[inp1] == 1)
                              			i_m = i;
                        		else
                              			i_m = (i%2==0)?(i+1):(i-1);
					for(int j=0; j<np[inp2]; ++j){
						if(np[inp2] == 1)
							j_m=j;
						else
							j_m = (j%2==0)?(j+1):(j-1);
					
						corr_dbarcp_gamma55[i0*n_rndvec_c+i1][inp1*number_of_mom+inp2][i*np[inp2]+j] = corr_cmbaru_gamma55[i1*n_rndvec_u+i0][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m].conjugate(); 
					}
				}
			}
		}

	time = clock() - time;
   printf("\t\t 2pt intermediate correlation function at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();


//two-point functions
        for(int i0=0; i0<n_rndvec_u; ++i0)
        for(int i1=0; i1<n_rndvec_c; ++i1){
                for(int inp=0; inp<number_of_mom; ++inp)
                for(int i=0; i<np[inp]; ++i){

			for(int gi=0; gi<4; gi++){
				if(gi==3){
	                		corr_dbarcp[inp][gi] += corr_dbarcp_gamma55[i0*n_rndvec_c+i1][inp*number_of_mom+inp][i*np[inp]+i];
	                		corr_cmbaru[inp][gi] += corr_cmbaru_gamma55[i1*n_rndvec_u+i0][inp*number_of_mom+inp][i*np[inp]+i];
				}
				else{
	                	corr_dbarcp[inp][gi] += corr_dbarcp_gammaii[i0*n_rndvec_c+i1][4*gi][inp*number_of_mom+inp][i*np[inp]+i];
	                	corr_cmbaru[inp][gi] += corr_cmbaru_gammaii[i1*n_rndvec_u+i0][4*gi][inp*number_of_mom+inp][i*np[inp]+i];
				}
			}
                }       
	 }

        for(int i0=0; i0<n_rndvec_c; ++i0)
        for(int i1=0; i1<n_rndvec_c; ++i1){
                for(int inp=0; inp<number_of_mom; ++inp)
                for(int i=0; i<np[inp]; ++i){

                        for(int gi=0; gi<4; gi++){
                                if(gi==3){
                                        corr_cmbarcp[inp][gi] += corr_cmbarcp_gamma55[i0*n_rndvec_c+i1][inp*number_of_mom+inp][i*np[inp]+i];
                                }
                                else{
                                corr_cmbarcp[inp][gi] += corr_cmbarcp_gammaii[i0*n_rndvec_c+i1][4*gi][inp*number_of_mom+inp][i*np[inp]+i];
                                }
                        }
                }
         }

        for(int i0=0; i0<n_rndvec_u; ++i0)
        for(int i1=0; i1<n_rndvec_u; ++i1){
                for(int inp=0; inp<number_of_mom; ++inp)
                for(int i=0; i<np[inp]; ++i){

                        for(int gi=0; gi<4; gi++){
                                if(gi==3){
                                        corr_dbaru[inp][gi] += corr_dbaru_gamma55[i0*n_rndvec_u+i1][inp*number_of_mom+inp][i*np[inp]+i];
                                }
                                else{
                                corr_dbaru[inp][gi] += corr_dbaru_gammaii[i0*n_rndvec_u+i1][4*gi][inp*number_of_mom+inp][i*np[inp]+i];
                                }
                        }
                }
         }


	time = clock() - time;
   printf("\t\t 2pt correlation function at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();

//four point functions
int k1, k2, k1m, k2m, k1i, k2i, k1mi, k2mi;
	for(int i0=0; i0<n_rndvec_u; ++i0)
	for(int i1=0; i1<n_rndvec_c; ++i1)
	for(int i2=0; i2<n_rndvec_c; ++i2)
	for(int i3=0; i3<n_rndvec_u; ++i3){
		if((i0==i3) || (i1==i2)) continue;

	//---------------TP0------------------
	        for(int inp1=0; inp1<TP0_dim; ++inp1)
                for(int inp2=0; inp2<TP0_dim; ++inp2){
                        for(int i=0; i<np[inp1]; ++i){
                        if(np[inp1] == 1)
                              i_m = i;
                        else
                              i_m = (i%2==0)?(i+1):(i-1);
                        for(int j=0; j<np[inp2]; ++j){
                        if(np[inp2] == 1)
                             j_m = j;
                        else
                             j_m = (j%2==0)?(j+1):(j-1);
				for(int gi=0; gi<3; ++gi){

                                        corr_DDstar_JPsiPi[inp1][inp2][i*np[inp2]+j][gi][0] -= Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(t_sink, &(SinkbarSink_dbarcp[i0*n_rndvec_c+i1][inp1][i][gi]), &(SinkbarSink_cmbaru[i2*n_rndvec_u+i3][inp1][i_m][3]), SourcebarSource_dbaru[i0*n_rndvec_u+i3][inp2][j][3], SourcebarSource_cmbarcp[i2*n_rndvec_c+i1][inp2][j_m][gi]);
                                        corr_DDstar_JPsiPi[inp1][inp2][i*np[inp2]+j][gi][1] -= Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(t_sink, &(SinkbarSink_cmbaru[i2*n_rndvec_u+i3][inp1][i][gi]), &(SinkbarSink_dbarcp[i0*n_rndvec_c+i1][inp1][i_m][3]), SourcebarSource_cmbarcp[i2*n_rndvec_c+i1][inp2][j][gi], SourcebarSource_dbaru[i0*n_rndvec_u+i3][inp2][j_m][3]);

					corr_DDstar_EtacRho[inp1][inp2][i*np[inp2]+j][gi][0] -= Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(t_sink, &(SinkbarSink_dbarcp[i0*n_rndvec_c+i1][inp1][i][gi]), &(SinkbarSink_cmbaru[i2*n_rndvec_u+i3][inp1][i_m][gi]), SourcebarSource_dbaru[i0*n_rndvec_u+i3][inp2][j][gi], SourcebarSource_cmbarcp[i2*n_rndvec_c+i1][inp2][j_m][3]);
                                        corr_DDstar_EtacRho[inp1][inp2][i*np[inp2]+j][gi][1] -= Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(t_sink, &(SinkbarSink_cmbaru[i2*n_rndvec_u+i3][inp1][i][gi]), &(SinkbarSink_dbarcp[i0*n_rndvec_u+i1][inp1][i_m][3]), SourcebarSource_cmbarcp[i2*n_rndvec_c+i1][inp2][j][3], SourcebarSource_dbaru[i0*n_rndvec_u+i3][inp2][j_m][gi]);	

					corr_JPsiPi_DDstar[inp1][inp2][i*np[inp2]+j][gi][0] -= Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(t_sink, &(SinkbarSink_cmbarcp[i2*n_rndvec_c+i1][inp1][i][gi]), &(SinkbarSink_dbaru[i0*n_rndvec_u+i3][inp1][i_m][3]), SourcebarSource_cmbaru[i2*n_rndvec_u+i3][inp2][j][3], SourcebarSource_dbarcp[i0*n_rndvec_c+i1][inp2][j_m][gi]);
                                        corr_JPsiPi_DDstar[inp1][inp2][i*np[inp2]+j][gi][1] -= Contractions_FourPoint_SiSiSiSi_SoSoSoSo_con(t_sink, &(SinkbarSink_cmbarcp[i2*n_rndvec_u+i1][inp1][i][gi]), &(SinkbarSink_dbaru[i0*n_rndvec_u+i3][inp1][i_m][3]), SourcebarSource_cmbaru[i2*n_rndvec_u+i3][inp2][j][gi], SourcebarSource_dbarcp[i0*n_rndvec_c+i1][inp2][j_m][3]);		
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti));

					for(int ti=0; ti<Lt; ++ti){
						corr_DstarDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][0](ti) += (corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti)
                                                                                                                + corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti)
                                                                                                                - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti)
                                                                                                                - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti));

                        			corr_DDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][0](ti) += ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][1]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti) 
													     - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][0]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti)); 
                        			corr_DDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][1](ti) += ii*(corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][0]][inp1*number_of_mom+inp2][i_m*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][1]][inp1*number_of_mom+inp2][i*np[inp2]+j_m](ti) 
													     - corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][1]][inp1*number_of_mom+inp2][i_m*np[inp2]+j](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][0]][inp1*number_of_mom+inp2][i*np[inp2]+j_m](ti));
 
                        			corr_DstarDstar_DDstar[inp1][inp2][i*np[inp2]+j][gi][0](ti) -= ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+gi][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][1]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti) 
													     - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+gi][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][0]][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti)); 
                        			corr_DstarDstar_DDstar[inp1][inp2][i*np[inp2]+j][gi][1](ti) -= ii*(corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][0]][inp1*number_of_mom+inp2][i*np[inp2]+j_m](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+gi][inp1*number_of_mom+inp2][i_m*np[inp2]+j](ti) 
													     - corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][1]][inp1*number_of_mom+inp2][i*np[inp2]+j_m](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+gi][inp1*number_of_mom+inp2][i_m*np[inp2]+j](ti));

						corr_JPsiPi_JPsiPi[inp1][inp2][i*np[inp2]+j][gi][0](ti) -= corr_cmbarcp_gammaii[i2*n_rndvec_c+i1][gi][inp1*number_of_mom+inp2][i*np[inp2]+j](ti)*corr_dbaru_gamma55[i0*n_rndvec_u+i3][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti);

						corr_EtacRho_EtacPho[inp1][inp2][i*np[inp2]+j][gi][0](ti) -= corr_cmbarcp_gamma55[i2*n_rndvec_c+i1][gi][inp1*number_of_mom+ip2][i*np[inp2]+j](ti)*corr_dbaru_gammaii[i0*n_rndvec_u+i3][inp1*number_of_mom+inp2][i_m*np[inp2]+j_m](ti);

						 
					}

				}
                        }
                        }
                }// TP0 inp1, inp2 loop end here

	time = clock() - time;
   printf("\t\t 4pt correlation function TP0 at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();
	//------------TP1----------------
        for(int n=0; n<np[1]/2; ++n)
                for(int inp1=0; inp1<TP1_dim; ++inp1)
                for(int inp2=0; inp2<TP1_dim; ++inp2){
                        for(int i=0; i<n_TP1[inp1]; ++i)
                        for(int j=0; j<n_TP1[inp2]; ++j){          
				k1=TP1[n][inp1][i][0][0];
				k1m=TP1[n][inp1][i][1][0];       
				k1i=TP1[n][inp1][i][0][1];
				k1mi=TP1[n][inp1][i][1][1]; 
				k2=TP1[n][inp2][j][0][0];
				k2m=TP1[n][inp2][j][1][0];     
				k2i=TP1[n][inp2][j][0][1];
				k2mi=TP1[n][inp2][j][1][1];
 
				for(int gi=0; gi<3; ++gi){
					for(int ti=0; ti<Lt; ++ti){
                        			corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][0](ti) -= corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma55[i2*n_rndvec_u+i3][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti);
                        			corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][1](ti) -= corr_dbarcp_gammai5[i0*n_rndvec_c+i1][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti);
	time = clock() - time;
   printf("\t\t 4pt correlation function TP0 at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();
	//------------TP1----------------
        for(int n=0; n<np[1]/2; ++n)
                for(int inp1=0; inp1<TP1_dim; ++inp1)
                for(int inp2=0; inp2<TP1_dim; ++inp2){
                        for(int i=0; i<n_TP1[inp1]; ++i)
                        for(int j=0; j<n_TP1[inp2]; ++j){          
				k1=TP1[n][inp1][i][0][0];
				k1m=TP1[n][inp1][i][1][0];       
				k1i=TP1[n][inp1][i][0][1];
				k1mi=TP1[n][inp1][i][1][1]; 
				k2=TP1[n][inp2][j][0][0];
				k2m=TP1[n][inp2][j][1][0];     
				k2i=TP1[n][inp2][j][0][1];
				k2mi=TP1[n][inp2][j][1][1];
 
				for(int gi=0; gi<3; ++gi){
					for(int ti=0; ti<Lt; ++ti){
                        			corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][0](ti) -= corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma55[i2*n_rndvec_u+i3][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti);
                        			corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][1](ti) -= corr_dbarcp_gammai5[i0*n_rndvec_c+i1][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti);
                        			corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][2](ti) -= corr_dbarcp_gamma55[i0*n_rndvec_c+i1][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti);
                        			corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][3](ti) -= corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti);

                        			corr_DstarDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][0](ti) += (corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														+ corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 

                        			corr_DDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][0](ti) += ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
													     	           - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 
                        			corr_DDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][1](ti) += ii*(corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][0]][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][1]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti) 
													     		   - corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][1]][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][0]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti));
 
                        			corr_DstarDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][0](ti) -= ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
													     - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 
                        			corr_DstarDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][1](ti) -= ii*(corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][0]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti) 
													     - corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][1]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti));
					}
				}
			}
                }// inp1, inp2 loop end here
	time = clock() - time;
   printf("\t\t 4pt correlation function TP1 at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();

	//------------TP2----------------
        for(int n=0; n<np[2]/2; ++n)
                for(int inp1=0; inp1<TP2_dim; ++inp1)
                for(int inp2=0; inp2<TP2_dim; ++inp2){
                        for(int i=0; i<n_TP2[inp1]; ++i)
                        for(int j=0; j<n_TP2[inp2]; ++j){          
				k1=TP2[n][inp1][i][0][0];
				k1m=TP2[n][inp1][i][1][0];       
				k1i=TP2[n][inp1][i][0][1];
				k1mi=TP2[n][inp1][i][1][1]; 
				k2=TP2[n][inp2][j][0][0];
				k2m=TP2[n][inp2][j][1][0];     
				k2i=TP2[n][inp2][j][0][1];
				k2mi=TP2[n][inp2][j][1][1];
 
				for(int gi=0; gi<3; ++gi){
					for(int ti=0; ti<Lt; ++ti){
                        			corr_DDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][0](ti) -= corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma55[i2*n_rndvec_u+i3][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti);
                        			corr_DDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][1](ti) -= corr_dbarcp_gammai5[i0*n_rndvec_c+i1][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti);
                        			corr_DDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][2](ti) -= corr_dbarcp_gamma55[i0*n_rndvec_c+i1][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti);
                        			corr_DDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][3](ti) -= corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti);

                        			corr_DstarDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][0](ti) += (corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														+ corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 

                        			corr_DDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][0](ti) += ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
													     	           - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 
                        			corr_DDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][1](ti) += ii*(corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][0]][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][1]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti) 
													     		   - corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][1]][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][0]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti));
 
                        			corr_DstarDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][0](ti) -= ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
													     - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 
                        			corr_DstarDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][1](ti) -= ii*(corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][0]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti) 
													     - corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][1]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti));
					}
				}
			}
                }// inp1, inp2 loop end here
	time = clock() - time;
   printf("\t\t 4pt correlation function TP2 at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();


	//------------TP3----------------
        for(int n=0; n<np[3]/2; ++n)
                for(int inp1=0; inp1<TP3_dim; ++inp1)
                for(int inp2=0; inp2<TP3_dim; ++inp2){
                        for(int i=0; i<n_TP3[inp1]; ++i)
                        for(int j=0; j<n_TP3[inp2]; ++j){          
				k1=TP3[n][inp1][i][0][0];
				k1m=TP3[n][inp1][i][1][0];       
				k1i=TP3[n][inp1][i][0][1];
				k1mi=TP3[n][inp1][i][1][1]; 
				k2=TP3[n][inp2][j][0][0];
				k2m=TP3[n][inp2][j][1][0];     
				k2i=TP3[n][inp2][j][0][1];
				k2mi=TP3[n][inp2][j][1][1];
 
				for(int gi=0; gi<3; ++gi){
					for(int ti=0; ti<Lt; ++ti){
                        			corr_DDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][0](ti) -= corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma55[i2*n_rndvec_u+i3][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti);
                        			corr_DDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][1](ti) -= corr_dbarcp_gammai5[i0*n_rndvec_c+i1][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti);
                        			corr_DDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][2](ti) -= corr_dbarcp_gamma55[i0*n_rndvec_c+i1][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti);
                        			corr_DDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][3](ti) -= corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][gi][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti);

                        			corr_DstarDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][0](ti) += (corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														+ corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
														- corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 

                        			corr_DDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][0](ti) += ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][0]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
													     	           - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][gi*3+vec_i[gi][1]][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gamma5i[i2*n_rndvec_u+i3][vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 
                        			corr_DDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][1](ti) += ii*(corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][0]][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][1]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti) 
													     		   - corr_dbarcp_gamma5i[i0*n_rndvec_c+i1][vec_i[gi][1]][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][gi*3+vec_i[gi][0]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti));
 
                        			corr_DstarDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][0](ti) -= ii*(corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][0]*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][1]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti) 
													     - corr_dbarcp_gammaii[i0*n_rndvec_c+i1][vec_i[gi][1]*3+gi][k1*number_of_mom+k2][k1i*np[k2]+k2i](ti)*corr_cmbaru_gammai5[i2*n_rndvec_u+i3][vec_i[gi][0]][k1m*number_of_mom+k2m][k1mi*np[k2m]+k2mi](ti)); 
                        			corr_DstarDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][1](ti) -= ii*(corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][0]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][1]*3+gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti) 
													     - corr_dbarcp_gammai5[i0*n_rndvec_c+i1][vec_i[gi][1]][k1*number_of_mom+k2m][k1i*np[k2m]+k2mi](ti)*corr_cmbaru_gammaii[i2*n_rndvec_u+i3][vec_i[gi][0]*3+gi][k1m*number_of_mom+k2][k1mi*np[k2]+k2i](ti));
					}
				}
			}
                }// inp1, inp2 loop end here
	time = clock() - time;
   printf("\t\t 4pt correlation function TP3 at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();



	} // i0, i1, i2, i3 loop ends here


	time = clock() - time;
   printf("\t\t 4pt correlation function at t=%d success - %.1f seconds \n ", t_sink, ((float) time)/CLOCKS_PER_SEC);
	time = clock();

	}//t_sink loop ends here

//normalization
	for(int inp=0; inp<number_of_mom; ++inp){
		for(int gi=0; gi<4; gi++){
			corr_dbarcp[inp][gi] /= (double)(norm_2pt_1*np[inp]*Lt);
			corr_cmbaru[inp][gi] /= (double)(norm_2pt_1*np[inp]*Lt);
		}
	}

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2)
        for(int i=0; i<np[inp1]*np[inp2]; ++i)
	for(int gi=0; gi<3; ++gi){
        for(int j=0; j<4; ++j){
                corr_DDstar_DDstar[inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        for(int j=0; j<2; ++j){
                corr_DDstar_DstarDstar[inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
                corr_DstarDstar_DDstar[inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        corr_DstarDstar_DstarDstar[inp1][inp2][i][gi][0] /= (double)(norm_4pt*Lt);
	}

        for(int n=0; n<np[1]/2; ++n)
	for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2)
        for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i)
	for(int gi=0; gi<3; ++gi){
        for(int j=0; j<4; ++j){
                corr_DDstar_DDstar_TP1[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        for(int j=0; j<2; ++j){
                corr_DDstar_DstarDstar_TP1[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
                corr_DstarDstar_DDstar_TP1[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        corr_DstarDstar_DstarDstar_TP1[n][inp1][inp2][i][gi][0] /= (double)(norm_4pt*Lt);
	}

        for(int n=0; n<np[2]/2; ++n)
	for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2)
        for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i)
	for(int gi=0; gi<3; ++gi){
        for(int j=0; j<4; ++j){
                corr_DDstar_DDstar_TP2[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        for(int j=0; j<2; ++j){
                corr_DDstar_DstarDstar_TP2[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
                corr_DstarDstar_DDstar_TP2[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        corr_DstarDstar_DstarDstar_TP2[n][inp1][inp2][i][gi][0] /= (double)(norm_4pt*Lt);
	}

        for(int n=0; n<np[3]/2; ++n)
	for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2)
        for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i)
	for(int gi=0; gi<3; ++gi){
        for(int j=0; j<4; ++j){
                corr_DDstar_DDstar_TP3[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        for(int j=0; j<2; ++j){
                corr_DDstar_DstarDstar_TP3[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
                corr_DstarDstar_DDstar_TP3[n][inp1][inp2][i][gi][j] /= (double)(norm_4pt*Lt);
        }
        corr_DstarDstar_DstarDstar_TP3[n][inp1][inp2][i][gi][0] /= (double)(norm_4pt*Lt);
	}

std::cout<<"normalization finished"<<std::endl;

	std::ofstream out;
        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<4; ++s){
        sprintf(outfile, "%s/DDstar_DDstar_corr%d_g%d_TP0_%d%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<np[inp1]; ++i)
        for(int j=0; j<np[inp2]; ++j){
        out<<p[inp1][i][0]<<" "<<p[inp1][i][1]<<" "<<p[inp1][i][2]<<std::endl;
        out<<p[inp2][j][0]<<" "<<p[inp2][j][1]<<" "<<p[inp2][j][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DDstar_DstarDstar_corr%d_g%d_TP0_%d%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<np[inp1]; ++i)
        for(int j=0; j<np[inp2]; ++j){
        out<<p[inp1][i][0]<<" "<<p[inp1][i][1]<<" "<<p[inp1][i][2]<<std::endl;
        out<<p[inp2][j][0]<<" "<<p[inp2][j][1]<<" "<<p[inp2][j][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DstarDstar_DDstar_corr%d_g%d_TP0_%d%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<np[inp1]; ++i)
        for(int j=0; j<np[inp2]; ++j){
        out<<p[inp1][i][0]<<" "<<p[inp1][i][1]<<" "<<p[inp1][i][2]<<std::endl;
        out<<p[inp2][j][0]<<" "<<p[inp2][j][1]<<" "<<p[inp2][j][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<1; ++s){
        sprintf(outfile, "%s/DstarDstar_DstarDstar_corr%d_g%d_TP0_%d%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<np[inp1]; ++i)
        for(int j=0; j<np[inp2]; ++j){
        out<<p[inp1][i][0]<<" "<<p[inp1][i][1]<<" "<<p[inp1][i][2]<<std::endl;
        out<<p[inp2][j][0]<<" "<<p[inp2][j][1]<<" "<<p[inp2][j][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DstarDstar[inp1][inp2][i*np[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

//TP1 output
	for(int n=0; n<np[1]/2; ++n)
        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<4; ++s){
        sprintf(outfile, "%s/DDstar_DDstar_corr%d_g%d_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP1[inp1]; ++i)
        for(int j=0; j<n_TP1[inp2]; ++j){
        out<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][0]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][1]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][0]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][1]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][0]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][1]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][0]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][1]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[1]/2; ++n)
        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DDstar_DstarDstar_corr%d_g%d_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP1[inp1]; ++i)
        for(int j=0; j<n_TP1[inp2]; ++j){
        out<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][0]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][1]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][0]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][1]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][0]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][1]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][0]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][1]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[1]/2; ++n)
        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DstarDstar_DDstar_corr%d_g%d_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP1[inp1]; ++i)
        for(int j=0; j<n_TP1[inp2]; ++j){
        out<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][0]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][1]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][0]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][1]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][0]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][1]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][0]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][1]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[1]/2; ++n)
        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<1; ++s){
        sprintf(outfile, "%s/DstarDstar_DstarDstar_corr%d_g%d_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP1[inp1]; ++i)
        for(int j=0; j<n_TP1[inp2]; ++j){
        out<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][0]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][1]<<" "<<p[TP1[n][inp1][i][0][0]][TP1[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][0]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][1]<<" "<<p[TP1[n][inp1][i][1][0]][TP1[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][0]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][1]<<" "<<p[TP1[n][inp2][j][0][0]][TP1[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][0]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][1]<<" "<<p[TP1[n][inp2][j][1][0]][TP1[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP1[n][inp1][inp2][i*n_TP1[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }
//TP2 output
	for(int n=0; n<np[2]/2; ++n)
        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<4; ++s){
        sprintf(outfile, "%s/DDstar_DDstar_corr%d_g%d_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP2[inp1]; ++i)
        for(int j=0; j<n_TP2[inp2]; ++j){
        out<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][0]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][1]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][0]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][1]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][0]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][1]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][0]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][1]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[2]/2; ++n)
        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DDstar_DstarDstar_corr%d_g%d_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP2[inp1]; ++i)
        for(int j=0; j<n_TP2[inp2]; ++j){
        out<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][0]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][1]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][0]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][1]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][0]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][1]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][0]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][1]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[2]/2; ++n)
        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DstarDstar_DDstar_corr%d_g%d_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP2[inp1]; ++i)
        for(int j=0; j<n_TP2[inp2]; ++j){
        out<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][0]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][1]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][0]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][1]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][0]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][1]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][0]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][1]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[2]/2; ++n)
        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<1; ++s){
        sprintf(outfile, "%s/DstarDstar_DstarDstar_corr%d_g%d_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP2[inp1]; ++i)
        for(int j=0; j<n_TP2[inp2]; ++j){
        out<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][0]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][1]<<" "<<p[TP2[n][inp1][i][0][0]][TP2[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][0]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][1]<<" "<<p[TP2[n][inp1][i][1][0]][TP2[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][0]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][1]<<" "<<p[TP2[n][inp2][j][0][0]][TP2[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][0]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][1]<<" "<<p[TP2[n][inp2][j][1][0]][TP2[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP2[n][inp1][inp2][i*n_TP2[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }



//TP3 output
	for(int n=0; n<np[3]/2; ++n)
        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<4; ++s){
        sprintf(outfile, "%s/DDstar_DDstar_corr%d_g%d_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP3[inp1]; ++i)
        for(int j=0; j<n_TP3[inp2]; ++j){
        out<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][0]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][1]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][0]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][1]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][0]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][1]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][0]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][1]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[3]/2; ++n)
        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DDstar_DstarDstar_corr%d_g%d_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP3[inp1]; ++i)
        for(int j=0; j<n_TP3[inp2]; ++j){
        out<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][0]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][1]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][0]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][1]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][0]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][1]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][0]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][1]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<" "<<imag(corr_DDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[3]/2; ++n)
        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<2; ++s){
        sprintf(outfile, "%s/DstarDstar_DDstar_corr%d_g%d_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP3[inp1]; ++i)
        for(int j=0; j<n_TP3[inp2]; ++j){
        out<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][0]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][1]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][0]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][1]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][0]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][1]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][0]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][1]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int n=0; n<np[3]/2; ++n)
        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
	for(int gi=0; gi<3; ++gi)
        for(int s=0; s<1; ++s){
        sprintf(outfile, "%s/DstarDstar_DstarDstar_corr%d_g%d_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), s+1, gi+1, inp1, inp2, n, config_i);
        out.open(outfile);
        if(out.is_open()){
        for(int i=0; i<n_TP3[inp1]; ++i)
        for(int j=0; j<n_TP3[inp2]; ++j){
        out<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][0]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][1]<<" "<<p[TP3[n][inp1][i][0][0]][TP3[n][inp1][i][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][0]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][1]<<" "<<p[TP3[n][inp1][i][1][0]][TP3[n][inp1][i][1][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][0]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][1]<<" "<<p[TP3[n][inp2][j][0][0]][TP3[n][inp2][j][0][1]][2]<<std::endl;
        out<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][0]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][1]<<" "<<p[TP3[n][inp2][j][1][0]][TP3[n][inp2][j][1][1]][2]<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP3[n][inp1][inp2][i*n_TP3[inp2]+j][gi][s](t))<<std::endl;
        }
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        }
        }

	for(int gi=0; gi<4; gi++)
        for(int inp=0; inp<number_of_mom; ++inp){
        sprintf(outfile, "%s/sum/dbarcp_gamma%d_corr_p%d.conf%d.dat", path_output.c_str(), (gi==3)?5:(gi+1), inp, config_i);
        out.open(outfile);
        if(out.is_open()){
        out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_dbarcp[inp][gi](t))<<" "<<imag(corr_dbarcp[inp][gi](t))<<std::endl;
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        } 

	for(int gi=0; gi<4; gi++)
        for(int inp=0; inp<number_of_mom; ++inp){
        sprintf(outfile, "%s/sum/cmbaru_gamma%d_corr_p%d.conf%d.dat", path_output.c_str(), (gi==3)?5:(gi+1), inp, config_i);
        out.open(outfile);
        if(out.is_open()){
        out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        for(int t=0; t<Lt; ++t)
                out<<std::setprecision(16)<<t<<" "<<real(corr_cmbaru[inp][gi](t))<<" "<<imag(corr_cmbaru[inp][gi](t))<<std::endl;
        }
        else {std::cout<<"error opening file"<<std::endl; exit(0);}
        out.close();
        } 
	
//do average
std::cout<<"start average"<<std::endl;
	for(int inp1=0; inp1<TP0_dim; inp1++)
		for(int inp2=0; inp2<TP0_dim; inp2++){
			for(int i=0; i<np[inp1]*np[inp2]; ++i)
			for(int k=0; k<3; ++k){
			for(int j=0; j<4; ++j){
				corr_DDstar_DDstar_sum[inp1][inp2] += corr_DDstar_DDstar[inp1][inp2][i][k][j];
			}
			for(int j=0; j<2; ++j){
				corr_DDstar_DstarDstar_sum[inp1][inp2] += corr_DDstar_DstarDstar[inp1][inp2][i][k][j];
				corr_DstarDstar_DDstar_sum[inp1][inp2] += corr_DstarDstar_DDstar[inp1][inp2][i][k][j];
			}
			corr_DstarDstar_DstarDstar_sum[inp1][inp2] += corr_DstarDstar_DstarDstar[inp1][inp2][i][k][0];
			}
		corr_DDstar_DDstar_sum[inp1][inp2] /= (double)(3.0*sqrt(np[inp1]*np[inp2]));
		corr_DDstar_DstarDstar_sum[inp1][inp2] /= (double)(3.0*sqrt(np[inp1]*np[inp2]));
		corr_DstarDstar_DDstar_sum[inp1][inp2] /= (double)(3.0*sqrt(np[inp1]*np[inp2]));
		corr_DstarDstar_DstarDstar_sum[inp1][inp2] /= (double)(3.0*sqrt(np[inp1]*np[inp2]));
	}

	for(int inp1=0; inp1<TP1_dim; ++inp1)
	for(int inp2=0; inp2<TP1_dim; ++inp2){
		for(int n=0; n<np[1]/2; ++n){
			for(int i=0; i<n_TP1[inp1]*n_TP1[inp2]; ++i)
			for(int k=0; k<3; ++k){
			for(int j=0; j<4; ++j){
				corr_DDstar_DDstar_TP1_sum[n][inp1][inp2] += corr_DDstar_DDstar_TP1[n][inp1][inp2][i][k][j];
			} 
			for(int j=0; j<2; ++j){
				corr_DDstar_DstarDstar_TP1_sum[n][inp1][inp2] += corr_DDstar_DstarDstar_TP1[n][inp1][inp2][i][k][j];
				corr_DstarDstar_DDstar_TP1_sum[n][inp1][inp2] += corr_DstarDstar_DDstar_TP1[n][inp1][inp2][i][k][j];
			} 
			corr_DstarDstar_DstarDstar_TP1_sum[n][inp1][inp2] += corr_DstarDstar_DstarDstar_TP1[n][inp1][inp2][i][k][0];
			}
			corr_DDstar_DDstar_TP1_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP1[inp1]*n_TP1[inp2]);
			corr_DDstar_DstarDstar_TP1_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP1[inp1]*n_TP1[inp2]);
			corr_DstarDstar_DDstar_TP1_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP1[inp1]*n_TP1[inp2]);
			corr_DstarDstar_DstarDstar_TP1_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP1[inp1]*n_TP1[inp2]);
			corr_DDstar_DDstar_TP1_sumsum[inp1][inp2] += corr_DDstar_DDstar_TP1_sum[n][inp1][inp2];
			corr_DDstar_DstarDstar_TP1_sumsum[inp1][inp2] += corr_DDstar_DstarDstar_TP1_sum[n][inp1][inp2];
			corr_DstarDstar_DDstar_TP1_sumsum[inp1][inp2] += corr_DstarDstar_DDstar_TP1_sum[n][inp1][inp2];
			corr_DstarDstar_DstarDstar_TP1_sumsum[inp1][inp2] += corr_DstarDstar_DstarDstar_TP1_sum[n][inp1][inp2];
		}
		corr_DDstar_DDstar_TP1_sumsum[inp1][inp2] /= (double)(np[1]/2);
		corr_DDstar_DstarDstar_TP1_sumsum[inp1][inp2] /= (double)(np[1]/2);
		corr_DstarDstar_DDstar_TP1_sumsum[inp1][inp2] /= (double)(np[1]/2);
		corr_DstarDstar_DstarDstar_TP1_sumsum[inp1][inp2] /= (double)(np[1]/2);
	}

	for(int inp1=0; inp1<TP2_dim; ++inp1)
	for(int inp2=0; inp2<TP2_dim; ++inp2){
		for(int n=0; n<np[2]/2; ++n){
			for(int i=0; i<n_TP2[inp1]*n_TP2[inp2]; ++i)
			for(int k=0; k<3; ++k){
			for(int j=0; j<4; ++j){
				corr_DDstar_DDstar_TP2_sum[n][inp1][inp2] += corr_DDstar_DDstar_TP2[n][inp1][inp2][i][k][j];
			} 
			for(int j=0; j<2; ++j){
				corr_DDstar_DstarDstar_TP2_sum[n][inp1][inp2] += corr_DDstar_DstarDstar_TP2[n][inp1][inp2][i][k][j];
				corr_DstarDstar_DDstar_TP2_sum[n][inp1][inp2] += corr_DstarDstar_DDstar_TP2[n][inp1][inp2][i][k][j];
			} 
			corr_DstarDstar_DstarDstar_TP2_sum[n][inp1][inp2] += corr_DstarDstar_DstarDstar_TP2[n][inp1][inp2][i][k][0];
			}
			corr_DDstar_DDstar_TP2_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP2[inp1]*n_TP2[inp2]);
			corr_DDstar_DstarDstar_TP2_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP2[inp1]*n_TP2[inp2]);
			corr_DstarDstar_DDstar_TP2_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP2[inp1]*n_TP2[inp2]);
			corr_DstarDstar_DstarDstar_TP2_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP2[inp1]*n_TP2[inp2]);
			corr_DDstar_DDstar_TP2_sumsum[inp1][inp2] += corr_DDstar_DDstar_TP2_sum[n][inp1][inp2];
			corr_DDstar_DstarDstar_TP2_sumsum[inp1][inp2] += corr_DDstar_DstarDstar_TP2_sum[n][inp1][inp2];
			corr_DstarDstar_DDstar_TP2_sumsum[inp1][inp2] += corr_DstarDstar_DDstar_TP2_sum[n][inp1][inp2];
			corr_DstarDstar_DstarDstar_TP2_sumsum[inp1][inp2] += corr_DstarDstar_DstarDstar_TP2_sum[n][inp1][inp2];
		}
		corr_DDstar_DDstar_TP2_sumsum[inp1][inp2] /= (double)(np[2]/2);
		corr_DDstar_DstarDstar_TP2_sumsum[inp1][inp2] /= (double)(np[2]/2);
		corr_DstarDstar_DDstar_TP2_sumsum[inp1][inp2] /= (double)(np[2]/2);
		corr_DstarDstar_DstarDstar_TP2_sumsum[inp1][inp2] /= (double)(np[2]/2);
	}

	for(int inp1=0; inp1<TP3_dim; ++inp1)
	for(int inp2=0; inp2<TP3_dim; ++inp2){
		for(int n=0; n<np[3]/2; ++n){
			for(int i=0; i<n_TP3[inp1]*n_TP3[inp2]; ++i)
			for(int k=0; k<3; ++k){
			for(int j=0; j<4; ++j){
				corr_DDstar_DDstar_TP3_sum[n][inp1][inp2] += corr_DDstar_DDstar_TP3[n][inp1][inp2][i][k][j];
			} 
			for(int j=0; j<2; ++j){
				corr_DDstar_DstarDstar_TP3_sum[n][inp1][inp2] += corr_DDstar_DstarDstar_TP3[n][inp1][inp2][i][k][j];
				corr_DstarDstar_DDstar_TP3_sum[n][inp1][inp2] += corr_DstarDstar_DDstar_TP3[n][inp1][inp2][i][k][j];
			} 
			corr_DstarDstar_DstarDstar_TP3_sum[n][inp1][inp2] += corr_DstarDstar_DstarDstar_TP3[n][inp1][inp2][i][k][0];
			}
			corr_DDstar_DDstar_TP3_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP3[inp1]*n_TP3[inp2]);
			corr_DDstar_DstarDstar_TP3_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP3[inp1]*n_TP3[inp2]);
			corr_DstarDstar_DDstar_TP3_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP3[inp1]*n_TP3[inp2]);
			corr_DstarDstar_DstarDstar_TP3_sum[n][inp1][inp2] /= 3.0*sqrt(n_TP3[inp1]*n_TP3[inp2]);
			corr_DDstar_DDstar_TP3_sumsum[inp1][inp2] += corr_DDstar_DDstar_TP3_sum[n][inp1][inp2];
			corr_DDstar_DstarDstar_TP3_sumsum[inp1][inp2] += corr_DDstar_DstarDstar_TP3_sum[n][inp1][inp2];
			corr_DstarDstar_DDstar_TP3_sumsum[inp1][inp2] += corr_DstarDstar_DDstar_TP3_sum[n][inp1][inp2];
			corr_DstarDstar_DstarDstar_TP3_sumsum[inp1][inp2] += corr_DstarDstar_DstarDstar_TP3_sum[n][inp1][inp2];
		}
		corr_DDstar_DDstar_TP3_sumsum[inp1][inp2] /= (double)(np[3]/2);
		corr_DstarDstar_DstarDstar_TP3_sumsum[inp1][inp2] /= (double)(np[3]/2);
	}


        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
        	sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP0_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_sum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_sum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
        	sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP0_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_sum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_sum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
        	sprintf(outfile, "%s/sum/DstarDstar_DDstar_corr_T1_TP0_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_sum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DDstar_sum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP0_dim; ++inp1)
        for(int inp2=0; inp2<TP0_dim; ++inp2){
        	sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP0_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_sum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_sum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
		for(int n=0; n<np[1]/2; ++n){
        		sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP1_sum[n][inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_TP1_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP1_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP1_sumsum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_TP1_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

	}

        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
		for(int n=0; n<np[1]/2; ++n){
        		sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP1_sum[n][inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_TP1_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP1_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP1_sumsum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_TP1_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

	}

        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
		for(int n=0; n<np[1]/2; ++n){
        		sprintf(outfile, "%s/sum/DstarDstar_DDstar_corr_T1_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP1_sum[n][inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DDstar_TP1_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DstarDstar_DDstar_corr_T1_TP1_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP1_sumsum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DDstar_TP1_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }


        for(int inp1=0; inp1<TP1_dim; ++inp1)
        for(int inp2=0; inp2<TP1_dim; ++inp2){
		for(int n=0; n<np[1]/2; ++n){
        		sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP1_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP1_sum[n][inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP1_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP1_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP1_sumsum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP1_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

 
        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
		for(int n=0; n<np[2]/2; ++n){
        		sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP2_sum[n][inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_TP2_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP2_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP2_sumsum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_TP2_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
		for(int n=0; n<np[2]/2; ++n){
        		sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP2_sum[n][inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_TP2_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP2_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP2_sumsum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_TP2_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
		for(int n=0; n<np[2]/2; ++n){
        		sprintf(outfile, "%s/sum/DstarDstar_DDstar_corr_T1_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP2_sum[n][inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DDstar_TP2_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DstarDstar_DDstar_corr_T1_TP2_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DDstar_TP2_sumsum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DDstar_TP2_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP2_dim; ++inp1)
        for(int inp2=0; inp2<TP2_dim; ++inp2){
		for(int n=0; n<np[2]/2; ++n){
        		sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP2_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP2_sum[n][inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP2_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP2_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP2_sumsum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP2_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
		for(int n=0; n<np[3]/2; ++n){
        		sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP3_sum[n][inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_TP3_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DDstar_DDstar_corr_T1_TP3_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DDstar_TP3_sumsum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DDstar_TP3_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
		for(int n=0; n<np[3]/2; ++n){
        		sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP3_sum[n][inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_TP3_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DDstar_DstarDstar_corr_T1_TP3_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DDstar_DstarDstar_TP3_sumsum[inp1][inp2](t))<<" "<<imag(corr_DDstar_DstarDstar_TP3_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

        for(int inp1=0; inp1<TP3_dim; ++inp1)
        for(int inp2=0; inp2<TP3_dim; ++inp2){
		for(int n=0; n<np[3]/2; ++n){
        		sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP3_%d%d.%d.conf%d.dat", path_output.c_str(), inp1, inp2, n, config_i);
        		out.open(outfile);
        		if(out.is_open()){
				out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        			for(int t=0; t<Lt; ++t)
                			out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP3_sum[n][inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP3_sum[n][inp1][inp2](t))<<std::endl;
        		}
        		else {std::cout<<"error opening file"<<std::endl; exit(0);}
        		out.close();
        	}

        	sprintf(outfile, "%s/sum/DstarDstar_DstarDstar_corr_T1_TP3_sum_%d%d.conf%d.dat", path_output.c_str(), inp1, inp2, config_i);
        	out.open(outfile);
        	if(out.is_open()){
			out<<1<<" "<<Lt<<" "<<1<<" "<<Lx<<" "<<1<<std::endl;
        		for(int t=0; t<Lt; ++t)
                		out<<std::setprecision(16)<<t<<" "<<real(corr_DstarDstar_DstarDstar_TP3_sumsum[inp1][inp2](t))<<" "<<imag(corr_DstarDstar_DstarDstar_TP3_sumsum[inp1][inp2](t))<<std::endl;
        		}
        	else {std::cout<<"error opening file"<<std::endl; exit(0);}
        	out.close();

        }

} // configuration loop ends here.


delete[] SourcebarSource_dbarcp;
delete[] SinkbarSink_dbarcp;
delete[] SourcebarSource_cmbaru;
delete[] SinkbarSink_cmbaru;
   time = clock() - time_start;
   printf("\t\t calculate correlation functions success - %.1f seconds \n ", ((float) time)/CLOCKS_PER_SEC);


    return 0;
}

