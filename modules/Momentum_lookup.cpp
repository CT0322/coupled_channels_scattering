
#include "Momentum_lookup.h"


void Momentum_lookup(Eigen::Vector3i TP, int mom1, int mom2, Eigen::Vector2i**& momentum_index, int& size){

        const int number_of_mom=5;
        const int np[5] = {1, 6, 12, 8, 6};
        int p0[1][3] = {{0, 0, 0}};
	int p1[6][3] = {{0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0}};
        int p2[12][3] = {{0, 1, 1}, {0, -1, -1}, {1, 1, 0}, {-1, -1, 0}, {1, 0, 1}, {-1, 0, -1}, {0, 1, -1}, {0, -1, 1}, {1, -1, 0}, {-1, 1, 0}, {1, 0, -1}, {-1, 0, 1}};
        int p3[8][3] = {{1, 1, 1}, {-1, -1, -1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, 1}, {-1, 1, -1}, {-1, 1, 1}, {1, -1, -1}};
        int p4[6][3] = {{0, 0, 2}, {0, 0, -2}, {0, 2, 0}, {0, -2, 0}, {2, 0, 0}, {-2, 0, 0}};
        Eigen::Vector3i** p = new Eigen::Vector3i*[number_of_mom];
        for(int i=0; i<number_of_mom; ++i){
                p[i] = new Eigen::Vector3i[np[i]];
                for(int j=0; j<np[i]; ++j){
                        p[i][j] = Eigen::Vector3i::Zero(3);
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

	int dim = 0;
	for(int ip1=0; ip1<np[mom1]; ip1++)
	for(int ip2=0; ip2<np[mom2]; ip2++){
		if(p[mom1][ip1] + p[mom2][ip2] == TP)
		dim++;
	}

//	std::cout<<dim<<std::endl;
	if(dim==0) {std::cout<<"no momentum combinations found"<<std::endl; exit(0);}

	Eigen::Vector2i** mom_index = new Eigen::Vector2i*[dim];
	for(int i=0; i<dim; i++)
		mom_index[i] = new Eigen::Vector2i[2];

	dim = 0;	
	for(int ip1=0; ip1<np[mom1]; ip1++)
	for(int ip2=0; ip2<np[mom2]; ip2++){
		if(p[mom1][ip1] + p[mom2][ip2] == TP){
		mom_index[dim][0] << mom1, ip1;
		mom_index[dim][1] << mom2, ip2;
		dim++;	
		}
	}


		momentum_index =  mom_index;
		size = dim;

}

