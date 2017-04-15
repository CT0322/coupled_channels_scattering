//everything to read and write from/to files
//readin the perambulators, rundom vectors and eigen vectors.
#include "ReadWrite.h"

// Definition of a pointer on global data
static GlobalData * const global_data = GlobalData::Instance();


/******************************************************************************/
/******************************************************************************/
// constructor ****************************************************************/
/******************************************************************************/
/******************************************************************************/

ReadWrite::ReadWrite () {
  try{
    const int Lt = global_data->get_Lt();
    const int dim_row = global_data->get_dim_row();
    const std::vector<quark> quarks = global_data->get_quarks();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
//    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;

    const int number_of_inversions = (quarks[0].number_of_dilution_T)
        * quarks[0].number_of_dilution_E * quarks[0].number_of_dilution_D;
	int number_of_flavor = quarks.size();
      V = Eigen::MatrixXcd::Zero(dim_row, number_of_eigen_vec);


    // memory for perambulator_c, random vector_c

	perambulator_c = NULL;
	perambulator_s = NULL;
	perambulator_u = NULL;
	rnd_vec_c = NULL;
	rnd_vec_s = NULL;
	rnd_vec_u = NULL;

	
	for(int j=0; j<number_of_flavor; ++j){
                int number_of_rnd_vec = quarks[j].number_of_rnd_vec;
	        if(quarks[j].type == "c"){
			perambulator_c = new Eigen::MatrixXcd[number_of_rnd_vec];


			rnd_vec_c = new Eigen::VectorXcd[number_of_rnd_vec];
			for(int i = 0; i < number_of_rnd_vec; ++i){
				perambulator_c[i] = Eigen::MatrixXcd::Zero(4 * Lt * number_of_eigen_vec, number_of_inversions);
				rnd_vec_c[i] = Eigen::VectorXcd::Zero(4 * Lt * number_of_eigen_vec);
			}

		}
		else if(quarks[j].type == "s"){
			perambulator_s = new Eigen::MatrixXcd[number_of_rnd_vec];
			rnd_vec_s = new Eigen::VectorXcd[number_of_rnd_vec];
			for(int i = 0; i < number_of_rnd_vec; ++i){
				perambulator_s[i] = Eigen::MatrixXcd::Zero(4 * Lt * number_of_eigen_vec, number_of_inversions);
				rnd_vec_s[i] = Eigen::VectorXcd::Zero(4 * Lt * number_of_eigen_vec);
			}
		}
		else if(quarks[j].type == "u"){
			perambulator_u = new Eigen::MatrixXcd[number_of_rnd_vec];
			rnd_vec_u = new Eigen::VectorXcd[number_of_rnd_vec];
			for(int i = 0; i < number_of_rnd_vec; ++i){
				perambulator_u[i] = Eigen::MatrixXcd::Zero(4 * Lt * number_of_eigen_vec, number_of_inversions);
				rnd_vec_u[i] = Eigen::VectorXcd::Zero(4 * Lt * number_of_eigen_vec);
			}
		}
		else {
			std::cout<<"wrong quark flavor"<<std::endl;
			exit(0);
		}
	}

			
std::cout<<"ReadWrite class costructor finished."<<std::endl;

    }//try ends here

  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::ReadWrite\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
// destructor *****************************************************************/
/******************************************************************************/
/******************************************************************************/

ReadWrite::~ReadWrite() {

    const std::vector<quark> quarks = global_data->get_quarks();
    int number_of_flavor = quarks.size();
  try{
	for(int i=0; i<number_of_flavor; ++i){
	
		if(quarks[i].type == "c"){
			delete[] perambulator_c;
			delete[] rnd_vec_c;
		}
		else if(quarks[i].type == "s"){
			delete[] perambulator_s;
			delete[] rnd_vec_s;
		}
		else if(quarks[i].type == "u"){
			delete[] perambulator_u;
			delete[] rnd_vec_u;
		}
		else {
			std::cout<<"wrong quark flavor"<<std::endl;
			exit(0);
		}
	}
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::~ReadWrite\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


void ReadWrite::read_eigenvectors_from_file (const int config_i, const int t) {

  try{
    clock_t time = clock();
//    const int Lt = global_data->get_Lt();
    const int dim_row = global_data->get_dim_row();
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();

    std::string filename = global_data->get_path_eigenvectors() + "/eigenvectors";
    //buffer for read in
    std::complex<double>* eigen_vec = new std::complex<double>[dim_row];

    if(verbose) printf("reading eigen vectors from files:\n");
    else printf("\treading eigenvectors:");
//    fflush(stdout);

      //setting up file
      char name[200];
      sprintf(name, "%s.%04d.%03d", filename.c_str(), config_i, t);
      if(verbose) std::cout << "Reading file: " << name << std::endl;
//      std::ifstream infile(name, std::ifstream::binary);
        std::ifstream infile;
        infile.open(name, std::ifstream::binary);
        if(infile.is_open()){
                 for (int nev = 0; nev < number_of_eigen_vec; ++nev) {
                        infile.read( (char*) eigen_vec, 2*dim_row*sizeof(double));
                        for(int nrow = 0; nrow < dim_row; ++nrow)
                        V(nrow, nev) = eigen_vec[nrow];
                 }
                infile.close();
      // small test of trace and sum over the eigen vector matrix!
                if(verbose){
                std::cout << "trace of V^d*V on t = " << t << ":\t"<< (V.adjoint() * V).trace() << std::endl;
                std::cout << "sum over all entries of V^d*V on t = " << t << ":\t"<< (V.adjoint() * V).sum() << std::endl;
                }
                delete[] eigen_vec;
                time = clock() - time;
                if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
        }
        else {std::cout<<"cann't open eigenvector file:"<<name<<std::endl; exit(0);}

  }


  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::read_eigenvectors_from_file\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::read_perambulators_from_file (const int config_i) {

  try{
    clock_t time = clock();
    char infile[400];
    FILE *fp = NULL;
    const int Lt = global_data->get_Lt();
    const int Lx = global_data->get_Lx();
    const int Ly = global_data->get_Ly();
    const int Lz = global_data->get_Lz();
    
    const int verbose = global_data->get_verbose();
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const std::vector<quark> quarks = global_data->get_quarks();
 //   const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_dilution_T = quarks[0].number_of_dilution_T;
    const int number_of_dilution_E = quarks[0].number_of_dilution_E;
    const int number_of_dilution_D = quarks[0].number_of_dilution_D;
    const int number_of_inversions = number_of_dilution_T * number_of_dilution_E * number_of_dilution_D;
    const int size_perambulator_entry = number_of_inversions * Lt * 4 * number_of_eigen_vec;

    // memory for reading perambulators
    std::complex<double>* perambulator_read = new std::complex<double>[size_perambulator_entry];

    if(verbose){
      printf("reading perambulators from files:\n");
    }
    else{
      printf("\treading perambulators:");
    }

//read perambulaors_c
    int number_of_flavor = quarks.size();
    for(int i=0; i<number_of_flavor; ++i){
        int number_of_rnd_vec = quarks[i].number_of_rnd_vec;
   
	if(quarks[i].type == "c"){
		for(int rnd_vec_0i = 0; rnd_vec_0i < number_of_rnd_vec; ++rnd_vec_0i){
		std::string filename = global_data->get_path_perambulators_c();

		sprintf(infile, "%s/cnfg%d/rnd_vec_%01d/perambulator.rndvecnb%02d.c.Tso%s%04d.Vso%s%04d.DsoF4.TsiF%04d.SsiF%d.DsiF4.CsiF3.smeared1.%05d", filename.c_str(),config_i, rnd_vec_0i, rnd_vec_0i, &quarks[i].dilution_T[1], number_of_dilution_T, &quarks[i].dilution_E[1], number_of_dilution_E, Lt, Lx*Ly*Lz, config_i);


		if((fp = fopen(infile, "rb")) == NULL){
			std::cout << "failed to open file: " << infile << "\n" << std::endl;
			exit(0);
		}
		if(verbose) printf("\tread file: %s\n", infile);

		fread(perambulator_read, sizeof(std::complex<double>), size_perambulator_entry, fp);

      // copy into matrix structure
      int col_i, row_i;
#pragma omp parallel for private(col_i, row_i) schedule(guided)
      for(int t1 = 0; t1 < Lt; ++t1)
        for(int ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
          for(int dirac1 = 0; dirac1 < 4; ++dirac1)
            for(int t2 = 0; t2 < number_of_dilution_T; ++t2)
              for(int ev2 = 0; ev2 < number_of_dilution_E; ++ev2)
                for(int dirac2 = 0; dirac2 < number_of_dilution_D; ++dirac2){
                  row_i = 4 * number_of_eigen_vec * t1 + 4 * ev1 + dirac1;
                  col_i = number_of_dilution_D * number_of_dilution_E * t2 + number_of_dilution_D * ev2 + dirac2;
                  perambulator_c[rnd_vec_0i](t1* 4 * number_of_eigen_vec +  dirac1 * number_of_eigen_vec + ev1, t2 * number_of_dilution_D * number_of_dilution_E + dirac2 * number_of_dilution_E + ev2) = perambulator_read[row_i * number_of_inversions + col_i];
            }
      fclose(fp);
		} // read perambulators_c ends here.

	}

// read perambulator_s

	else if(quarks[i].type == "s"){
      for(int rnd_vec_0i = 0; rnd_vec_0i < number_of_rnd_vec; ++rnd_vec_0i){
      std::string filename = global_data->get_path_perambulators_s();


      sprintf(infile,
          "%s/cnfg%04d/rnd_vec_%02d/perambulator.rndvecnb%02d.u.Tso%s%04d.Vso%s%04d.DsoF4.TsiF%04d."
          "SsiF%05d.DsiF4.CsiF3.smeared1.%05d", filename.c_str(),config_i, rnd_vec_0i, rnd_vec_0i, &quarks[i].dilution_T[1],
          number_of_dilution_T, &quarks[i].dilution_E[1], number_of_dilution_E, Lt, Lx*Ly*Lz, config_i);

      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      if(verbose) printf("\tread file: %s\n", infile);

      fread(perambulator_read, sizeof(std::complex<double>),
          size_perambulator_entry, fp);

      // copy into matrix structure
      int col_i, row_i;
#pragma omp parallel for private(col_i, row_i) schedule(guided)
      for(int t1 = 0; t1 < Lt; ++t1)
        for(int ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
          for(int dirac1 = 0; dirac1 < 4; ++dirac1)
            for(int t2 = 0; t2 < number_of_dilution_T; ++t2)
              for(int ev2 = 0; ev2 < number_of_dilution_E; ++ev2)
                for(int dirac2 = 0; dirac2 < number_of_dilution_D; ++dirac2){
                  row_i = 4 * number_of_eigen_vec * t1 + 4 * ev1 + dirac1;
                  col_i = number_of_dilution_D * number_of_dilution_E * t2 + number_of_dilution_D * ev2 + dirac2;
                  perambulator_s[rnd_vec_0i](t1 * 4 * number_of_eigen_vec + dirac1 * number_of_eigen_vec + ev1, t2 * number_of_dilution_D * number_of_dilution_E  + dirac2 * number_of_dilution_E + ev2) = perambulator_read[row_i * number_of_inversions + col_i];
            }
        fclose(fp);
      } //read perambulators_s ends here

	}

// read perambulaotrs_u
	else if(quarks[i].type == "u"){

      for(int rnd_vec_0i = 0; rnd_vec_0i < number_of_rnd_vec; ++rnd_vec_0i){
      std::string filename = global_data->get_path_perambulators_u();

      sprintf(infile,
          "%s/cnfg%d/rnd_vec_%01d/perambulator.rndvecnb%02d.u.Tso%s%04d.Vso%s%04d.DsoF4.TsiF%04d."
          "SsiF%d.DsiF4.CsiF3.smeared1.%05d", filename.c_str(),config_i, rnd_vec_0i, rnd_vec_0i, &quarks[i].dilution_T[1],
          number_of_dilution_T, &quarks[i].dilution_E[1], number_of_dilution_E, Lt, Lx*Ly*Lz, config_i);


      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      if(verbose) printf("\tread file: %s\n", infile);

      fread(perambulator_read, sizeof(std::complex<double>),
          size_perambulator_entry, fp);

      // copy into matrix structure
      int col_i, row_i;
#pragma omp parallel for private(col_i, row_i) schedule(guided)
      for(int t1 = 0; t1 < Lt; ++t1)
        for(int ev1 = 0; ev1 < number_of_eigen_vec; ++ev1)
          for(int dirac1 = 0; dirac1 < 4; ++dirac1)
            for(int t2 = 0; t2 < number_of_dilution_T; ++t2)
              for(int ev2 = 0; ev2 < number_of_dilution_E; ++ev2)
                for(int dirac2 = 0; dirac2 < number_of_dilution_D; ++dirac2){
                  row_i = 4 * number_of_eigen_vec * t1 + 4 * ev1 + dirac1;
                  col_i = number_of_dilution_D * number_of_dilution_E * t2 + number_of_dilution_D * ev2 + dirac2;
                  perambulator_u[rnd_vec_0i](t1 * 4 * number_of_eigen_vec + dirac1 * number_of_eigen_vec + ev1, t2 * number_of_dilution_D * number_of_dilution_E  + dirac2 * number_of_dilution_E + ev2) = perambulator_read[row_i * number_of_inversions + col_i];
            }
        fclose(fp);
      } //read perambulator_u ends here

	}

	else{
		std::cout<<"wrong quark flavor"<<std::endl;
		exit(0);
	}

	} // number of quarks flavor loop ends here

    delete[] perambulator_read;
    time = clock() - time;
    if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
  }
  catch(std::exception& e){
    std::cout << e.what()
        << "in: ReadWrite::read_perambulators_from_file\n";
    exit(0);
  }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void ReadWrite::read_rnd_vectors_from_file (const int config_i) {

  try{
    clock_t time = clock();
    char infile[400];
    FILE *fp = NULL;
    const int Lt = global_data->get_Lt();
    const int verbose = global_data->get_verbose();
    const std::vector<quark> quarks = global_data->get_quarks();
//    const int number_of_rnd_vec = quarks[0].number_of_rnd_vec;
    const int number_of_eigen_vec = global_data->get_number_of_eigen_vec();
    const int rnd_vec_length = Lt * number_of_eigen_vec * 4;
    // memory for reading random vectors
    std::complex<double>* rnd_vec_read =
        new std::complex<double>[rnd_vec_length];

    if(verbose){
      printf("reading random vectors from files:\n");
    }
    else{
      printf("\treading random vectors:");
    }


    int number_of_flavor = quarks.size();

	for(int i=0; i<number_of_flavor; ++i){
            int number_of_rnd_vec = quarks[i].number_of_rnd_vec;

//read rnd_vec_c
	if(quarks[i].type == "c"){

    for(int rnd_vec_0i = 0; rnd_vec_0i < number_of_rnd_vec; ++rnd_vec_0i){
      std::string filename = global_data->get_path_perambulators_c();
      sprintf(infile, "%s/cnfg%d/rnd_vec_%01d/randomvector.rndvecnb%02d.c.nbev%04d.%04d", 
          filename.c_str(), config_i, rnd_vec_0i, rnd_vec_0i, number_of_eigen_vec, config_i);

      if(verbose) printf("\tread file: %s\n", infile);
      if((fp = fopen(infile, "rb")) == NULL){
        std::cout << "failed to open file: " << infile << "\n" << std::endl;
        exit(0);
      }
      fread(rnd_vec_read, sizeof(std::complex<double>), rnd_vec_length, fp);

      // copy into matrix structure
      for(int t = 0; t < Lt; ++t)
        for(int nev = 0; nev < number_of_eigen_vec; ++nev)
          for(int dirac = 0; dirac < 4; ++dirac){
        rnd_vec_c[rnd_vec_0i](t*4*number_of_eigen_vec + dirac*number_of_eigen_vec + nev) = rnd_vec_read[t*number_of_eigen_vec*4 + nev*4 + dirac];
      }
      fclose(fp);
    } // read rnd_vec_c ends here
	
	}

//read rnd_vec_s
	else if(quarks[i].type == "s"){
       for(int rnd_vec_0i = 0; rnd_vec_0i < number_of_rnd_vec; ++rnd_vec_0i){
       std::string filename = global_data->get_path_perambulators_s();

       sprintf(infile, "%s/cnfg%d/rnd_vec_%01d/randomvector.rndvecnb%02d.u.nbev%04d.%04d",
          filename.c_str(), config_i, rnd_vec_0i, rnd_vec_0i, number_of_eigen_vec, config_i);
       if(verbose) printf("\tread file: %s\n", infile);
       if((fp = fopen(infile, "rb")) == NULL){ 
         std::cout << "failed to open file: " << infile << "\n" << std::endl;
         exit(0);
       }
       fread(rnd_vec_read, sizeof(std::complex<double>), rnd_vec_length, fp);

      // copy into matrix structure
      for(int t = 0; t < Lt; ++t)
        for(int nev = 0; nev < number_of_eigen_vec; ++nev)
          for(int dirac = 0; dirac < 4; ++dirac){
             rnd_vec_s[rnd_vec_0i](t*4*number_of_eigen_vec + dirac*number_of_eigen_vec + nev) = rnd_vec_read[t*number_of_eigen_vec*4 + nev*4 + dirac];
          }
       fclose(fp);
      } //read rnd_vec_s ends here
	
	}

//read rnd_vec_u
	else if(quarks[i].type == "u"){
       for(int rnd_vec_0i = 0; rnd_vec_0i < number_of_rnd_vec; ++rnd_vec_0i){
       std::string filename = global_data->get_path_perambulators_u();
       sprintf(infile, "%s/cnfg%d/rnd_vec_%01d/randomvector.rndvecnb%02d.u.nbev%04d.%04d",
          filename.c_str(), config_i, rnd_vec_0i, rnd_vec_0i, number_of_eigen_vec, config_i);

       if(verbose) printf("\tread file: %s\n", infile);
       if((fp = fopen(infile, "rb")) == NULL){ 
         std::cout << "failed to open file: " << infile << "\n" << std::endl;
         exit(0);
       }
       fread(rnd_vec_read, sizeof(std::complex<double>), rnd_vec_length, fp);

      // copy into matrix structure
      for(int t = 0; t < Lt; ++t)
        for(int nev = 0; nev < number_of_eigen_vec; ++nev)
          for(int dirac = 0; dirac < 4; ++dirac){
             rnd_vec_u[rnd_vec_0i](t*4*number_of_eigen_vec + dirac*number_of_eigen_vec + nev) = rnd_vec_read[t*number_of_eigen_vec*4 + nev*4 + dirac];
          }
       fclose(fp);
      } // read rnd_vec_u ends here
	
	}
	
	else{
		std::cout<<"wrong quark flavor"<<std::endl;
		exit(0);
	}

	} // quark flavor loop ends here
    delete[] rnd_vec_read;
    time = clock() - time;
    if(!verbose) printf("\t\tSUCCESS - %.1f seconds\n", ((float) time)/CLOCKS_PER_SEC);
  }
  catch(std::exception& e){
    std::cout << e.what() << "in: ReadWrite::read_eigenvectors_from_file\n";
    exit(0);
  }
}
