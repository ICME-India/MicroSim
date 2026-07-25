#include "utilities/IO.hpp"
#include <gsl/gsl_errno.h>

#include <hdf5.h>
#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cstring>

IO::Read_input::Read_input(const std::string Infile, const std::string Filling_file, microsim::Info *Parameters, sycl::queue &dev, std::string &filename_, Logger &Logger)
    :infile(Infile), filling_file(Filling_file), queue(dev), filename_(filename_), logger(Logger), global_param(Parameters){}

void IO::Read_input::Read_infile(){
    std::ifstream Infile(infile);

    if(!Infile){
        logger.log("Unable to find infile", LogLevel::ERROR);
    }

    std::string line;
    while (std::getline(Infile, line)) {

        // Skip comments and empty lines
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//") {
            continue;
        }

        std::string Key, Value;
        const char* param;
        vector<string> tokens;
        Parse_scalars(line, Key, Value);
        //logger.log("Key: " + Key, LogLevel::DEBUG);
        get_tokens(line, Key, Value, tokens);
        
        //logger.log("Processing Key: " + Key, LogLevel::DEBUG); // Added debug print

        param = Value.c_str();

        Input_parameters *INP = (global_param->parameters);
        Settings *flags = (global_param->Flags);
        Domain *domain = (global_param->domainInfo);

        if (Key == "DIMENSION"){
            domain->DIMENSION = std::stoi(Value);
            if(domain->DIMENSION != 2 && domain->DIMENSION != 3){
                logger.log("Dimention must be either 2 or 3",LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
        else if (Key == "MESH_X") { domain->MESH_X = std::stoi(param);}
        else if (Key == "MESH_Y") { domain->MESH_Y = std::stoi(param);}
        else if (Key == "MESH_Z") { domain->MESH_Z = (domain->DIMENSION == 3) ? std::stod(param) : 1;}
        else if (Key == "DELTA_X"){ domain->DELTA_X = std::stod(param);}
        else if (Key == "DELTA_Y"){ domain->DELTA_Y = std::stod(param);}
        else if (Key == "DELTA_Z"){ domain->DELTA_Z = (domain->DIMENSION == 3) ? std::stod(param) : 1;}
        else if (Key == "DELTA_t"){ domain->DELTA_T = std::stod(param);}
        
        // Reading Input_parameters
        else if (Key == "NSMOOTH"){  INP->Nsmooth =    std::stoi(param);}
        else if (Key == "STARTTIME"){INP->start_time = std::stoi(param);}
        else if (Key == "NUMPHASES"){ 
            INP->NUM_PHASES = std::stoi(param);
            nPhase = std::atoi(param);
            if(INP->NUM_PHASES <=1){
                logger.log("Number of phases should be >= 2",LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }

        else if (Key == "NUMCOMPONENTS"){ 
            INP->NUM_COMPONENTS = atoi(param);
            nComp = atoi(param);
            if(INP->NUM_COMPONENTS<=1){
                logger.log("Number of components should be >= 2",LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
            flags->BINARY  = INP->NUM_COMPONENTS == 2;
            flags->TERNARY = INP->NUM_COMPONENTS == 3;
        }

    }
    
    nPhase = global_param->parameters->NUM_PHASES;
    nComp  = global_param->parameters->NUM_COMPONENTS;  

    Initialize_PF_matrix(nPhase, nComp);
    Init_default_conditions(global_param->bcs);

    Infile.clear();
    Infile.seekg(0, std::ios::beg);

    while (std::getline(Infile, line)) {
        
        // Skip comments and empty lines
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//") {
            continue;
        }

        std::string Key, Value;
        const char* param;
        vector<string> tokens;
        vector<double> tuples;
        Parse_scalars(line, Key, Value);
        get_tokens(line, Key, Value, tokens);
        
        //logger.log("Processing Key: " + Key, LogLevel::DEBUG); // Added debug print

        param = Value.c_str();

        Input_parameters *INP = (global_param->parameters);
        Settings *flags = (global_param->Flags);
        phases_components *ref = (global_param->reference);
        Phasefield_matraix *pfmatrix = (global_param->PFMAT);
        DomainBCs *BCs = (global_param->bcs);

        if(Key == "Restart")        flags->restart   = (std::atoi(param) == 1) ? true : false;
        else if(Key == "NTIMESTEPS")INP->ntime_steps = std::atoi(param);
        else if(Key == "SAVET")     INP->saveT       = std::atoi(param);

        // setting flag
        else if(Key == "WRITE_FORMAT")flags->ASCII      = (!strcmp(param, "ASCII")) ? true : false;
        else if(Key == "RESTART")     flags->restart    = (std::atoi(param) == 1) ? true : false;
        else if(Key == "HDF5")        flags->HDF5       = (std::atoi(param) == 1) ? true : false;
        else if(Key == "ISOTHERMAL")  flags->ISOTHERMAL = (std::atoi(param) == 1) ? true : false;
        else if(Key == "DILUTE")      flags->Dilute     = (std::atoi(param) == 1) ? true : false;
        else if(Key == "SHIFT")       flags->shift      = (std::atoi(param) == 1) ? true : false;
        else if(Key == "NOISE")       flags->Noise      = (std::atoi(param) == 1) ? true : false;
        else if(Key == "ELASTICITY")  flags->elasticity = (std::atoi(param) == 1) ? true : false;

        // Fill parameters
        else if(Key == "COMPONENTS")         ref->componenets         = tokens;
        else if(Key == "PHASES")             ref->phases              = tokens;
        else if(Key == "num_thermo_phase")   INP->Num_thermo_phase    = std::stoi(param);
        else if(Key == "tdb_phases")         ref->phases_tdb          = tokens;
        else if(Key == "phase_map")          ref->phases_map          = tokens;
        else if(Key == "epsilon")            INP->epsilon             = std::stod(Value);
        else if(Key == "tau")                INP->tau                 = std::stod(Value);
        else if(Key == "R")                  INP->R                   = std::stod(Value);
        else if(Key == "V")                  INP->Volm                = std::stod(Value);
        else if(Key == "OBSTACLE")           INP->Obstacle            = std::stoi(Value);
        else if(Key == "Function_W")         INP->Function_W          = std::stoi(Value);
        else if(Key == "Function_anisotropy")INP->function_anisotropy = std::stoi(Value);
        else if(Key == "Anisotropy_type")    INP->FOLD                = std::stoi(Value);
        else if(Key == "Function_F")         flags->Function_F        = std::stoi(Value);
        
        // Parameters specific to KKS model
        else if(Key == "refD")           INP->refD            = std::stod(Value);
        else if(Key == "Er")             INP->Er              = std::stod(Value);
        else if(Key == "epsc")           INP->epsc            = std::stod(Value);
        else if(Key == "epsm")           INP->epsm            = std::stod(Value);
        else if(Key == "W")              INP->W               = std::stod(Value);
        else if(Key == "ee")             INP->ee              = std::stod(Value);
        else if(Key == "IntMobInv")      INP->IntMobInv       = std::stod(Value);
        else if(Key == "a2")             INP->a2              = std::stod(Value);

        else if(Key == "GAMMA"){
            int expected_size = nPhase * (nPhase - 1) * 0.5;
            get_doubles(line, Key, Value, tuples);

            if(expected_size != tuples.size()){
                logger.log("Size of token is not equals the expected size", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_symmetric_matrix_ab(pfmatrix->gamma, nPhase, tuples);

        }

        else if(Key == "Tau"){
            int expected_size = nPhase * (nPhase - 1) * 0.5;
            get_doubles(line, Key, Value, tuples);

            if(expected_size != tuples.size()){
                logger.log("Size of token is not equals the expected size", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_symmetric_matrix_ab(pfmatrix->tau_ab, nPhase, tuples);
        }

        else if(Key == "dab"){
            logger.log("Found dab", LogLevel::DEBUG);
            int expected_size = nPhase * (nPhase - 1) * 0.5;
            get_doubles(line, Key, Value, tuples);

            if(expected_size != tuples.size()){
                logger.log("Size of token is not equals the expected size", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_symmetric_matrix_ab(pfmatrix->dab, nPhase, tuples);
            logger.log("Populated dab", LogLevel::DEBUG);
        }

        else if(Key == "fab"){
            int expected_size = nPhase * (nPhase - 1) * 0.5;
            get_doubles(line, Key, Value, tuples);

            if(expected_size != tuples.size()){
                logger.log("Size of token is not equals the expected size", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_symmetric_matrix_ab(pfmatrix->fab, nPhase, tuples);
        }

        else if(Key == "Gamma_abc" && !flags->BINARY){
            int expected_size = nPhase * (nPhase - 1) * (nPhase - 2) / 6;
            get_doubles(line, Key, Value, tuples);

            if(expected_size != tuples.size()){
                logger.log("Size of token is not equals the expected size", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_symmetric_matrix_abc(pfmatrix->gamma_abc, nPhase, tuples);
        }

        else if(Key == "A"){
            int expected_size = 1 + 0.5 * nComp * (nComp-1);
            get_doubles(line, Key, Value, tuples);

            if(expected_size != tuples.size()){
                logger.log("Size of A matrix is not equals the expected size", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_matrix_A(pfmatrix->A, nPhase, nComp, tuples);
        }

        else if(Key == "Equilibrium_temperature"){
            INP->T_eq = stod(Value);
        }

        else if(Key == "Filling_temperature"){
            INP->Filling_temp = stod(Value);
        }
        
        else if(Key == "T"){
            INP->T = stod(Value);
        }

        else if(Key == "DIFFUSIVITY"){
            int expected_size = 2 + 0.5 * nComp * (nComp-1);
            get_doubles(line, Key, Value, tuples);

            //if(expected_size != tuples.size()){
            //    logger.log("Size of Diffusivity matrix is not equals the expected size", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            populate_diffusivity(pfmatrix->Diffusivity, pfmatrix->Inv_Diffusivity, nPhase, nComp, tuples);
            std::cout<<"populated diffusivity"<<std::endl;
        }

        else if(Key == "ceq"){
            int expected_size = 2 + (nComp-1);
            get_doubles(line, Key, Value, tuples);

            //if(expected_size != tuples.size()){
            //    logger.log("Size of ceq matrix is not equals the expected size", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            populate_thermodynamic_matrix(pfmatrix->ceq, nPhase, nComp, tuples);
        }

        else if(Key == "cfill"){
            int expected_size =  2 + (nComp-1);
            get_doubles(line, Key, Value, tuples);

            //if(expected_size != tuples.size()){
            //    logger.log("Size of cfill matrix is not equals the expected size", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            logger.log("Populating cfill for " + std::to_string(tuples[0]) + " " + std::to_string(tuples[1]) + " val: " + std::to_string(tuples[2]), LogLevel::INFO);

            populate_thermodynamic_matrix(pfmatrix->cfill, nPhase, nComp, tuples);

        }


        else if(Key == "c_guess"){
            int expected_size =  2 + (nComp-1);
            get_doubles(line, Key, Value, tuples);

            //if(expected_size != tuples.size()){
            //    logger.log("Size of c_guess matrix is not equals the expected size", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            populate_thermodynamic_matrix(pfmatrix->c_guess, nPhase, nComp, tuples);

        }

        else if(Key == "slopes"){
            int expected_size =  2 + (nComp-1);
            get_doubles(line, Key, Value, tuples);

            //if(expected_size != tuples.size()){
            //    logger.log("Size of slopes matrix is not equals the expected size", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            populate_thermodynamic_matrix(pfmatrix->slopes, nPhase, nComp, tuples);

        }

        else if(flags->anisotropy && Key == "Rotation_matrix"){
            logger.log("Found Rotation_matrix", LogLevel::DEBUG);
            int expected_size = 5;
            get_doubles(line, Key, Value, tuples);

            //if(expected_size != tuples.size()){
            //    logger.log("Size of token is not equals the expected size", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            populate_rotation_matrix(pfmatrix->Rotation_matrix, pfmatrix->Inv_Rotation_matrix, nPhase, tuples);
            logger.log("Populated Rotation_matrix", LogLevel::DEBUG);
        }


        else if(Key == "Tempgrady"){
            int exp_size = 5;
            get_doubles(line, Key, Value, tuples);

            //if(tuples.size() < exp_size){
            //    logger.log("Temperature gradient parameters are not enough, expected size 5", LogLevel::ERROR);
            //    MPI_Abort(MPI_COMM_WORLD, 0);
            //}

            global_param->Temp->base_Temp       = stod(tokens[0]);
            global_param->Temp->DeltaT          = stod(tokens[1]);
            global_param->Temp->Distance        = stod(tokens[2]);
            global_param->Temp->gradient_offset = stod(tokens[3]);
            global_param->Temp->velocity        = stod(tokens[4]);
        }

        else if(Key == "BOUNDARY"){
            std::string fd;
            std::vector<int> bc_vals;
            get_boundary(line, Key, Value, fd, bc_vals);

            if(bc_vals.size() != 6){
                logger.log("Boundary conditions are not enough to apply for all phases", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_boundary_conditions(BCs, bc_vals, fd);
        }

        else if(Key == "BOUNDARY_VALUE"){
            std::string fd;
            std::vector<int> bc_vals;
            get_boundary(line, Key, Value, fd, bc_vals);

            if(bc_vals.size() != 6){
                logger.log("Boundary values are not enough to apply for all conditions", LogLevel::ERROR);
                MPI_Abort(MPI_COMM_WORLD, 0);
            }

            populate_boundary_values(BCs, bc_vals, fd);
        }

        else if(Key == "PERIODICITY"){
            std::string fd;
            std::vector<bool> perd; 
            get_bools(line, Key, Value, perd);

            assert(perd.size() == 3 && "Periodicity vector must have three boolean values.");
            for (int i = 0; i < 3; ++i) {
                flags->periodic[i] = perd[i];
            }
        }
    }

    // Check for missing parameters and set defaults
    Input_parameters *INP = global_param->parameters;
    if (INP->epsilon <= 0.0) {
        INP->epsilon = 4.0 * global_param->domainInfo->DELTA_X;
        std::cout << "WARNING: epsilon not set or invalid. Using default: " << INP->epsilon << std::endl;
        logger.log("WARNING: epsilon not set or invalid. Using default: " + std::to_string(INP->epsilon), LogLevel::WARNING);
    }
    if (INP->tau <= 0.0) {
        // Calculate stable tau based on explicit Euler stability criterion:
        // D_phi * dt / dx^2 <= 0.25 (for 2D/3D safety)
        // D_phi = gamma / tau
        // => (gamma / tau) * dt / dx^2 <= 0.25
        // => tau >= 4 * gamma * dt / dx^2
        
        // Find max gamma
        double max_gamma = 0.0;
        int nP = global_param->parameters->NUM_PHASES;
        for(int i=0; i<nP*nP; ++i) {
            if(global_param->PFMAT->gamma[i] > max_gamma) max_gamma = global_param->PFMAT->gamma[i];
        }
        if(max_gamma == 0.0) max_gamma = 1.0; // Fallback

        double dx = global_param->domainInfo->DELTA_X;
        double dt = global_param->domainInfo->DELTA_T;
        
        // Use safety factor of 1.1
        INP->tau = 1.1 * 4.0 * max_gamma * dt / (dx * dx);
        
        std::cout << "WARNING: tau not set. Calculated stable value: " << INP->tau << std::endl;
        logger.log("WARNING: tau not set. Calculated stable value: " + std::to_string(INP->tau), LogLevel::WARNING);
    }

    Infile.close();

    logger.log("Completed reading the infile", LogLevel::DEBUG);
}



void IO::Read_input::populate_symmetric_matrix_ab(double *Mat, int nPhase, std::vector<double> tokens){

    double **temp = malloc_2d_host(nPhase, nPhase);

    size_t k = 0;
    for(int i=0; i<nPhase; ++i){
        for(int j=i+1; j<nPhase; ++j){
            temp[i][i] = 0.0;
            temp[i][j] = tokens[k++];
            temp[j][i] = temp[i][j];
        }
    }

    for(int i=0; i<nPhase; ++i){
        for(int j=0; j<nPhase; ++j){
            Mat[i*nPhase + j] = temp[i][j];
        }
    }

    free_2d_host(temp, nPhase);
}

void IO::Read_input::populate_symmetric_matrix_abc(double *Mat, int nPhase, std::vector<double> tokens){

    double ***temp;
    temp = malloc_3d_host(nPhase, nPhase, nPhase);
    
    size_t l = 0;
    for(int i=0; i<nPhase; ++i){
        for(int j=i+1; j<nPhase; ++j){
            for(int k= j+1; k<nPhase; ++k){
                temp[i][i][i] = 0.0;
                temp[i][j][j] = 0.0;
                temp[i][k][k] = 0.0;
                temp[i][j][k] = tokens[l++];
                temp[i][k][j] = temp[i][j][k];
                temp[j][i][k] = temp[i][j][k];
                temp[j][k][i] = temp[i][j][k];
                temp[k][i][j] = temp[i][j][k];
                temp[k][j][i] = temp[i][j][k];
            }
        }
    }

    for(int i=0; i<nPhase; ++i){
        for(int j=0; j<nPhase; ++j){
            for(int k=0; k<nPhase; ++k){
                Mat[i*nPhase*nPhase + j*nPhase + k] = temp[i][j][k];
            }
        }
    }

    free_3d_host(temp, nPhase, nPhase);
}



void IO::Read_input::populate_matrix_A(double *Mat, int nPhase, int nComp, std::vector<double> tokens){
    int index_i, l = 0;
    int phase    = static_cast<int>(tokens[l++]);
    
    double*** temp = malloc_3d_host(nPhase, nComp-1, nComp-1);

    for(int i=0; i < nComp-1; ++i){
        temp[phase][i][i] = tokens[l++];
    }

    for(int i=0; i < nComp-1; ++i){
        for(int j=i+1; j < nComp-1; ++j){
            temp[phase][i][j] = tokens[l++];
            temp[phase][j][i] = temp[phase][i][j];
        }
    }

    auto cu_index = [nComp](int a, int b, int c){
        return a * (nComp-1) * (nComp-1) + b * (nComp-1) + c;
    };

    for(int i=0; i<nComp-1; ++i){
        for(int j=0; j<nComp-1; ++j){
            index_i = cu_index(phase, i, j);
            Mat[index_i] = temp[phase][i][j];
        }
    }

    free_3d_host(temp, nPhase, nComp-1);
}

void IO::Read_input::populate_diffusivity(double *Mat, double *InvMat, int nPhase, int nComp, std::vector<double> tokens){
    int l = 0;
    const int diag_mat = static_cast<int>(tokens[l++]);
    const int phase    = static_cast<int>(tokens[l++]);

    const int N = nComp - 1;

    if (phase < 0 || phase >= nPhase) {
        logger.log("Invalid phase index in DIFFUSIVITY", LogLevel::ERROR);
        return;
    }

    std::vector<double> D(N * N, 0.0);
    constexpr double Dmin = 1.0e-12;

    for (int i = 0; i < N; ++i) {
        double val = tokens[l++];
        if (val <= 0.0) {
            logger.log("Zero/negative diffusivity detected; clamping", LogLevel::WARNING);
            val = Dmin;
        }
        D[i * N + i] = val;
    }

    // Full matrix override
    if (!diag_mat) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                D[i * N + j] = tokens[l++];
            }
        }
    }


    std::vector<double> Dinv(N * N, 0.0);
    bool invert_ok = false;

    if (N == 1) {
        // Analytical inversion (binary system)
        if (D[0] > 0.0) {
            Dinv[0] = 1.0 / D[0];
            invert_ok = true;
        }
    } else {
        invert_ok = matrix_inversion_GSL(D.data(), Dinv.data(), N);
    }

    if (!invert_ok) {
        logger.log("Diffusivity matrix is singular; inverse set to zero",
                   LogLevel::WARNING);
    }

    const int offset = phase * N * N;

    for (int i = 0; i < N * N; ++i) {
        Mat[offset + i]    = D[i];
        InvMat[offset + i] = Dinv[i];
    }

    logger.log("Diffusivity matrix populated for phase " + std::to_string(phase),
               LogLevel::DEBUG);
}

void IO::Read_input::populate_thermodynamic_matrix(double *Mat, int nPhase, int nComp, std::vector<double> tokens){
    int l=0;
    int phase1 = tokens[l++];
    int phase2 = tokens[l++];

    double ***temp = malloc_3d_host(nPhase, nPhase, nComp-1);

    for(int i=0; i<nComp-1; ++i){
        temp[phase1][phase2][i] = tokens[l++];
    }

    auto cu_index = [nPhase, nComp](int a, int b, int c){
        return a * nPhase * (nComp-1) + b * (nComp-1) + c;
    };

    for(int i=0; i<nComp-1; ++i){
        int index_i = cu_index(phase1, phase2, i);
        Mat[index_i] = temp[phase1][phase2][i];
    }

    free_3d_host(temp, nPhase, nPhase);
}

double** IO::Read_input::malloc_2d_host(size_t m, size_t n){
    double **Mat = nullptr;
    Mat = new double*[m];
    for(int i=0; i<m; ++i){
        Mat[i] = new double [n]; 
    }
    return Mat;
}

double*** IO::Read_input::malloc_3d_host(size_t m , size_t n, size_t t){
    double ***Mat = nullptr;
    Mat = new double**[m];
    for(int i=0; i<m; ++i){
        Mat[i] = new double*[n];
        for(int j=0; j<n; ++j){
            Mat[i][j] = new double[t];
        }
    }
    return Mat;
}

void IO::Read_input::free_2d_host(double** Mat, size_t m){
    for(int i=0; i<m; ++i){
        delete[] Mat[i];
    }
    delete[] Mat;
}

void IO::Read_input::free_3d_host(double*** Mat, size_t m, size_t n){
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j){
            delete[] Mat[i][j];
        }
        delete[] Mat[i];
    }
    delete[] Mat;
}

bool IO::Read_input::matrix_inversion_GSL(double *Mat, double *InvMat, const int N){
    int s;

    // Allocate memory
    gsl_matrix *matrix      = gsl_matrix_alloc(N, N);
    gsl_matrix *Inv_matrix  = gsl_matrix_alloc(N, N);
    gsl_permutation *perm = gsl_permutation_alloc(N);

    gsl_matrix_set_all(matrix, 0.0);
    for(int i=0; i<N; ++i){
        for(int j=0; j<N; ++j){
            gsl_matrix_set(matrix, i, j, Mat[i*N+j]);
        }
    }

    // Disable default error handler to prevent abort on singular matrix
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();

    int status = gsl_linalg_LU_decomp(matrix, perm, &s);
    
    if(status != GSL_SUCCESS){
        logger.log("GSL matrix inversion failed (LU decomp)", LogLevel::ERROR);
        // Set InvMat to zero
        for(int i=0; i<N*N; ++i) InvMat[i] = 0.0;
    }else{
        status = gsl_linalg_LU_invert(matrix, perm, Inv_matrix);
        if(status != GSL_SUCCESS){
             logger.log("GSL matrix inversion failed (LU invert)", LogLevel::ERROR);
             for(int i=0; i<N*N; ++i) InvMat[i] = 0.0;
        } else {
             logger.log("GSL matrix inversion completed", LogLevel::DEBUG);
             for(int i=0; i<N; ++i){
                for(int j=0; j<N; ++j){
                    InvMat[i*N+j] = gsl_matrix_get(Inv_matrix, i, j);
                }
            }
        }
    }

    // Restore error handler
    gsl_set_error_handler(old_handler);

    // Free the allocated memory
    gsl_permutation_free(perm);
    gsl_matrix_free(matrix);
    gsl_matrix_free(Inv_matrix);
    return status == GSL_SUCCESS;
}

void IO::Read_input::get_rotation_matrix(double **R, double theta, int axis){
    double costheta, sintheta;

    costheta = cos(theta*M_PI/180.0);
    sintheta = sin(theta*M_PI/180.0);

    if (axis == 0) {
        R[0][0] = 1.0;
        R[0][1] = 0.0;
        R[0][2] = 0.0;
        R[1][0] = 0.0;
        R[2][0] = 0.0;
    
        R[1][1] = costheta;
        R[1][2] = -sintheta;
        R[2][1] = sintheta;
        R[2][2] = costheta;
    }

    if (axis == 1) {
        R[0][0] = costheta;
        R[0][1] = 0.0;
        R[0][2] = sintheta;
        R[1][0] = 0.0;
        R[2][0] = -sintheta;

        R[1][1] = 1.0;
        R[1][2] = 0.0;
        R[2][1] = 0.0;
        R[2][2] = costheta;
    }

    if (axis == 2) {
        R[0][0] = costheta;
        R[0][1] = -sintheta;
        R[0][2] = 0.0;
        R[1][0] = sintheta;
        R[2][0] = 0.0;
    
        R[1][1] = costheta;
        R[1][2] = 0.0;
        R[2][1] = 0.0;
        R[2][2] = 1.0;
    }
}

void IO::Read_input::populate_rotation_matrix(double *Mat, double *InvMat, int nPhase, std::vector<double> tokens){
    int a = static_cast<int>(tokens[0]);
    int b = static_cast<int>(tokens[1]);

    int phase1 = fmax(a, b);
    int phase2 = fmin(a, b); 

    double thetax = tokens[2];
    double thetay = tokens[3];
    double thetaz = tokens[4];

    double **RX    = malloc_2d_host(3, 3);
    double **RY    = malloc_2d_host(3, 3);
    double **RZ    = malloc_2d_host(3, 3);
    double **temp1 = malloc_2d_host(3, 3);
    double **temp2 = malloc_2d_host(3, 3);

    get_rotation_matrix(RX, thetax, 0);
    get_rotation_matrix(RY, thetay, 1);
    get_rotation_matrix(RZ, thetaz, 2);

    multiply_2d(RX, RY, temp1, 3);
    multiply_2d(temp1, RZ, temp2, 3);

    double* res_1 = new double[9];
    double* res_2 = new double[9];

    // linearise the rotation matrix
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            res_1[i*3+j] = temp2[i][j];
        }
    }

    matrix_inversion_GSL(res_1, res_2, 3);

    auto cu_index = [nPhase](int p1, int p2, int a, int b){
        return p1 * nPhase * 3 * 3 + p2 * 3 * 3 + a * 3 + b;
    };

    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            int index_i = cu_index(phase1, phase2, i, j);
            Mat[index_i]    = res_1[i*3 + j];
            InvMat[index_i] = res_2[i*3 + j];
        }
    }

    // deallocate the memory
    for(int i=0; i<3; ++i){
        delete[] RX[i];
        delete[] RY[i];
        delete[] RZ[i];
        delete[] temp1[i];
        delete[] temp2[i];
    }
    delete[] RX;
    delete[] RY;
    delete[] RZ;
    delete[] temp1;
    delete[] temp2;
    delete[] res_1;
    delete[] res_2;
}



void IO::Read_input::Initialize_PF_matrix(int nPhase, int nComp){
    Phasefield_matraix *pfmatrix = (global_param->PFMAT);


    pfmatrix->Rotation_matrix     = sycl::malloc_shared<double>(nPhase * nPhase * 3 * 3, queue);
    pfmatrix->Inv_Rotation_matrix = sycl::malloc_shared<double>(nPhase * nPhase * 3 * 3, queue);

    pfmatrix->cmu = sycl::malloc_shared<double>(nPhase * (nComp-1) * (nComp-1), queue);
    pfmatrix->muc = sycl::malloc_shared<double>(nPhase * (nComp-1) * (nComp-1), queue);

    pfmatrix->Diffusivity = sycl::malloc_shared<double>( nPhase * (nComp-1)    * (nComp-1), queue);
    pfmatrix->Inv_Diffusivity = sycl::malloc_shared<double>( nPhase * (nComp-1)    * (nComp-1), queue);

    pfmatrix->A           = sycl::malloc_shared<double>( nPhase * (nComp-1)    * (nComp-1), queue);
    pfmatrix->ceq         = sycl::malloc_shared<double>( nPhase * nPhase * (nComp-1), queue);
    pfmatrix->cfill       = sycl::malloc_shared<double>( nPhase * nPhase * (nComp-1), queue);
    for(size_t i=0; i<nPhase * nPhase * (nComp-1); ++i) pfmatrix->cfill[i] = 0.0;
    // Initialize cfill to zero
    for(int i=0; i<nPhase * nPhase * (nComp-1); ++i) pfmatrix->cfill[i] = 0.0;

    pfmatrix->c_guess     = sycl::malloc_shared<double>( nPhase * nPhase * (nComp-1), queue);
    pfmatrix->slopes      = sycl::malloc_shared<double>( nPhase * nPhase * (nComp-1), queue);
    pfmatrix->dcbdT       = sycl::malloc_shared<double>( nPhase * nPhase * (nComp-1), queue);

    pfmatrix->ceq_coeffs  = sycl::malloc_shared<double>(nPhase * (nComp-1) * 4, queue);
    pfmatrix->gamma_abc   = sycl::malloc_shared<double>(nPhase * nPhase * nPhase, queue);

    pfmatrix->B           = sycl::malloc_shared<double>(nPhase * (nComp-1), queue);
    pfmatrix->Beq         = sycl::malloc_shared<double>(nPhase * (nComp-1), queue);
    pfmatrix->dBbdT       = sycl::malloc_shared<double>(nPhase * (nComp-1), queue);
    pfmatrix->DELTA_T     = sycl::malloc_shared<double>(nPhase *  nPhase  , queue);
    pfmatrix->DELTA_C     = sycl::malloc_shared<double>(nPhase * (nComp-1), queue);
    pfmatrix->dcbdT_phase = sycl::malloc_shared<double>(nPhase * (nComp-1), queue);

    pfmatrix->C           = sycl::malloc_shared<double>(nPhase, queue);
    pfmatrix->Rotated_qab = sycl::malloc_shared<double>(3, queue);

    pfmatrix->gamma  = sycl::malloc_shared<double>(nPhase * nPhase, queue);
    pfmatrix->tau_ab = sycl::malloc_shared<double>(nPhase * nPhase, queue);
    pfmatrix->dab    = sycl::malloc_shared<double>(nPhase * nPhase, queue);
    pfmatrix->fab    = sycl::malloc_shared<double>(nPhase * nPhase, queue);

    logger.log("Initialized the PFMAT matrix", LogLevel::DEBUG);
}


void IO::Read_input::Read_filling_file(){
    std::ifstream Filling_file(filling_file);

    if(!Filling_file){
        logger.log("Unable to find filling file ",LogLevel::ERROR);
    }

    fillParameters *fill_info = global_param->Filling;
    
    fill_info->countFill = 0;

    string line;
    while (std::getline(Filling_file, line)) {

        if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//") {
            continue;
        }else fill_info->countFill++;

    }
    Filling_file.close();

    fill_info->fillType     = sycl::malloc_shared<fill>(fill_info->countFill, queue);
    fill_info->xC           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->yC           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->zC           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->xS           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->xE           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->yS           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->yE           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->zS           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->zE           = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->radius       = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->phase        = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->major_axis   = sycl::malloc_shared<double>(fill_info->countFill, queue); 
    fill_info->eccentricity = sycl::malloc_shared<double>(fill_info->countFill, queue);
    fill_info->rot_angle    = sycl::malloc_shared<double>(fill_info->countFill, queue);
    fill_info->seed         = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->volFrac      = sycl::malloc_shared<double>(fill_info->countFill, queue);
    fill_info->shieldDist   = sycl::malloc_shared<int>(fill_info->countFill, queue);
    fill_info->shiftFrac    = sycl::malloc_shared<int>(fill_info->countFill, queue); 
    fill_info->radVar       = sycl::malloc_shared<double>(fill_info->countFill, queue);


    Filling_file.open(filling_file);

    int count =0;

    do {
        std::getline(Filling_file, line);
        if (line.empty() || line[0] == '#' || line.substr(0, 2) == "//") {
            continue;
        }

        string Key, Value;
        vector<string> tokens;

        get_tokens(line, Key, Value, tokens);

        if(Key == "FILLCYLINDER"){
            fill_info->fillType[count] = FILL_CYLINDER;
            fill_info->phase[count]    = stoi(tokens[0]);
            fill_info->xC[count]       = stoi(tokens[1]);
            fill_info->yC[count]       = stoi(tokens[2]);
            fill_info->zS[count]       = stoi(tokens[3]);
            fill_info->zE[count]       = stoi(tokens[4]);
            fill_info->radius[count]   = stoi(tokens[5]);
            count++;
        }else if(Key == "FILLSPHERE"){
            fill_info->fillType[count] = FILL_SPHERE;
            fill_info->phase[count]    = stoi(tokens[0]);
            fill_info->xC[count]       = stoi(tokens[1]);
            fill_info->yC[count]       = stoi(tokens[2]);
            fill_info->zC[count]       = stoi(tokens[3]);
            fill_info->radius[count]   = stoi(tokens[4]);
            count++;
        }else if(Key == "FILLCUBE"){
            fill_info->fillType[count] = FILL_CUBE;
            fill_info->phase[count]    = stoi(tokens[0]);
            fill_info->xS[count]       = stoi(tokens[1]);
            fill_info->yS[count]       = stoi(tokens[2]);
            fill_info->zS[count]       = stoi(tokens[3]);
            fill_info->xE[count]       = stoi(tokens[4]);
            fill_info->yE[count]       = stoi(tokens[5]);
            fill_info->zE[count]       = stoi(tokens[6]);
            count++;
        }else if(Key == "FILLELLIPSE"){
            fill_info->fillType[count] = FILL_ELLIPSE;
            fill_info->phase[count]        = stoi(tokens[0]);
            fill_info->xC[count]           = stoi(tokens[1]);
            fill_info->yC[count]           = stoi(tokens[2]);
            fill_info->zC[count]           = stoi(tokens[3]);
            fill_info->major_axis[count]   = stof(tokens[5]);
            fill_info->eccentricity[count] = stof(tokens[6]);
            fill_info->rot_angle[count]    = stof(tokens[7]);
            count++;
        }

    }while (count < fill_info->countFill);

    Filling_file.close();

    logger.log("Filling parameters read successfully", LogLevel::DEBUG);

}


IO::Read_input::~Read_input(){}

void IO::Read_input::Parse_scalars(std::string &line, std::string &key, std::string &value){
    const char* whitespace = " \t\n\r\f\v;";
    
    size_t pos = line.find('=');
    if (pos != string::npos) {

        key = line.substr(0, pos);
        value = line.substr(pos + 1);

        // Remove any leading or trailing whitespace
        key.erase(0, key.find_first_not_of(whitespace));
        key.erase(key.find_last_not_of(whitespace) + 1);
        value.erase(0, value.find_first_not_of(whitespace));
        value.erase(value.find_last_not_of(whitespace) + 1);
    }
}

void IO::Read_input::get_tokens(std::string &line, std::string &key, std::string &value, std::vector<std::string> &tokens){
    Parse_scalars(line, key, value);
    std::string trimmed = value.substr(1, value.size() - 2);
    std::stringstream ss(trimmed);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        tokens.push_back(item);
    }
}

void IO::Read_input::get_doubles(std::string &line, std::string &key, std::string &value, std::vector<double> &tokens){
    Parse_scalars(line, key, value);
    std::string trimmed = value.substr(1, value.size() - 2);
    std::stringstream ss(trimmed);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        tokens.push_back(std::stod(item));
    }
}

void IO::Read_input::get_bools(std::string &line, std::string &key, std::string &value, std::vector<bool> &tokens){
    Parse_scalars(line, key, value);
    std::string trimmed = value.substr(1, value.size() - 2);
    std::stringstream ss(trimmed);
    std::string item;
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        
        std::string lower = item;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

        tokens.push_back(lower == "true");
    }
}

void IO::Read_input::get_boundary(string &line, string &key, string &value, string &field, vector<int> &tokens){
    Parse_scalars(line, key, value);
    std::string trimmed = value.substr(1, value.size() - 2);
    std::stringstream ss(trimmed);
    std::string item;
    int count = 0;
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        if(count == 0)field = item;else tokens.push_back(std::stoi(item));
        count++;
    }
}

void IO::Read_input::populate_boundary_conditions(DomainBCs *bcs, std::vector<int> tokens, const std::string fd){
    logger.log("populate_boundary_conditions called for " + fd, LogLevel::INFO);
    int nPhase = global_param->parameters->NUM_PHASES;
    int nComp  = global_param->parameters->NUM_COMPONENTS;

    field val = string_to_field(fd);

    auto set_phase_bc = [nPhase](FaceBCs& f, BCType t){
        for(int i=0; i<nPhase; ++i) f.phi[i].type = t;
    };

    auto set_comp_bc = [nComp](FaceBCs& f, BCType t){
        for(int i=0; i<nComp-1; ++i) { f.comp[i].type = t; f.mu[i].type = t; }
    };

    auto set_mu_bc = [nComp](FaceBCs& f, BCType t){
        for(int i=0; i<nComp-1; ++i) f.mu[i].type = t;
    };

    // Apply boundary conditions based on the field type
    if(val == field::Order_Parameter){
        set_phase_bc(bcs->left,   static_cast<BCType>(tokens[0]));
        set_phase_bc(bcs->right,  static_cast<BCType>(tokens[1]));
        set_phase_bc(bcs->top,    static_cast<BCType>(tokens[2]));
        set_phase_bc(bcs->bottom, static_cast<BCType>(tokens[3]));
        set_phase_bc(bcs->front,  static_cast<BCType>(tokens[4]));
        set_phase_bc(bcs->back,   static_cast<BCType>(tokens[5]));
    }else if(val == field::Composition){
        set_comp_bc(bcs->left,   static_cast<BCType>(tokens[0]));
        set_comp_bc(bcs->right,  static_cast<BCType>(tokens[1]));
        set_comp_bc(bcs->top,    static_cast<BCType>(tokens[2]));
        set_comp_bc(bcs->bottom, static_cast<BCType>(tokens[3]));
        set_comp_bc(bcs->front,  static_cast<BCType>(tokens[4]));
        set_comp_bc(bcs->back,   static_cast<BCType>(tokens[5]));
    }else if(val == field::Chemical_potential){
        set_mu_bc(bcs->left,   static_cast<BCType>(tokens[0]));
        set_mu_bc(bcs->right,  static_cast<BCType>(tokens[1]));
        set_mu_bc(bcs->top,    static_cast<BCType>(tokens[2]));
        set_mu_bc(bcs->bottom, static_cast<BCType>(tokens[3]));
        set_mu_bc(bcs->front,  static_cast<BCType>(tokens[4]));
        set_mu_bc(bcs->back,   static_cast<BCType>(tokens[5]));
    }else if(val == field::Temperature){
        bcs->left.temp->type   = static_cast<BCType>(tokens[0]);
        bcs->right.temp->type  = static_cast<BCType>(tokens[1]);
        bcs->top.temp->type    = static_cast<BCType>(tokens[2]);
        bcs->bottom.temp->type = static_cast<BCType>(tokens[3]);
        bcs->front.temp->type  = static_cast<BCType>(tokens[4]);
        bcs->back.temp->type   = static_cast<BCType>(tokens[5]);
    }
    
}

void IO::Read_input::populate_boundary_values(DomainBCs *bcs, std::vector<int> tokens, const std::string fd){
    int nPhase = global_param->parameters->NUM_PHASES;
    int nComp  = global_param->parameters->NUM_COMPONENTS;

    field val = string_to_field(fd);

    auto set_phase_bc_value = [nPhase](FaceBCs& f, double v){
        for(int i=0; i<nPhase; ++i) f.phi[i].value = v;
    };

    auto set_comp_bc_value = [nComp](FaceBCs& f, double v){
        for(int i=0; i<nComp-1; ++i) { f.comp[i].value = v; f.mu[i].value = v; }
    };

    auto set_mu_bc_value = [nComp](FaceBCs& f, double v){
        for(int i=0; i<nComp-1; ++i) f.mu[i].value = v;
    };

    // Apply boundary values based on the field type
    if(val == field::Order_Parameter){
        set_phase_bc_value(bcs->left,   tokens[0]);
        set_phase_bc_value(bcs->right,  tokens[1]);
        set_phase_bc_value(bcs->top,    tokens[2]);
        set_phase_bc_value(bcs->bottom, tokens[3]);
        set_phase_bc_value(bcs->front,  tokens[4]);
        set_phase_bc_value(bcs->back,   tokens[5]);
    }else if(val == field::Composition){
        set_comp_bc_value(bcs->left,   tokens[0]);
        set_comp_bc_value(bcs->right,  tokens[1]);
        set_comp_bc_value(bcs->top,    tokens[2]);
        set_comp_bc_value(bcs->bottom, tokens[3]);
        set_comp_bc_value(bcs->front,  tokens[4]);
        set_comp_bc_value(bcs->back,   tokens[5]);
    }else if(val == field::Chemical_potential){
        set_mu_bc_value(bcs->left,   tokens[0]);
        set_mu_bc_value(bcs->right,  tokens[1]);
        set_mu_bc_value(bcs->top,    tokens[2]);
        set_mu_bc_value(bcs->bottom, tokens[3]);
        set_mu_bc_value(bcs->front,  tokens[4]);
        set_mu_bc_value(bcs->back,   tokens[5]);
    }else if(val == field::Temperature){
        bcs->left.temp->value   = tokens[0];
        bcs->right.temp->value  = tokens[1];
        bcs->top.temp->value    = tokens[2];
        bcs->bottom.temp->value = tokens[3];
        bcs->front.temp->value  = tokens[4];
        bcs->back.temp->value   = tokens[5];
    }
}


field IO::Read_input::string_to_field(const std::string &str){
    auto toLower = [](std::string s) {
        std::transform(s.begin(), s.end(), s.begin(),
            [](unsigned char c){ return std::tolower(c); });
        return s;
    };

    static const std::map<std::string, field> fieldMap = {
        {"phi", field::Order_Parameter},
        {"c", field::Composition},
        {"mu", field::Chemical_potential},
        {"t", field::Temperature}
    };

    std::string lower_str = toLower(str);

    auto it = fieldMap.find(lower_str);

    if(it != fieldMap.end()){
        return it->second;
    }else{
        logger.log("Unknown boundary condition applied : "+str , LogLevel::ERROR);
        MPI_Abort(MPI_COMM_WORLD, 0);
        return field::Unknown;
    }

}


void IO::Read_input::multiply_2d(double **m1,double **m2,double **prod,long size){
    double sum;
    for(int k=0;k<size;k++)
    {
        for(int i=0;i<size;i++)
        {
            sum=0;
            for(int j=0;j<size;j++)
                sum=sum+m1[k][j]*m2[j][i];
            prod[k][i]=sum;
        }
    }
}


// Initialize default boundary conditions
void IO::Read_input::Init_default_conditions(DomainBCs *bcs){
    int nPhase = global_param->parameters->NUM_PHASES;
    int nComp  = global_param->parameters->NUM_COMPONENTS;

    sycl::queue queue = (this->queue);

    auto set_face_default = [nPhase, nComp](FaceBCs& f, BCType t, double v){
        f.temp->type = t; f.temp->value = v;
        for(int i=0; i<nPhase; ++i) f.phi[i] = {t, v};
        for(int i=0; i<nComp-1; ++i) { f.comp[i] = {t, v}; f.mu[i] = {t, v}; }
    };

    auto alloc_face = [nPhase, nComp, queue](FaceBCs& f){
        f.temp = sycl::malloc_shared<BCValue>(1, queue);
        f.phi  = sycl::malloc_shared<BCValue>(nPhase, queue);
        f.comp = sycl::malloc_shared<BCValue>(nComp-1, queue);
        f.mu   = sycl::malloc_shared<BCValue>(nComp-1, queue);
    };

    // Allocate each face
    alloc_face(bcs->left);  
    alloc_face(bcs->right);
    alloc_face(bcs->top);
    alloc_face(bcs->bottom);
    alloc_face(bcs->front);
    alloc_face(bcs->back);

    // 1. Initialize all to Periodic (Standard case)
    logger.log("Initializing default BCs to PERIODIC (0)", LogLevel::INFO);
    set_face_default(bcs->left  , BC_PERIODIC, 0);
    set_face_default(bcs->right , BC_PERIODIC, 0);
    set_face_default(bcs->top   , BC_PERIODIC, 0);
    set_face_default(bcs->bottom, BC_PERIODIC, 0);
    set_face_default(bcs->front , BC_PERIODIC, 0);
    set_face_default(bcs->back  , BC_PERIODIC, 0);
}


template<int Dim>
void IO::Read_input::write_hdf5(int step, int rank, MPI_Comm comm){
    MPI_Barrier(comm);

    Domain *domain = global_param->domainInfo;
    workers *worker = global_param->workers_mpi;
    DomainDecomp *decomp = global_param->decomp;
    OFFSETS *offsets = global_param->offsets;
    Input_parameters *param = global_param->parameters;

    size_t vol = worker->local_d * worker->local_h * worker->local_w;

    //if(ht.size() != worker->num_pts * offsets->NUM_FIELDS)
    //    ht.resize(worker->num_pts * offsets->NUM_FIELDS);
    //
    //if(hi.size() != vol * offsets->NUM_FIELDS)
    //    hi.resize(vol * offsets->NUM_FIELDS);
    std::vector<double> ht(worker->num_pts * offsets->NUM_FIELDS);
    std::vector<double> hi(vol * offsets->NUM_FIELDS);

    queue.memcpy(ht.data(), global_param->grid_data, worker->num_pts * offsets->NUM_FIELDS * sizeof(double)).wait();
    
    // Extract only the local domain data
    for(int z=0; z<worker->local_d; ++z){
        for(int y=0; y<worker->local_h; ++y){
            for(int x=0; x<worker->local_w; ++x) {
                int s = idx(z+(Dim==3?1:0), y+1, x+1, worker->ghost_h, worker->ghost_w);
                int d = idx(z, y, x, worker->local_h, worker->local_w);
                std::memcpy(&hi[d*offsets->NUM_FIELDS], &ht[s*offsets->NUM_FIELDS], offsets->NUM_FIELDS*sizeof(double));
            }
        }
    }

    std::ostringstream ss;
    ss << filename_<< "_" << std::setw(5) << std::setfill('0') << step;
    std::string h5_name = ss.str() + ".h5"; std::string xmf_name = "DATA/" + ss.str() + ".xmf";

    if(rank == 0 && step == 0){
        std::remove(h5_name.c_str());
        std::remove(xmf_name.c_str());
    }

    MPI_Barrier(comm);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, comm, MPI_INFO_NULL);
    hid_t file = H5Fcreate(("DATA/" + h5_name).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    H5Pclose(plist);

    int g_d = (Dim==2) ? 1 : domain->MESH_Z;
    int g_h = domain->MESH_Y;
    int g_w = domain->MESH_X;
    hsize_t g_dims[3] = { (hsize_t)g_d, (hsize_t)g_h, (hsize_t)g_w};
    hid_t fspace = H5Screate_simple(3, g_dims, NULL);

    hsize_t offset[3];
    if constexpr (Dim == 3) {
        offset[0] = decomp->coords[0] * worker->local_d;
        offset[1] = decomp->coords[1] * worker->local_h;
        offset[2] = decomp->coords[2] * worker->local_w;
    }else {
        offset[0]=0;
        offset[1]=decomp->coords[0] * worker->local_h;
        offset[2]=decomp->coords[1] * worker->local_w;
    }

    int loc_d = worker->local_d;
    int loc_h = worker->local_h;
    int loc_w = worker->local_w;

    hsize_t count[3] = { (hsize_t)loc_d, (hsize_t)loc_h, (hsize_t)loc_w };
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    hid_t mspace = H5Screate_simple(3, count, NULL);
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

    std::vector<double> buf(vol);
    auto write_var = [&](std::string n, int off) {
        for(size_t i=0; i<vol; ++i) buf[i] = hi[i*offsets->NUM_FIELDS + off];
        hid_t dset = H5Dcreate2(file, n.c_str(), H5T_NATIVE_DOUBLE, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, mspace, fspace, dxpl, buf.data());
        H5Dclose(dset);
    };

    int nPhase = param->NUM_PHASES;
    int nComp  = param->NUM_COMPONENTS;

    int off_phi  = offsets->OFF_PHI;
    int off_comp = offsets->OFF_COMP;

    #ifndef __FID
    int off_mu   = offsets->OFF_MU;
    #elif defined(__FID)
    int off_mus = offsets->OFF_MUS;
    int off_mul = offsets->OFF_MUL;
    #endif


    int off_temp   = offsets->OFF_TEMP;

    for(int p=0; p<nPhase; ++p){
        write_var("/phi_"+std::to_string(p), off_phi+p);
    }
    for(int c=0; c<nComp-1; ++c){
        write_var("/comp_"+std::to_string(c), off_comp+c);
    }

    #ifndef __FID
    for(int c=0; c<nComp-1; ++c){
        write_var("/mu_"+std::to_string(c), off_mu+c);
    }
    #elif defined(__FID)
    for(int c=0; c<nComp-1; ++c){
        write_var("/mu_s_"+std::to_string(c), off_mus+c);
        write_var("/mu_l_"+std::to_string(c), off_mul+c);
    }
    #endif

    write_var("/temperature", off_temp);

    H5Pclose(dxpl);
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Fclose(file);


    // writing the XMF file
    if (rank == 0) {
        std::ofstream xmf(xmf_name);
        xmf << "<?xml version=\"1.0\" ?>\n<Xdmf Version=\"2.0\">\n <Domain>\n";
        xmf << "   <Grid Name=\"Mesh\" GridType=\"Uniform\">\n";
        xmf << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << g_d << " " << g_h << " " << g_w << "\"/>\n";
        xmf << "     <Geometry GeometryType=\"Origin_DxDyDz\">\n";
        xmf << "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Format=\"XML\">0 0 0</DataItem>\n";
        xmf << "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Format=\"XML\">1 1 1</DataItem>\n";
        xmf << "     </Geometry>\n";
        auto add = [&](std::string n, std::string p) {
            xmf << "     <Attribute Name=\"" << n << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
            xmf << "       <DataItem Dimensions=\"" << g_d << " " << g_h << " " << g_w << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
            xmf << "        " << h5_name << ":" << p << "\n";
            xmf << "       </DataItem>\n     </Attribute>\n";
        };

        for(int p=0; p<nPhase; ++p){
            add("Phi_"+std::to_string(p), "/phi_"+std::to_string(p));
        }
        for(int c=0; c<nComp-1; ++c){
            add("Comp_"+std::to_string(c), "/comp_"+std::to_string(c));
        }

        #ifndef __FID
        for(int c=0; c<nComp-1; ++c){
            add("Mu_"+std::to_string(c), "/mu_"+std::to_string(c));
        }
        #elif defined(__FID)
        for(int c=0; c<nComp-1; ++c){
            add("Mu_s_"+std::to_string(c), "/mu_s_"+std::to_string(c));
            add("Mu_l_"+std::to_string(c), "/mu_l_"+std::to_string(c));
        }
        #endif

        add("Temperature", "/temperature");
        xmf << "   </Grid>\n </Domain>\n</Xdmf>\n";
        xmf.flush(); xmf.close();
        std::cout << "Saved " << xmf_name << std::endl;
    }
    MPI_Barrier(comm);
}


template<int Dim>
void IO::Read_input::write_vtk(int step, int rank, MPI_Comm comm){

    mkdir(("DATA/processor_" + std::to_string(rank)).c_str(), 0777);

    MPI_Barrier(comm);
    workers *worker = global_param->workers_mpi;
    DomainDecomp *decomp = global_param->decomp;
    OFFSETS *offsets = global_param->offsets;
    Input_parameters *param = global_param->parameters;
    
    size_t vol = worker->local_d * worker->local_h * worker->local_w;
    
    if(ht.size() != worker->num_pts * offsets->NUM_FIELDS)
        ht.resize(worker->num_pts * offsets->NUM_FIELDS);
    
    if(hi.size() != vol * offsets->NUM_FIELDS)
        hi.resize(vol * offsets->NUM_FIELDS);
    
    queue.memcpy(ht.data(), global_param->grid_data, 
                 worker->num_pts * offsets->NUM_FIELDS * sizeof(double)).wait();
    
    // Extract only the local domain data
    for(int z=0; z<worker->local_d; ++z){
        for(int y=0; y<worker->local_h; ++y){
            for(int x=0; x<worker->local_w; ++x) {
                int s = idx(z+(Dim==3?1:0), y+1, x+1, worker->ghost_h, worker->ghost_w);
                int d = idx(z, y, x, worker->local_h, worker->local_w);
                std::memcpy(&hi[d*offsets->NUM_FIELDS], &ht[s*offsets->NUM_FIELDS], offsets->NUM_FIELDS*sizeof(double));
            }
        }
    }
    
    // Generate filename
    std::ostringstream ss;
    ss << filename_ << "_" << std::setw(5) << std::setfill('0') << step;
    std::string vtk_name = "DATA/processor_" + std::to_string(rank) + "/" + ss.str() + ".vtk";
    
    // Open file for writing
    std::ofstream outfile(vtk_name);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open file " << vtk_name << std::endl;
        return;
    }
    
    int loc_d = (Dim == 2) ? 1 : worker->local_d;
    int loc_h = worker->local_h;
    int loc_w = worker->local_w;
    
    // Write VTK header
    outfile << "# vtk DataFile Version 3.0\n";
    outfile << "Phase field simulation data\n";
    outfile << "ASCII\n";
    outfile << "DATASET STRUCTURED_POINTS\n";
    outfile << "DIMENSIONS " << loc_w << " " << loc_h << " " << loc_d << "\n";
    
    // Calculate origin based on rank position
    double origin_x = 0.0, origin_y = 0.0, origin_z = 0.0;
    if constexpr (Dim == 3) {
        origin_z = decomp->coords[0] * worker->local_d;
        origin_y = decomp->coords[1] * worker->local_h;
        origin_x = decomp->coords[2] * worker->local_w;
    } else {
        origin_y = decomp->coords[0] * worker->local_h;
        origin_x = decomp->coords[1] * worker->local_w;
    }
    
    outfile << "ORIGIN " << origin_x << " " << origin_y << " " << origin_z << "\n";
    outfile << "SPACING 1.0 1.0 1.0\n";
    outfile << "POINT_DATA " << vol << "\n";
    
    int nPhase = param->NUM_PHASES;
    int nComp  = param->NUM_COMPONENTS;
    int off_phi  = offsets->OFF_PHI;
    int off_comp = offsets->OFF_COMP;
    
#ifndef __FID
    int off_mu   = offsets->OFF_MU;
#elif defined(__FID)
    int off_mus = offsets->OFF_MUS;
    int off_mul = offsets->OFF_MUL;
#endif
    int off_temp = offsets->OFF_TEMP;
    
    // Lambda to write scalar field
    auto write_scalar = [&](std::string name, int offset) {
        outfile << "SCALARS " << name << " double 1\n";
        outfile << "LOOKUP_TABLE default\n";
        //for(size_t i=0; i<vol; ++i) {
        //    outfile << hi[i*offsets->NUM_FIELDS + offset] << "\n";
        //}

        for(size_t i=0; i<worker->local_d; ++i){
            for(size_t j=0; j<worker->local_h; ++j){
                for(size_t k=0; k<worker->local_w; ++k){
                    int idx = i * worker->local_h * worker->local_w + j * worker->local_w + k;
                    outfile << hi[idx * offsets->NUM_FIELDS + offset] << "\t";
                }
                outfile << "\n";  // Newline after each row
            }
            outfile << "\n";  // Newline after each slice
        }
    };
    
    // Write phase fields
    for(int p=0; p<nPhase; ++p){
        write_scalar("phi_" + std::to_string(p), off_phi + p);
    }
    
    // Write composition fields
    for(int c=0; c<nComp-1; ++c){
        write_scalar("comp_" + std::to_string(c), off_comp + c);
    }
    
    // Write chemical potential fields
#ifndef __FID
    for(int c=0; c<nComp-1; ++c){
        write_scalar("mu_" + std::to_string(c), off_mu + c);
    }
#elif defined(__FID)
    for(int c=0; c<nComp-1; ++c){
        write_scalar("mu_s_" + std::to_string(c), off_mus + c);
        write_scalar("mu_l_" + std::to_string(c), off_mul + c);
    }
#endif
    
    // Write temperature field
    write_scalar("temperature", off_temp);
    
    outfile.close();
    MPI_Barrier(comm);
}

template <typename T>
void IO::Read_input::PrintI(const string &name, T value){
    logger.log(name+" : "+std::to_string(value));
}

void IO::Read_input::PrintS(const string &name, const std::vector<std::string> tokens){
    if(tokens.empty()){return;}
    std::string temp = name;
    for(const auto &s : tokens){
        temp += (", "+s);
    }
    logger.log(temp);
}

void IO::Read_input::print_diffusivity_matrix(){
    std::ostringstream dfm;
    dfm<<"Diffusivity matrix: \n";

    for(int a=0; a<nPhase; ++a){
        dfm<<"\nPhase "<<a<<":\n";
        for(int j=0; j<nComp-1; ++j){
            for(int k=0; k<nComp-1; ++k){
                dfm << std::setw(width)
                    <<std::scientific
                    <<std::setprecision(precision)

                    <<global_param->PFMAT->Diffusivity[a*(nComp-1)*(nComp-1)+j*(nComp-1) +k]<<"\t";
            }
            dfm<<"\n";
        }
    }
    logger.log(dfm.str());
}

/// @brief  Print 3D matrix of nPhase * nComp-1 *nComp-1
/// @param Mat 
void IO::Read_input::print_3D_matrix(double *Mat){
    std::ostringstream oss;

    for(int a=0; a<nPhase; ++a){
        oss<<"\n Phase "<<a<<":\n";
        for(int j=0; j<nComp-1; ++j){
            for(int k=0; k<nComp-1; ++k){
                oss<<std::setw(width)
                    <<std::fixed
                    <<std::setprecision(precision)

                    <<Mat[a*(nComp-1)*(nComp-1)+j*(nComp-1) +k];
            }
            oss<<"\n";
        }
    }
    logger.log(oss.str());
}

/// @brief Print 2D matrix of nPhase * nPhase
/// @param Mat 
void IO::Read_input::print_2D_matrix(double *Mat){
    std::ostringstream oss;

    for (int i = 0; i < nPhase; ++i) {
        oss << "\nPhase " << i << ":\n";
        for (int j = 0; j < nPhase; ++j) {
            oss << std::setw(width)
                << std::fixed << std::setprecision(precision)
                << Mat[i*nPhase + j];
        }
        oss << "\n";  // newline after each row
    }

    logger.log(oss.str());
}

/// @brief print the composition information in a nPhase * nPhase * nComp-1 matrix
/// @param Mat 
void IO::Read_input::print_composition_info(double *Mat){
    std::ostringstream oss;

    const int width = 10;
    const int precision = 4;

    for (int a = 0; a < nPhase; ++a) {
        oss << "\n=== Phase1 " << a << " ===\n";

        // Print component headers
        oss << std::setw(8) << "Phase2";
        for (int k = 0; k < nComp - 1; ++k) {
            oss << std::setw(width) << "Comp" + std::to_string(k);
        }
        oss << "\n";

        for (int j = 0; j < nPhase; ++j) {
            oss << std::setw(8) << j;  // Phase2 index
            for (int k = 0; k < nComp - 1; ++k) {
                oss << std::setw(width) 
                    << std::fixed << std::setprecision(precision) 
                    << Mat[a*nPhase*(nComp-1) + j*(nComp-1) + k];
            }
            oss << "\n";
        }
    }
    logger.log(oss.str());
}

const char* IO::Read_input::bc_type_to_str(BCType t) {
    switch (t) {
        case BC_PERIODIC:  return "PER";
        case BC_NEUMANN:   return "NEU";
        case BC_DIRICHLET: return "DIR";
        case BC_FLUX:      return "FLX";
        default:           return "UNK";
    }
}

const FaceBCs& IO::Read_input::get_face_bcs(const DomainBCs& bcs, int face) {
    switch (face) {
        case 0: return bcs.left;
        case 1: return bcs.right;
        case 2: return bcs.top;
        case 3: return bcs.bottom;
        case 4: return bcs.front;
        default:return bcs.back;
    }
}


void IO::Read_input::print_boundary_condition(){
    const DomainBCs &Bcs = *(global_param->bcs);

    const int width = 18;
    const char* faces[6] = {"RIGHT", "LEFT", "FRONT", "BACK", "TOP", "BOTTOM"};

    std::ostringstream out;

    int nface = global_param->domainInfo->DIMENSION == 2 ? 4 : 6;

    // ---------- HEADER ----------
    out << "\n" << std::setw(width) << "Field/Face";
    for (int f = 0; f < nface; ++f)
        out << std::setw(width) << faces[f];
    out << "\n" << std::string(width * 7, '-') << "\n";

    // ---------- PHASE FIELDS ----------
    for (int p = 0; p < nPhase; ++p) {
        out << std::setw(width) << ("phi[" + std::to_string(p) + "]");

        for (int f = 0; f < nface; ++f) {
            const FaceBCs& face = get_face_bcs(Bcs, f);
            const BCValue& bc = face.phi[p];

            std::ostringstream cell;
            cell << bc_type_to_str(bc.type);
            if (bc.type == BC_DIRICHLET || bc.type == BC_FLUX)
                cell << ":" << bc.value;

            out << std::setw(width) << cell.str();
        }
        out << "\n";
    }

    // ---------- COMPOSITION ----------
    for (int c = 0; c < nComp - 1; ++c) {
        out << std::setw(width) << ("comp[" + std::to_string(c) + "]");

        for (int f = 0; f < nface; ++f) {
            const FaceBCs& face = get_face_bcs(Bcs, f);
            const BCValue& bc = face.comp[c];

            std::ostringstream cell;
            cell << bc_type_to_str(bc.type);
            if (bc.type == BC_DIRICHLET || bc.type == BC_FLUX)
                cell << ":" << bc.value;

            out << std::setw(width) << cell.str();
        }
        out << "\n";
    }

    // ---------- CHEMICAL POTENTIAL ----------
    for (int c = 0; c < nComp - 1; ++c) {
        out << std::setw(width) << ("mu[" + std::to_string(c) + "]");

        for (int f = 0; f < nface; ++f) {
            const FaceBCs& face = get_face_bcs(Bcs, f);
            const BCValue& bc = face.mu[c];

            std::ostringstream cell;
            cell << bc_type_to_str(bc.type);
            if (bc.type == BC_DIRICHLET || bc.type == BC_FLUX)
                cell << ":" << bc.value;

            out << std::setw(width) << cell.str();
        }
        out << "\n";
    }

    // ---------- TEMPERATURE ----------
    out << std::setw(width) << "temp";

    for (int f = 0; f < nface; ++f) {
        const FaceBCs& face = get_face_bcs(Bcs, f);
        const BCValue& bc = face.temp[0];

        std::ostringstream cell;
        cell << bc_type_to_str(bc.type);
        if (bc.type == BC_DIRICHLET || bc.type == BC_FLUX)
            cell << ":" << bc.value;

        out << std::setw(width) << cell.str();
    }

    out << "\n";

    logger.log(out.str());
}


void IO::Read_input::print_infile(){
    
    Settings *flags = global_param->Flags;
    Domain   *domain = global_param->domainInfo;
    Input_parameters *INP = global_param->parameters;
    Temp_gradient *TMPG = global_param->Temp;
    Phasefield_matraix *pfmatrix = global_param->PFMAT;
    
    auto skipline = [&](){
        logger.log("\n");
    };

    auto skiplines = [&](){
        logger.log("\n\n\n");
    };


    skiplines();

    logger.log("--------------------------------------------------------------");
    logger.log("                           MODULES                              ");
    logger.log("--------------------------------------------------------------");

    logger.log(std::string("BINARY      : ") + (flags->BINARY     ? "true"    : "false"   ));
    logger.log(std::string("TERNARY     : ") + (flags->TERNARY    ? "true"    : "false"   ));
    logger.log(std::string("RESTART     : ") + (flags->restart    ? "ENABLED" : "DISABLED"));
    logger.log(std::string("ASCII       : ") + (flags->ASCII      ? "ENABLED" : "DISABLED"));
    logger.log(std::string("HDF5        : ") + (flags->HDF5       ? "ENABLED" : "DISABLED"));
    logger.log(std::string("ELASTICITY  : ") + (flags->elasticity ? "ENABLED" : "DISABLED"));
    logger.log(std::string("ISOTHERMAL  : ") + (flags->ISOTHERMAL ? "ENABLED" : "DISABLED"));
    logger.log(std::string("DILUTE      : ") + (flags->Dilute     ? "ENABLED" : "DISABLED"));
    logger.log(std::string("SHIFT       : ") + (flags->shift      ? "ENABLED" : "DISABLED"));
    logger.log(std::string("NOISE       : ") + (flags->Noise      ? "ENABLED" : "DISABLED"));
    logger.log(std::string("ANISOTROPY  : ") + (flags->anisotropy ? "ENABLED" : "DISABLED"));
    logger.log(std::string("FUNCTION_F  : ") + (std::to_string(flags->Function_F)));

    skiplines();

    logger.log("--------------------------------------------------------------");
    logger.log("                      Domain Information                      ");
    logger.log("--------------------------------------------------------------");


    PrintI(std::string("DIMENSION         "), domain->DIMENSION);
    PrintI(std::string("MESH_X            "), domain->MESH_X   );
    PrintI(std::string("MESH_Y            "), domain->MESH_Y   );
    PrintI(std::string("MESH_Z            "), domain->MESH_Z   );
    PrintI(std::string("DELTA_X           "), domain->DELTA_X  );
    PrintI(std::string("DELTA_Y           "), domain->DELTA_Y  );
    PrintI(std::string("DELTA_Z           "), domain->DELTA_Z  );
    PrintI(std::string("DELTA_T           "), domain->DELTA_T  );

    skiplines();

    logger.log("--------------------------------------------------------------");
    logger.log("                      Input Parameters                        ");
    logger.log("--------------------------------------------------------------");

    PrintS("Components          ", global_param->reference->componenets);
    PrintS("Phases              ", global_param->reference->phases     );
    PrintS("Phases TDB          ", global_param->reference->phases_tdb );
    PrintS("Phases map          ", global_param->reference->phases_map );

    skipline();

    PrintI("Equilibrium temp    ", INP->T_eq);
    PrintI("Filling temp        ", INP->Filling_temp);
    PrintI("T                   ", INP->T);
    PrintI("FOLD                ", INP->FOLD);
    PrintI("Function anisotropy ", INP->function_anisotropy);
    PrintI("Noise phase field   ", INP->Noise_phase_field);
    PrintI("Epsilon             ", INP->epsilon);
    PrintI("Tau                 ", INP->tau);
    PrintI("R                   ", INP->R);
    PrintI("V                   ", INP->Volm);
    PrintI("Obstacle            ", INP->Obstacle);
    PrintI("Function W          ", INP->Function_W);

    skiplines();

    if(!flags->ISOTHERMAL){
        logger.log("--------------------------------------------------------------");
        logger.log("                      Temperature Gradient                    ");
        logger.log("--------------------------------------------------------------");
        PrintI("Base temperature     ", TMPG->base_Temp);
        PrintI("Distance             ", TMPG->Distance);
        PrintI("Gradient offset      ", TMPG->gradient_offset);
        PrintI("Delta T              ", TMPG->DeltaT);
        PrintI("Offset               ", TMPG->offset);
        PrintI("Velocity             ", TMPG->velocity);
        
        skiplines();
    }

    logger.log("--------------------------------------------------------------");
    logger.log("                      Phasefield Matraix                      ");
    logger.log("--------------------------------------------------------------");

    skipline();

    print_diffusivity_matrix();

    skipline();

    logger.log("Gamma Matrix");
    print_2D_matrix(pfmatrix->gamma);

    skipline();

    logger.log("Tau matrix");
    print_2D_matrix(pfmatrix->tau_ab);

    skipline();

    logger.log("fab matrix");
    print_2D_matrix(pfmatrix->fab);

    if(global_param->Flags->Function_F !=4){
        skipline();

        logger.log("A matrix");
        print_3D_matrix(pfmatrix->A);
    }

    skipline();
    logger.log("Equivalent composition (Ceq)");
    print_composition_info(pfmatrix->ceq);

    skipline();
    logger.log("Filling composition (C_fill)");
    print_composition_info(pfmatrix->cfill);

    if(global_param->Flags->Function_F == 1){
        skipline();
        logger.log("Slope of slovus line");
        print_composition_info(pfmatrix->slopes);
    }

    skiplines();

    logger.log("-------------------------------------------------------------");
    logger.log("                    Boundary conditions                      ");
    logger.log("-------------------------------------------------------------");

    print_boundary_condition();

    skiplines();

}


void IO::Read_input::print_filling_file(){
    fillParameters *fill_info = global_param->Filling;

    auto skipline = [&](){
        logger.log("\n");
    };

    logger.log("--------------------------------------------------------------");
    logger.log("                         FILLING                              ");
    logger.log("--------------------------------------------------------------");

    logger.log("Number of seeds : "+std::to_string(fill_info->countFill));

    skipline();

    auto printSeed = [fill_info, this](int i) {
        int type = fill_info->fillType[i];
        std::ostringstream oss;
            if(type == FILL_CUBE){
            oss << "Filling pattern : CUBE\n"
                << "Xs : " << fill_info->xS[i] << "\n"
                << "Xe : " << fill_info->xE[i] << "\n"
                << "Ys : " << fill_info->yS[i] << "\n"
                << "Ye : " << fill_info->yE[i] << "\n"
                << "Zs : " << fill_info->zS[i] << "\n"
                << "Ze : " << fill_info->zE[i];
        }else if(type == FILL_SPHERE){
            oss << "Filling pattern : SPHERE\n"
                << "Xc : " << fill_info->xC[i] << "\n"
                << "Yc : " << fill_info->yC[i] << "\n"
                << "Zc : " << fill_info->zC[i] << "\n"
                << "Radius : " << fill_info->radius[i];
        }

        this->logger.log(oss.str());
    };


    for(int i=0; i<fill_info->countFill; ++i){
        printSeed(i);
        skipline();
    }
}
template void IO::Read_input::write_hdf5<2>(int, int, MPI_Comm);
template void IO::Read_input::write_hdf5<3>(int, int, MPI_Comm);
template void IO::Read_input::write_vtk<2>(int, int, MPI_Comm);
template void IO::Read_input::write_vtk<3>(int, int, MPI_Comm);