#pragma once

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <string>

#include <map>
#include <algorithm>

#include "utilities/Input_parameters.hpp"
#include "utilities/Logger.hpp"
#include "utilities/microsim.hpp"
#include "utilities/helper_functions.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>

namespace IO{

    class Read_input{
        private:
            std::string infile;

            std::string filling_file;

            sycl::queue &queue;

            std::string &filename_;

            Logger &logger;

            std::vector<double> ht;
            std::vector<double> hi;

            microsim::Info *global_param;

            int nPhase;
            int nComp;

            const int width = 10;
            const int precision = 4;

            void Parse_scalars(std::string &line, std::string &key, std::string &value);
            
            void get_tokens(std::string &line, std::string &key, std::string &value, std::vector<std::string> &tokens);
            
            void get_doubles(std::string &line, std::string &key, std::string &value, std::vector<double> &tokens);

            void get_bools(std::string &line, std::string &key, std::string &value, std::vector<bool> &tokens);
            
            void get_boundary(std::string &line, std::string &key, std::string &value, std::string &field, std::vector<int> &tokens);
        
            void Initialize_PF_matrix(int nPhase, int nComp);

            void populate_boundary_conditions(DomainBCs *bcs, std::vector<int> tokens, const std::string fd);
        
            void populate_boundary_values(DomainBCs *bcs, std::vector<int> tokens, const std::string fd);

            void populate_symmetric_matrix_ab(double *Mat, int nPhase, std::vector<double> tokens);

            void populate_symmetric_matrix_abc(double *Mat, int nPhase, std::vector<double> tokens);

            void populate_matrix_A(double *Mat, int nPhase, int nComp, std::vector<double> tokens);

            void populate_diffusivity(double *Mat, double *InvMat, int nPhase, int nComp, std::vector<double> tokens);

            void populate_thermodynamic_matrix(double *Mat, int nPhase, int nComp, std::vector<double> tokens);

            void populate_rotation_matrix(double *Mat, double *InvMat, int nPhase, std::vector<double> tokens);

            double** malloc_2d_host(size_t m, size_t n);

            double*** malloc_3d_host(size_t m , size_t n, size_t k);

            void free_2d_host(double** Mat, size_t m);

            void free_3d_host(double*** Mat, size_t m, size_t n);

            void multiply_2d(double **m1,double **m2,double **prod,long size);

            bool matrix_inversion_GSL(double *Mat, double *InvMat, const int N);

            void get_rotation_matrix(double **R, double theta, int axis);

            field string_to_field(const std::string &str);

            void Init_default_conditions(DomainBCs *bcs);

            // ---- Function to print the info to log file ---- //
            template <typename T>
            void PrintI(const string &name, T value);

            void PrintS(const string &name, const std::vector<std::string> tokens);

            void print_diffusivity_matrix();

            void print_3D_matrix(double *Mat);

            void print_2D_matrix(double *Mat);

            void print_composition_info(double *Mat);

            static const FaceBCs& get_face_bcs(const DomainBCs& bcs, int face);

            static const char* bc_type_to_str(BCType t);

            //void print_boundary_condition();

        public:
            Read_input(std::string Infile, std::string Filling_file, microsim::Info *Parameters, sycl::queue &dev, std::string &filename, Logger &Logger);
            
            void Read_infile();

            void Read_filling_file();

            void print_infile();

            void print_filling_file();

            void print_boundary_condition();

            void print_mpi_info();

            template<int Dim>
            void write_hdf5(int step, int rank, MPI_Comm comm);

            template<int Dim>
            void write_vtk(int step, int rank, MPI_Comm comm);

            ~Read_input();
    };
}