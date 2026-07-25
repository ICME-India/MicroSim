#pragma once

#include <sycl/sycl.hpp>
#include <vector>

#include "utilities/Input_parameters.hpp"


namespace microsim{
    class Info{
        public:
            Domain *domainInfo;
            Input_parameters *parameters;
            Settings *Flags;
            phases_components *reference;
            Temp_gradient *Temp;
            fillParameters *Filling;
            Phasefield_matraix *PFMAT;
            DomainBCs *bcs;
            workers *workers_mpi;

            DomainDecomp *decomp;
            OFFSETS *offsets;

            fields *gridfields;
            fields *gridfields_new;
            double *grid_data;
            double *grid_data_new;

            double *d_send;
            FieldRange *d_active_ranges;
            
            double *d_recv_l;
            double *d_recv_r;
            double *d_recv_t;
            double *d_recv_b;
            double *d_recv_f;
            double *d_recv_k;

            std::vector<double> h_recv_l;
            std::vector<double> h_recv_r;
            std::vector<double> h_recv_t;
            std::vector<double> h_recv_b;
            std::vector<double> h_recv_f;
            std::vector<double> h_recv_k;

            sycl::queue dev_q;

            Info(sycl::queue &q);
            void Allocate_field();

            ~Info();

        private:

            
    };
}
