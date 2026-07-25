#include "utilities/microsim.hpp"
#include <new>

namespace microsim{

    Info::Info(sycl::queue &q):dev_q(q){
        domainInfo  = sycl::malloc_shared<Domain>(1, dev_q);
        parameters  = sycl::malloc_shared<Input_parameters>(1,dev_q);
        Flags       = sycl::malloc_shared<Settings>(1, dev_q);
        new (Flags) Settings();

        reference   = new struct phases_components;
        Temp        = sycl::malloc_shared<Temp_gradient>(1,dev_q);
        Filling     = sycl::malloc_shared<fillParameters>(1,dev_q);
        PFMAT       = sycl::malloc_shared<Phasefield_matraix>(1, dev_q);
        workers_mpi = sycl::malloc_shared<workers>(1, dev_q);
        bcs         = sycl::malloc_shared<DomainBCs>(1, dev_q);
        offsets     = sycl::malloc_shared<OFFSETS>(1, dev_q);
        decomp      = sycl::malloc_shared<DomainDecomp>(1, dev_q);
        new (decomp) DomainDecomp();
    }

    void Info::Allocate_field(){
        // Calculate the number of fields based on the physics model
        int nPhase = parameters->NUM_PHASES;
        int nComp  = parameters->NUM_COMPONENTS;

        offsets->OFF_PHI   = 0;
        
        #ifdef __KKS
        offsets->OFF_COMP     = offsets->OFF_PHI      + nPhase;
        offsets->OFF_C_PHASE  = offsets->OFF_COMP     + (nComp -1);
        offsets->OFF_MU       = offsets->OFF_C_PHASE  + nPhase * (nComp -1);
        offsets->OFF_DELTAPHI = offsets->OFF_MU       + (nComp -1);
        offsets->OFF_TEMP     = offsets->OFF_DELTAPHI + nPhase;
        offsets->NUM_FIELDS   = offsets->OFF_TEMP     + 1;

        #elif defined(__FID)
        offsets->OFF_COMP     = offsets->OFF_PHI  +  nPhase;
        offsets->OFF_CS       = offsets->OFF_COMP + (nComp -1);
        offsets->OFF_CL       = offsets->OFF_CS   + (nComp -1);
        offsets->OFF_MUS      = offsets->OFF_CL   + (nComp -1);
        offsets->OFF_MUL      = offsets->OFF_MUS  + (nComp -1);
        offsets->OFF_DELTAPHI = offsets->OFF_MUL  + (nComp -1);
        offsets->OFF_TEMP     = offsets->OFF_MUL  + nPhase;
        offsets->NUM_FIELDS = offsets->OFF_TEMP + 1;

        #else
        offsets->OFF_COMP  = offsets->OFF_PHI   + nPhase;
        offsets->OFF_MU    = offsets->OFF_COMP  + (nComp -1);
        offsets->OFF_TERM1   = offsets->OFF_MU    + (nComp -1);
        offsets->OFF_TERM2   = offsets->OFF_TERM1   + nPhase;
        offsets->OFF_TERM3   = offsets->OFF_TERM2   + nPhase;
        offsets->OFF_DELTAPHI = offsets->OFF_TERM3   + nPhase;
        offsets->OFF_TEMP  = offsets->OFF_DELTAPHI   + nPhase;
        offsets->NUM_FIELDS = offsets->OFF_TEMP + 1;

        #endif

        gridfields             = sycl::malloc_device<fields>(workers_mpi->num_pts, dev_q);
        gridfields_new         = sycl::malloc_device<fields>(workers_mpi->num_pts, dev_q);
        grid_data              = sycl::malloc_device<double>(workers_mpi->num_pts * offsets->NUM_FIELDS, dev_q);
        grid_data_new          = sycl::malloc_device<double>(workers_mpi->num_pts * offsets->NUM_FIELDS, dev_q);

        d_active_ranges = sycl::malloc_device<FieldRange>(offsets->NUM_FIELDS, dev_q);

        #ifdef CPU_H
            // need to commend later for device aware mpi - later in the mpi_info.cpp line 116 this dsend will be allocated inside the mpi shared memory
            #ifndef ENABLE_CPU_SHARED_MPI
            d_send = sycl::malloc_shared<double>(workers_mpi->total_send_bytes/sizeof(double), dev_q);
            #endif

            d_recv_l = sycl::malloc_shared<double>(workers_mpi->sz_yz, dev_q);
            d_recv_r = sycl::malloc_shared<double>(workers_mpi->sz_yz, dev_q);
            d_recv_t = sycl::malloc_shared<double>(workers_mpi->sz_xz, dev_q);
            d_recv_b = sycl::malloc_shared<double>(workers_mpi->sz_xz, dev_q);
            d_recv_f = sycl::malloc_shared<double>(workers_mpi->sz_xy, dev_q);
            d_recv_k = sycl::malloc_shared<double>(workers_mpi->sz_xy, dev_q);
        #endif

        #ifdef GPU
            d_send = sycl::malloc_device<double>(workers_mpi->total_send_bytes/sizeof(double), dev_q);

            d_recv_l = sycl::malloc_device<double>(workers_mpi->sz_yz, dev_q);
            d_recv_r = sycl::malloc_device<double>(workers_mpi->sz_yz, dev_q);
            d_recv_t = sycl::malloc_device<double>(workers_mpi->sz_xz, dev_q);
            d_recv_b = sycl::malloc_device<double>(workers_mpi->sz_xz, dev_q);
            d_recv_f = sycl::malloc_device<double>(workers_mpi->sz_xy, dev_q);
            d_recv_k = sycl::malloc_device<double>(workers_mpi->sz_xy, dev_q);
        #endif
        
        #ifndef ENABLE_GPU_AWARE_MPI
        #ifndef ENABLE_CPU_SHARED_MPI
            // Host buffers for MPI communication
            h_recv_l.resize(workers_mpi->sz_yz);
            h_recv_r.resize(workers_mpi->sz_yz);
            h_recv_t.resize(workers_mpi->sz_xz);
            h_recv_b.resize(workers_mpi->sz_xz);
            h_recv_f.resize(workers_mpi->sz_xy);
            h_recv_k.resize(workers_mpi->sz_xy);
        #endif
        #endif

        auto init_ptrs = [&](double* dat, fields* fs) {
            OFFSETS offsets_l = *offsets;
            dev_q.submit([&](sycl::handler& h) {
                h.parallel_for(sycl::range<1>(workers_mpi->num_pts), [=](sycl::id<1> i) {
                    long id = i[0];
                    double* b = dat + id * offsets_l.NUM_FIELDS;

                    #ifdef __KKS
                    fs[id] = {b+offsets_l.OFF_PHI, b+offsets_l.OFF_COMP, b+offsets_l.OFF_C_PHASE, b+offsets_l.OFF_MU, b+offsets_l.OFF_DELTAPHI, b+offsets_l.OFF_TEMP};

                    #elif defined(__FID)
                    fs[id] = {b+offsets_l.OFF_PHI, b+offsets_l.OFF_COMP, b+offsets_l.OFF_CS, b+offsets_l.OFF_CL, b+offsets_l.OFF_MUS, b+offsets_l.OFF_MUL, b+offsets_l.OFF_DELTAPHI, b+offsets_l.OFF_TEMP};

                    #else
                    fs[id] = {b+offsets_l.OFF_PHI, b+offsets_l.OFF_COMP, b+offsets_l.OFF_MU, b+offsets_l.OFF_TERM1, b+offsets_l.OFF_TERM2, b+offsets_l.OFF_TERM3, b+offsets_l.OFF_DELTAPHI, b+offsets_l.OFF_TEMP};
                    #endif

                    
                });
            }).wait();
        };
        init_ptrs(grid_data, gridfields);
        init_ptrs(grid_data_new, gridfields_new);
    }


    Info::~Info(){

        std::cout<<"Freeing microsim::Info allocated memory\n";
        sycl::free(domainInfo, dev_q);
        sycl::free(parameters, dev_q);
        sycl::free(Flags, dev_q);
        delete reference;
        sycl::free(Temp, dev_q);
        sycl::free(Filling, dev_q);
        sycl::free(PFMAT, dev_q);
        sycl::free(workers_mpi, dev_q);
        sycl::free(bcs, dev_q);
        sycl::free(offsets, dev_q);
        sycl::free(decomp, dev_q);

        sycl::free(gridfields, dev_q);
        sycl::free(gridfields_new, dev_q);
        sycl::free(grid_data, dev_q);
        sycl::free(grid_data_new, dev_q);

        sycl::free(d_send, dev_q);
        sycl::free(d_active_ranges, dev_q);

        sycl::free(d_recv_l, dev_q);
        sycl::free(d_recv_r, dev_q);
        sycl::free(d_recv_t, dev_q);
        sycl::free(d_recv_b, dev_q);
        sycl::free(d_recv_f, dev_q);
        sycl::free(d_recv_k, dev_q);
    }
}