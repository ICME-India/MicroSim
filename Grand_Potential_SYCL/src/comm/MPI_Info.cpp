#include "comm/MPI_Info.hpp"
#include "utilities/helper_functions.hpp"

#include <vector>
#include <cstring> // For std::memcpy


// Constructor 
mpi::mpi_comm::mpi_comm(int num_arg, char** args, microsim ::Info *g_parameters, sycl::queue &q, IO::Logger &log):
                argc(num_arg), argv(args), global_param(g_parameters), dev_q(q), logger(log) {
    
    comm_world = MPI_COMM_WORLD;
    Settings *flags = global_param->Flags;
    workers *worker_info = global_param->workers_mpi;
    Domain *domain = global_param->domainInfo;
    DomainDecomp *decomp = global_param->decomp;

    int nPhase = global_param->parameters->NUM_PHASES;
    int nComp  = global_param->parameters->NUM_COMPONENTS;

    int NUM_FIELDS;
    #ifdef __KKS
        NUM_FIELDS = nPhase + 4 * (nComp - 1) + 1;
    #elif defined(__FID)
        NUM_FIELDS = nPhase + 5 * (nComp - 1) + 1;
    #else
        NUM_FIELDS = nPhase + 2 * (nComp - 1) + 1;
    #endif

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &size);

    bool isValid = false;
    int dimention = global_param->domainInfo->DIMENSION;

    if((dimention == 2 && argc == 6) || (dimention == 3 && argc ==7)){
        decomp->Dims[0] = std::atoi(argv[4]);
        decomp->Dims[1] = atoi(argv[5]);
        decomp->Dims[2] = (dimention == 2) ? 1 : atoi(argv[6]);

        if (dimention == 2) {
            isValid = (size == decomp->Dims[0] * decomp->Dims[1]);
        } else {
            isValid = (size == decomp->Dims[0] * decomp->Dims[1] * decomp->Dims[2]);
        }
    }

    if(!isValid){
        logger.log("Invalid domain decomposition",LogLevel::ERROR);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    int periodic[3] = {1, 1, 1};
    for (int i = 0; i < dimention; ++i){
        std::cout<<"Periodicity in dimension " << i << " : " << flags->periodic[i] << std::endl;
        decomp->periods[i] = 1.0;//flags->periodic[i] ? 1 : 0;
    }

    if(dimention == 2){periodic[0] = decomp->periods[1]; periodic[1] = decomp->periods[0];} // Mapping 2D Y/X
    else{periodic[0] = decomp->periods[2]; periodic[1] = decomp->periods[1]; periodic[2] = decomp->periods[0];} // Mapping 3D Z/Y/X

    MPI_Cart_create(MPI_COMM_WORLD, dimention,  decomp->Dims, periodic, 1, &comm_cart);
    MPI_Cart_coords(comm_cart, rank, dimention, decomp->coords);
    
    worker_info->front=MPI_PROC_NULL;
    worker_info->back=MPI_PROC_NULL;

    if(dimention == 2) {
        MPI_Cart_shift(comm_cart, 0, 1, &worker_info->up,   &worker_info->down);
        MPI_Cart_shift(comm_cart, 1, 1, &worker_info->left, &worker_info->right);
    } else {
        MPI_Cart_shift(comm_cart, 0, 1, &worker_info->front, &worker_info->back);
        MPI_Cart_shift(comm_cart, 1, 1, &worker_info->up,    &worker_info->down);
        MPI_Cart_shift(comm_cart, 2, 1, &worker_info->left,  &worker_info->right);
    }

    int idx_x = (dimention == 3) ? 2 : 1;
    int idx_y = (dimention == 3) ? 1 : 0;
    int idx_z = 0; // Only used if dimention == 3

    // Determine physical edges based on coordinates
    worker_info->is_left_edge  = (decomp->coords[idx_x] == 0);
    worker_info->is_right_edge = (decomp->coords[idx_x] == decomp->Dims[idx_x] - 1);
    
    worker_info->is_top_edge   = (decomp->coords[idx_y] == 0);
    worker_info->is_bot_edge   = (decomp->coords[idx_y] == decomp->Dims[idx_y] - 1);

    // Front/Back only exist in 3D
    worker_info->is_front_edge = (dimention == 3) && (decomp->coords[idx_z] == 0);
    worker_info->is_back_edge  = (dimention == 3) && (decomp->coords[idx_z] == decomp->Dims[idx_z] - 1);

    MPI_Comm_split_type(comm_world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shm);
    MPI_Comm_rank(comm_shm, &shm_rank);
    MPI_Comm_size(comm_shm, &shm_size);

    worker_info->local_d = (dimention==3) ? domain->MESH_Z / decomp->Dims[0] : 1;
    worker_info->local_h = domain->MESH_Y / decomp->Dims[(dimention==3)?1:0];
    worker_info->local_w = domain->MESH_X / decomp->Dims[(dimention==3)?2:1];
    
    worker_info->ghost_d = worker_info->local_d + 2;
    worker_info->ghost_h = worker_info->local_h + 2;
    worker_info->ghost_w = worker_info->local_w + 2;
    worker_info->num_pts = worker_info->ghost_d * worker_info->ghost_h * worker_info->ghost_w;

    worker_info->sz_yz = (size_t)worker_info->local_d * worker_info->local_h * NUM_FIELDS;
    worker_info->sz_xz = (size_t)worker_info->local_d * worker_info->local_w * NUM_FIELDS;
    worker_info->sz_xy = (size_t)worker_info->local_h * worker_info->local_w * NUM_FIELDS;
    
    worker_info->off_l=0;
    worker_info->off_r=worker_info->off_l+worker_info->sz_yz;
    worker_info->off_t=worker_info->off_r+worker_info->sz_yz;
    worker_info->off_b=worker_info->off_t+worker_info->sz_xz;
    worker_info->off_f=worker_info->off_b+worker_info->sz_xz;
    worker_info->off_bk=worker_info->off_f+worker_info->sz_xy;

    worker_info->total_send_bytes = (dimention == 2) ? (worker_info->off_b + worker_info->sz_xz) * sizeof(double) : (worker_info->off_bk + worker_info->sz_xy) * sizeof(double);

    MPI_Win_allocate_shared(worker_info->total_send_bytes, sizeof(double), MPI_INFO_NULL, comm_shm, &my_send_buf_base, &win_send);
    neighbor_ptr_cache.resize(size, nullptr);

    // This will allocate the send buffer inside mpi shared memory (Win_allocate_shared)
    #ifdef ENABLE_CPU_SHARED_MPI
        global_param->d_send = my_send_buf_base;
    #endif

    MPI_Group wg, sg;
    MPI_Comm_group(comm_world, &wg);
    MPI_Comm_group(comm_shm, &sg);
    std::vector<int> sr(shm_size), wr(shm_size);

    for(int i=0; i<shm_size; ++i) sr[i] = i;   
    MPI_Group_translate_ranks(sg, shm_size, sr.data(), wg, wr.data());

    for(int i=0; i<shm_size; ++i) {
        MPI_Aint sz;
        int unit;
        double* base;
        
        MPI_Win_shared_query(win_send, i, &sz, &unit, &base);
        neighbor_ptr_cache[wr[i]] = base;
    }
    
    MPI_Group_free(&wg);
    MPI_Group_free(&sg);
}


// Distructor
mpi::mpi_comm::~mpi_comm() {
    if(win_send != MPI_WIN_NULL) MPI_Win_free(&win_send);
    if(comm_shm != MPI_COMM_NULL)MPI_Comm_free(&comm_shm);
}

bool mpi::mpi_comm::is_local(int t){
    return (t >= 0 && t < size) && neighbor_ptr_cache[t];
}

double* mpi::mpi_comm::get_my_send_ptr(size_t offset){
    return my_send_buf_base + offset;
}

void mpi::mpi_comm::copy_from_neighbor(int t, size_t off, double* d, size_t c){
    std::memcpy(d, neighbor_ptr_cache[t] + off, c * sizeof(double));
}

void mpi::mpi_comm::node_barrier(){
    MPI_Barrier(comm_shm);
}

template<int Dim>
void mpi::mpi_comm::apply_physical_bcs(double *data, int offset_fields, int num_fields){
    
    workers *worker   = global_param->workers_mpi;
    OFFSETS *offset_l = global_param->offsets;
    DomainBCs *bc     = global_param->bcs;
    Domain *domain    = global_param->domainInfo;

    //FieldRange *d_active_range = global_param->d_active_ranges;
    //int num_ranges = ranges.size();
    
    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    dev_q.submit([&](sycl::handler &h){
        sycl::stream out(1024, 256, h);
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int z=c.z, y=c.y, x=c.x;

            bool is_ghost = (x==0 || x==worker->local_w+1 || y==0 || y==worker->local_h+1);
            if constexpr (Dim==3) is_ghost |= (z==0 || z==worker->local_d+1);
            if (!is_ghost) return;

            int current_idx = idx(z, y, x, worker->ghost_h, worker->ghost_w);
            int inner_idx = -1;

            const FaceBCs* active_bc = nullptr;

            if(worker->is_left_edge && x == 0){
                active_bc = &bc->left;
                inner_idx = idx(z, y, 1, worker->ghost_h, worker->ghost_w);
            }
            if(worker->is_right_edge && x == worker->local_w+1){
                active_bc = &bc->right;
                inner_idx = idx(z, y, worker->local_w, worker->ghost_h, worker->ghost_w);
            }
            if(worker->is_top_edge && y == 0){
                active_bc = &bc->top;
                inner_idx = idx(z, 1, x, worker->ghost_h, worker->ghost_w);
            }
            if(worker->is_bot_edge && y == worker->local_h+1){
                active_bc = &bc->bottom;
                inner_idx = idx(z, worker->local_h, x, worker->ghost_h, worker->ghost_w);
            }
            if constexpr(Dim == 3){
                if(worker->is_front_edge && z == 0){
                    active_bc = &bc->front;
                    inner_idx = idx(1, y, x, worker->ghost_h, worker->ghost_w);
                }
                if(worker->is_back_edge && z == worker->local_d+1){
                    active_bc = &bc->back;
                    inner_idx = idx(worker->local_d, y, x, worker->ghost_h, worker->ghost_w);
                }
            }

            if(active_bc && inner_idx != -1){
                for(int k=offset_fields; k<offset_fields + num_fields; ++k){

                    // RESOLVE FIELD BC
                    BCValue conf = {BC_NEUMANN, 0.0};

                    #ifdef __KKS
                    if (k >= offset_l->OFF_PHI && k < offset_l->OFF_COMP) {
                        conf = active_bc->phi[k - offset_l->OFF_PHI];
                    } else if (k >= offset_l->OFF_COMP && k < offset_l->OFF_MU) {
                        conf = active_bc->comp[k - offset_l->OFF_COMP];
                    } else if (k >= offset_l->OFF_MU && k < offset_l->OFF_TEMP) {
                        conf = active_bc->mu[k- offset_l->OFF_MU];
                    } else if (k == offset_l->OFF_TEMP) {
                        conf = active_bc->temp[0];
                    }
                    
                    #elif defined(__FID)
                    if (k >= offset_l->OFF_PHI && k < offset_l->OFF_COMP) {
                        conf = active_bc->phi[k - offset_l->OFF_PHI];
                    } else if (k >= offset_l->OFF_COMP && k < offset_l->OFF_MUS) {
                        conf = active_bc->comp[k - offset_l->OFF_COMP];
                    } else if (k >= offset_l->OFF_MUS && k < offset_l->OFF_TEMP) {
                        conf = active_bc->mu[k - offset_l->OFF_MUS];
                    } else if (k == offset_l->OFF_TEMP) {
                        conf = active_bc->temp[0];
                    }
                    
                    #else
                    if (k >= offset_l->OFF_PHI && k < offset_l->OFF_COMP) {
                        //out<<"Appling BC for phi field "<< k - offset_l->OFF_PHI << sycl::endl;
                        conf = active_bc->phi[k - offset_l->OFF_PHI];
                    } else if (k >= offset_l->OFF_COMP && k < offset_l->OFF_MU) {
                        //out<<"Appling BC for comp field "<< k - offset_l->OFF_COMP << sycl::endl;
                        conf = active_bc->comp[k - offset_l->OFF_COMP];
                    } else if (k >= offset_l->OFF_MU && k < offset_l->OFF_TERM1) {
                        conf = active_bc->mu[k - offset_l->OFF_MU];
                    } else if (k == offset_l->OFF_TEMP) {
                        conf = active_bc->temp[0];
                    }
                    #endif

                    // APPLY
                    if (conf.type == BC_PERIODIC) continue;

                    double val_inner = data[inner_idx*offset_l->NUM_FIELDS + k];
                    double res = val_inner;

                    if (conf.type == BC_DIRICHLET) {
                        res = conf.value;
                    } else if (conf.type == BC_FLUX) {
                        double del = 0.0;

                        if(active_bc == &bc->left || active_bc == &bc->right)
                            del = domain->DELTA_X;
                        else if(active_bc == &bc->top || active_bc == &bc->bottom)
                            del = domain->DELTA_Y;
                        else if constexpr(Dim == 3){
                            del = domain->DELTA_Z;
                        }

                        res = val_inner - conf.value * del; 
                    } else {
                        if(id[0]==1&&id[1]==0){
                        }
                        res = val_inner;
                    }
                    data[current_idx*offset_l->NUM_FIELDS + k] = res;
                }
            }
        });
    }).wait();

    MPI_Barrier(comm_world);

}


template<int Dim>
void mpi::mpi_comm::exchange_boundary(std::vector<FieldRange> ranges){
    FieldRange *d_active_range = global_param->d_active_ranges;
    workers *worker   = global_param->workers_mpi;
    OFFSETS *offset_l = global_param->offsets;
    double *dsent     = global_param->d_send;
    double *d_data    = global_param->grid_data;

    double *d_recv_l = global_param->d_recv_l;
    double *d_recv_r = global_param->d_recv_r;
    double *d_recv_t = global_param->d_recv_t;
    double *d_recv_b = global_param->d_recv_b;
    double *d_recv_f = global_param->d_recv_f;
    double *d_recv_k = global_param->d_recv_k;
    

    int num_ranges = ranges.size();
    dev_q.memcpy(d_active_range, ranges.data(), num_ranges * sizeof(FieldRange)).wait();
    
    size_t active_doubles = 0; 
    for(auto &r : ranges){
        active_doubles += r.count;
    }

    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    
    // PACKING KERNEL
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            // Bounds check - ensure valid threads
            if(c.x < 1 || c.x > worker->local_w || c.y < 1 || c.y > worker->local_h) return;
            if constexpr (Dim==3) { if(c.z < 1 || c.z > worker->local_d) return; }

            auto pack = [&](int db, int gz, int gy, int gx, int pi) {
                int s_idx = idx(gz, gy, gx, worker->ghost_h, worker->ghost_w);
                int cur_off = 0;
                for(int r=0; r<num_ranges; ++r) {
                    FieldRange rng = d_active_range[r];
                    for(int k=0; k<rng.count; ++k){
                        dsent[db + pi*active_doubles + cur_off + k] = d_data[s_idx*offset_l->NUM_FIELDS + rng.start_idx + k];
                    }
                    cur_off += rng.count;
                }
            };

            int z=c.z, y=c.y, x=c.x;
            int z0 = (Dim==3) ? (z-1) : 0; int y0 = y-1; int x0 = x-1;

            if(x==1){
                pack(worker->off_l, z, y, 1, z0*worker->local_h+y0);
            }
            if(x==worker->local_w){
                pack(worker->off_r, z, y, worker->local_w, z0*worker->local_h+y0);
            }
            if(y==1){
                pack(worker->off_t, z, 1, x, z0*worker->local_w+x0);
            }
            if(y==worker->local_h){
                pack(worker->off_b, z, worker->local_h, x, z0*worker->local_w+x0);
            }
            if constexpr (Dim==3) {
                if(z==1){
                    pack(worker->off_f, 1, y, x, y0*worker->local_w+x0);
                }
                if(z==worker->local_d){
                    pack(worker->off_bk, worker->local_d, y, x, y0*worker->local_w+x0);
                }
            }
        });
    }).wait();

    // -----------------------------------------------------------
    // IMPORTANT: Barrier/Wait to ensure GPU has finished packing
    // before MPI (Host) tries to read 'dsent'.
    // -----------------------------------------------------------
    dev_q.wait(); 

    #ifdef ENABLE_GPU_AWARE_MPI
    // === GPU AWARE MPI (RDMA) ===
    node_barrier();
    
    std::vector<MPI_Request> rq;
    
    // Helper with TAG support to prevent crosstalk
    auto gpu_ex = [&](int n, size_t om, double* dst_ptr, size_t pts, int tag){
        if(n == MPI_PROC_NULL) return; 
        size_t cnt = pts * active_doubles;
        
        MPI_Request r1, r2;
        MPI_Isend(dsent + om, cnt, MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &r1);
        MPI_Irecv(dst_ptr,    cnt, MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &r2);
        rq.push_back(r1); rq.push_back(r2);
    };

    size_t nyz=worker->local_d*worker->local_h;
    size_t nxz=worker->local_d*worker->local_w;
    size_t nxy=worker->local_h*worker->local_w;

    // Use DISTINCT tags for different directions
    // Tag 0: Left/Right
    gpu_ex(worker->left,  worker->off_l, d_recv_l, nyz, 0);
    gpu_ex(worker->right, worker->off_r, d_recv_r, nyz, 0);
    
    // Tag 1: Up/Down
    gpu_ex(worker->up,    worker->off_t, d_recv_t, nxz, 1);
    gpu_ex(worker->down,  worker->off_b, d_recv_b, nxz, 1);
    
    // Tag 2: Front/Back
    if constexpr (Dim==3){
        gpu_ex(worker->front, worker->off_f,  d_recv_f, nxy, 2);
        gpu_ex(worker->back,  worker->off_bk, d_recv_k, nxy, 2);
    }
    
    if(!rq.empty()){
        MPI_Waitall(rq.size(), rq.data(), MPI_STATUS_IGNORE);
    }

    #elif defined ENABLE_CPU_SHARED_MPI
    // === CPU SHARED MEMORY PATH ===
    node_barrier();
    std::vector<MPI_Request> rq;
    
    auto cpu_ex = [&](int n, size_t om, size_t on, double* dst_ptr, size_t pts, int tag){
        if(n == MPI_PROC_NULL) return; 
        size_t cnt = pts * active_doubles;
        
        if(is_local(n)){
            copy_from_neighbor(n, on, dst_ptr, cnt);
            return;
        }
        
        MPI_Request r1, r2;
        MPI_Isend(dsent + om, cnt, MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &r1);
        MPI_Irecv(dst_ptr,    cnt, MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &r2);
        rq.push_back(r1); rq.push_back(r2);
    };

    size_t nyz=worker->local_d*worker->local_h;
    size_t nxz=worker->local_d*worker->local_w;
    size_t nxy=worker->local_h*worker->local_w;

    cpu_ex(worker->left,  worker->off_l, worker->off_r, d_recv_l, nyz, 0);
    cpu_ex(worker->right, worker->off_r, worker->off_l, d_recv_r, nyz, 0);
    cpu_ex(worker->up,    worker->off_t, worker->off_b, d_recv_t, nxz, 1);
    cpu_ex(worker->down,  worker->off_b, worker->off_t, d_recv_b, nxz, 1);
    if constexpr (Dim==3){
        cpu_ex(worker->front, worker->off_f,  worker->off_bk, d_recv_f, nxy, 2);
        cpu_ex(worker->back,  worker->off_bk, worker->off_f, d_recv_k, nxy, 2);
    }
    
    if(!rq.empty()){
        MPI_Waitall(rq.size(), rq.data(), MPI_STATUS_IGNORE);
    }

    #else

    // === STANDARD MPI (HOST BUFFER) PATH ===
    dev_q.wait();
    dev_q.memcpy(get_my_send_ptr(0), dsent, worker->total_send_bytes).wait();
    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<MPI_Request> rq;
    auto ex = [&](int n, size_t om, size_t on, std::vector<double>& hr, size_t pts, int tag) {
        if(n == MPI_PROC_NULL) return; 
        size_t cnt = pts * active_doubles;
        
        // Safety: Resize buffer if needed
        if(hr.size() < cnt) hr.resize(cnt);

        if(is_local(n)){
            copy_from_neighbor(n, on, hr.data(), cnt);
        }
        else {
            MPI_Request r1, r2;
            MPI_Isend(get_my_send_ptr(om), cnt, MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &r1);
            MPI_Irecv(hr.data(),           cnt, MPI_DOUBLE, n, tag, MPI_COMM_WORLD, &r2);
            rq.push_back(r1); rq.push_back(r2);
        }
    };

    size_t nyz = worker->local_d * worker->local_h;
    size_t nxz = worker->local_d * worker->local_w;
    size_t nxy = worker->local_h * worker->local_w;

    ex(worker->left,  worker->off_l, worker->off_r, global_param->h_recv_l, nyz, 0);
    ex(worker->right, worker->off_r, worker->off_l, global_param->h_recv_r, nyz, 0);
    ex(worker->up,    worker->off_t, worker->off_b, global_param->h_recv_t, nxz, 1);
    ex(worker->down,  worker->off_b, worker->off_t, global_param->h_recv_b, nxz, 1);


    if constexpr(Dim ==3){
        ex(worker->front, worker->off_f,  worker->off_bk, global_param->h_recv_f, nxy, 2);
        ex(worker->back,  worker->off_bk, worker->off_f,  global_param->h_recv_k, nxy, 2);
    }
    
    if(!rq.empty()) MPI_Waitall(rq.size(), rq.data(), MPI_STATUS_IGNORE);

    size_t sd = sizeof(double);
    // Copy Host -> Device
    if(worker->left!=MPI_PROC_NULL){
        dev_q.memcpy(d_recv_l, global_param->h_recv_l.data(), nyz*active_doubles*sd);
    }
    if(worker->right!=MPI_PROC_NULL){
        dev_q.memcpy(d_recv_r, global_param->h_recv_r.data(), nyz*active_doubles*sd);
    }
    if(worker->up!=MPI_PROC_NULL){
        dev_q.memcpy(d_recv_t, global_param->h_recv_t.data(), nxz*active_doubles*sd);
    }
    if(worker->down!=MPI_PROC_NULL){
        dev_q.memcpy(d_recv_b, global_param->h_recv_b.data(), nxz*active_doubles*sd);
    }
    if constexpr (Dim==3) {
        if(worker->front!=MPI_PROC_NULL){
            dev_q.memcpy(d_recv_f, global_param->h_recv_f.data(), nxy*active_doubles*sd);
        }
        if(worker->back!=MPI_PROC_NULL){
            dev_q.memcpy(d_recv_k, global_param->h_recv_k.data(), nxy*active_doubles*sd);
        }
    }
    #endif

    // -----------------------------------------------------------
    // IMPORTANT: Barrier/Wait to ensure Recv Buffers are ready 
    // on Device before Unpacking starts.
    // -----------------------------------------------------------
    dev_q.wait();

    // UNPACKING KERNEL
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int z=c.z, y=c.y, x=c.x;
            int z0 = (Dim==3) ? (z-1) : 0; int y0 = y-1; int x0 = x-1;
            
            auto unp = [&](double* s, int bi, int gz, int gy, int gx) {
                int d_idx = idx(gz, gy, gx, worker->ghost_h, worker->ghost_w);
                int cur_off = 0;
                for(int r=0; r<num_ranges; ++r) {
                    FieldRange rng = d_active_range[r];
                    for(int k=0; k<rng.count; ++k){
                        d_data[d_idx*offset_l->NUM_FIELDS + rng.start_idx + k] = s[bi*active_doubles + cur_off + k];
                    }
                    cur_off += rng.count;
                }
            };

            bool in_y = (y >= 1 && y <= worker->local_h);
            bool in_z = (Dim==3) ? (z >= 1 && z <= worker->local_d) : true;
            bool in_x = (x >= 1 && x <= worker->local_w);

            if(in_y && in_z) {
                if(x==1 && worker->left!=MPI_PROC_NULL){
                    unp(d_recv_l, z0*worker->local_h+y0, z, y, 0);
                }
                if(x==worker->local_w && worker->right!=MPI_PROC_NULL){
                    unp(d_recv_r, z0*worker->local_h+y0, z, y, worker->local_w+1);
                }
            }
            if(in_x && in_z) {
                if(y==1 && worker->up!=MPI_PROC_NULL){
                    unp(d_recv_t, z0*worker->local_w+x0, z, 0, x);
                }
                if(y==worker->local_h && worker->down!=MPI_PROC_NULL){
                    unp(d_recv_b, z0*worker->local_w+x0, z, worker->local_h+1, x);
                }
            }
            if constexpr (Dim==3) {
                if(in_x && in_y) {
                    if(z==1 && worker->front!=MPI_PROC_NULL){
                        unp(d_recv_f, y0*worker->local_w+x0, 0, y, x);
                    }
                    if(z==worker->local_d && worker->back!=MPI_PROC_NULL){
                        unp(d_recv_k, y0*worker->local_w+x0, worker->local_d+1, y, x);
                    }
                }
            }

        });
    }).wait();


    //apply_bcs bcs_ptr[4] = {nullptr, nullptr, &mpi::mpi_comm::apply_physical_bcs<2>, &mpi::mpi_comm::apply_physical_bcs<3>};
    
    //(this->*bcs_ptr[domain->DIMENSION])(d_data);
    for(int r=0; r<num_ranges; ++r){
        FieldRange rng = ranges[r];
        //(this->*bcs_ptr[domain->DIMENSION])(d_data, rng.start_idx, rng.count);
        this->apply_physical_bcs<Dim>(d_data, rng.start_idx, rng.count);
    }

}

template void mpi::mpi_comm::exchange_boundary<2>(std::vector<FieldRange> ranges);
template void mpi::mpi_comm::exchange_boundary<3>(std::vector<FieldRange> ranges);
template void mpi::mpi_comm::apply_physical_bcs<2>(double *data, int offset_fields, int num_fields);
template void mpi::mpi_comm::apply_physical_bcs<3>(double *data, int offset_fields, int num_fields);
//using apply_bcs = void (mpi::mpi_comm::mpi_comm::*)(double*, int, int);