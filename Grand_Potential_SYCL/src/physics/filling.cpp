#include "physics/filling.hpp"
#include "utilities/helper_functions.hpp"
#include "thermo/Function_F.hpp"
#include "thermo/Function_F4.hpp"

#include <sycl/sycl.hpp>

using null_func = void (microsim::filling::*)();
using fill_func = void (microsim::filling::*)(int, int);

microsim::filling::filling(microsim::Info *global_param, sycl::queue &q, IO::Logger &logg):
    global_fd(global_param), logger(logg), dev_q(q) {
    
    nPhase = global_fd->parameters->NUM_PHASES;
    nComp  = global_fd->parameters->NUM_COMPONENTS;
    Dimension = global_fd->domainInfo->DIMENSION;
}

void microsim::filling::fill_domain(){
    fillParameters *Fill = this->global_fd->Filling;

    null_func init_zero[] = {nullptr, nullptr, &microsim::filling::Initialize_zero<2>, &microsim::filling::Initialize_zero<3>};
    (this->*init_zero[Dimension])();

    for(int count =0; count< Fill->countFill; count++){
        int phase = Fill->phase[count];
        
        fill_func func_ptr[4]{nullptr, nullptr,nullptr, nullptr};

        if(Fill->fillType[count] == fill::FILL_CUBE){
            func_ptr[2] = &microsim::filling::fill_cube<2>;
            func_ptr[3] = &microsim::filling::fill_cube<3>;;
        }else if(Fill->fillType[count] == fill::FILL_SPHERE){
            func_ptr[2] = &microsim::filling::fill_sphere<2>;
            func_ptr[3] = &microsim::filling::fill_sphere<3>;
        }

        (this->*func_ptr[Dimension])(phase, count);
        (this->*func_ptr[Dimension])(nPhase-1, count);
    }

    std::cout<<"Filling phase fields completed. \n";

    null_func fill_temp[] = {nullptr, nullptr, &microsim::filling::fill_temperature<2>, &microsim::filling::fill_temperature<3>};
    (this->*fill_temp[Dimension])();

    std::cout<<"Filling temperature field completed. \n";

    null_func fill_comp[] = {nullptr, nullptr, &microsim::filling::fill_composition<2>, &microsim::filling::fill_composition<3>};
    (this->*fill_comp[Dimension])();

    std::cout<<"Filling composition field completed. \n";

    dev_q.memcpy(global_fd->grid_data_new, global_fd->grid_data, sizeof(double) * global_fd->workers_mpi->num_pts * global_fd->offsets->NUM_FIELDS).wait();

    logger.log("completed filling opetations", LogLevel::INFO);
}

template<int Dim>
void microsim::filling::Initialize_zero(){
    workers *worker   = global_fd->workers_mpi;
    OFFSETS *offset_l = global_fd->offsets;
    double *grid_data = global_fd->grid_data;
    fields *gridfield = global_fd->gridfields;
    int nPhase = this->nPhase;

    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

            for(int f =0; f< offset_l->NUM_FIELDS; f++){
                grid_data[idx_g * offset_l->NUM_FIELDS + f] = 0.0;    
            }

            gridfield[idx_g].phia[nPhase-1] = 1.0;
        });
    }).wait();
}

template<int Dim>
void microsim::filling::fill_cube(int phase, int count){
    workers *worker      = global_fd->workers_mpi;
    fillParameters *Fill = global_fd->Filling;
    fields *gridfield    = global_fd->gridfields;
    DomainDecomp *Decomp = global_fd->decomp;
    int nPhase = this->nPhase;

    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

            int gz = (Dim==3) ? (c.z-1) + Decomp->coords[0]*worker->local_d : 0;
            int gy = (c.y-1) + Decomp->coords[(Dim==3)?1:0]*worker->local_h;
            int gx = (c.x-1) + Decomp->coords[(Dim==3)?2:1]*worker->local_w;

            double xs = Fill->xS[count], xe = Fill->xE[count];
            double ys = Fill->yS[count], ye = Fill->yE[count];
            double zs = 0.0f, ze = 0.0f;
            if constexpr(Dim == 3){
                zs = Fill->zS[count], ze = Fill->zE[count];
            }

            fields& f = gridfield[idx_g];

            if(phase < (nPhase-1)){
                if((gx >= xs && gy >= ys && gx <= xe && gy <= ye) 
                && (Dim == 3 ? (gz >= zs && gz <= ze): true )){
                    f.phia[phase] = 1.0;

                    for(int i=0; i<nPhase; ++i){
                        if(phase != i){
                            f.phia[i] = 0.0;
                        }
                    }
                }else{
                    if(f.phia[phase] != 1.0){
                        f.phia[phase] = 0.0;
                    }
                }
            }else{
                double sum = 0.0f;

                for(int i=0; i<nPhase; ++i){
                    sum += f.phia[i];
                }

                if(sum > 1.0){
                    f.phia[phase] = 1.0;
                    for(int i=0; i<nPhase; ++i){
                        if(i != phase) f.phia[i] = 0.0;
                    }
                }else{
                    f.phia[nPhase] = 1.0 - sum;
                }
            }
            
        });
    }).wait();
}


template<int Dim>
void microsim::filling::fill_sphere(int phase, int count){
    workers *worker      = global_fd->workers_mpi;
    fillParameters *Fill = global_fd->Filling;
    fields *gridfield    = global_fd->gridfields;
    DomainDecomp *Decomp = global_fd->decomp;
    int nPhase = this->nPhase;

    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

            int gz = (Dim==3) ? (c.z-1) + Decomp->coords[0]*worker->local_d : 0;
            int gy = (c.y-1) + Decomp->coords[(Dim==3)?1:0]*worker->local_h;
            int gx = (c.x-1) + Decomp->coords[(Dim==3)?2:1]*worker->local_w;

            double xc = Fill->xC[count];
            double yc = Fill->yC[count];
            double zc = 0.0f;
            if constexpr(Dim == 3){
                zc = Fill->zC[count];
            }

            double rad = Fill->radius[count];

            fields& f = gridfield[idx_g];

            double dist_sq = (gx-xc) * (gx-xc) + (gy-yc) * (gy-yc);
            if(Dim == 3) dist_sq += (gz-zc) * (gz-zc);

            //double dist = sycl::sqrt(dist_sq);

            //double profile = 0.5 * (1.0 - std::tanh((dist - rad)/(4.0)));

            if(phase < (nPhase-1)){
                  if(dist_sq <= rad * rad){
                    f.phia[phase] = 1.0;

                    for(int i=0; i<nPhase; ++i){
                        if(phase != i){
                            f.phia[i] = 0.0;
                        }
                    }
                }else{
                    if(f.phia[phase] != 1.0){
                        f.phia[phase] = 0.0;
                    }
                }

                //f.phia[phase] = profile;
            }else{
                double sum = 0.0f;

                for(int i=0; i<nPhase-1; ++i){
                    sum += f.phia[i];
                }

                if(sum > 1.0){
                    f.phia[phase] = 1.0;
                    for(int i=0; i<nPhase-1; ++i){
                        if(i != phase) f.phia[i] = 0.0;
                    }
                }else{
                    f.phia[phase] = 1.0 - sum;
                }
            }
            
        });
    }).wait();
}

template<int Dim>
void microsim::filling::fill_composition(){
    workers *worker      = global_fd->workers_mpi;
    Phasefield_matraix *pfmatrix = global_fd->PFMAT;
    fields *gridfield    = global_fd->gridfields;
    int nPhase = this->nPhase;
    int nComp  = this->nComp;
    int function_F = global_fd->Flags->Function_F;

    #ifndef __KKS
    if(function_F == 1){
        std::cout<<"Using Thermo-Calc Function_F = 1 for filling composition. \n";
        microsim::ThermoDynamics thermo_F(global_fd, true);
        thermo_F.Initialize_host();
        
        auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
        dev_q.submit([&](sycl::handler &h){
            h.parallel_for(loop_range, [=](sycl::id<Dim> id){
                auto c = Coords<Dim>::from_id(id);
                int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

                bool phasefield = false;
                int active_phase = nPhase - 1;

                for(int b=0; b<nPhase-1; ++b){
                    if(gridfield[idx_g].phia[b] == 1.0){
                        phasefield = true;
                        active_phase = b;
                        for(int k=0; k<nComp-1; ++k){
                            gridfield[idx_g].C[k] = pfmatrix->ceq[b*nPhase*(nComp-1)+b*(nComp-1)+k];
                        }
                        break;
                    }
                }

                if(!phasefield){
                    for(int k=0; k<nComp-1; ++k){
                        gridfield[idx_g].C[k] = pfmatrix->cfill[(nPhase-1)*nPhase*(nComp-1)+(nPhase-1)*(nComp-1)+k];
                    }
                }

                thermo_F.Function_Mu(gridfield[idx_g].C, gridfield[idx_g].temperature[0], active_phase, gridfield[idx_g].mu);
            });
        }).wait();

    }else if(function_F ==4){
        std::cout<<"Using Thermo-Calc Function_F = 4 for filling composition. \n";
        microsim::ThermoDynamicsF4 thermo_F(global_fd, true);
        thermo_F.queue_ptr = &dev_q;
        thermo_F.Initialize_host(global_fd->reference);

        auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
        dev_q.submit([&](sycl::handler &h){
            h.parallel_for(loop_range, [=](sycl::id<Dim> id){
                auto c = Coords<Dim>::from_id(id);
                int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

                bool phasefield = false;
                int active_phase = nPhase - 1;

                for(int b=0; b<nPhase-1; ++b){
                    if(gridfield[idx_g].phia[b] == 1.0){
                        phasefield = true;
                        active_phase = b;
                        for(int k=0; k<nComp-1; ++k){
                            gridfield[idx_g].C[k] = pfmatrix->ceq[b*nPhase*(nComp-1)+b*(nComp-1)+k];
                        }
                        break;
                    }
                }

                if(!phasefield){
                    for(int k=0; k<nComp-1; ++k){
                        gridfield[idx_g].C[k] = pfmatrix->cfill[(nPhase-1)*nPhase*(nComp-1)+(nPhase-1)*(nComp-1)+k];
                    }
                }

                thermo_F.Function_Mu(gridfield[idx_g].C, gridfield[idx_g].temperature[0], active_phase, gridfield[idx_g].mu);
            });
        }).wait();
    }else{
        throw std::runtime_error("Unsupported Function_F type in filling::fill_composition");
    }
    #else
    // KKS Implementation (No Thermo needed for filling composition)
    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

            bool phasefield = false;
            // int active_phase = nPhase - 1;

            for(int b=0; b<nPhase-1; ++b){
                if(gridfield[idx_g].phia[b] == 1.0){
                    phasefield = true;
                    // active_phase = b;
                    for(int k=0; k<nComp-1; ++k){
                        gridfield[idx_g].C[k] = pfmatrix->ceq[b*nPhase*(nComp-1)+b*(nComp-1)+k];
                    }
                    break;
                }
            }

            if(!phasefield){
                for(int k=0; k<nComp-1; ++k){
                    gridfield[idx_g].C[k] = pfmatrix->cfill[(nPhase-1)*nPhase*(nComp-1)+(nPhase-1)*(nComp-1)+k];
                }
            }
        });
    }).wait();
    #endif

}


template<int Dim>
void microsim::filling::fill_temperature(){
    fields *gridfield    = global_fd->gridfields;
    workers *worker      = global_fd->workers_mpi;
    Settings *flags      = global_fd->Flags;
    Domain *dm         = global_fd->domainInfo;
    Input_parameters *inp = global_fd->parameters;
    Temp_gradient *tg = global_fd->Temp;
    DomainDecomp *Decomp = global_fd->decomp;

    double deltay = dm->DELTA_Y;
    double deltat = dm->DELTA_T;
    long t = inp->start_time;
    double shift_OFFSET = tg->offset;

    double BASE_POS = (tg->gradient_offset/deltay) - ((tg->velocity/deltay)*(t*deltat));
    double GRADIENT = (tg->DeltaT)*deltay/(tg->Distance);
    double base_Temp = tg->base_Temp;
    double T_iso = inp->T;
    bool is_isothermal = flags->ISOTHERMAL;

    auto loop_range = Coords<Dim>::to_range(worker->ghost_d, worker->ghost_h, worker->ghost_w);
    dev_q.submit([&](sycl::handler &h){
        h.parallel_for(loop_range, [=](sycl::id<Dim> id){
            auto c = Coords<Dim>::from_id(id);
            int idx_g = idx(c.z, c.y, c.x, worker->ghost_h, worker->ghost_w);

            if(is_isothermal){
                gridfield[idx_g].temperature[0] = T_iso;
            }else{
                int gy = (c.y-1) + Decomp->coords[(Dim==3)?1:0]*worker->local_h;
                gridfield[idx_g].temperature[0] = base_Temp + (gy - BASE_POS) * GRADIENT;
            }
        });
    }).wait();

}

microsim::filling::~filling(){
    std::cout<<"Freeing microsim::filling allocated memory\n";
}