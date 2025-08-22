void calculate_distribution_function_2D(long x) {
  long dim, k, l, a, gidy;
  double feq;
  double magu;
  double c_u;
  double GD[2];
  double GB[2];
  double GS[2];
  double G1[2];
//   double nu=0.07;
  double Gfdis=0.0;
//   double W_0 = 0.8;
//   double g_y = 0.01;
//   double beta_c = 0.2;
//   double delta_rho = 0.0001;
  double df;
  double gradphi[2];
  double n[2];
  double rho;
  
  for (gidy=0; gidy <= (workers_mpi.end[Y]+2); gidy++) {
    center        =  gidy   + (x)*workers_mpi.layer_size;
    front         =  gidy   + (x+1)*workers_mpi.layer_size;
    back          =  gidy   + (x-1)*workers_mpi.layer_size;
    right         =  center + 1;
    left          =  center - 1;
    
    magu = 0.0;
    
    if ((x <= workers_mpi.end[X]+2) && (x>0)) {
      gradphi[0] = 0.5*(gridinfo_w[front].phia[NUMPHASES-1] - gridinfo_w[back].phia[NUMPHASES-1])/deltax;
    } else {
      gradphi[0] = 0.0;
    }
  
    if (gidy > 0) {
      gradphi[1] = 0.5*(gridinfo_w[right].phia[NUMPHASES-1] - gridinfo_w[left].phia[NUMPHASES-1])/deltay;
    } else {
      gradphi[1] = 0.0;
    }
    
    if ((fabs(gradphi[0]) > 1e-12) || (fabs(gradphi[1]) > 1e-12)) {
      n[0] = gradphi[0]/((gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1] + 1e-6));
      n[1] = gradphi[1]/((gradphi[0]*gradphi[0] + gradphi[1]*gradphi[1] + 1e-6));
    }
    
    
    for (dim=0; dim<DIMENSION; dim++) {
      magu += lbm_gridinfo_w[center].u[dim]*lbm_gridinfo_w[center].u[dim];
    }
    
    for(dim=0; dim < DIMENSION; dim++) {
      GD[dim] = -0.5/(W_0*W_0)*nu*10*(1.0-gridinfo_w[center].phia[NUMPHASES-1])*(1.0-gridinfo_w[center].phia[NUMPHASES-1])*lbm_gridinfo_w[center].u[dim];
      GB[dim] = 0.0;
      GS[dim] = 0.0;
    }
    for (k=0; k < NUMCOMPONENTS-1; k++) {
      GB[0] += -0.50*lbm_gridinfo_w[center].rho*g_y*beta_c[k]*(gridinfo_w[center].composition[k]-ceq[NUMPHASES-1][NUMPHASES-1][k])*(gridinfo_w[center].phia[NUMPHASES-1]);
      GB[1] = 0.0;
    }
    
    for (dim=0; dim < DIMENSION; dim++) {
      for (a=0; a < NUMPHASES-1; a++) {
        GS[dim] += -((rho_LBM[a]-rho_LBM[NUMPHASES-1])/rho_LBM[NUMPHASES-1])*n[dim]*(gridinfo_w[center].deltaphi[a])/deltat;
      }
    }
      
    for (k=0; k < 9; k++) {
      c_u  = 0.0;
      for (dim=0; dim < DIMENSION; dim++) {
        c_u += lbm_gridinfo_w[center].u[dim]*vgrid[dim][k];
      }
      Gfdis = 0.0;
      for (dim=0; dim<DIMENSION; dim++) {
        G1[dim]  = w_2D[k]*(3.0*(vgrid[dim][k]-lbm_gridinfo_w[center].u[dim]) + 9*c_u*vgrid[dim][k]);
        Gfdis   += G1[dim]*(GD[dim]+GB[dim]+GS[dim]);
      }
      feq = lbm_gridinfo_w[center].rho*w_2D[k]*(1.0 + 3.0*c_u + 4.5*c_u*c_u - 1.50*magu);
      lbm_gridinfo_w[center].fdis[k] = (1.0-omega)*lbm_gridinfo_w[center].fdis[k] + omega*feq + Gfdis*dt;
    }
  }
}

void calculate_streaming_row_forward_2D(long x) {
  long gidy;
  for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
    center        =  gidy   + (x)*workers_mpi.layer_size;
    back          =  gidy   + (x-1)*workers_mpi.layer_size;
    
    lbm_gridinfo_w[back].fdis[2]   =  lbm_gridinfo_w[center].fdis[2];
    lbm_gridinfo_w[back].fdis[5]   =  lbm_gridinfo_w[center].fdis[5];
    lbm_gridinfo_w[back].fdis[6]   =  lbm_gridinfo_w[center].fdis[6];
  }
}

void calculate_streaming_row_back_2D(long x) {
  long gidy;
  for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
    center        =  gidy   + (x)*workers_mpi.layer_size;
    back          =  gidy   + (x-1)*workers_mpi.layer_size;
    lbm_gridinfo_w[center].fdis[0] =  lbm_gridinfo_w[back].fdis[0];
    lbm_gridinfo_w[center].fdis[4] =  lbm_gridinfo_w[back].fdis[4];
    lbm_gridinfo_w[center].fdis[7] =  lbm_gridinfo_w[back].fdis[7];
  }
}

void calculate_streaming_col_forward_2D(long x) {
  long gidy;
  long left;
  long center;
  for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
    center        =  gidy   + (x)*workers_mpi.layer_size;

    left          =  center - 1;
    lbm_gridinfo_w[left].fdis[7]   =  lbm_gridinfo_w[center].fdis[7];
    lbm_gridinfo_w[left].fdis[6]   =  lbm_gridinfo_w[center].fdis[6];
    lbm_gridinfo_w[left].fdis[3]   =  lbm_gridinfo_w[center].fdis[3];
  }
}
void calculate_streaming_col_back_2D(long x) {
  long gidy;
  
  for (gidy=workers_mpi.end[Y]+2; gidy > 1; gidy--) {
    center        =  gidy   + (x)*workers_mpi.layer_size;

    left          =  center - 1;

    lbm_gridinfo_w[center].fdis[5] =  lbm_gridinfo_w[left].fdis[5];
    lbm_gridinfo_w[center].fdis[4] =  lbm_gridinfo_w[left].fdis[4];
    lbm_gridinfo_w[center].fdis[1] =  lbm_gridinfo_w[left].fdis[1];
  }
}
void calculate_streaming(long *start, long *end) {
  long x, gidy, k;
  long center;
  long back, left;
  for (x=start[X]-2; x<=end[X]+2; x++) {
    for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
      center        =  gidy   + (x)*workers_mpi.layer_size;
      for (k=0; k < 9; k++) {
        temp_fdis[center][k] = lbm_gridinfo_w[center].fdis[k]; 
      }
    }
  }
  for (x=start[X]-2; x<=end[X]+2; x++) {
    for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
      center        =  gidy   + (x)*workers_mpi.layer_size;
      back          =  gidy   + (x-1)*workers_mpi.layer_size;
      lbm_gridinfo_w[back].fdis[2]   =  temp_fdis[center][2];
      lbm_gridinfo_w[back].fdis[5]   =  temp_fdis[center][5];
      lbm_gridinfo_w[back].fdis[6]   =  temp_fdis[center][6];
      
      lbm_gridinfo_w[center].fdis[0] =  temp_fdis[back][0];
      lbm_gridinfo_w[center].fdis[4] =  temp_fdis[back][4];
      lbm_gridinfo_w[center].fdis[7] =  temp_fdis[back][7];
    }
  }
  for (x=start[X]-2; x<=end[X]+2; x++) {
    for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
      center        =  gidy   + (x)*workers_mpi.layer_size;
      for (k=0; k < 9; k++) {
        temp_fdis[center][k] = lbm_gridinfo_w[center].fdis[k]; 
      }
    }
  }
  for (x=start[X]-2; x<=end[X]+2; x++) {
    for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {
      center        =  gidy   + (x)*workers_mpi.layer_size;
      left          =  center - 1;
      lbm_gridinfo_w[left].fdis[7]   =  temp_fdis[center][7];
      lbm_gridinfo_w[left].fdis[6]   =  temp_fdis[center][6];
      lbm_gridinfo_w[left].fdis[3]   =  temp_fdis[center][3];
      
      lbm_gridinfo_w[center].fdis[5] =  temp_fdis[left][5];
      lbm_gridinfo_w[center].fdis[4] =  temp_fdis[left][4];
      lbm_gridinfo_w[center].fdis[1] =  temp_fdis[left][1];
    }
  }
}


void compute_field_variables_LBM_2D(long x) {
  long k,gidy,dim;
  double rho=0.0, u[3];
  
  for (gidy=2; gidy <= (workers_mpi.end[Y]); gidy++) {
    center        =  gidy   + (x)*workers_mpi.layer_size;
    rho = 0.0;
    for(dim=0; dim < DIMENSION; dim++) {
      u[dim] = 0.0;
    }
    for (k=0; k < 9; k++) {
      rho += lbm_gridinfo_w[center].fdis[k];
      for (dim=0; dim < DIMENSION; dim++) {
        u[dim] += lbm_gridinfo_w[center].fdis[k]*vgrid[dim][k];
      }
    }
    lbm_gridinfo_w[center].rho = rho;
    for(dim=0; dim < DIMENSION; dim++) {
//       if (rho !=0.0) {
        lbm_gridinfo_w[center].u[dim] = u[dim]/rho;
//       }
    }
  }
}


