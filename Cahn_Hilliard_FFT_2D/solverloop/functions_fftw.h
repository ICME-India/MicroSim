#include<complex.h>
#include <gsl/gsl_rng.h>
#include "fftw3.h"


void evolve_fftw() {
  calculate_df();
  call_fftwF();
  update_Fourier_com_phi();
  call_fftwB();
  Fouriertoreal();
}

void prepfftw() {

  long index;

  delta_kx = (2.0*M_PI)/(MESH_X*deltax);
  delta_ky = (2.0*M_PI)/(MESH_Y*deltay);

  kx = (double *) malloc(MESH_X*sizeof(double));
  ky = (double *) malloc(MESH_Y*sizeof(double));

  half_MESH_X = (long) (MESH_X/2);
  half_MESH_Y = (long) (MESH_Y/2);

  for (i = 0; i < MESH_X; i++) {
    if (i < half_MESH_X) {
      kx[i] = i*delta_kx;
    }
    else {
      kx[i] = (i-MESH_X)*delta_kx;
    }
  }

  for (j = 0; j < MESH_Y; j++) {
    if (j < half_MESH_Y) {
      ky[j] = j*delta_ky;
    }
    else {
      ky[j] = (j-MESH_Y)*delta_ky;
    }
  }

  P = (double *) malloc((NUMPHASES-1)*sizeof(double)); 
  for (b=0; b < NUMPHASES-1; b++) {
    P[b] = 1.0;
  }

  W = (double *) malloc((NUMPHASES-1)*sizeof(double));


  phi = fftw_malloc((NUMPHASES-1)*sizeof(**phi));
  dfdphi = fftw_malloc((NUMPHASES-1)*sizeof(**dfdphi));
  
  for (b=0; b < NUMPHASES-1; b++) {

    phi[b] = fftw_malloc(index_count*sizeof(fftw_complex));

    dfdphi[b] = fftw_malloc(index_count*sizeof(fftw_complex));

  }

  com = fftw_malloc((NUMCOMPONENTS-1)*sizeof(**com));
  dfdc = fftw_malloc((NUMCOMPONENTS-1)*sizeof(**dfdc));
  for (a=0; a < NUMCOMPONENTS-1; a++) {
    
    com[a] = fftw_malloc(index_count*sizeof(fftw_complex));

    dfdc[a] = fftw_malloc(index_count*sizeof(fftw_complex));

  }

  planF = fftw_plan_dft_2d(MESH_X, MESH_Y, phi[0], phi[0],  FFTW_FORWARD, FFTW_ESTIMATE);
  planB = fftw_plan_dft_2d(MESH_X, MESH_Y, phi[0], phi[0], FFTW_BACKWARD, FFTW_ESTIMATE);

  realtoFourier();

  if (NUMCOMPONENTS > 2 || NUMPHASES > 2) {
    printf("************\n");
    printf("Current version works for NUMPHASES = 2 and NUMCOMPONENTS = 2\n");
    printf("Otherwise results may be wrong\n");
    printf("************\n");
    //exit(1);
  }

}


void realtoFourier() {

  long index;

  for (x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        index = x*layer_size + z*rows_y + y;
        for (a=0; a < NUMCOMPONENTS-1; a++) {
          com[a][index] = gridinfo[index].compi[a] + 0.0*_Complex_I;
        }
        for (b=0; b < NUMPHASES-1; b++) {
          phi[b][index] = gridinfo[index].phia[b] + 0.0*_Complex_I;
        }
      }
    }
  }

}

void Fouriertoreal() {

  long index;

  inv_denom = 1.0 / (MESH_X*MESH_Y);

  for (b=0; b < NUMPHASES-1; b++) {
    global_max_min.rel_change_phi[b] = 0.0;
  }
  for (a=0; a < NUMCOMPONENTS-1; a++) {
    global_max_min.rel_change_com[a] = 0.0;
  }

  for (x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        index = x*layer_size + z*rows_y + y;
        for (a=0; a < NUMCOMPONENTS-1; a++) {
          com[a][index] = creal(com[a][index])*inv_denom + 0.0*_Complex_I;
          global_max_min.rel_change_com[a] += (creal(com[a][index])-gridinfo[index].compi[a])*(creal(com[a][index])-gridinfo[index].compi[a]);
          gridinfo[index].compi[a] = creal(com[a][index]);
          if (gridinfo[index].compi[a] > global_max_min.com_max[a]) { 
            global_max_min.com_max[a] = gridinfo[index].compi[a];
          }
          if (gridinfo[index].compi[a] < global_max_min.com_min[a]) {
            global_max_min.com_min[a] = gridinfo[index].compi[a];
          }
        }
        for (b=0; b < NUMPHASES-1; b++) {
          phi[b][index] = creal(phi[b][index])*inv_denom + 0.0*_Complex_I;
          global_max_min.rel_change_phi[b] += (creal(phi[b][index])-gridinfo[index].phia[b])*(creal(phi[b][index])-gridinfo[index].phia[b]);
          gridinfo[index].phia[b] = creal(phi[b][index]);
          if (gridinfo[index].phia[b] > global_max_min.phi_max[b]) { 
            global_max_min.phi_max[b] = gridinfo[index].phia[b];
          }
          if (gridinfo[index].phia[b] < global_max_min.phi_min[b]) { 
            global_max_min.phi_min[b] = gridinfo[index].phia[b];
          }
        }
      }
    }
  }

}

void calculate_df() {

  long index;

  for (x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {

        index = x*layer_size + z*rows_y + y;

        for (b=0; b < NUMPHASES-1; b++) {

          if (gridinfo[index].phia[b] < 0) {
            W[b] = 0.0;
          }
          else if (gridinfo[index].phia[b] > 1) {
            W[b] = 1.0;
          }
          else {
            W[b] = gridinfo[index].phia[b]*gridinfo[index].phia[b]*gridinfo[index].phia[b]*(10.0-15.0*gridinfo[index].phia[b]+6.0*gridinfo[index].phia[b]*gridinfo[index].phia[b]);
          }

        }

        sum1 = 0;
        for (b=0; b < NUMPHASES-1; b++) {
          sum1 += W[b];
        }

        a=0; //For 1 solute or 2 components only 
        dfdc[a][index] = 2.0*A_fm[0]*gridinfo[index].compi[a]*(1.0-sum1) - 2.0*B_fp[0]*(1.0-gridinfo[index].compi[a])*sum1;

        for (b=0; b < NUMPHASES-1; b++) {

          if (gridinfo[index].phia[b] < 0 || gridinfo[index].phia[b] > 1) {
            dWdphi = 0.0;
          }
          else {
            dWdphi = 30.0*gridinfo[index].phia[b]*gridinfo[index].phia[b]*(1.0-gridinfo[index].phia[b])*(1.0-gridinfo[index].phia[b]);
          }
          a = 0; //For 1 solute or 2 components
          tmp1 = dWdphi*(-A_fm[0]*gridinfo[index].compi[a]*gridinfo[index].compi[a] + B_fp[0]*(1.0-gridinfo[index].compi[a])*(1.0-gridinfo[index].compi[a]));
          tmp2 = 2.0*P[0]*gridinfo[index].phia[b]*(1.0-gridinfo[index].phia[b])*(1.0-2.0*gridinfo[index].phia[b]); 
          dfdphi[b][index] = tmp1 + tmp2;
        }
      }
    }
  }
}

void call_fftwF() {

  long index;
  
  for (a=0; a < NUMCOMPONENTS-1; a++) {
    fftw_execute_dft(planF,dfdc[a],dfdc[a]);
  }
  
  for (b=0; b < NUMPHASES-1; b++) {
    fftw_execute_dft(planF,dfdphi[b],dfdphi[b]);
  }
  
  for (b=0; b < NUMPHASES-1; b++) {
    fftw_execute_dft(planF,phi[b],phi[b]);
  }
  
  for (a=0; a < NUMCOMPONENTS-1; a++) {
    fftw_execute_dft(planF,com[a],com[a]);
  }

}

void update_Fourier_com_phi() {

  long index;
  
  for (x=0;x < rows_x; x++) {
    for (z=0; z < rows_z; z++) {
      for (y=0; y < rows_y; y++) {
        
        index = x*layer_size + z*rows_y + y;
        
        k2 = kx[x]*kx[x] + ky[y]*ky[y]; 
        // For 2 components and 2 phases only
        for (b=0; b < NUMPHASES-1; b++) { 
          inv_denom = 1.0 / (1.0+2.0*L_phi[0][1]*Kappa_phi[0][1]*k2*deltat);
          phi[b][index] = inv_denom*(phi[b][index]-dfdphi[b][index]*deltat*L_phi[0][1]);
        }

        for (a=0; a < NUMCOMPONENTS-1; a++) {
          inv_denom = 1.0 / (1.0+2.0*AtomicMobility[0][0][0]*Kappa_c[0][1]*k2*k2*deltat); 
          com[a][index] = inv_denom*(com[a][index]-dfdc[a][index]*k2*deltat*AtomicMobility[0][0][0]);
        }
      }
    }
  }
}

void call_fftwB() {

  long index;
  
  for (b=0; b < NUMPHASES-1; b++) {
    fftw_execute_dft(planB,phi[b],phi[b]);
  }
  
  for (a=0; a < NUMCOMPONENTS-1; a++) {
    fftw_execute_dft(planB,com[a],com[a]);
  }
  
}
