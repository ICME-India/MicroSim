#ifndef function_A_00_H_
#define function_A_00_H_

double function_A_00_dAdphi(double *phi, struct gradlayer **gradient, long gidy, long a) {
  long b;
  double scalprod=0.0;
  double sum=0.0;
  long dir;
  
  struct gradlayer *grad1;
  grad1 = &gradient[1][gidy];
  
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      scalprod = 0.0;
      for (dir=0; dir < DIMENSION; dir++)  {
	      scalprod += (phi[a]*grad1->gradphi_c[dir][b] - phi[b]*grad1->gradphi_c[dir][a])*grad1->gradphi_c[dir][b];
      }
      scalprod *= Gamma[a][b];
      sum += scalprod;
    }
  }
  return 2.0*sum;
}
double function_A_00_divdAdgradphi(struct gradlayer **gradient, long index, long gidy, long a) {
  long b;
  double divergence=0.0;
  double sum=0.0;
  
  struct gradlayer *grad1;
  struct gradlayer *grad1_left;
  struct gradlayer *grad1_back;
  struct gradlayer *grad1_bottom;
  
  grad1      =  &gradient[1][gidy];
  grad1_left =  &gradient[1][gidy-1];
  grad1_back =  &gradient[0][gidy];
  
  if (DIMENSION !=2){
    grad1_bottom =  &gradient[1][gidy-workers_mpi.rows_y];
  }
  
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      divergence = ((grad1->phistagg[X][a]*grad1->gradphi[X][b] - grad1->phistagg[X][b]*grad1->gradphi[X][a])*(-grad1->phistagg[X][b]) - ((grad1_back->phistagg[X][a]*grad1_back->gradphi[X][b] - grad1_back->phistagg[X][b]*grad1_back->gradphi[X][a])*(-grad1_back->phistagg[X][b])))/deltax;
      
      divergence += ((grad1->phistagg[Y][a]*grad1->gradphi[Y][b] - grad1->phistagg[Y][b]*grad1->gradphi[Y][a])*(-grad1->phistagg[Y][b]) - ((grad1_left->phistagg[Y][a]*grad1_left->gradphi[Y][b] - grad1_left->phistagg[Y][b]*grad1_left->gradphi[Y][a])*(-grad1_left->phistagg[Y][b])))/deltay;
      
      if (DIMENSION !=2) {
         divergence += (((grad1->phistagg[Z][a]*grad1->gradphi[Z][b] - grad1->phistagg[Z][b]*grad1->gradphi[Z][a])*(-grad1->phistagg[Z][b])) - ((grad1_bottom->phistagg[Z][a]*grad1_bottom->gradphi[Z][b] - grad1_bottom->phistagg[Z][b]*grad1_bottom->gradphi[Z][a])*(-grad1_bottom->phistagg[Z][b])))/deltay;
      }
      
      divergence *= Gamma[a][b];
      sum        += divergence;
    }
  }
  return 2.0*sum;
}

#endif
