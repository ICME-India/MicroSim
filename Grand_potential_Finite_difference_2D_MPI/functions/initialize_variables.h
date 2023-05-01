#ifndef INITIALIAZE_VARIABLES_H_
#define INITIALIAZE_VARIABLES_H_

void initialize_variables() {
  long a, i, j;
  for (a=0;a<NUMPHASES;a++) {
    for (i=0;i<NUMCOMPONENTS-1;i++) {
      for (j=0;j<NUMCOMPONENTS-1;j++) {
        if (i==j) {
          muc[a][i][j]=2.0*A[a][i][j];
        } else {
          muc[a][i][j]=A[a][i][j];
        }
      }
    }
    matinvnew(muc[a], cmu[a], NUMCOMPONENTS-1);
  }
  
// #ifndef ISOTHERMAL
  for (a=0;a<NUMPHASES;a++) {
    for (i=0; i < NUMCOMPONENTS-1; i++) {      
      dcbdT_phase[a][i] = 0.0;
      for (j=0; j < NUMCOMPONENTS-1; j++) {
        dcbdT_phase[a][i] += cmu[a][i][j]*(-dBbdT[a][j]);
//           dcbdT_phase[a][i] += dc_dmu(a,i,j)*(-dBbdT[a][j]);
      }
    }
  }
// #endif
  dcdmu       = MallocM((NUMCOMPONENTS-1),(NUMCOMPONENTS-1));
  dcdmu_phase = Malloc3M(NUMPHASES, NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  Ddcdmu      = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  inv_dcdmu   = MallocM(NUMCOMPONENTS-1,NUMCOMPONENTS-1);
  deltamu     = MallocV(NUMCOMPONENTS-1);
  deltac      = MallocV(NUMCOMPONENTS-1);
  sum         = MallocV(NUMCOMPONENTS-1);
  divphi      = MallocV(NUMPHASES);
  lambda_phi  = MallocV(NUMPHASES);
  divflux     = MallocV(NUMCOMPONENTS-1);
  c_old       = MallocV(NUMCOMPONENTS-1);
  c_new       = MallocV(NUMCOMPONENTS-1);
  c_tdt       = MallocV(NUMCOMPONENTS-1);
  divjat      = MallocV(NUMCOMPONENTS-1);
  deltac      = MallocV(NUMCOMPONENTS-1);
  
  gridinfo_instance = (struct fields *)malloc(sizeof(*gridinfo_instance));
  allocate_memory_fields(gridinfo_instance);
}

#endif
