#ifndef FUNCTION_Q_H_
#define FUNCTION_Q_H_

void q_divx (struct gradlayer *grad1, struct gradlayer *grad1_front, long a, long b, double *qab) {
  qab[X] =       (grad1->phistagg[X][a]*grad1->gradphi[X][b]    - grad1->phistagg[X][b]*grad1->gradphi[X][a]);
  qab[Y] =   0.5*(grad1->phistagg[X][a]*(grad1->gradphi_c[Y][b] + grad1_front->gradphi_c[Y][b]) - grad1->phistagg[X][b]*(grad1->gradphi_c[Y][a] + grad1_front->gradphi_c[Y][a]));
  if (DIMENSION != 2) {
    qab[Z] = 0.5*(grad1->phistagg[X][a]*(grad1->gradphi_c[Z][b] + grad1_front->gradphi_c[Z][b]) - grad1->phistagg[X][b]*(grad1->gradphi_c[Z][a] + grad1_front->gradphi_c[Z][a]));
  } else {
    qab[Z] = 0.0;
  }
//   qab[Y]  = (grad1->phi_center[a]*(grad1->gradphi_c[Y][b])             - grad1->phi_center[b]*(grad1->gradphi_c[Y][a]));
//   qab[Y] += (grad1_front->phi_center[a]*(grad1_front->gradphi_c[Y][b]) - grad1_front->phi_center[b]*(grad1_front->gradphi_c[Y][a]));
//   qab[Y] *= 0.5;
  
  multiply(Rotation_matrix[a][b], qab, Rotated_qab, DIMENSION);
  long dir;
  for(dir=0; dir<DIMENSION; dir++) {
    qab[dir] = Rotated_qab[dir];
  }
}
void q_divy (struct gradlayer *grad1, struct gradlayer *grad1_right, long a, long b, double *qab) {
  qab[X] =   0.5*(grad1->phistagg[Y][a]*(grad1->gradphi_c[X][b] + grad1_right->gradphi_c[X][b]) - grad1->phistagg[Y][b]*(grad1->gradphi_c[X][a]+grad1_right->gradphi_c[X][a]));
//   qab[X]  = (grad1->phi_center[a]*(grad1->gradphi_c[X][b])             - grad1->phi_center[b]*(grad1->gradphi_c[X][a]));
//   qab[X] += (grad1_right->phi_center[a]*(grad1_right->gradphi_c[X][b]) - grad1_right->phi_center[b]*(grad1_right->gradphi_c[X][a]));
//   qab[X] *= 0.5;
  qab[Y] =   (grad1->phistagg[Y][a]*grad1->gradphi[Y][b]    - grad1->phistagg[Y][b]*grad1->gradphi[Y][a]);
  
  if (DIMENSION != 2) {
    qab[Z] =   0.5*(grad1->phistagg[Y][a]*(grad1->gradphi_c[Z][b] + grad1_right->gradphi_c[Z][b]) - grad1->phistagg[Y][b]*(grad1->gradphi_c[Z][a] + grad1_right->gradphi_c[Z][a]));
  } else {
    qab[Z] = 0.0;
  }
  
  multiply(Rotation_matrix[a][b], qab, Rotated_qab, DIMENSION);
  long dir;
  for(dir=0; dir<DIMENSION; dir++) {
    qab[dir] = Rotated_qab[dir];
  }
}
void q_divz (struct gradlayer *grad1, struct gradlayer *grad1_top, long a, long b, double *qab) {
  qab[X]  =   0.5*(grad1->phistagg[Z][a]*(grad1->gradphi_c[X][b] + grad1_top->gradphi_c[X][b]) - grad1->phistagg[Z][b]*(grad1->gradphi_c[X][a] + grad1_top->gradphi_c[X][a]));
  qab[Y]  =   0.5*(grad1->phistagg[Z][a]*(grad1->gradphi_c[Y][b] + grad1_top->gradphi_c[Y][b]) - grad1->phistagg[Z][b]*(grad1->gradphi_c[Y][a] + grad1_top->gradphi_c[Y][a]));
  
  qab[Z]  =  (grad1->phistagg[Z][a]*grad1->gradphi[Z][b]    - grad1->phistagg[Z][b]*grad1->gradphi[Z][a]);
  
  multiply(Rotation_matrix[a][b], qab, Rotated_qab, DIMENSION);
  long dir;
  for(dir=0; dir<DIMENSION; dir++) {
    qab[dir] = Rotated_qab[dir];
  }
}


void q_dadphi (double *phi, struct gradlayer *grad1, long a,  long b , double *qab) {
  qab[X] =  (phi[a]*grad1->gradphi_c[X][b] - phi[b]*grad1->gradphi_c[X][a]);
  qab[Y] =  (phi[a]*grad1->gradphi_c[Y][b] - phi[b]*grad1->gradphi_c[Y][a]);
  if (DIMENSION !=2) {
    qab[Z] =  (phi[a]*grad1->gradphi_c[Z][b] - phi[b]*grad1->gradphi_c[Z][a]);
  } else {
    qab[Z] = 0.0;
  }
  multiply(Rotation_matrix[a][b], qab, Rotated_qab, DIMENSION);
  long dir;
  for(dir=0; dir<DIMENSION; dir++) {
    qab[dir] = Rotated_qab[dir];
  }
}
void rotate_vector(double *q, double *rot_q, long a, long b) {
  multiply(Inv_Rotation_matrix[a][b], q, rot_q, DIMENSION);
}

#endif