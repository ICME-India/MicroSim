#ifndef FUNCTION_A_01_H_
#define FUNCTION_A_01_H_

double function_A_01_dAdphi(double *phi, struct gradlayer **gradient, long gidy, long a) {
  long b;
  double scalprod1=0.0, scalprod2=0.0;
  double sum1=0.0, sum2=0.0,sum3=0.0;
  long dir;
  double dadq[3];
  double qab[3];
  double ac1, ac2, qab2;
  double dqdphi[3], Rotated_dqdphi[3];
  
  struct gradlayer *grad1,*grad1_back,*grad1_front,*grad1_right,*grad1_left, *grad1_top, *grad1_bottom;
  
  grad1       =  &gradient[1][gidy];
  grad1_left  =  &gradient[1][gidy-1];
  grad1_back  =  &gradient[0][gidy];
  grad1_right =  &gradient[1][gidy+1];
  grad1_front =  &gradient[2][gidy];
  
  if(DIMENSION ==3) {
    grad1_top    =  &gradient[1][gidy+rows_y];
    grad1_bottom =  &gradient[1][gidy-rows_y];
  }
  
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      q_dadphi(phi, grad1, a,  b , qab);    
      dAdq(qab, dadq, a, b);
      qab2 = (qab[X]*qab[X] + qab[Y]*qab[Y] + qab[Z]*qab[Z]);
      
      ac1  = function_ac(qab, a, b);
      ac2  = ac1*ac1;
      
      scalprod1 = 0.0;
      scalprod2 = 0.0;
      
      dqdphi[X] = grad1->gradphi_c[X][b];
      dqdphi[Y] = grad1->gradphi_c[Y][b];
      if (DIMENSION == 3) {
        dqdphi[Z] = grad1->gradphi_c[Z][b];
      } else {
        dqdphi[Z] = 0.0;
      }
      
      multiply(Rotation_matrix[a][b], dqdphi, Rotated_dqdphi, DIMENSION);
      
      for(dir=0; dir<DIMENSION; dir++) {
	dqdphi[dir] = Rotated_dqdphi[dir];
      }
      
      for (dir=0; dir < DIMENSION; dir++)  {
// 	scalprod1 += (phi[a]*grad1->gradphi_c[dir][b] - phi[b]*grad1->gradphi_c[dir][a])*grad1->gradphi_c[dir][b];
        scalprod1 += qab[dir]*dqdphi[dir];
      }

      for (dir=0; dir < DIMENSION; dir++)  {
	scalprod2 += (dadq[dir])*dqdphi[dir];
      }
       
      scalprod1 *= Gamma[a][b]*ac2;
      
      scalprod2 *= Gamma[a][b]*ac1*qab2;
      
      sum1 += scalprod1;
      sum2 += scalprod2;
    }
  }
  return 2.0*(sum1 + sum2);
}

double function_A_01_dAdphi_smooth(double *phi, struct gradlayer **gradient, long gidy, long a) {
  long b;
  double scalprod1=0.0, scalprod2=0.0;
  double sum1=0.0, sum2=0.0;
  long dir;
  double dadq[3];
  double qab[3];
  double ac1, ac2, qab2;
  double dqdphi[3], Rotated_dqdphi[3];

  struct gradlayer *grad1;
  grad1 = &gradient[1][gidy];
   
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      q_dadphi(phi, grad1, a,  b , qab);
      scalprod1 = 0.0;
      
      for (dir=0; dir < DIMENSION; dir++)  {
        scalprod1 += (phi[a]*grad1->gradphi_c[dir][b] - phi[b]*grad1->gradphi_c[dir][a])*grad1->gradphi_c[dir][b];
      }
      scalprod1 *= Gamma[a][b];
      
      sum1 += scalprod1;
    }
  }
  return 2.0*(sum1 + sum2);
}
double function_A_01_divdAdgradphi(struct gradlayer **gradient, long index, long gidy, long a) {
  long b;
  double sum1 = 0.0, sum2 = 0.0;
  double divergence1=0.0, divergence2=0.0, divergence3=0.0;
  double dadq_front[3],                 dadq_back[3],         dadq_left[3],         dadq_right[3],         dadq_top[3],         dadq_bottom[3];
  double rotated_dadq_front[3], rotated_dadq_back[3], rotated_dadq_left[3], rotated_dadq_right[3], rotated_dadq_top[3], rotated_dadq_bottom[3];
  double qab_front[3],  qab_back[3],  qab_left[3],  qab_right[3], qab_top[3], qab_bottom[3];
  
  double ac_front,  ac_back,   ac_right,  ac_left, ac_top, ac_bottom;
  double ac_front2, ac_back2,  ac_right2, ac_left2, ac_top2, ac_bottom2;
  double qab_front2, qab_back2, qab_right2, qab_left2, qab_top2, qab_bottom2;
   
  struct gradlayer *grad1;
  struct gradlayer *grad1_left;
  struct gradlayer *grad1_back;
  struct gradlayer *grad1_right;
  struct gradlayer *grad1_front;
  struct gradlayer *grad1_top;
  struct gradlayer *grad1_bottom;
  
  
  grad1       =  &gradient[1][gidy];
  grad1_left  =  &gradient[1][gidy-1];
  grad1_back  =  &gradient[0][gidy];
  grad1_right =  &gradient[1][gidy+1];
  grad1_front =  &gradient[2][gidy];
  
  if (DIMENSION !=2) {
    grad1_top    =  &gradient[1][gidy+rows_y];
    grad1_bottom =  &gradient[1][gidy-rows_y];
  }
  
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      q_divx (grad1, grad1_front, a,  b, qab_front);
      ac_front = function_ac(qab_front, a, b);
      ac_front2 = ac_front*ac_front;
      
      q_divx (grad1_back, grad1, a,  b, qab_back);
      ac_back   = function_ac(qab_back, a, b);
      ac_back2  = ac_back*ac_back;
      
      q_divy (grad1, grad1_right, a,  b, qab_right);
      ac_right   = function_ac(qab_right, a, b);
      ac_right2  = ac_right*ac_right;
      
      
      q_divy (grad1_left, grad1, a,  b, qab_left);
      ac_left = function_ac(qab_left, a, b);
      ac_left2  = ac_left*ac_left;
      
      
      if(DIMENSION != 2) {
        q_divz (grad1, grad1_top, a, b, qab_top);
        ac_top  = function_ac(qab_top, a, b);
        ac_top2 = ac_top*ac_top; 
        
        q_divz (grad1_bottom, grad1, a, b, qab_bottom);
        ac_bottom  = function_ac(qab_bottom, a, b);
        ac_bottom2 = ac_bottom*ac_bottom; 
      }
      
      divergence1  = ((ac_front2*(grad1->phistagg[X][a]*grad1->gradphi[X][b] - grad1->phistagg[X][b]*grad1->gradphi[X][a])*(-grad1->phistagg[X][b])) - (ac_back2*(grad1_back->phistagg[X][a]*grad1_back->gradphi[X][b] - grad1_back->phistagg[X][b]*grad1_back->gradphi[X][a])*(-grad1_back->phistagg[X][b])))/deltax;
      
      divergence1 += ((ac_right2*(grad1->phistagg[Y][a]*grad1->gradphi[Y][b] - grad1->phistagg[Y][b]*grad1->gradphi[Y][a])*(-grad1->phistagg[Y][b])) - (ac_left2*(grad1_left->phistagg[Y][a]*grad1_left->gradphi[Y][b] - grad1_left->phistagg[Y][b]*grad1_left->gradphi[Y][a])*(-grad1_left->phistagg[Y][b])))/deltay;
      
      if (DIMENSION !=2) {
        divergence1 += ((ac_top2*(grad1->phistagg[Z][a]*grad1->gradphi[Z][b] - grad1->phistagg[Z][b]*grad1->gradphi[Z][a])*(-grad1->phistagg[Z][b])) - (ac_bottom2*(grad1_bottom->phistagg[Z][a]*grad1_bottom->gradphi[Z][b] - grad1_bottom->phistagg[Z][b]*grad1_bottom->gradphi[Z][a])*(-grad1_bottom->phistagg[Z][b])))/deltaz;
      }
      
      divergence1 *= Gamma[a][b];
      sum1        += divergence1;
      
      q_divx (grad1, grad1_front, a,  b, qab_front);
      dAdq(qab_front, dadq_front, a, b);
      rotate_vector(dadq_front, rotated_dadq_front, a, b);
      qab_front2 = qab_front[X]*qab_front[X] + qab_front[Y]*qab_front[Y] + qab_front[Z]*qab_front[Z];
      
      q_divx (grad1_back, grad1, a,  b, qab_back);
      dAdq(qab_back, dadq_back, a, b);
      rotate_vector(dadq_back, rotated_dadq_back, a, b);
      qab_back2 = qab_back[X]*qab_back[X] + qab_back[Y]*qab_back[Y] + qab_back[Z]*qab_back[Z];
      
      q_divy (grad1, grad1_right, a,  b, qab_right);
      dAdq(qab_right, dadq_right, a, b);
      rotate_vector(dadq_right, rotated_dadq_right, a, b);
      qab_right2 = qab_right[X]*qab_right[X] + qab_right[Y]*qab_right[Y] + qab_right[Z]*qab_right[Z];
      
      q_divy (grad1_left, grad1, a,  b, qab_left);
      dAdq(qab_left, dadq_left, a, b);
      rotate_vector(dadq_left, rotated_dadq_left, a, b);
      qab_left2 = qab_left[X]*qab_left[X] + qab_left[Y]*qab_left[Y] + qab_left[Z]*qab_left[Z];
      
      if (DIMENSION != 2) {
        q_divz (grad1, grad1_top, a,  b, qab_top);
        dAdq(qab_top, dadq_top, a, b);
        rotate_vector(dadq_top, rotated_dadq_top, a, b);
        qab_top2 = qab_top[X]*qab_top[X] + qab_top[Y]*qab_top[Y] + qab_top[Z]*qab_top[Z];
        
        q_divz (grad1_bottom, grad1, a,  b, qab_bottom);
        dAdq(qab_bottom, dadq_bottom, a, b);
        rotate_vector(dadq_bottom, rotated_dadq_bottom, a, b);
        qab_bottom2 = qab_bottom[X]*qab_bottom[X] + qab_bottom[Y]*qab_bottom[Y] + qab_bottom[Z]*qab_bottom[Z];
      }
      
      divergence2  =  (ac_front*qab_front2*rotated_dadq_front[X]*(-grad1->phistagg[X][b])   - ac_back*qab_back2*rotated_dadq_back[X]*(-grad1_back->phistagg[X][b]))/deltax;
      divergence2 +=  (ac_right*qab_right2*rotated_dadq_right[Y]*(-grad1->phistagg[Y][b])   - ac_left*qab_left2*rotated_dadq_left[Y]*(-grad1_left->phistagg[Y][b]))/deltay;
      
      if (DIMENSION !=2) {
        divergence2 +=  (ac_top*qab_top2*rotated_dadq_top[Z]*(-grad1->phistagg[Z][b])       - ac_bottom*qab_bottom2*rotated_dadq_bottom[Z]*(-grad1_bottom->phistagg[Z][b]))/deltaz;
      }
      
      
      divergence2 *= Gamma[a][b];
      sum2        += divergence2;
    }
  }
  return (2.0*(sum1 + sum2) + divergence3);
}
double function_A_01_divdAdgradphi_smooth(struct gradlayer **gradient, long index, long gidy, long a) {
  long b;
  double sum1 = 0.0, sum2 = 0.0;
  double divergence1=0.0, divergence2=0.0;
  
  struct gradlayer *grad1;
  struct gradlayer *grad1_left;
  struct gradlayer *grad1_back;
  struct gradlayer *grad1_right;
  struct gradlayer *grad1_front;
  struct gradlayer *grad1_top;
  struct gradlayer *grad1_bottom;
  
  grad1       =  &gradient[1][gidy];
  grad1_left  =  &gradient[1][gidy-1];
  grad1_back  =  &gradient[0][gidy];
  grad1_right =  &gradient[1][gidy+1];
  grad1_front =  &gradient[2][gidy];
  
  if (DIMENSION !=2) {
    grad1_top    =  &gradient[1][gidy+rows_y];
    grad1_bottom =  &gradient[1][gidy-rows_y];
  }
  
  for (b=0; b < NUMPHASES; b++) {
    if (b!=a) {
      divergence1  = (((grad1->phistagg[X][a]*grad1->gradphi[X][b] - grad1->phistagg[X][b]*grad1->gradphi[X][a])*(-grad1->phistagg[X][b])) - ((grad1_back->phistagg[X][a]*grad1_back->gradphi[X][b] - grad1_back->phistagg[X][b]*grad1_back->gradphi[X][a])*(-grad1_back->phistagg[X][b])))/deltax;
      
      divergence1 += (((grad1->phistagg[Y][a]*grad1->gradphi[Y][b] - grad1->phistagg[Y][b]*grad1->gradphi[Y][a])*(-grad1->phistagg[Y][b])) - ((grad1_left->phistagg[Y][a]*grad1_left->gradphi[Y][b] - grad1_left->phistagg[Y][b]*grad1_left->gradphi[Y][a])*(-grad1_left->phistagg[Y][b])))/deltay;
      
      if (DIMENSION !=2) {
         divergence1 += (((grad1->phistagg[Z][a]*grad1->gradphi[Z][b] - grad1->phistagg[Z][b]*grad1->gradphi[Z][a])*(-grad1->phistagg[Z][b])) - ((grad1_bottom->phistagg[Z][a]*grad1_bottom->gradphi[Z][b] - grad1_bottom->phistagg[Z][b]*grad1_bottom->gradphi[Z][a])*(-grad1_bottom->phistagg[Z][b])))/deltay;
      }
      
      divergence1 *= Gamma[a][b];
      sum1        += divergence1;
    }
  }
  return 2.0*(sum1 + sum2);
}

#endif
