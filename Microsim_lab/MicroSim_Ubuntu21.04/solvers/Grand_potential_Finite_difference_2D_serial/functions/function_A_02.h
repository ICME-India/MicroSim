#ifdef ANISOTROPY_POTENTIAL 
      scalprod1 = 0.0;
      
      for (dir=0; dir < 2; dir++)  {
	scalprod1 += (phi[a]*grad1->gradphi_c[dir][b] - phi[b]*grad1->gradphi_c[dir][a])*grad1->gradphi_c[dir][b];
      }
      scalprod1 *= Gamma[a][b];
      
      sum1 += scalprod1;
#endif
#ifdef ANISOTROPY_GRADIENT_POTENTIAL
      dAdq(qab, dadq, a, b);
      qab2 = (qab[X]*qab[X] + qab[Y]*qab[Y]);
      ac1  = function_ac(qab, a, b);
      
      scalprod1 =0.0;
      scalprod2 =0.0;
      
      for (dir=0; dir < 2; dir++) {
	scalprod1 += (phi[a]*grad1->gradphi_c[dir][b] - phi[b]*grad1->gradphi_c[dir][a])*grad1->gradphi_c[dir][b];
      }
      
      scalprod1 *= Gamma[a][b]*ac1;
      
      dqdphi[X] = grad1->gradphi_c[X][b];
      dqdphi[Y] = grad1->gradphi_c[Y][b];
      
      multiply(Rotation_matrix[a][b], dqdphi, Rotated_dqdphi, DIMENSION);
      
      for(dir=0; dir<DIMENSION; dir++) {
	dqdphi[dir] = Rotated_dqdphi[dir];
      }
     
      for (dir=0; dir < 2; dir++)  {
	scalprod2 += (dadq[dir])*dqdphi[dir];
      }
     
      scalprod2 *= Gamma[a][b]*qab2;
      
      sum1 += scalprod1;
      sum2 += 0.5*scalprod2;   
      
#endif