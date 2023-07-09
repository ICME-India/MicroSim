
void CL_initialize_variables() {

  printf("In  CL_initialize_variables\n");
  int ip, is, is1, is2, ip1, ip2, ip3, i1, i2;

  pfmdat.nproc              = numtasks;
  pfmdat.ntimesteps         = ntimesteps;
  pfmdat.savetime           = saveT;
  pfmdat.Tr                 = TLiquidus;
  pfmdat.interfaceDownlimit = 1e-3;
  pfmdat.interfaceUplimit   = 1.0 - pfmdat.interfaceDownlimit;


  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      pfmdat.gamma[ip1*NUMPHASES+ip2] = Gamma[ip1][ip2];
      //printf("gamma: ip1 = %d, ip2 = %d, %le, %le\n", ip1, ip2, pfmdat.gamma[ip1][ip2], Gamma[ip1][ip2]);
    }
  }

  pfmdat.Vm                 = V;
  
  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    for ( is1 = 0; is1 < (NUMCOMPONENTS-1); is1++ ) { 
      for ( is2 =0; is2 < (NUMCOMPONENTS-1); is2++ ) { 
        pfmdat.D[ip][is1*(NUMCOMPONENTS-1)+is2] = Diffusivity[ip][is1][is2];
        //printf("D: ip = %d, is1 = %d, is2 = %d, %le, %le\n", ip, is1, is2, pfmdat.D[ip][is1][is2], Diffusivity[ip][is1][is2]);
      }
    }
  }
  //DiffusivityInv = Malloc3M(NUMPHASES, NUMCOMPONENTS-1, NUMCOMPONENTS-1);  // Done in reading_input_parameters.h
  //matinvnew(Diffusivity,DiffusivityInv,NUMCOMPONENTS-1);                   // Done in reading_input_parameters.h
  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    for ( is1 = 0; is1 < (NUMCOMPONENTS-1); is1++ ) { 
      for ( is2 =0; is2 < (NUMCOMPONENTS-1); is2++ ) { 
        pfmdat.DInv[ip][is1*(NUMCOMPONENTS-1)+is2] = DiffusivityInv[ip][is1][is2];
        //printf("DInv: ip = %d, is1 = %d, is2 = %d, %le, %le\n", ip, is1, is2, pfmdat.DInv[ip][is1][is2], DiffusivityInv[ip][is1][is2]);
        //printf("%le\n", pfmdat.DInv[ip][is1][is2]);
      }
    }
  }

  pfmdat.phisolid           = 1.0;
  pfmdat.philiquid          = 0.0;
  pfmdat.Rg                 = R;
  
  pfmdat.T0                 = T;
  
  pfmdat.Teq                 = Teq;
  pfmdat.Tfill               = Tfill;
  pfmdat.ISOTHERMAL          = ISOTHERMAL;
  
  pfmdat.lrep               = 1.0;//(pfmdat.gamma*pfmdat.Vm)/(TLiquidus*R);//1.0;//CharLength*1e2;
  
  pfmdat.a2                 = 47.0/60.0; // phi^3(10-15phi+6phi^2)
  //pfmdat.a2                 = 5.0/6.0; // Abhik hphi = 3 phi phi - 2 phi phi phi +..
  //pfmdat.a2                 = 1.0; // hphi = phi
  //pfmdat.a2                 = 0.785398; //Nele  phi phi / (sum phi phi)

  pfmdat.rad                = fill_cylinder_parameters.radius;
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      if ( FUNCTION_ANISOTROPY ) { 
        pfmdat.epsc[ip1*NUMPHASES+ip2] = dab[ip1][ip2]; 
        //printf("epsc: ip1 = %d, ip2 = %d, %le, %le\n", ip1, ip2, pfmdat.epsc[ip1][ip2], dab[ip1][ip2]);
      }
      else { 
        pfmdat.epsc[ip1*NUMPHASES+ip2] = 0.0;
      }
    }
  }
  
  pfmdat.intwidth           = epsilon;
  pfmdat.dxIntfcPoints      = 1;//NoOfPointsInInterface;//Not Used in NSM
  pfmdat.dtParam            = 1;//dtDenomParam; //Not Used in NSM
  
  pfmdat.TLiquidus          = TLiquidus;
  pfmdat.Toffset            = temperature_gradientY.base_temp;
  //pfmdat.PosOffset          = PositionOffset;
  pfmdat.TPosOffset         = temperature_gradientY.gradient_OFFSET/deltay;
  //pfmdat.TG                 = TemperatureGradient;
  pfmdat.TGRADIENT          = temperature_gradientY.DeltaT/temperature_gradientY.Distance;
  //pfmdat.Vp                 = PullingVelocity;
  pfmdat.velocity           = temperature_gradientY.velocity;
  pfmdat.NoiseFac           = AMP_NOISE_PHASE;
  pfmdat.atr                = atr;
  pfmdat.Function_anisotropy = FUNCTION_ANISOTROPY;
  //pfmdat.anisotropy_type     = Anisotropy_type;
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      for ( i1 = 0; i1 < 3; i1++ ) { 
        if ( FUNCTION_ANISOTROPY ) { 
          pfmdat.angle[(ip1*NUMPHASES+ip2)*NUMPHASES+i1] = RotAngles[ip1][ip2][i1]*M_PI/180.0;
          //printf("angle: ip = %d, is1 = %d, is2 = %d, %le, %le\n", ip1, ip2, i1, pfmdat.angle[ip1][ip2][i1]/(M_PI/180.0), RotAngles[ip1][ip2][i1]);
        }
        else { 
          pfmdat.angle[(ip1*NUMPHASES+ip2)*NUMPHASES+i1] = 0.0;
        }
      }
    }
  }
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      for ( i1 = 0; i1 < 3; i1++ ) { 
        for ( i2 = 0; i2 < 3; i2++ ) { 
          if ( FUNCTION_ANISOTROPY ) { 
            pfmdat.Rotation_matrix[ip1][ip2][i1][i2] = Rotation_matrix[ip1][ip2][i1][i2]; 
            pfmdat.Inv_Rotation_matrix[ip1][ip2][i1][i2] = Inv_Rotation_matrix[ip1][ip2][i1][i2]; 
          }
          else { 
            pfmdat.Rotation_matrix[ip1][ip2][i1][i2] = 0.0; 
            pfmdat.Inv_Rotation_matrix[ip1][ip2][i1][i2] = 0.0; 
          }
        }
      }
    }
  }

  //pfmdat.angle              = RotAngles[0]*M_PI/180.0;
  pfmdat.shift_OFFSET       = shift_OFFSET; /* */

  //pfmdat.RefD = pfmdat.D[1][0][0];
  pfmdat.RefD = pfmdat.D[0][0];
  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    for ( is1 = 0; is1 < (NUMCOMPONENTS-1); is1++ ) { 
      for ( is2 =0; is2 < (NUMCOMPONENTS-1); is2++ ) { 
        if ( pfmdat.RefD < pfmdat.D[ip][is1*(NUMCOMPONENTS-1)+is2])  { 
          pfmdat.RefD = pfmdat.D[ip][is1*(NUMCOMPONENTS-1)+is2]; 
        }
      }
    }
  }
  pfmdat.RefD = 1.0;
  //printf("RefD: %le\n", pfmdat.RefD);

  pfmvar.E0                 = 1.0;//8.314*pfmdat.Tr/pfmdat.Vm;
  pfmvar.deltax             = deltax; 
  pfmvar.deltay             = deltay;
  pfmvar.deltaz             = deltaz;
  pfmvar.deltat             = deltat; 
  pfmvar.Er                 = 8.314*pfmdat.Tr;
  //pfmvar.dx                 = pfmvar.deltax*pfmdat.lrep; // m
  //pfmvar.dy                 = pfmvar.deltay*pfmdat.lrep; // m
  //pfmvar.dt                 = pfmvar.deltat*pfmdat.lrep*pfmdat.lrep/(pfmdat.RefD);

  //printf("deltax = %le, deltat = %le \n", pfmvar.deltax, pfmvar.deltat);
  //printf("dx = %le, dt = %le \n", pfmvar.dx, pfmvar.dt);
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2 ++ ) {  
      if ( ip1 != ip2 ) { 
      pfmvar.surfTen[ip1*NUMPHASES+ip2] = pfmdat.gamma[ip1*NUMPHASES+ip2];
      
      //pfmvar.w[ip1*NUMPHASES+ip2]       = (3.0*pfmvar.surfTen[ip1*NUMPHASES+ip2]*4.39445)/pfmdat.intwidth;
      pfmvar.w[ip1*NUMPHASES+ip2]       = (6.0*pfmvar.surfTen[ip1*NUMPHASES+ip2]*2.94)/pfmdat.intwidth;
      pfmvar.wh[ip1*NUMPHASES+ip2]  = pfmvar.w[ip1*NUMPHASES+ip2];
      
      //pfmvar.ee[ip1*NUMPHASES+ip2]      = 3.0*pfmvar.surfTen[ip1*NUMPHASES+ip2]*pfmdat.intwidth/(4.39445);
      pfmvar.ee[ip1*NUMPHASES+ip2]      = 3.0*pfmvar.surfTen[ip1*NUMPHASES+ip2]*pfmdat.intwidth/(2.0*2.94);
      pfmvar.eesqrt[ip1*NUMPHASES+ip2]  = sqrt(pfmvar.ee[ip1*NUMPHASES+ip2]);
      
      printf("ee: ip1 = %d, ip2 = %d, %le\n", ip1, ip2, pfmvar.ee[ip1*NUMPHASES+ip2]);
      }
      else {
      pfmvar.surfTen[ip1*NUMPHASES+ip2] = 0.0;
      pfmvar.w[ip1*NUMPHASES+ip2]       = 0.0;
      pfmvar.wh[ip1*NUMPHASES+ip2]  = pfmvar.w[ip1*NUMPHASES+ip2];
      pfmvar.ee[ip1*NUMPHASES+ip2]      = 0.0;
      pfmvar.eesqrt[ip1*NUMPHASES+ip2]  = 0.0;

      }
    }
  }
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      for ( ip3 = 0; ip3 < NUMPHASES; ip3++ ) { 
        if ( npha > 2 ) {
          pfmvar.w_abc[(ip1*NUMPHASES+ip2)*NUMPHASES+ip3] = Gamma_abc[ip1][ip2][ip3]; 
          //printf("Gammaabc %le, %le\n", pfmvar.w_abc[(ip1*NUMPHASES+ip2)*NUMPHASES+ip3], Gamma_abc[ip1][ip2][ip3]);
        }
        else { 
          pfmvar.w_abc[(ip1*NUMPHASES+ip2)*NUMPHASES+ip3] = 0.0; 
        }
      }
    }
  }

  
  printf("****************************\n");

  for ( ip = 0; ip < NUMPHASES; ip++ ) { 
    pfmdat.thermophase[ip] = thermo_phase[ip];
    printf("ip = %d, thermo_phase[%d]\n", ip, ip);
  }
  
  for ( ip1 = 0; ip1 < NUMPHASES; ip1++ ) { 
    for ( ip2 = 0; ip2 < NUMPHASES; ip2++ ) { 
      for ( is = 0; is < (NUMCOMPONENTS-1); is++ ) { 
        pfmdat.cguess[ip1*NUMPHASES+ip2][is] = c_guess[ip1][ip2][is];
          pfmdat.c_eq[ip1*NUMPHASES+ip2][is] =     ceq[ip1][ip2][is];
        pfmdat.c_fill[ip1*NUMPHASES+ip2][is] =   cfill[ip1][ip2][is];
      }
    }
  }

  pfmdat.ELASTICITY = ELASTICITY;
  pfmdat.deltat_e = deltat_e;
  pfmdat.rho = rho;
  pfmdat.damping_factor = damping_factor;
  pfmdat.DIMENSION = DIMENSION;

  nx = mpiparam.rows_x;
  ny = mpiparam.rows_y;
  nz = mpiparam.rows_z;

  pfmdat.Nx = nx;
  pfmdat.Ny = ny;
  pfmdat.Nz = nz;

  printf("rank = %d, nx = %d, ny = %d, nz = %d\n", rank, nx, ny, nz);

  //nxtotny = nx*totny;
  nynz = ny*nz;
  nxnynz = nx*ny*nz;
  
  globaldim[0] = ny;
  globaldim[1] = nz;
  globaldim[2] = nx;

  //globaldim[0] = globaldim0;
  //globaldim[1] = globaldim1;

  // gridNew  = (struct grid*)malloc(nxny*sizeof(struct grid));
  // gridOld  = (struct grid*)malloc(nxny*sizeof(struct grid));
  cscl     = (struct csle*)malloc(nxnynz*sizeof(struct csle));
  //temp     = (double*)malloc(nx*sizeof(double)); //Changed to nx, According to MESH_Y****
  tstep    = (long*)malloc(sizeof(long));
  propf4spline     = (struct propmatf4spline*)malloc(ny*sizeof(struct propmatf4spline));
  propf4spline1     = (struct propmatf4spline*)malloc(ny*sizeof(struct propmatf4spline));

  pfmdat.myrank = rank;

  t = STARTTIME;//1;
  tstart = STARTTIME;//1;
  tstep[0] = STARTTIME;//1;

  printf("Out CL_initialize_variables\n");
}
