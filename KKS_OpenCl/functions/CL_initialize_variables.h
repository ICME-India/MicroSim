void CL_initialize_variables() {

  pfmdat.nproc              = 1;
  pfmdat.jDimX              = MESH_X;
  pfmdat.iDimY              = MESH_Y;
  pfmdat.ntimesteps         = ntimesteps;
  pfmdat.savetime           = saveT;
  pfmdat.Tr                 = TLiquidus;
  pfmdat.sigma              = Gamma[0][1];
  pfmdat.Vm                 = V;
  pfmdat.D11s               = Diffusivity[0][0][0];
  pfmdat.D11l               = Diffusivity[1][0][0];

  pfmdat.phisolid           = 1.0;
  pfmdat.philiquid          = 0.0;
  pfmdat.Rg                 = R;
  pfmdat.T0                 = T;
  pfmdat.lrep               = (pfmdat.sigma*pfmdat.Vm)/(TLiquidus*R);//1.0;//CharLength*1e2;
  pfmdat.c1s_Initial        = cfill[0][0][0];
  pfmdat.c1l_Initial        = cfill[1][0][0];
  pfmdat.c1s_1stguess       = ceq[0][0][0];
  pfmdat.c1l_1stguess       = ceq[1][0][0];

  pfmdat.a2                 = 47.0/60.0;
  pfmdat.rad                = fill_cylinder_parameters.radius;
  pfmdat.epsc               = dab[0][1];
  pfmdat.intwidth           = epsilon;
  pfmdat.dxIntfcPoints      = 1;//NoOfPointsInInterface;//Not Used in NSM
  pfmdat.dtParam            = 1;//dtDenomParam; //Not Used in NSM
  pfmdat.InterfaceMobility  = -1.0;
  pfmdat.epsm               = 0.0;
  pfmdat.TLiquidus          = TLiquidus;
  pfmdat.Toffset            = T_Offset;
  pfmdat.PosOffset          = PositionOffset;
  pfmdat.TG                 = TemperatureGradient;
  pfmdat.Vp                 = PullingVelocity;
  pfmdat.NoiseFac           = AMP_NOISE_PHASE;
  pfmdat.angle              = RotAngles[0]*M_PI/180.0;

  pfmdat.RefD = pfmdat.D11l;

  pfmvar.E0                 = 8.314*pfmdat.Tr/pfmdat.Vm;
  pfmvar.deltax             = deltax; 
  pfmvar.deltay             = deltay;
  pfmvar.deltat             = deltat; 
  pfmvar.Er                 = 8.314*pfmdat.Tr;
  pfmvar.dx                 = pfmvar.deltax*pfmdat.lrep; // m
  pfmvar.dy                 = pfmvar.deltay*pfmdat.lrep; // m
  pfmvar.dt                 = pfmvar.deltat*pfmdat.lrep*pfmdat.lrep/(pfmdat.RefD);

  printf("deltax = %le, deltat = %le \n", pfmvar.deltax, pfmvar.deltat);
  printf("dx = %le, dt = %le \n", pfmvar.dx, pfmvar.dt);
    
  pfmvar.surfTen            = pfmdat.sigma/(pfmdat.lrep*pfmvar.E0);//Nothing but 1
  pfmvar.ee                 = (6.0*pfmvar.surfTen*pfmdat.intwidth)/4.39445;
  pfmvar.w                  = (3.0*pfmvar.surfTen*4.39445)/pfmdat.intwidth;
  pfmvar.eesqrt             = sqrt(pfmvar.ee);
  pfmvar.IntMob             = pfmdat.InterfaceMobility*pfmdat.lrep*pfmvar.E0/pfmdat.RefD;

  printf("surfTen = %le, ee = %le, w = %le \n", pfmvar.surfTen, pfmvar.ee, pfmvar.w);
  printf("sigma = %le, intwidth = %le \n", pfmdat.sigma, pfmdat.intwidth);

  printf("sigma=%le,lrep=%le,E0=%le,intwidth=%le,RefD=%le\n", pfmdat.sigma, pfmdat.lrep, pfmvar.E0, pfmdat.intwidth,pfmdat.RefD);

  if ( pfmdat.InterfaceMobility < 0.0 ) { 
    pfmvar.IntMobInv = 0.0;
  }
  else {
    pfmvar.IntMobInv = 1.0/pfmvar.IntMob;
  }

  NX = pfmdat.jDimX;
  NY = pfmdat.iDimY;

  NXNY = NX*NY;

  totny = NY + 2*pfmdat.nproc;
  totsliceny = (NY/pfmdat.nproc)+2;

  nx = NX+2;
  ny = totsliceny;

  pfmdat.jNx = nx;
  pfmdat.iNy = ny;

  nxtotny = nx*totny;
  nxny = nx*ny;

  globaldim0 = nx;
  globaldim1 = ny;

  globaldim[0] = globaldim0;
  globaldim[1] = globaldim1;

  gridNew  = (struct grid*)malloc(nxny*sizeof(struct grid));
  gridOld  = (struct grid*)malloc(nxny*sizeof(struct grid));
  cscl     = (struct csle*)malloc(nxny*sizeof(struct csle));
  temp     = (double*)malloc(ny*sizeof(double));
  tstep    = (int*)malloc(sizeof(int));

  //pfmdat.myrank = rank;
  pfmdat.myrank = 0;
  istart = 0;
  iend = ny;

  t = 1;
  tstart = 1;
  tstep[0] = 1;

  //deltax = pfmvar.deltax;
  //deltay = deltax;
  //deltaz = deltax;
  //deltat = pfmvar.deltat;

}
