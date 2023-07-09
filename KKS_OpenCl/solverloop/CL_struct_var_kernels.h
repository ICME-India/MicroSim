
struct fields {
  double phi[npha];
  double mu[nsol];
  double com[nsol];
  double temperature;
};
struct csle {
  double comie[npha][nsol];
};
struct pfmpar {
  double surfTen[npha*npha];
  double ee[npha*npha];
  double w[npha*npha];
  double wh[npha*npha];
  double eesqrt[npha*npha];
  double w_abc[npha*npha*npha];
  double deltax;
  double deltay;
  double deltaz;
  double deltat;
  double Er;
  double E0;
};
struct pfmval {
  double Rotation_matrix[npha][npha][3][3];
  double Inv_Rotation_matrix[npha][npha][3][3];
  double D[npha][nsol*nsol];
  double DInv[npha][nsol*nsol];
  double cguess[npha*npha][nsol];
  double c_eq[npha*npha][nsol];
  double c_fill[npha*npha][nsol];
  double angle[npha*npha*3];
  double gamma[npha*npha];
  double epsc[npha*npha];
  double Tr;
  double Vm;
  double phisolid;
  double philiquid;
  double Rg;
  double T0;
  double Teq;
  double Tfill;
  double lrep;
  double com_1_Initial;
  double com_2_Initial;
  double com_1_1stguess;
  double com_2_1stguess;
  double a2;
  double rad;
  double intwidth;
  double dxIntfcPoints;
  double dtParam;
  double RefD;
  double TLiquidus;
  double Toffset;
  double TPosOffset;
  double TGRADIENT;
  double velocity;
  double NoiseFac;
  double interfaceUplimit;
  double interfaceDownlimit;
  double deltat_e;
  double rho;
  double damping_factor;
  long shift_OFFSET;
  int thermophase[npha];
  int   nproc;
  int   Nx;
  int   Ny;
  int   Nz;
  int   ntimesteps;
  int   savetime;
  int   myrank;
  int   ISOTHERMAL;
  int atr;
  int Noise_phasefield;
  int Function_anisotropy;
  int anisotropy_type;
  int ELASTICITY;
  int DIMENSION;
};
struct propmatf3 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
  double A[npha][nsol][nsol];
  double DELTA_T[npha][npha];
  double DELTA_C[npha][nsol];
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
};
struct propmatf4 {
  double ceq[npha][npha][nsol];
  double cfill[npha][npha][nsol];
  double slopes[npha][npha][nsol];
  double dcbdT[npha][npha][nsol];
  double cmu[npha][nsol][nsol];
  double muc[npha][nsol][nsol];
  double A[npha][nsol][nsol];
  double DELTA_T[npha][npha];
  double DELTA_C[npha][nsol];
  double dcbdT_phase[npha][nsol];
  double B[npha][nsol];
  double Beq[npha][nsol];
  double dBbdT[npha][nsol];
  double C[npha];
};
struct propmatf4spline {
  double A[npha][nsol][nsol];
  double B[npha][nsol];
  double C[npha];
};
struct symmetric_tensor {
  double xx;
  double yy;
  double zz;
  double yz;
  double xz;
  double xy;
};

struct Stiffness_cubic {
  double C11;
  double C12;
  double C44;
};
struct iter_variables {
 double disp[3][3];
};

