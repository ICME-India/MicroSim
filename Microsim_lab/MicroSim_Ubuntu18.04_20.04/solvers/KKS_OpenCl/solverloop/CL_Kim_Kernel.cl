#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "solverloop/GibbsEnergyData.h"

/*
      0          ---- j ----->           nx (X)
    
 0    #  #  #  #  #  #  #  #  #  #  #  #  # 
               
      #  #  #  #  #  #  #  #  #  #  #  #  #
               
 |    #  #  #  #  #  #  #  #  #  #  #  #  # 
 |             
 |    #  #  #  #  0  1  2  #  #  #  #  #  # 
               
 i    #  #  #  #  3  4  5  #  #  #  #  #  # 
               
 |    #  #  #  #  6  7  8  #  #  #  #  #  # 
 |             
 V    #  #  #  #  #  #  #  #  #  #  #  #  # 
               
      #  #  #  #  #  #  #  #  #  #  #  #  # 
               
 ny   #  #  #  #  #  #  #  #  #  #  #  #  # 
 
 (Y)

*/

struct grid {
    double phi;
    double c1;
};
struct csle {
    double c1l;
    double c1s;
};
struct pfmpar {
    double E0;
    double surfTen;
    double ee;
    double w;
    double deltax;
    double deltay;
    double deltat;
    //double dx;
    //double dy;
    //double dt;
    double eesqrt;
    double Er;
    double IntMob;
    double IntMobInv;
};
struct pfmval {
    double Tr;
    double sigma;
    double Vm;
    double D11l;
    double D11s;
    double phisolid;
    double philiquid;
    double Rg;
    double T0;
    double lrep;
    double c1l_Initial;
    double c1s_Initial;
    double c1l_1stguess;
    double c1s_1stguess;
    double a2;
    double rad;
    double epsc;
    double intwidth;
    double dxIntfcPoints;
    double dtParam;
    double epsm;
    double InterfaceMobility;
    double RefD;
    double angle;
  double TLiquidus;
  double Toffset;
  //double PosOffset;
  //double TG;
  //double Vp;
  double TPosOffset; 
  double TGRADIENT;
  double velocity;
  double NoiseFac;
  long shift_OFFSET;
    int   nproc;
    int   jNx;
    int   iNy;
    int   jDimX;
    int   iDimY;
    int   ntimesteps;
    int   savetime;
    int   myrank;
};

__kernel void SolverCsClEq(__global struct grid *gridOld, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp) {
    
    int i;
    int j;
    int nx;
    int ny;
    int index;
    int count;

    double dgldx1l;
    double dgsdx1s;
    double ddgldx1ldx1l;
    double ddgsdx1sdx1s;
    
    double hphi;
    double x1, x1l, x1s;
    double detJac;
    double fi;
 
    double ca[2];
    double c_new[2];
    double dumcn[2];
    double fun[2];
    double Jac[4];
    double JacInv[4];
    
    double Ti;

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);
    
    Ti = temp[j];

 
    fi = gridOld[index].phi;
    x1 = gridOld[index].c1;

    x1l = cscl[index].c1l;
    x1s = cscl[index].c1s;
    
    ca[0] = x1l;
    ca[1] = x1s;

    c_new[0] = ca[0];
    c_new[1] = ca[1];
    
    dumcn[0] = ca[0];
    dumcn[1] = ca[1];
    
    hphi = fi;
    
    count=0;
    do {

        count=count+1;

        dgldx1l = (dGLIQdX1L(Ti, x1l))/(pfmvar->Er);
        dgsdx1s = (dGSOLdX1S(Ti, x1s))/(pfmvar->Er);

        ddgldx1ldx1l = (ddGLIQdX1LdX1L(Ti, x1l))/(pfmvar->Er);

        ddgsdx1sdx1s = (ddGSOLdX1SdX1S(Ti, x1s))/(pfmvar->Er);

        fun[0] = hphi*x1s + (1.0-hphi)*x1l - x1;
        fun[1] = dgldx1l - dgsdx1s;

        Jac[0] = 1.0-hphi;
        Jac[1] = hphi;
        Jac[2] = ddgldx1ldx1l;
        Jac[3] = -(ddgsdx1sdx1s);

        detJac = (Jac[0]*Jac[3]-Jac[1]*Jac[2]);

        JacInv[0] = Jac[3]/detJac;
        JacInv[1] = -Jac[1]/detJac;
        JacInv[2] = -Jac[2]/detJac;
        JacInv[3] = Jac[0]/detJac;

        dumcn[0] = c_new[0];
        dumcn[1] = c_new[1];

        c_new[0] = ca[0] - ( fun[0]*JacInv[0] + fun[1]*JacInv[1] );
        c_new[1] = ca[1] - ( fun[0]*JacInv[2] + fun[1]*JacInv[3] );

        ca[0] = c_new[0];
        ca[1] = c_new[1];

        x1l = ca[0];  
        x1s = ca[1];
		
		if (count > 1000) {
		    printf("Too many iterations ~ %d\n", count);
			printf("Exited the loop. Results may not be correct\n");
			printf("Stop the code and choose correct parameters\n");
		    break;
		}


	} while ( (fabs(dumcn[0]-ca[0]) > 1e-6) && (fabs(dumcn[1]-ca[1]) > 1e-6) );
	
	

    cscl[index].c1l = ca[0];
    cscl[index].c1s = ca[1];

    //printf("%le,%le\n",cscl[index].c1l,cscl[index].c1s);
 
}

__kernel void SolverCWoatr(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp) { 

    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;

    double hphi[5];
    double term1lx;
    double term1sx;
    double term1ly;
    double term1sy;
    double term1l;
    double term1s;
    double c1dot;

    struct grid stgridO[9];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        hphi[0] = stgridO[1].phi;
        hphi[1] = stgridO[3].phi;
        hphi[4] = stgridO[4].phi;
        hphi[2] = stgridO[5].phi;
        hphi[3] = stgridO[7].phi;

        term1lx = ( ((1.0-hphi[2])+(1.0-hphi[4]))*(stcscl[2].c1l-stcscl[4].c1l) - ((1.0-hphi[4])+(1.0-hphi[1]))*(stcscl[4].c1l-stcscl[1].c1l) )/(2.0*(pfmvar->deltax)*(pfmvar->deltax));
	    term1ly = ( ((1.0-hphi[3])+(1.0-hphi[4]))*(stcscl[3].c1l-stcscl[4].c1l) - ((1.0-hphi[4])+(1.0-hphi[0]))*(stcscl[4].c1l-stcscl[0].c1l) )/(2.0*(pfmvar->deltay)*(pfmvar->deltay));


        term1sx = ( ((hphi[2])+(hphi[4]))*(stcscl[2].c1s-stcscl[4].c1s) - ((hphi[4])+(hphi[1]))*(stcscl[4].c1s-stcscl[1].c1s) )/(2.0*(pfmvar->deltax)*(pfmvar->deltax));
        term1sy = ( ((hphi[3])+(hphi[4]))*(stcscl[3].c1s-stcscl[4].c1s) - ((hphi[4])+(hphi[0]))*(stcscl[4].c1s-stcscl[0].c1s) )/(2.0*(pfmvar->deltay)*(pfmvar->deltay));

        term1l = term1lx + term1ly;
        term1s = term1sx + term1sy;

        c1dot = ( pfmdat->D11l*term1l + pfmdat->D11s*term1s ) / pfmdat->RefD;

        gridNew[( i )*nx + ( j )].c1 = stgridO[4].c1 + (pfmvar->deltat)*( c1dot );
    }
}

__kernel void SolverPhiIso(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp) {   
  
    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double div_phi;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_iso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Ti;

    double tmp1;
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        //printf("%le,%le,%le,%le\n",cliq1, csol1, cscl[( i )*nx + ( j )].c1l, cscl[( i )*nx + ( j )].c1s);

        D11invL = 1.0/pfmdat->D11l;

        gl = (GLIQ(Ti, cliq1))/(pfmvar->Er);
        gs = (GSOL(Ti, csol1))/(pfmvar->Er);

        dgldc1l = (dGLIQdX1L(Ti, cliq1))/(pfmvar->Er);

        ddgldc1ldc1l = (ddGLIQdX1LdX1L(Ti, cliq1))/(pfmvar->Er);

        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l);

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;

        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

        

        Mphi = (pfmvar->w)/(3.0*(pfmvar->ee)*(pfmdat->a2)*zeta);
        
        ////For FiniteAnisoPFMobility
        //Mphi_FiniteMobility = sqrt(pfmvar->w)/( (3.0*sqrt(2.0*pfmvar->ee))*( (pfmvar->IntMobInv) + ((pfmdat->a2*sqrt(pfmvar->ee)*zeta)/(sqrt(2.0*pfmvar->w))) ) );
        //Mphi_FiniteMobilityAniso = Mphi_FiniteMobility;
        //Mphi = Mphi_FiniteMobilityAniso;

        phidot_iso =  Mphi*((pfmvar->ee)*div_phi - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_iso;
        
        //printf("%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", Mphi, pfmvar->w, pfmvar->ee, pfmdat->a2, zeta, cliq1, csol1, ddgldc1ldc1l);
    }
}

__kernel void SolverPhi(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar, __global double *temp) {   
  
    int i;
    int j;
    int k;
    int ii;    
    int nx;
    int ny;
    int nz;
    int index;

    double gphi;
    double dgphidphi;
    double dhphidphi; 
    double DetDL;
    double D11invL;
    double gl;
    double gs;
    double dgldc1l;
    double ddgldc1ldc1l;
    double zeta;
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double div_phi;
    double angleig;
    double Nnx;
    double Nny;
    double neta[5];
    double etam[5];
    double B[5];
    double neta_x;
    double neta_y;
    double alp1;
    double eep_x;
    double eep_y;
    double alp2;
    double alp3;
    double alp;
    double Mphi;
    double Mphi_aniso;
    double Mphi_FiniteMobility;
    double Mphi_FiniteMobilityAniso;
    double phidot_aniso;
    
    double cliq1;
    double csol1;
    double cliq2;
    double csol2;

    double Ti;

    double tmp1;
    
    struct grid stgridO[9];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        Ti = temp[j];

        index = (i)*nx + (j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        cliq1 = cscl[( i )*nx + ( j )].c1l;
        csol1 = cscl[( i )*nx + ( j )].c1s;

        D11invL = 1.0/pfmdat->D11l;

        gl = (GLIQ(Ti, cliq1))/(pfmvar->Er);
        gs = (GSOL(Ti, csol1))/(pfmvar->Er);

        dgldc1l = (dGLIQdX1L(Ti, cliq1))/(pfmvar->Er);

        ddgldc1ldc1l = (ddGLIQdX1LdX1L(Ti, cliq1))/(pfmvar->Er);

        zeta = ((cliq1-csol1)*(cliq1-csol1)*ddgldc1ldc1l);

        gphi = stgridO[4].phi*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-stgridO[4].phi);
        dgphidphi = 2.0*stgridO[4].phi*(1.0-stgridO[4].phi)*(1.0-2.0*stgridO[4].phi);
        dhphidphi = 30.0*gphi;

        div_phi = ( 4.0*(stgridO[1].phi+stgridO[3].phi+stgridO[7].phi+stgridO[5].phi ) + stgridO[0].phi+stgridO[6].phi+stgridO[8].phi+stgridO[2].phi - 20.0*stgridO[4].phi ) / ( 6.0*(pfmvar->deltax)*(pfmvar->deltax) );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        for (ii = 0; ii <5; ii++) {
            if ( ( phix[ii] == 0.0 ) && ( phiy[ii] == 0.0 ) ) {

                angleig = (pfmdat->angle);
                Nnx = cos(angleig);
                Nny = sin(angleig);

                tmp1 = ( Nnx*Nnx*Nnx*Nnx + Nny*Nny*Nny*Nny );

                neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) ) * ( 1.0 + ( (4.0*(pfmdat->epsc))/(1.0-3.0*(pfmdat->epsc)) )*tmp1 );

                B[ii] = neta[ii] * 16.0 * (pfmdat->epsc) *  Nnx*Nny * ( Nny*Nny - Nnx*Nnx );

                etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) ) * ( 1.0 + ( (4.0*(pfmdat->epsm))/(1.0-3.0*(pfmdat->epsm)) )*tmp1 );

            }
            else if ( phix[ii] == 0.0 ) {

                angleig = asin(1.0) + (pfmdat->angle);
                Nnx = cos(angleig);
                Nny = sin(angleig);

                tmp1 = ( Nnx*Nnx*Nnx*Nnx + Nny*Nny*Nny*Nny );

                neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) ) * ( 1.0 + ( (4.0*(pfmdat->epsc))/(1.0-3.0*(pfmdat->epsc)) )*tmp1 );

                B[ii] = neta[ii] * 16.0 * (pfmdat->epsc) *  Nnx*Nny * ( Nny*Nny - Nnx*Nnx );

                etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) ) * ( 1.0 + ( (4.0*(pfmdat->epsm))/(1.0-3.0*(pfmdat->epsm)) )*tmp1 );
            }
            else {

                angleig = atan( phiy[ii] / phix[ii] ) + (pfmdat->angle);
                Nnx = cos(angleig);
                Nny = sin(angleig);

                tmp1 = ( Nnx*Nnx*Nnx*Nnx + Nny*Nny*Nny*Nny );

                neta[ii] = ( 1.0 - 3.0*(pfmdat->epsc) ) * ( 1.0 + ( (4.0*(pfmdat->epsc))/(1.0-3.0*(pfmdat->epsc)) )*tmp1 );

                B[ii] = neta[ii] * 16.0 * (pfmdat->epsc) *  Nnx*Nny * ( Nny*Nny - Nnx*Nnx );

                etam[ii] = ( 1.0 - 3.0*(pfmdat->epsm) ) * ( 1.0 + ( (4.0*(pfmdat->epsm))/(1.0-3.0*(pfmdat->epsm)) )*tmp1 );
            
            }
        }

        neta_x = ( neta[1] - neta[2] ) / ( (pfmvar->deltax) );
        neta_y = ( neta[3] - neta[4] ) / ( (pfmvar->deltay) );

        alp1 = neta[0]*neta[0]*div_phi + 2.0*neta[0] * ( neta_x * gradx_phi[0] + neta_y * grady_phi[0] );

        eep_x = ( B[1] - B[2] ) / ( (pfmvar->deltax) );
        eep_y = ( B[3] - B[4] ) / ( (pfmvar->deltay) );

        alp2 = eep_x * grady_phi[0];
        alp3 = eep_y * gradx_phi[0];

        alp = alp1 - alp2 + alp3;

        Mphi = (pfmvar->w)/(3.0*neta[0]*neta[0]*(pfmvar->ee)*(pfmdat->a2)*zeta);
	    
	    ////For FiniteAnisoPFMobility
	    //Mphi_FiniteMobility = sqrt(pfmvar->w)/( (3.0*sqrt(2.0*pfmvar->ee)*neta[0])*( (pfmvar->IntMobInv) + ((pfmdat->a2*sqrt(pfmvar->ee)*neta[0]*zeta)/(sqrt(2.0*pfmvar->w))) ) );
        //Mphi_FiniteMobilityAniso = etam[0]*Mphi_FiniteMobility;
        //Mphi = Mphi_FiniteMobilityAniso;

        phidot_aniso =  Mphi*((pfmvar->ee)*alp - (pfmvar->w)*dgphidphi - dhphidphi*(gs-gl-(csol1-cliq1)*dgldc1l));

        gridNew[ ( ( i )*nx + ( j ) ) ].phi = stgridO[4].phi + (pfmvar->deltat)*phidot_aniso;
		//printf("%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le,\t%le\n", Mphi, pfmvar->w, pfmvar->ee, pfmdat->a2, zeta, cliq1, csol1, ddgldc1ldc1l);
		

    }
}

__kernel void SolverCatr(__global struct grid *gridOld, __global struct grid *gridNew, __global struct csle *cscl, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;
    int j;
    int k;
    int jj;
    int nx;
    int ny;
    int nz;
    int index;

    double hphi[5];
    double A1;
    double alpha1[5];
    double dphidt[5];
    double gradx_phi[3];
    double grady_phi[3];
    double phix[5];
    double phiy[5];
    double modgradphi[5];
    double alphidot1[4];
    double jat1[4];
    double c1jatx;
    double c1jaty;
    double c1jat;
    double term1lx;
    double term1ly;
    double term1l;
    double c1dot;
    
    struct grid stgridO[9];
    struct grid stgridN[5];
    struct csle stcscl[5];

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i != 0 && i != ny-1 && j != 0 && j != nx-1 ) {

        index = (i*nx + j);

        stgridO[0] = gridOld[(i-1)*nx + (j-1)];
        stgridO[1] = gridOld[(i-1)*nx + ( j )];
        stgridO[2] = gridOld[(i-1)*nx + (j+1)];
        stgridO[3] = gridOld[( i )*nx + (j-1)];
        stgridO[4] = gridOld[( i )*nx + ( j )];
        stgridO[5] = gridOld[( i )*nx + (j+1)];
        stgridO[6] = gridOld[(i+1)*nx + (j-1)];
        stgridO[7] = gridOld[(i+1)*nx + ( j )];
        stgridO[8] = gridOld[(i+1)*nx + (j+1)];

        hphi[0] = stgridO[1].phi;
        hphi[1] = stgridO[3].phi;
        hphi[4] = stgridO[4].phi;
        hphi[2] = stgridO[5].phi;
        hphi[3] = stgridO[7].phi;

        stcscl[0] = cscl[(i-1)*nx + ( j )];
        stcscl[1] = cscl[( i )*nx + (j-1)];
        stcscl[4] = cscl[( i )*nx + ( j )];
        stcscl[2] = cscl[( i )*nx + (j+1)];
        stcscl[3] = cscl[(i+1)*nx + ( j )];

        stgridN[0].phi = gridNew[(i-1)*nx + ( j )].phi;
	    stgridN[1].phi = gridNew[( i )*nx + (j-1)].phi;
	    stgridN[4].phi = gridNew[( i )*nx + ( j )].phi;
	    stgridN[2].phi = gridNew[( i )*nx + (j+1)].phi;
	    stgridN[3].phi = gridNew[(i+1)*nx + ( j )].phi;

        dphidt[0] = ( stgridN[4].phi - stgridO[4].phi ) / ( (pfmvar->deltat) );
	    dphidt[1] = ( stgridN[2].phi - stgridO[5].phi ) / ( (pfmvar->deltat) );
	    dphidt[2] = ( stgridN[1].phi - stgridO[3].phi ) / ( (pfmvar->deltat) );
	    dphidt[3] = ( stgridN[3].phi - stgridO[7].phi ) / ( (pfmvar->deltat) );
	    dphidt[4] = ( stgridN[0].phi - stgridO[1].phi ) / ( (pfmvar->deltat) );

        A1 = (sqrt((pfmvar->ee)))/(sqrt(2.0*(pfmvar->w)));

        alpha1[0] = (A1)*( stcscl[4].c1l - stcscl[4].c1s );
        alpha1[1] = (A1)*( stcscl[2].c1l - stcscl[2].c1s );
        alpha1[2] = (A1)*( stcscl[1].c1l - stcscl[1].c1s );
        alpha1[3] = (A1)*( stcscl[3].c1l - stcscl[3].c1s );
        alpha1[4] = (A1)*( stcscl[0].c1l - stcscl[0].c1s );

        gradx_phi[0] = ( stgridO[5].phi - stgridO[3].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[1] = ( stgridO[8].phi - stgridO[6].phi ) / ( 2.0*(pfmvar->deltax) );
        gradx_phi[2] = ( stgridO[2].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltax) );

        grady_phi[0] = ( stgridO[7].phi - stgridO[1].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[1] = ( stgridO[8].phi - stgridO[2].phi ) / ( 2.0*(pfmvar->deltay) );
        grady_phi[2] = ( stgridO[6].phi - stgridO[0].phi ) / ( 2.0*(pfmvar->deltay) );

        phix[0] = gradx_phi[0];
        phix[1] = ( stgridO[5].phi - stgridO[4].phi ) / ( (pfmvar->deltax) );
        phix[2] = ( stgridO[4].phi - stgridO[3].phi ) / ( (pfmvar->deltax) );
        phix[3] = ( gradx_phi[0] + gradx_phi[1] ) / ( 2.0 );
        phix[4] = ( gradx_phi[0] + gradx_phi[2] ) / ( 2.0 );

        phiy[0] = grady_phi[0];
        phiy[1] = ( grady_phi[0] + grady_phi[1] ) / ( 2.0 );
        phiy[2] = ( grady_phi[0] + grady_phi[2] ) / ( 2.0 );
        phiy[3] = ( stgridO[7].phi - stgridO[4].phi ) / ( (pfmvar->deltay) );
        phiy[4] = ( stgridO[4].phi - stgridO[1].phi ) / ( (pfmvar->deltay) );

        alphidot1[0] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[1] * dphidt[1] ) ) / ( 2.0 );
        alphidot1[1] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[2] * dphidt[2] ) ) / ( 2.0 );
        alphidot1[2] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[3] * dphidt[3] ) ) / ( 2.0 );
        alphidot1[3] = ( ( alpha1[0] * dphidt[0] ) + ( alpha1[4] * dphidt[4] ) ) / ( 2.0 );

        for (jj=0; jj<5; jj++) {
            modgradphi[jj] = sqrt( phix[jj]*phix[jj] + phiy[jj]*phiy[jj] );
        }

        jat1[0] = ( alphidot1[0] * phix[1] ) / ( modgradphi[1] );
        jat1[1] = ( alphidot1[1] * phix[2] ) / ( modgradphi[2] );
        jat1[2] = ( alphidot1[2] * phiy[3] ) / ( modgradphi[3] );
        jat1[3] = ( alphidot1[3] * phiy[4] ) / ( modgradphi[4] );

        for (jj=0; jj<4; jj++) {
            if ( modgradphi[jj+1] == 0.0 ) {
                jat1[jj] = 0.0;
            }
        }

        c1jatx = ( jat1[0] - jat1[1] ) / ( (pfmvar->deltax) );
        c1jaty = ( jat1[2] - jat1[3] ) / ( (pfmvar->deltay) );

        c1jat = c1jatx + c1jaty;

        term1lx = ( ((1.0-hphi[2])+(1.0-hphi[4]))*(stcscl[2].c1l-stcscl[4].c1l) - ((1.0-hphi[4])+(1.0-hphi[1]))*(stcscl[4].c1l-stcscl[1].c1l) )/(2.0*(pfmvar->deltax)*(pfmvar->deltax));
        term1ly = ( ((1.0-hphi[3])+(1.0-hphi[4]))*(stcscl[3].c1l-stcscl[4].c1l) - ((1.0-hphi[4])+(1.0-hphi[0]))*(stcscl[4].c1l-stcscl[0].c1l) )/(2.0*(pfmvar->deltay)*(pfmvar->deltay));

        term1l = term1lx + term1ly;

        c1dot = ((pfmdat->D11l)/(pfmdat->RefD))*(term1l);

        gridNew[( i )*nx + ( j )].c1 = stgridO[4].c1 + (pfmvar->deltat)*( c1dot + c1jat );

    }
}

__kernel void apply_BC(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    
    int index;
    int indexa;
    
    int me;
    
    int istart;
    int iend;
    int totny;
    int nxtotny;
    int totsliceny;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    istart = 0;
    iend = ny;
    
    me = pfmdat->myrank;
    
    index = (i*nx + j);
    
      if ( i == 0 ) {
	gridNew[(i*nx + j)] = gridNew[((i+1)*nx + j)];
      }  
    
      if ( i == (iend-1) ) {
	gridNew[(i*nx + j)] = gridNew[((i-1)*nx + j)];
      }  
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j+1))];
    }
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j-1))];
    }
    
}

__kernel void apply_BC_phi(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    
    int index;
    int indexa;
    
    int me;
    
    int istart;
    int iend;
    int totny;
    int nxtotny;
    int totsliceny;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    istart = 0;
    iend = ny;
    
    me = pfmdat->myrank;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
	gridNew[(i*nx + j)].phi = gridNew[((i+1)*nx + j)].phi;
    }
    
    if ( i == (iend-1) ) {
	gridNew[(i*nx + j)].phi = gridNew[((i-1)*nx + j)].phi;
    }
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j+1))].phi;
    }
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j-1))].phi;
    }
    
}

__kernel void apply_BC_com(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    
    int index;
    int indexa;
    
    int me;
    
    int istart;
    int iend;
    int totny;
    int nxtotny;
    int totsliceny;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    istart = 0;
    iend = ny;
    
    me = pfmdat->myrank;
    
    index = (i*nx + j);
    
      if ( i == 0 ) {
	gridNew[(i*nx + j)].c1 = gridNew[((i+1)*nx + j)].c1;
      }  
    
      if ( i == (iend-1) ) {
	gridNew[(i*nx + j)].c1 = gridNew[((i-1)*nx + j)].c1;
      }  
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j+1))].c1;
    }
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j-1))].c1;
    }
    
}

__kernel void apply_BC_it_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)] = gridNew[((i+1)*nx + j)];
    }
    
}

__kernel void apply_BC_ib_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)] = gridNew[((i-1)*nx + j)];  
    }
    
}

__kernel void apply_BC_it_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)] = gridNew[((ny-2)*nx + j)];
    }
    
}

__kernel void apply_BC_ib_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)] = gridNew[((1)*nx + j)];  
    }
    
}

__kernel void apply_BC_jl_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j+1))];
    }
    
}

__kernel void apply_BC_jr_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (j-1))];
    }
    
}

__kernel void apply_BC_jl_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (nx-2))];
    }
    
}

__kernel void apply_BC_jr_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)] = gridNew[(i*nx + (1))];
    }
    
}

__kernel void apply_BC_phi_it_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
	gridNew[(i*nx + j)].phi = gridNew[((i+1)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_ib_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
	gridNew[(i*nx + j)].phi = gridNew[((i-1)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_it_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
	gridNew[(i*nx + j)].phi = gridNew[((ny-2)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_ib_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
	gridNew[(i*nx + j)].phi = gridNew[((1)*nx + j)].phi;
    }
    
}

__kernel void apply_BC_phi_jl_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j+1))].phi;
    }
    
}

__kernel void apply_BC_phi_jr_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (j-1))].phi;
    }
    
}

__kernel void apply_BC_phi_jl_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (nx-2))].phi;
    }
    
}

__kernel void apply_BC_phi_jr_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].phi = gridNew[(i*nx + (1))].phi;
    }
    
}

__kernel void apply_BC_com_it_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)].c1 = gridNew[((i+1)*nx + j)].c1;
    }
    
}

__kernel void apply_BC_com_ib_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)].c1 = gridNew[((i-1)*nx + j)].c1;
    }
    
}

__kernel void apply_BC_com_it_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == 0 ) {
      gridNew[(i*nx + j)].c1 = gridNew[((ny-2)*nx + j)].c1;
    }
    
}

__kernel void apply_BC_com_ib_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( i == (ny-1) ) {
      gridNew[(i*nx + j)].c1 = gridNew[((1)*nx + j)].c1;
    }
    
}

__kernel void apply_BC_com_jl_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j+1))].c1;
    }
    
}

__kernel void apply_BC_com_jr_noflux(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (j-1))].c1;
    }
    
}

__kernel void apply_BC_com_jl_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == 0 ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (nx-2))].c1;
    }
    
}

__kernel void apply_BC_com_jr_periodic(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {
    int i, j, k;
    int nx;
    int ny;
    int nz;
    int index;
    
    j = get_global_id(0);
    i = get_global_id(1);
    
    nx = pfmdat->jNx;
    ny = pfmdat->iNy;
    
    index = (i*nx + j);
    
    if ( j == (nx-1) ) {
        gridNew[(i*nx + j)].c1 = gridNew[(i*nx + (1))].c1;
    }
    
}

__kernel void addNoise(__global struct grid *gridNew, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);

    double minr = -1.0;
    double maxr = 1.0;
    double rt, noise;

    if ( (gridNew[index].phi > 0.1) && (gridNew[index].phi < 0.5) ) {

        unsigned long x = get_local_id(0);
        unsigned long y = get_local_id(1);
        unsigned long z = get_global_id(0);
        unsigned long w = get_global_id(1);

        ulong t = x ^ (x << 11);
        x = y; 
        y = z; 
        z = w;
        w=index;
        w = w ^ (w >> 51) ^ (t ^ (t >> 8));
        rt = w%2222;
        noise = minr + rt/3333 * (maxr - minr + 1);

        gridNew[index].phi = gridNew[index].phi - ( pfmdat->NoiseFac * noise );

    }

    if (gridNew[index].phi < 0.0) {
        gridNew[index].phi = 0.0;
    }
    if (gridNew[index].phi > 1.0) {
        gridNew[index].phi = 1.0;
    }

}

__kernel void copy_New_To_Old(__global struct grid *gridOld, __global struct grid *gridNew, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    index = (i*nx + j);

    gridOld[index] = gridNew[index];

}

__kernel void update_temp_UC(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    temp[j] = pfmdat->T0;
    

}

__kernel void update_temp_DS(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;
    int i0;
    int j0;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    //i0 = i + (pfmdat->myrank*(ny-2))-1; 
    //temp[i] = pfmdat->Toffset + pfmdat->TG*((i0-pfmdat->PosOffset)*pfmvar->dx - pfmdat->Vp*tstep[0]*pfmvar->dt);

    //i0 = i + (pfmdat->myrank*(ny-2));
    //temp[i] = pfmdat->Toffset + pfmdat->TGRADIENT*((i0-pfmdat->TPosOffset+pfmdat->shift_OFFSET)*pfmvar->deltax-(pfmdat->velocity*tstep[0]*pfmvar->deltat));

    j0 = j + (pfmdat->myrank*(nx-2));
    temp[j] = pfmdat->Toffset + pfmdat->TGRADIENT*((j0-pfmdat->TPosOffset+pfmdat->shift_OFFSET)*pfmvar->deltax-(pfmdat->velocity*tstep[0]*pfmvar->deltat));

}
  
  
__kernel void update_temp_CR(__global double *temp, __global int *tstep, __constant struct pfmval *pfmdat, __constant struct pfmpar *pfmvar) { 
    
    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    //	temp[i] = pfmdat->Tempstart - pfmdat->CR*tstep[0]*pfmvar->dt;

}
  
__kernel void apply_BC_temp_it_noflux(__global double *temp, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;
    int i0;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i == 0 ) { 
        temp[i] = temp[i+1];
    }

}
  
__kernel void apply_BC_temp_ib_noflux(__global double *temp, __constant struct pfmval *pfmdat) {

    int i;
    int j;
    int k;
    int nx;
    int ny;
    int nz;
    int index;
    int i0;

    j = get_global_id(0);
    i = get_global_id(1);

    nx = pfmdat->jNx;
    ny = pfmdat->iNy;

    if ( i == ny-1 ) { 
        temp[i] = temp[i-1];
    }

}
