#ifndef CALCULATE_DIVERGENCE_STRESS_H_
#define CALCULATE_DIVERGENCE_STRESS_H_

void calculate_divergence_stress_2D(long x, struct gradlayer **gradient) {
 int i,j;
 double trace_strain_right;
 double trace_strain_left;
 double trace_strain_back;
 double trace_strain_front;
 
 double sigma_xx_front;
 double sigma_xx_back;
 double sigma_yx_right;
 double sigma_yx_left;
 double sigma_xy_front;
 double sigma_xy_back;
 double sigma_yy_right;
 double sigma_yy_left;
 double forceX;
 double forceY;
 
 double div_phi_front;
 double div_phi_back;
 double div_phi_right;
 double div_phi_left;
  
 double strain_normal[3];
 double eigen_strain_normal[3];
 double stress_normal[3];
 double sigma[3];
 double strain_yy;
 double strain_xx;
//  double stiffness_1(long index, long k);

 for (gidy=1; gidy <= (workers_mpi.end[Y]+2); gidy++) {			//gidy=1  end[Y]+2
  center        =  gidy  + x*(workers_mpi.rows_y);		               //center =  gidy + x*rows_y
  
  grad_front    =  &gradient[2][gidy];
  grad_back     =  &gradient[0][gidy];
  grad          =  &gradient[1][gidy];
  
  grad_right    =  &gradient[1][gidy+1];
  grad_left     =  &gradient[1][gidy-1];
    
//   for(i=0; i<DIMENSION; i++) {
    trace_strain_front      = 0.0;
    trace_strain_back       = 0.0;
    trace_strain_right      = 0.0;
    trace_strain_left       = 0.0;
//   }
  
//   for(j=0; j<DIMENSION; j++) {
//   for(i=0; i<DIMENSION; i++) {
   trace_strain_front = grad_front->strain[X].xx + grad_front->strain[X].yy;
   trace_strain_back  = grad_back->strain[X].xx  + grad_back->strain[X].yy;
   trace_strain_right = grad_right->strain[Y].xx + grad_right->strain[Y].yy;
   trace_strain_left  = grad_left->strain[Y].xx  + grad_left->strain[Y].yy;
//   }

// #ifdef TENSORIAL  
//   if (x <= (end[X]+1)) {
//    div_phi_front = (gridinfo_w[center+2*rows_y].phia[0] -2.0*gridinfo_w[center+rows_y].phia[0] + gridinfo_w[center].phia[0])/(deltax*deltax);
//   } else {
//    div_phi_front = 0.0;
//   }
//   if (x >= 2) {
//     div_phi_back  = (gridinfo1[center-2*rows_y].phia[0]   -2.0*gridinfo1[center-rows_y].phia[0] + gridinfo1[center].phia[0])/(deltax*deltax);
//   } else {
//     div_phi_back  = 0.0;
//   }
//   div_phi_right = (gridinfo1[center+2].phia[0]          -2.0*gridinfo1[center+1].phia[0]      + gridinfo1[center].phia[0])/(deltax*deltax);
//   div_phi_left  = (gridinfo1[center-2].phia[0]          -2.0*gridinfo1[center-1].phia[0]      + gridinfo1[center].phia[0])/(deltax*deltax);
//   
//   
// //   }
//   
// //   sigma_xx_front  = grad_front->stiffness[X][0]*trace_strain_front  + 2.0*grad_front->stiffness[X][1]*grad_front->strain[X][X][X];		//isotropic
// //   sigma_xx_back   = grad_back->stiffness[X][0]*trace_strain_back    + 2.0*grad_back->stiffness[X][1]*grad_back->strain[X][X][X];
// 
//   //Functions for the tensorial interpolation
//   //Convention (XX, YY, XY), (NN, TT, NT)
//    
//   if (fabs(div_phi_front) > 0.0) {
//     strain_yy = grad_front->strain[Y][Y][X] + grad_front->eigen_strain[Y][Y][X];
//     strain_xx = grad_front->strain[X][X][X] + grad_front->eigen_strain[X][X][X];
//     get_normal_strain(gidy,   x+1,  strain_xx,                           strain_yy , grad_front->strain[X][Y][X], strain_normal);
//     get_normal_strain(gidy,   x+1,  eigen_strain_const[0],   eigen_strain_const[0] , 0.0,                   eigen_strain_normal);
//     get_normal_stress(gidy,   x+1,  strain_normal,           eigen_strain_normal, stress_normal);
//     transform_stress(gidy,    x+1,  stress_normal, sigma);
//     sigma_xx_front = sigma[0];
//     sigma_xy_front = sigma[2];
//   } else {
//     sigma_xx_front = grad_front->stiffness[X][0]*trace_strain_front  + 2.0*grad_front->stiffness[X][1]*grad_front->strain[X][X][X]		
//                   + grad_front->stiffness[X][2]*grad_front->strain[X][X][X];
//     sigma_xy_front = 2.0*grad_front->stiffness[Y][1]*grad_front->strain[X][Y][Y];
//   }
//   
//   if (fabs(div_phi_back) > 0.0) {
//     strain_yy = grad_back->strain[Y][Y][X] + grad_back->eigen_strain[Y][Y][X];
//     strain_xx = grad_back->strain[X][X][X] + grad_back->eigen_strain[X][X][X];
//     get_normal_strain(gidy,   x-1,  strain_xx,           strain_yy,  grad_back->strain[X][Y][X],  strain_normal);
//     get_normal_strain(gidy,   x-1,  eigen_strain_const[0],   eigen_strain_const[0] , 0.0,                   eigen_strain_normal);
//     get_normal_stress(gidy,   x-1, strain_normal,  eigen_strain_normal, stress_normal);
//     transform_stress(gidy,   x-1,  stress_normal,  sigma);
//     sigma_xx_back  = sigma[0];
//     sigma_xy_back  = sigma[2];
//   } else {
//     sigma_xx_back  = grad_back->stiffness[X][0]*trace_strain_back    + 2.0*grad_back->stiffness[X][1]*grad_back->strain[X][X][X]
// 		   + grad_back->stiffness[X][2]*grad_back->strain[X][X][X];
//     sigma_xy_back   = 2.0*grad_back->stiffness[Y][1]*grad_back->strain[X][Y][Y];
//   }
//   
//   if (fabs(div_phi_right) > 0.0) {
//     strain_yy = grad_right->strain[Y][Y][Y] + grad_right->eigen_strain[Y][Y][Y];
//     strain_xx = grad_right->strain[X][X][Y] + grad_right->eigen_strain[X][X][Y];
//     get_normal_strain(gidy+1, x  ,  strain_xx,          strain_yy, grad_right->strain[X][Y][Y], strain_normal);
//     get_normal_strain(gidy+1, x  ,  eigen_strain_const[0],   eigen_strain_const[0] , 0.0,                   eigen_strain_normal);
//     get_normal_stress(gidy+1, x  , strain_normal, eigen_strain_normal, stress_normal);
//     transform_stress(gidy+1, x,   stress_normal, sigma);
//     sigma_yy_right = sigma[1];
//     sigma_yx_right = sigma[2];
//   } else {
//     sigma_yx_right  = 2.0*grad_right->stiffness[X][1]*grad_right->strain[Y][X][X];
//     sigma_yy_right  = grad_right->stiffness[Y][0]*trace_strain_right + 2.0*grad_right->stiffness[Y][1]*grad_right->strain[Y][Y][Y];
// 		   + grad_right->stiffness[Y][2]*grad_right->strain[Y][Y][Y];
//   }
//   
//   if (fabs(div_phi_left) > 0.0) {
//     strain_yy = grad_left->strain[Y][Y][Y] + grad_left->eigen_strain[Y][Y][Y];
//     strain_xx = grad_left->strain[X][X][Y] + grad_left->eigen_strain[X][X][Y];
//     get_normal_strain(gidy-1, x  ,  strain_xx,           strain_yy,  grad_left->strain[X][Y][Y],  strain_normal);
//     get_normal_strain(gidy-1, x  ,  eigen_strain_const[0],   eigen_strain_const[0] , 0.0,         eigen_strain_normal);
//     get_normal_stress(gidy-1, x  , strain_normal,  eigen_strain_normal, stress_normal);
//     transform_stress(gidy-1, x,    stress_normal,  sigma);
//     sigma_yx_left  = sigma[2];
//     sigma_yy_left  = sigma[1];
//   } else {
//     sigma_yx_left   = 2.0*grad_left->stiffness[X][1]*grad_left->strain[Y][X][X];
//     sigma_yy_left   = grad_left->stiffness[Y][0]*trace_strain_left   + 2.0*grad_left->stiffness[Y][1]*grad_left->strain[Y][Y][Y];
// 		   + grad_left->stiffness[Y][2]*grad_left->strain[Y][Y][Y];
//   }
// #endif
// KHACHATURIYAN
// #ifdef KHACHATURIYAN

  double lambda_front   = grad_front->stiffness_c[X].C12;  //C12
  double lambda_back    = grad_back->stiffness_c[X].C12;   //C12
  
  double mu_front       = grad_front->stiffness_c[X].C44;  //C44
  double mu_back        = grad_back->stiffness_c[X].C44;   //C44
  
  double mu_prime_front = grad_front->stiffness_c[X].C11 - grad_front->stiffness_c[X].C12 - 2.0*grad_front->stiffness_c[X].C44;
  double mu_prime_back  = grad_back->stiffness_c[X].C11  - grad_back->stiffness_c[X].C12  - 2.0*grad_back->stiffness_c[X].C44;
  
  double lambda_right   = grad_right->stiffness_c[Y].C12;  //C12
  double lambda_left    = grad_left->stiffness_c[Y].C12;   //C12
  
  double mu_right       = grad_right->stiffness_c[Y].C44;  //C44
  double mu_left        = grad_left->stiffness_c[Y].C44;   //C44
  
  double mu_prime_right = grad_right->stiffness_c[Y].C11 - grad_right->stiffness_c[Y].C12 - 2.0*grad_right->stiffness_c[Y].C44;
  double mu_prime_left  = grad_left->stiffness_c[Y].C11  - grad_left->stiffness_c[Y].C12  - 2.0*grad_left->stiffness_c[Y].C44;
  
//   printf("mu_prime_right=%le, mu_prime_left=%le\n", mu_prime_right, mu_prime_left);
  
  sigma_xx_front  = lambda_front*trace_strain_front  + 2.0*mu_front*grad_front->strain[X].xx		//Cubic
		              + mu_prime_front*grad_front->strain[X].xx;
  sigma_xx_back   = lambda_back*trace_strain_back  + 2.0*mu_back*grad_back->strain[X].xx
                  + mu_prime_back*grad_back->strain[X].xx;
  
  //sigma_yx_right  = 2.0*grad_right->stiffness[X][1]*grad_right->strain[Y][X][X];
  //sigma_yx_left   = 2.0*grad_left->stiffness[X][1]*grad_left->strain[Y][X][X];
  
  sigma_yx_right  = 2.0*mu_right*grad_right->strain[X].xy;
  sigma_yx_left   = 2.0*mu_left*grad_left->strain[X].xy;    
//   printf("%le %le\n", sigma_xx_front,sigma_xx_back);
  
  //sigma_xy_front  = 2.0*grad_front->stiffness[Y][1]*grad_front->strain[X][Y][Y];
  //sigma_xy_back   = 2.0*grad_back->stiffness[Y][1]*grad_back->strain[X][Y][Y];
  
  sigma_xy_front  = 2.0*mu_front*grad_front->strain[Y].xy;
  sigma_xy_back   = 2.0*mu_back*grad_back->strain[Y].xy;
  
//   printf("sigma_xx_front=%le, sigma_yx_right=%le,  sigma_xy_front=%le\n", sigma_xx_front, sigma_yx_right, sigma_xy_front);
     
//   sigma_yy_right  = grad_right->stiffness[Y][0]*trace_strain_right + 2.0*grad_right->stiffness[Y][1]*grad_right->strain[Y][Y][Y]
// 		              + grad_right->stiffness[Y][2]*grad_right->strain[Y][Y][Y];
//   sigma_yy_left   = grad_left->stiffness[Y][0]*trace_strain_left   + 2.0*grad_left->stiffness[Y][1]*grad_left->strain[Y][Y][Y]
// 		              + grad_left->stiffness[Y][2]*grad_left->strain[Y][Y][Y];
  
  sigma_yy_right  = lambda_right*trace_strain_right + 2.0*mu_right*grad_right->strain[Y].yy
		              + mu_prime_right*grad_right->strain[Y].yy;
                  
  sigma_yy_left   = lambda_left*trace_strain_left   + 2.0*mu_left*grad_left->strain[Y].yy
		              + mu_prime_left*grad_left->strain[Y].yy;
// #endif
             
//----------------------------------------------------------------------------------------------------------------------------------------------                  
      /*      sigma_xx_front  = grad_front->stiffness[X][0]*trace_strain_front  + 2.0*grad_front->stiffness[X][1]*grad_front->strain[X][X][X]		//Cubic
		   + grad_front->stiffness[X][2]*grad_front->strain[X][X][X];
  sigma_xx_back   = grad_back->stiffness[X][0]*trace_strain_back    + 2.0*grad_back->stiffness[X][1]*grad_back->strain[X][X][X]
		   + grad_back->stiffness[X][2]*grad_back->strain[X][X][X];
  
  sigma_yx_right  = 2.0*grad_right->stiffness[X][1]*grad_right->strain[Y][X][X];
  sigma_yx_left   = 2.0*grad_left->stiffness[X][1]*grad_left->strain[Y][X][X];
  
//   printf("%le %le\n", sigma_xx_front,sigma_xx_back);
  
  sigma_xy_front  = 2.0*grad_front->stiffness[Y][1]*grad_front->strain[X][Y][Y];
  sigma_xy_back   = 2.0*grad_back->stiffness[Y][1]*grad_back->strain[X][Y][Y];
     
  sigma_yy_right  = grad_right->stiffness[Y][0]*trace_strain_right + 2.0*grad_right->stiffness[Y][1]*grad_right->strain[Y][Y][Y]
		   + grad_right->stiffness[Y][2]*grad_right->strain[Y][Y][Y];
  sigma_yy_left   = grad_left->stiffness[Y][0]*trace_strain_left   + 2.0*grad_left->stiffness[Y][1]*grad_left->strain[Y][Y][Y]
		   + grad_left->stiffness[Y][2]*grad_left->strain[Y][Y][Y];      */ 
                  
//-------------------------------------------------------------------------------------------------------------------------------------------------                  

//   forceX          = (sigma_xx_front-sigma_xx_back)/deltax           +  (sigma_yx_right-sigma_yx_left)/deltay;
  forceX          = (sigma_xx_front-sigma_xx_back)                 +  (sigma_yx_right-sigma_yx_left);               
  forceX         *= 0.5;
       
//   forceY          = (sigma_xy_front - sigma_xy_back)/deltax        + (sigma_yy_right-sigma_yy_left)/deltay;
  forceY          = (sigma_xy_front - sigma_xy_back)               + (sigma_yy_right-sigma_yy_left);
  forceY         *= 0.5;
//   printf("%le %le\n",forceX,forceY);
  
  iter_gridinfo_w[center].disp[X][0] = iter_gridinfo_w[center].disp[X][1];
  iter_gridinfo_w[center].disp[Y][0] = iter_gridinfo_w[center].disp[Y][1];
  iter_gridinfo_w[center].disp[X][1] = iter_gridinfo_w[center].disp[X][2];
  iter_gridinfo_w[center].disp[Y][1] = iter_gridinfo_w[center].disp[Y][2];
  
//   iter_gridinfo1[center].disp[X][2] = (((deltat*deltat)/rho)*forceX - iter_gridinfo1[center].disp[X][0] + (2+damping_factor*deltat)*iter_gridinfo1[center].disp[X][1])/(1.0 + damping_factor*deltat);
//   
//   iter_gridinfo1[center].disp[Y][2] = (((deltat*deltat)/rho)*forceY - iter_gridinfo1[center].disp[Y][0] + (2+damping_factor*deltat)*iter_gridinfo1[center].disp[Y][1])/(1.0 + damping_factor*deltat); 
  
  
  iter_gridinfo_w[center].disp[X][2] = (((deltat_e*deltat_e)/rho)*forceX - (1 - damping_factor*deltat_e)*iter_gridinfo_w[center].disp[X][0] + 2*iter_gridinfo_w[center].disp[X][1])/(1.0 + damping_factor*deltat_e);

  iter_gridinfo_w[center].disp[Y][2] = (((deltat_e*deltat_e)/rho)*forceY - (1 - damping_factor*deltat_e)*iter_gridinfo_w[center].disp[Y][0] + 2*iter_gridinfo_w[center].disp[Y][1])/(1.0 + damping_factor*deltat_e); 
   
//     iter_gridinfo_w[center].disp[X][2] = (((0.001*0.001)/rho)*forceX - (1 - damping_factor*0.001)*iter_gridinfo_w[center].disp[X][0] + 2*iter_gridinfo_w[center].disp[X][1])/(1.0 + damping_factor*0.001);
//   
//    iter_gridinfo_w[center].disp[Y][2] = (((0.001*0.001)/rho)*forceY - (1 - damping_factor*0.001)*iter_gridinfo_w[center].disp[Y][0] + 2*iter_gridinfo_w[center].disp[Y][1])/(1.0 + damping_factor*0.001); 
  
//    	printf("%le\n",(iter_gridinfo1[center].disp[X][2]-iter_gridinfo1[center].disp[X][1]));
//       	printf("%le\n",(iter_gridinfo1[center].disp[X][2]));
   
//   iter_gridinfo_w[center].disp[X][0] = iter_gridinfo_w[center].disp[X][1];
//   iter_gridinfo_w[center].disp[Y][0] = iter_gridinfo_w[center].disp[Y][1];
//   iter_gridinfo_w[center].disp[X][1] = iter_gridinfo_w[center].disp[X][2];
//   iter_gridinfo_w[center].disp[Y][1] = iter_gridinfo_w[center].disp[Y][2];

  }
 
}
void calculate_divergence_stress_3D(long x, struct gradlayer **gradient) {
 int i,j;
 double trace_strain_right;
 double trace_strain_left;
 double trace_strain_back;
 double trace_strain_front;
 double trace_strain_top;
 double trace_strain_bottom;
 
 struct symmetric_tensor sigma_front;
 struct symmetric_tensor sigma_back;
 struct symmetric_tensor sigma_right;
 struct symmetric_tensor sigma_left;
 struct symmetric_tensor sigma_top;
 struct symmetric_tensor sigma_bottom;
 
 double forceX;
 double forceY;
 double forceZ;
 
 double div_phi_front;
 double div_phi_back;
 double div_phi_right;
 double div_phi_left;
 double div_phi_top;
 double div_phi_bottom;
  
 for (gidy=workers_mpi.rows_y+1; gidy <= (workers_mpi.layer_size - workers_mpi.rows_y-2); gidy++) {			//gidy=1  end[Y]+2
  center        =  gidy  + x*(workers_mpi.layer_size);		               //center =  gidy + x*rows_y
  
  grad_front    =  &gradient[2][gidy];
  grad_back     =  &gradient[0][gidy];
  grad          =  &gradient[1][gidy];
  
  grad_right    =  &gradient[1][gidy+1];
  grad_left     =  &gradient[1][gidy-1];
  grad_top      =  &gradient[1][gidy+workers_mpi.rows_y];
  grad_bottom   =  &gradient[1][gidy-workers_mpi.rows_y];
  
//   trace_strain_front      = 0.0;
//   trace_strain_back       = 0.0;
//   trace_strain_right      = 0.0;
//   trace_strain_left       = 0.0;
//   trace_strain_top        = 0.0;
//   trace_strain_bottom     = 0.0;

  trace_strain_front     = grad_front->strain[X].xx   + grad_front->strain[X].yy   + grad_front->strain[X].zz;
  trace_strain_back      = grad_back->strain[X].xx    + grad_back->strain[X].yy    + grad_back->strain[X].zz;
  trace_strain_right     = grad_right->strain[Y].xx   + grad_right->strain[Y].yy   + grad_right->strain[Y].zz;
  trace_strain_left      = grad_left->strain[Y].xx    + grad_left->strain[Y].yy    + grad_left->strain[Y].zz;
  trace_strain_top       = grad_top->strain[Z].xx     + grad_top->strain[Z].yy     + grad_top->strain[Z].zz;
  trace_strain_bottom    = grad_bottom->strain[Z].xx  + grad_bottom->strain[Z].yy  + grad_bottom->strain[Z].zz;

  double lambda_front    = grad_front->stiffness_c[X].C12;  //C12
  double lambda_back     = grad_back->stiffness_c[X].C12;   //C12
  
  double mu_front        = grad_front->stiffness_c[X].C44;  //C44
  double mu_back         = grad_back->stiffness_c[X].C44;   //C44
  
  double mu_prime_front  = grad_front->stiffness_c[X].C11 - grad_front->stiffness_c[X].C12 - 2.0*grad_front->stiffness_c[X].C44;
  double mu_prime_back   = grad_back->stiffness_c[X].C11  - grad_back->stiffness_c[X].C12  - 2.0*grad_back->stiffness_c[X].C44;
  
  double lambda_right    = grad_right->stiffness_c[Y].C12;  //C12
  double lambda_left     = grad_left->stiffness_c[Y].C12;   //C12
  
  double mu_right        = grad_right->stiffness_c[Y].C44;  //C44
  double mu_left         = grad_left->stiffness_c[Y].C44;   //C44
  
  double mu_prime_right  = grad_right->stiffness_c[Y].C11 - grad_right->stiffness_c[Y].C12 - 2.0*grad_right->stiffness_c[Y].C44;
  double mu_prime_left   = grad_left->stiffness_c[Y].C11  - grad_left->stiffness_c[Y].C12  - 2.0*grad_left->stiffness_c[Y].C44;
  
  double lambda_top      = grad_top->stiffness_c[Z].C12;  //C12
  double lambda_bottom   = grad_bottom->stiffness_c[Z].C12;   //C12
  
  double mu_top          = grad_top->stiffness_c[Z].C44;  //C44
  double mu_bottom       = grad_bottom->stiffness_c[Z].C44;   //C44
  
  double mu_prime_top    = grad_top->stiffness_c[Z].C11    - grad_top->stiffness_c[Z].C12    - 2.0*grad_top->stiffness_c[Z].C44;
  double mu_prime_bottom = grad_bottom->stiffness_c[Z].C11 - grad_bottom->stiffness_c[Z].C12 - 2.0*grad_bottom->stiffness_c[Z].C44;
  
  sigma_front.xx  = lambda_front*trace_strain_front  + 2.0*mu_front*grad_front->strain[X].xx		//Cubic
		              + mu_prime_front*grad_front->strain[X].xx;
  sigma_back.xx   = lambda_back*trace_strain_back  + 2.0*mu_back*grad_back->strain[X].xx
                  + mu_prime_back*grad_back->strain[X].xx;
    
  sigma_right.xy  = 2.0*mu_right*grad_right->strain[X].xy;
  sigma_left.xy   = 2.0*mu_left*grad_left->strain[X].xy;    
  
  sigma_right.yz  = 2.0*mu_right*grad_right->strain[Z].yz;
  sigma_left.yz   = 2.0*mu_left*grad_left->strain[Z].yz; 
  
  sigma_top.xz    = 2.0*mu_top*grad_top->strain[X].xz;
  sigma_bottom.xz = 2.0*mu_bottom*grad_bottom->strain[X].xz;
  
  sigma_top.yz    = 2.0*mu_top*grad_top->strain[Y].yz;
  sigma_bottom.yz = 2.0*mu_bottom*grad_bottom->strain[Y].yz;
  
  sigma_front.xy  = 2.0*mu_front*grad_front->strain[Y].xy;
  sigma_back.xy   = 2.0*mu_back*grad_back->strain[Y].xy;
   
  sigma_front.xz  = 2.0*mu_front*grad_front->strain[Z].xz;
  sigma_back.xz   = 2.0*mu_back*grad_back->strain[Z].xz;
  
  sigma_right.yy  = lambda_right*trace_strain_right   + 2.0*mu_right*grad_right->strain[Y].yy
		              + mu_prime_right*grad_right->strain[Y].yy;
                  
  sigma_left.yy   = lambda_left*trace_strain_left     + 2.0*mu_left*grad_left->strain[Y].yy
		              + mu_prime_left*grad_left->strain[Y].yy;
                  
  sigma_top.zz    = lambda_top*trace_strain_top       + 2.0*mu_top*grad_top->strain[Z].zz
		              + mu_prime_top*grad_top->strain[Z].zz;
                  
  sigma_bottom.zz = lambda_bottom*trace_strain_bottom + 2.0*mu_bottom*grad_bottom->strain[Z].zz
		              + mu_prime_bottom*grad_bottom->strain[Z].zz;                

  forceX          = (sigma_front.xx - sigma_back.xx)  + (sigma_right.xy - sigma_left.xy) + (sigma_top.xz   - sigma_bottom.xz);
  forceX         *= 0.5;
       
  forceY          = (sigma_front.xy - sigma_back.xy)  + (sigma_top.yz - sigma_bottom.yz) + (sigma_right.yy - sigma_left.yy);
  forceY         *= 0.5;
  
  forceZ          = (sigma_front.xz - sigma_back.xz)  + (sigma_right.yz - sigma_left.yz) + (sigma_top.zz   - sigma_bottom.zz);
  forceZ         *= 0.5;
  
  
  iter_gridinfo_w[center].disp[X][0] = iter_gridinfo_w[center].disp[X][1];
  iter_gridinfo_w[center].disp[Y][0] = iter_gridinfo_w[center].disp[Y][1];
  iter_gridinfo_w[center].disp[Z][0] = iter_gridinfo_w[center].disp[Z][1];
  iter_gridinfo_w[center].disp[X][1] = iter_gridinfo_w[center].disp[X][2];
  iter_gridinfo_w[center].disp[Y][1] = iter_gridinfo_w[center].disp[Y][2];
  iter_gridinfo_w[center].disp[Z][1] = iter_gridinfo_w[center].disp[Z][2];
  
//   iter_gridinfo1[center].disp[X][2] = (((deltat*deltat)/rho)*forceX - iter_gridinfo1[center].disp[X][0] + (2+damping_factor*deltat)*iter_gridinfo1[center].disp[X][1])/(1.0 + damping_factor*deltat);
//   
//   iter_gridinfo1[center].disp[Y][2] = (((deltat*deltat)/rho)*forceY - iter_gridinfo1[center].disp[Y][0] + (2+damping_factor*deltat)*iter_gridinfo1[center].disp[Y][1])/(1.0 + damping_factor*deltat); 
  
//    double deltat_e = 0.1;
  iter_gridinfo_w[center].disp[X][2] = (((deltat_e*deltat_e)/rho)*forceX - (1 - damping_factor*deltat_e)*iter_gridinfo_w[center].disp[X][0] + 2*iter_gridinfo_w[center].disp[X][1])/(1.0 + damping_factor*deltat_e);

  iter_gridinfo_w[center].disp[Y][2] = (((deltat_e*deltat_e)/rho)*forceY - (1 - damping_factor*deltat_e)*iter_gridinfo_w[center].disp[Y][0] + 2*iter_gridinfo_w[center].disp[Y][1])/(1.0 + damping_factor*deltat_e); 
  
  iter_gridinfo_w[center].disp[Z][2] = (((deltat_e*deltat_e)/rho)*forceZ - (1 - damping_factor*deltat_e)*iter_gridinfo_w[center].disp[Z][0] + 2*iter_gridinfo_w[center].disp[Z][1])/(1.0 + damping_factor*deltat_e); 
  
//    	printf("%le\n",(iter_gridinfo1[center].disp[X][2]-iter_gridinfo1[center].disp[X][1]));
//       	printf("%le\n",(iter_gridinfo1[center].disp[X][2]));

  }
 
}

//computing error
void compute_error_2D(long x, double *error) {
 for (gidy=workers_mpi.start[Y]; gidy <= (workers_mpi.end[Y]); gidy++) {			//gidy=1  end[Y]+2
   center        = gidy  + x*(workers_mpi.rows_y);
   (*error) += ((iter_gridinfo_w[center].disp[X][2]-iter_gridinfo_w[center].disp[X][1])*(iter_gridinfo_w[center].disp[X][2]-iter_gridinfo_w[center].disp[X][1]) 
	      + (iter_gridinfo_w[center].disp[Y][2]-iter_gridinfo_w[center].disp[Y][1])*(iter_gridinfo_w[center].disp[Y][2]-iter_gridinfo_w[center].disp[Y][1]));
  }
}

void compute_error_3D(long x, double *error) {
 for (gidy=3; gidy <= (workers_mpi.end[Y]); gidy++) {			//gidy=1  end[Y]+2
   center        = gidy  + x*(workers_mpi.rows_y);
   (*error) += ((iter_gridinfo_w[center].disp[X][2]-iter_gridinfo_w[center].disp[X][1])*(iter_gridinfo_w[center].disp[X][2]-iter_gridinfo_w[center].disp[X][1]) 
	      + (iter_gridinfo_w[center].disp[Y][2]-iter_gridinfo_w[center].disp[Y][1])*(iter_gridinfo_w[center].disp[Y][2]-iter_gridinfo_w[center].disp[Y][1])
        + (iter_gridinfo_w[center].disp[Z][2]-iter_gridinfo_w[center].disp[Z][1])*(iter_gridinfo_w[center].disp[Z][2]-iter_gridinfo_w[center].disp[Z][1])     
        );
  }
}
#endif

