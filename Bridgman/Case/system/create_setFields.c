//Create random filling of precipitates in the setFields
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "ran2.c"
int main(int argc, char *argv[]) {
 FILE *fp;
 static long A = 1000;
 long limit_x, limit_y;
 long xstart, xend;
 long ystart, yend;
 long zstart, zend;
 
 long rand_x, rand_y;
 long rand_x_, rand_y_;
 
 long nppt=0;
 
 long x,y;
 
 double rcrit;
 
 double volfrac;
 
 double phi_value;
 double mu_1_value;
 double mu_2_value;
 
 char name[10000];
 
 fp = fopen("Infile.dat","r");
 
 fscanf(fp,"%s %ld",  name, &xstart);
 fscanf(fp,"%s %ld",  name, &xend);
 fscanf(fp,"%s %ld",  name, &ystart);
 fscanf(fp,"%s %ld",  name, &yend);
 fscanf(fp,"%s %ld",  name, &zstart);
 fscanf(fp,"%s %ld",  name, &zend);
 fscanf(fp,"%s %le",  name, &volfrac);
 fscanf(fp,"%s %le",  name, &rcrit);
 fscanf(fp,"%s %le",  name, &phi_value);
 fscanf(fp,"%s %le",  name, &mu_1_value);
 fscanf(fp,"%s %le",  name, &mu_2_value);
 
 fclose(fp);
 
 
 
//  xstart = atol(argv[1]);
//  xend   = atol(argv[2]);
//  ystart = atol(argv[3]);
//  yend   = atol(argv[4]);
//  zstart = atol(argv[5]);
//  zend   = atol(argv[6]);
//  
//  volfrac = atof(argv[7]);
//  rcrit   = atof(argv[8]);
 
 limit_x  =  xend - xstart;
 limit_y  =  yend - ystart;
 
 int FLAG[limit_x][limit_y];
 int FILLED=0;
 
 nppt = ceil(((limit_x-2.0*rcrit)*(limit_y-2.0*rcrit)*volfrac)/(M_PI*rcrit*rcrit));

 long xpoints[nppt];
 long ypoints[nppt];
 
 long ncount=0;
 long iter=0;
 long i;
 double vol_ppt=0.0;
 
 for(x=xstart; x<=xend; x++) {
   for(x=xstart; x<=xend; x++) {
     FLAG[x][y]=0;
   }
 }
 
 
 while((ncount<nppt) && (iter<1000*nppt)) {
  
  long num = rand()%20000;
  rand_x = (xstart+rcrit) + ran2(&num)*(limit_x-2.0*rcrit);
  rand_y = (ystart+rcrit) + ran2(&num)*(limit_y-2.0*rcrit);
  
  if ((rand_x!=rand_x_) && (rand_y!=rand_y_)) {
    xpoints[ncount] = rand_x;
    ypoints[ncount] = rand_y;
    
    ncount++;
    vol_ppt += M_PI*rcrit*rcrit;
  }
  rand_x_ = rand_x;  
  rand_y_ = rand_y;
  iter++;
  printf("ncount=%ld, rand_x=%ld, rand_y=%ld\n", ncount, rand_x, rand_y);
 }
  
//   for(x=(rand_x-rcrit); x<=(rand_x+rcrit); x++) {
//      for(y=(rand_y-rcrit); y<=(rand_y+rcrit); y++) {
//        
//      }
//   }
  
  
//   for(x=(rand_x-rcrit); x<=(rand_x+rcrit); x++) {
//     for(y=(rand_y-rcrit); y<=(rand_y+rcrit); y++) {
//       if (FLAG[x][y] == 1) {
//         FILLED=1;
//         break;
//       }
//     }
//     if(FILLED==1) {
//       break;
//     }
//   }
//   
//   if(FILLED !=1) {
//     for(x=(rand_x-rcrit); x<=(rand_x+rcrit); x++) {
//      for(y=(rand_y-rcrit); y<=(rand_y+rcrit); y++) {
//        if ((x-rand_x)*(x-rand_x) + (y-rand_y)*(y-rand_y) <= rcrit*rcrit) {
//          FLAG[x][y]=1;
//        }
//      }
//     }
//     xpoints[ncount] = rand_x;
//     ypoints[ncount] = rand_y;
//     ncount++;
//     FILLED=0;
//     vol_ppt += M_PI*rcrit*rcrit;
//   }
//   iter++;
//   printf("ncount=%ld\n", ncount);
//  }
 printf("Volume-fraction of precipitates, vol_ppt, domain_volume, nppt = %le %le %le %ld\n",vol_ppt/((double)(limit_x-2.0*rcrit)*(limit_y-2.0*rcrit)), vol_ppt, (limit_x-2.0*rcrit)*(limit_y-2.0*rcrit), nppt);
 

  fp = fopen("random_filling_setfields.dat","w");

  fprintf(fp,"regions (\n");

  for(i=0; i<ncount; i++){
    fprintf(fp,"cylinderToCell {\n");
    fprintf(fp,"p1 (%ld %ld %ld);\n", xpoints[i], ypoints[i], zstart);
    fprintf(fp,"p2 (%ld %ld %ld);\n", xpoints[i], ypoints[i], zend);
    fprintf(fp,"radius %le;\n",rcrit);
    fprintf(fp,"fieldValues (\n");
    fprintf(fp,"volScalarFieldValue phi   %le\n", phi_value);
    fprintf(fp,"volScalarFieldValue mu_1  %le\n", mu_1_value);
    fprintf(fp,"volScalarFieldValue mu_2  %le\n", mu_2_value);
    fprintf(fp,");\n");
    fprintf(fp,"}\n");
  }
  fprintf(fp,");");
  
  fclose(fp);
 
}