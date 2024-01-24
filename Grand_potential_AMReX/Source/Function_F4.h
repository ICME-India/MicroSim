#ifndef _FUNCTION_F4_H_
#define _FUNCTION_F4_H_

#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "Variables.h"
#include "Filling_fv.h"
#include "Chkpnt.h"
#include "GP_Utility.h"

using namespace amrex;
using namespace std;



void readC()
{       
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> cval;

        for(int a=0; a<nump-1;a++){
            int title = 1;
            fstream fout;
            //fout.open("constant/Composition_"+ phasemap[a] +".csv");
            fout.open("tdbs_encrypted/Composition_"+ phasemap[a] +".csv");
            
            if(title==1) 
            {
                getline(fout, line);
                title = 0;
            }
            
            while(!fout.eof())
            {
                getline(fout, line);
                stringstream s(line);

                while(getline(s, value, ','))
                {
                    data.push_back(stod(value));
                }

                s.str("");
                s.clear();
                cval.push_back(data);
                data.clear();
            }

            
            conc_Sol.resize(nump-1);
            for (int i = 0; i < nump-1; i++)
            {
                conc_Sol[i].resize(cval.size()-1);
                for (int j = 0; j < cval.size()-1; j++)
                {
                conc_Sol[i][j].resize(numcom-1);
                }
            }

            conc_Liq.resize(nump-1);
            for (int i = 0; i < nump-1; i++)
            {
                conc_Liq[i].resize(cval.size()-1);
                for (int j = 0; j < cval.size()-1; j++)
                {
                conc_Liq[i][j].resize(numcom-1);
                }
            }


            temprt.resize(nump-1);
            temprt[a].resize(cval.size()-1);
            
         
            for(int r=0; r<numcom-1;r++){
                for(int b=0; b < cval.size() - 1; b++)
                {   
                   
                    conc_Sol[a][b][r] = cval[b][1+r];
                    
                    conc_Liq[a][b][r] = cval[b][numcom+r];
                   
                    if(r==0){
                    temprt[a][b] = cval[b][0];
                    }
                }
                
            }
         

            cval.clear();
            
            fout.close();

        }

}

double findC(int phase, double temp, int i)
{   
    if(phase!=(nump-1)){
    
    A_accel_ptr = gsl_interp_accel_alloc();
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Sol[phase].size());

    
    double x_array[conc_Sol[phase].size()];
    double y_array[conc_Sol[phase].size()];

    for (int m=0; m < conc_Sol[phase].size(); m++)
    {
        x_array[m] = temprt[phase][m];
        y_array[m] = conc_Sol[phase][m][i];

    }
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Sol[phase].size());
    }
    
    else if(phase==(nump-1)){
    A_accel_ptr = gsl_interp_accel_alloc();
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Liq[0].size());

    
    double x_array[conc_Liq[0].size()];
    double y_array[conc_Liq[0].size()];

    for (int m=0; m < conc_Liq[0].size(); m++)
    {
        x_array[m] = temprt[0][m];
        y_array[m] = conc_Liq[0][m][i];

    }
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Liq[0].size());
    }

    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y;

}

void getc(){
    
    readC();

    conc = Vector<Vector<Real>>(nump,Vector<Real>(numcom-1,0.0));
    conceq = Vector<Vector<Real>>(nump,Vector<Real>(numcom-1,0.0));

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
                
            conc[a][i] = findC(a,T,i);
            conceq[a][i] = findC(a,Teq,i);

        } 
    }

}

void readA()
{       
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> Aval;

        for(int a=0; a<nump; a++)
        {   
            int title = 1;
            fstream fout;
            //fout.open("constant/HSN_"+ phasemap[a] +".csv");
            fout.open("tdbs_encrypted/HSN_"+ phasemap[a] +".csv");


            if(title==1) 
            {
                getline(fout, line);
                title = 0;
            }
            
            while(!fout.eof())
            {
                getline(fout, line);
                stringstream s(line);

                while(getline(s, value, ','))
                {
                    data.push_back(stod(value));
                }

                s.str("");
                s.clear();
                Aval.push_back(data);
                data.clear();
            }


            A_values.resize(nump);
            for (int i = 0; i < nump; i++)
            {
                A_values[i].resize(Aval.size()-1);
                for (int j = 0; j < Aval.size()-1; j++)
                {
                A_values[i][j].resize(Aval[0].size()-1);
                }
            }

            A_temp.resize(nump);
            A_temp[a].resize(Aval.size()-1);


            for(int i=0; i < Aval.size() - 1; i++)
            {   for(int j=0; j<Aval[0].size(); j++){
                    if(j==0){
                        A_temp[a][i] = Aval[i][j];
                    }
                    else{
                        A_values[a][i][j-1] = Aval[i][j];
                    }
                }
            }

            Aval.clear();
            
            fout.close();
        
        }

}

double findA(int phase, double temp, int o, int r)
 {
     A_accel_ptr = gsl_interp_accel_alloc();
    
     A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, A_values[phase].size());

    double x_array[A_values[phase].size()];
    double y_array[A_values[phase].size()];

    int u=0;
    Vector<Vector<Real>>L(numcom-1,Vector<Real>(numcom-1,0.0));
    
    for(int p=0; p<numcom-1; p++){
        for(int m=p+1; m<numcom-1; m++){
            L[p][m] = numcom-1+u;
            L[m][p] = numcom-1+u;
            u++;
        }
    }

    if(o==r){
        for (int m=0; m < A_values[phase].size(); m++)
        {
            x_array[m] = A_temp[phase][m];
            y_array[m] = A_values[phase][m][o];
        }
    }
    else{
        for (int m=0; m < A_values[phase].size(); m++)
        {
            x_array[m] = A_temp[phase][m];
            y_array[m] = A_values[phase][m][L[o][r]];
        }
    }

    gsl_spline_init(A_spline_ptr, x_array, y_array, A_values[phase].size());
    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y;


 }

void function_F_04_function_A(){

	BL_PROFILE("function_F_04_function_A()");

    readA();
    
    A = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    Aeq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));

    for (int a=0; a < nump; a++) {
        for(int i=0; i<numcom-1; i++) {
            for(int j=0; j<numcom-1; j++) {
                if (i==j) {
                    A[a][i][j] = 0.5*findA(a,T,i,j);
                    Aeq[a][i][j] = 0.5*findA(a,Teq,i,j);
                } else {
                    A[a][i][j] = findA(a,T,i,j);
                    Aeq[a][i][j] = findA(a,Teq,i,j);
                }
            }
        }
    }

}

void function_F_04_function_B(){
    
    getc();
    
    B = Vector<Vector<Real>>(nump,Vector<Real>(numcom-1,0.0));

    double sum_c = 0.0;
    for(int a=0; a< nump; a++){
        if (a != (nump-1)) {
        for(int i = 0; i< numcom-1; i++){    
            for (int k=0; k < numcom-1; k++) {
                if (k!=i) {
                    sum_c += A[nump-1][k][i]*conc[nump-1][k] - A[a][k][i]*conc[a][k];
            }
            }
            B[a][i] = (2.0*(A[nump-1][i][i]*conc[nump-1][i] - A[a][i][i]*conc[a][i]) + sum_c);
            sum_c=0.0;
        }
        }
    }

}

void function_F_04_function_C(){

    C=Vector<Real>(nump,0.0);

    for(int a =0; a<nump;a++){
        if (a != (nump-1)) {
        for(int i=0; i<numcom-1; i++){
            for(int j=0; j<numcom-1; j++){
                if (i <= j) {
                    C[a] += (A[a][i][j]*conc[a][i]*conc[a][j] - A[nump-1][i][j]*conc[nump-1][i]*conc[nump-1][j]);
                }
            }
        }
        }
    }
    
}

void function_F_04_Mu(MultiFab& mu_new){

    init_mu(mu_new);

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void c_mu(int i, int j, int k, amrex::Array4<Real const> const& mu, Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> &c, Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> BB, Array3D<Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA, int numcomp, int a){

            Array2D <Real,0,compcount-2,0,compcount-2> muc{};
            Array2D <Real,0,compcount-2,0,compcount-2> cmu{};

            
                for (int l=0;l<numcomp-1;l++) {
                    for (int m=0;m<numcomp-1;m++) {
                        if (l==m) {
                            muc(l,m)=2.0*AA(a,l,m);
                        
                        } else {
                            muc(l,m)=AA(a,l,m);
                            
                    }
                }
                }

                if(numcomp==3){
                    double det = muc(0,0)*muc(1,1) - muc(0,1)*muc(1,0);
                    cmu(0,0) = muc(1,1)/det;
                    cmu(0,1) = -muc(0,1)/det;
                    cmu(1,0) = -muc(1,0)/det;
                    cmu(1,1) = muc(0,0)/det;

                }

                if(numcomp==2){
                    cmu(0,0) = 1/muc(0,0);
                }
                

                double sum = 0.0;
                
                for(int l=0; l < numcomp-1; l++) {
                    
                    for (int m=0; m < numcomp-1; m++) {
                        sum += cmu(l,m)*(mu(i,j,k,m)-BB(a,m)); 
                    }
                    c(a,l) = sum;
                    sum=0.0;
                }
}


void function_F_04_dc_dmu(){
    
    dcdmu = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    dcdmu_eq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));

     for(int a=0; a<nump; a++){

        Array2D <double,0,compcount-2,0,compcount-2> muc{};
        Array2D <double,0,compcount-2,0,compcount-2> muc_eq{};

            
                for (int l=0;l<numcom-1;l++) {
                    for (int m=0;m<numcom-1;m++) {
                        if (l==m) {
                            muc(l,m)=2.0*A[a][l][m];
                            muc_eq(l,m)=2.0*Aeq[a][l][m];
                        
                        } else {
                            muc(l,m)=A[a][l][m];
                            muc_eq(l,m)=Aeq[a][l][m];         
                    }
                }
                }

                if(numcom==3){
                    double det = muc(0,0)*muc(1,1) - muc(0,1)*muc(1,0);
                    dcdmu[a][0][0] = muc(1,1)/det;
                    dcdmu[a][0][1] = -muc(0,1)/det;
                    dcdmu[a][1][0] = -muc(1,0)/det;
                    dcdmu[a][1][1] = muc(0,0)/det;

                    double det_eq = muc_eq(0,0)*muc_eq(1,1) - muc_eq(0,1)*muc_eq(1,0);
                    dcdmu_eq[a][0][0] = muc_eq(1,1)/det;
                    dcdmu_eq[a][0][1] = -muc_eq(0,1)/det;
                    dcdmu_eq[a][1][0] = -muc_eq(1,0)/det;
                    dcdmu_eq[a][1][1] = muc_eq(0,0)/det;

                }

                if(numcom==2){
                    dcdmu[a][0][0] = 1/muc(0,0);
                    dcdmu_eq[a][0][0] = 1/muc_eq(0,0);
                }
     }   
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void free_energy(int i, int j, int k, int numphase, Array3D<Real,0,phasecount-1,0,compcount-2,0,compcount-2, Order::C> AA, Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> BB, Array1D <Real,0,phasecount-1> CC, Array1D<Real,0,phasecount-1> &fe ,Array2D<Real,0,phasecount-1,0,compcount-2, Order::C> &c, int numcomp, int a){

    for(int a=0; a<numphase; a++){
    double sum=0.0;
        for (int l=0;l<numcomp-1;l++) {
            for (int m=0;m<numcomp-1;m++) {
                if (l<=m) {
                    sum += AA(a,l,m)*c(a,l)*c(a,m);
                }
            }
            sum += BB(a,l)*c(a,l);
        }
    sum += CC(a);
    fe(a) = sum;
    }
}


#endif
