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
#include "Variables.H"
#include "Initialize.H"

using namespace amrex;
using namespace std;


void readA()
{
    if(ntp != 2)
    {
        Print()<<"This code supports two phases only\n";
    }

    else
    {   
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> Aval;

        for(int i=0; i<ntp; i++)
        {
            int title = 1;
            //Print() << i <<"\n";
            fstream fout;
            fout.open("constant/HSN_"+ tdb_phase[i] +".csv");
            
            if(title==1) 
            {
                getline(fout, line);
                title = 0;
            }
            
            while(!fout.eof())
            {
                getline(fout, line);
                //Print() <<line<<"\n";
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

            //Print() << Aval[77][0]<<"\n";

            A_values.resize(ntp);
            A_values[i].resize(Aval.size()-1);
            A_temp.resize(ntp);
            A_temp[i].resize(Aval.size()-1);

            for(int j=0; j < Aval.size() - 1; j++)
            {
                //Print()<<j<<"\n";
                //Print() << A_values[i].size()<<"\n";
                A_values[i][j] = Aval[j][1];
                //Print()<<A_values[i][j]<<"\n";
                A_temp[i][j] = Aval[j][0];

                //Print()<<A_temp[i][j]<<","<<A_values[i][j]<<"\n";
            }

            Aval.clear();
            
            fout.close();
        
        }
        
    }

}

double findA(int phase, double temp)
{
    A_accel_ptr = gsl_interp_accel_alloc();
    //Print() <<"Size of A_values: " <<A_values[phase].size() <<"\n";
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, A_values[phase].size());

    
    double x_array[A_values[phase].size()];
    double y_array[A_values[phase].size()];

    for (int i=0; i < A_values[phase].size(); i++)
    {
        x_array[i] = A_temp[phase][i];
        y_array[i] = A_values[phase][i];
    }
    
    gsl_spline_init(A_spline_ptr, x_array, y_array, A_values[phase].size());
    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y/2.0;
}


void function_F_04_function_A(){
    readA();
    //Print()<<"In A"<<"\n";
    
    A = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));
    Aeq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                
                A[a][i][j] = findA(a,T);

           } 
        }
    }

    

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                
                Aeq[a][i][j] = findA(a,Teq);

           } 
        }
    }

}

void function_F_04_function_B(){
    B = Vector <Real> (nump,0);
    
   //BB = 2.0*A1[1][0]*cguess[1][2] - 2.0*A1[0][0]*cguess[0][2];

    for(int a=0; a<nump-1; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
    
            B[a] = 2.0*A[nump-1][i][j]*c_guess[nump-1][nump-1][a] - 2.0*A[a][i][j]*c_guess[i][j][a];

           }
        }
    }
    
}

void function_F_04_function_D(){
    D = Vector <Real> (nump,0);

    //DD = A1[0][0]*cguess[0][2]*cguess[0][2] - A1[1][0]*cguess[1][2]*cguess[1][2];

    for(int a=0; a<nump-1; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
    
            D[a] = A[a][i][j]*c_guess[i][j][a]*c_guess[i][j][a] - A[nump-1][i][j]*c_guess[nump-1][nump-1][a]*c_guess[nump-1][nump-1][a];
           
           }
        }
    }
}

void function_F_04_Mu(MultiFab& mu_new){
    
    init_mu(mu_new);
    
}

void function_F_04_c_mu(MultiFab& mu_new){}

void function_F_04_dc_dmu(){
    
     for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                dcdmu.push_back(1.0/(2.0*A[a][i][j]));
           }
        }
     }
    
}


#endif