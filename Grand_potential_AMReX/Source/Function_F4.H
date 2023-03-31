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
#include "Chkpnt.H"

using namespace amrex;
using namespace std;


void readc(){
    if(ntp != 2)
    {
        Print()<<"This code supports two phases only\n";
    }

    else
    {   
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> cval;

        
            int title = 1;
            //Print() << i <<"\n";
            fstream fout;
            //fout.open("constant/HSN_"+ tdb_phase[i] +".csv");
            fout.open("constant/Composition_"+ tdb_phase[0] +".csv");
            
            if(title==1) 
            {
                getline(fout, line);
                //Print()<<line<<"\n";
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
                cval.push_back(data);
                //Print()<<data[2]<<"\n";
                data.clear();
            }

            //Print() << Aval[77][0]<<"\n";

            conc_Sol.resize(ntp-1);
            conc_Sol[0].resize(cval.size()-1);
            conc_Liq.resize(ntp-1);
            conc_Liq[0].resize(cval.size()-1);
            temprt.resize(ntp-1);
            temprt[0].resize(cval.size()-1);

            for(int j=0; j < cval.size() - 1; j++)
            {
                //Print()<<j<<"\n";
                //Print() << conc_Al[i].size()<<"\n";
                conc_Sol[0][j] = cval[j][1];
                conc_Liq[0][j] = cval[j][2];
                //Print()<<conc_Al[i][j]<<"\n";
                temprt[0][j] = cval[j][0];


                //Print()<<conc_Li[i][j]<<","<<conc_Al[i][j]<<"\n";
            }

            cval.clear();
            
            fout.close();
        
        
        
    }
}

double findc(int phase, double temp)
{
    A_accel_ptr = gsl_interp_accel_alloc();
    //Print() <<"Size of conc_Al: " <<conc_Al[phase].size() <<"\n";
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Sol[0].size());

    
    double x_array[conc_Sol[0].size()];
    double y_array[conc_Sol[0].size()];

    for (int i=0; i < conc_Sol[0].size(); i++)
    {
        x_array[i] = temprt[0][i];
        if(phase == 0){
        y_array[i] = conc_Sol[0][i];
        }
        else{
        y_array[i] = conc_Liq[0][i];
        }
    }
    
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Sol[0].size());
    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y;
}

void getc(){
    readc();

    conc = Vector <Vector<Real>>(nump,Vector<Real>(nump,0));
    conc[0][0] = findc(0,T);
    conc[1][1] = findc(1,T);
    
}

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

	BL_PROFILE("function_F_04_function_A()");

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
    
    getc();

   //BB = 2.0*A1[1][0]*cguess[1][2] - 2.0*A1[0][0]*cguess[0][2];

    for(int a=0; a<nump-1; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
    
            //B[a] = 2.0*A[nump-1][i][j]*c_guess[nump-1][nump-1][a] - 2.0*A[a][i][j]*c_guess[i][j][a];
            B[a] = 2.0*A[nump-1][i][j]*conc[nump-1][nump-1] - 2.0*A[a][i][j]*conc[i][j];

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
    
            //D[a] = A[a][i][j]*c_guess[i][j][a]*c_guess[i][j][a] - A[nump-1][i][j]*c_guess[nump-1][nump-1][a]*c_guess[nump-1][nump-1][a];
            D[a] = A[a][i][j]*conc[i][j]*conc[i][j] - A[nump-1][i][j]*conc[nump-1][nump-1]*conc[nump-1][nump-1];
           
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
