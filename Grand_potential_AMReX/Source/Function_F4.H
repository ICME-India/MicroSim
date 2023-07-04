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
    
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> cval;

        for(int i=0; i<nump-1;i++){
            int title = 1;
            //Print() << i <<"\n";
            fstream fout;
            //fout.open("constant/HSN_"+ tdb_phase[i] +".csv");
            fout.open("constant/Composition_"+ phasemap[i] +".csv");
            
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

            conc_Sol.resize(nump-1);
            conc_Sol[i].resize(cval.size()-1);
            conc_Liq.resize(nump-1);
            conc_Liq[i].resize(cval.size()-1);
            temprt.resize(nump-1);
            temprt[i].resize(cval.size()-1);

            for(int j=0; j < cval.size() - 1; j++)
            {
                //Print()<<j<<"\n";
                //Print() << conc_Al[i].size()<<"\n";
                conc_Sol[i][j] = cval[j][1];
                conc_Liq[i][j] = cval[j][2];
                //Print()<<conc_Al[i][j]<<"\n";
                temprt[i][j] = cval[j][0];


                //Print()<<conc_Li[i][j]<<","<<conc_Al[i][j]<<"\n";
            }

            cval.clear();
            
            fout.close();
        }

        //Print()<<"Temp  "<<temprt[0].size()<<"\n";
       // Print()<<"Al  "<<conc_Sol[0].size()<<"\n";
        //Print()<<"Li  "<<conc_Liq[0].size()<<"\n";
        
}

double findc(int phase, double temp)
{
    if(phase!=nump-1){
    
    A_accel_ptr = gsl_interp_accel_alloc();
    //Print() <<"Size of conc_Al: " <<conc_Al[phase].size() <<"\n";
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Sol[phase].size());

    
    double x_array[conc_Sol[phase].size()];
    double y_array[conc_Sol[phase].size()];

    for (int i=0; i < conc_Sol[phase].size(); i++)
    {
        x_array[i] = temprt[phase][i];
        y_array[i] = conc_Sol[phase][i];

    }
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Sol[phase].size());
    }
    
    else{
    A_accel_ptr = gsl_interp_accel_alloc();
    //Print() <<"Size of conc_Al: " <<conc_Al[phase].size() <<"\n";
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, conc_Liq[0].size());

    
    double x_array[conc_Liq[0].size()];
    double y_array[conc_Liq[0].size()];

    for (int i=0; i < conc_Liq[0].size(); i++)
    {
        x_array[i] = temprt[0][i];
        y_array[i] = conc_Liq[0][i];

    }
    gsl_spline_init(A_spline_ptr, x_array, y_array, conc_Liq[0].size());
    }

    double y = gsl_spline_eval(A_spline_ptr, temp, A_accel_ptr);
    return y;
}

void getc(){
    Print()<<"ValC1in"<<"\n";
    readc();
    Print()<<"ValC2in"<<"\n";

    conc = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));
    conceq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0)));

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                
                conc[a][i][j] = findc(a,T);

           } 
        }
    }

    Print()<<"ValC3in"<<"\n";

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                
                conceq[a][i][j] = findc(a,Teq);

           } 
        }
    }

    Print()<<"ValC4in"<<"\n";

    //Print()<<"conc"<<conc[0][0][0]<<"  "<<conc[1][0][0]<<"\n";
   // Print()<<"conceq"<<conceq[0][0][0]<<"  "<<conceq[1][0][0]<<"\n";
    
}

void readA()
{  
        string line, value;
        Vector<double> data;
        Vector<Vector<double>> Aval;

        for(int i=0; i<nump; i++)
        {
            int title = 1;
            //Print() << i <<"\n";
            fstream fout;
            fout.open("constant/HSN_"+ phasemap[i] +".csv");
            
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

            A_values.resize(nump);
            A_values[i].resize(Aval.size()-1);
            A_temp.resize(nump);
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

    Print()<<"ValA1in"<<"\n";

    readA();
    //Print()<<"In A"<<"\n";
    
    A = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    Aeq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                
                A[a][i][j] = findA(a,T);

                Print()<<"A: "<< A[a][i][j]<<"\n";

           } 
        }
    }

    Print()<<"ValA2in"<<"\n";

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
                
                Aeq[a][i][j] = findA(a,Teq);

                 Print()<<"Aeq: "<< Aeq[a][i][j]<<"\n";

           } 
        }
    }
    //Print()<<"A"<<A[0][0][0]<<"  "<<A[1][0][0]<<"\n";
    //Print()<<"Aeq"<<Aeq[0][0][0]<<"  "<<Aeq[1][0][0]<<"\n";
}

void function_F_04_function_B(){
    Print()<<"ValB1in"<<"\n";
    getc();
    Print()<<"ValB2in"<<"\n";
    
    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
    
            B.push_back(2.0*A[nump-1][i][j]*conc[nump-1][i][j] - 2.0*A[a][i][j]*conc[a][i][j]);

             Print()<<"B: "<< B[a]<<"\n";

           }
        }
    }
    Print()<<"ValB3in"<<"\n";
}


void function_F_04_function_D(){

    for(int a=0; a<nump; a++){
        for(int i=0; i<numcom-1; i++){
           for(int j=0; j<numcom-1; j++){
    
            D.push_back(A[a][i][j]*conc[a][i][j]*conc[a][i][j] - A[nump-1][i][j]*conc[nump-1][i][j]*conc[nump-1][i][j]);

            Print()<<"D " <<a<<": "<<D[a]<<"\n";

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
                dcdmu_eq.push_back(1.0/(2.0*Aeq[a][i][j]));
           }
        }
     }
    
}

void function_F_04_free_energy(){

    //     fe=Vector<Real>(nump,0);
    //  for(int a=0; a<nump; a++){
    //     for(int i=0; i<numcom-1; i++){
    //        for(int j=0; j<numcom-1; j++){
    //            fe[a] = A[a][i][j]*
    //        }
    //     }
    //  }
    
}


#endif
