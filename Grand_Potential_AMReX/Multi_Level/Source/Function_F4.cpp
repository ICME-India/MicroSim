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
#include "AmrCoreGP.H"
using namespace amrex;
using namespace std;

//####################################################################################################################
void
AmrCoreGP::readC()
{       
    string line, value;
    Vector<double> data;
    Vector<Vector<double>> cval;

    for(int a=0; a<nump-1;a++)
    {
        int title = 1;
        fstream fout;
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
        for(int r=0; r<numcom-1;r++)
        {
            for(int b=0; b < cval.size() - 1; b++)
            {   
                conc_Sol[a][b][r] = cval[b][1+r];
                conc_Liq[a][b][r] = cval[b][numcom+r];
                if(r==0)
                {
                    temprt[a][b] = cval[b][0];
                }
            }
                
        }
        cval.clear();
        fout.close();
    }
}

//####################################################################################################################
double 
AmrCoreGP::findC(int phase, double temp, int i)
{   
    if(phase!=(nump-1))
    {
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
    
    else if(phase==(nump-1))
    {
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

//####################################################################################################################
void 
AmrCoreGP::getc()
{
    readC();
    conc = Vector<Vector<Real>>(nump,Vector<Real>(numcom-1,0.0));
    conceq = Vector<Vector<Real>>(nump,Vector<Real>(numcom-1,0.0));
    for(int a=0; a<nump; a++)
    {
        for(int i=0; i<numcom-1; i++)
        {
            conc[a][i] = findC(a,T,i);
            conceq[a][i] = findC(a,Teq,i);
        } 
    }
}

//####################################################################################################################
void
AmrCoreGP::readA()
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

//####################################################################################################################
double
AmrCoreGP::findA(int phase, double temp, int o, int r)
 {
    A_accel_ptr = gsl_interp_accel_alloc();
    A_spline_ptr = gsl_spline_alloc(gsl_interp_cspline, A_values[phase].size());
    double x_array[A_values[phase].size()];
    double y_array[A_values[phase].size()];
    int u=0;
    Vector<Vector<Real>>L(numcom-1,Vector<Real>(numcom-1,0.0));
    for(int p=0; p<numcom-1; p++)
    {
        for(int m=p+1; m<numcom-1; m++)
        {
            L[p][m] = numcom-1+u;
            L[m][p] = numcom-1+u;
            u++;
        }
    }
    if(o==r)
    {
        for (int m=0; m < A_values[phase].size(); m++)
        {
            x_array[m] = A_temp[phase][m];
            y_array[m] = A_values[phase][m][o];
        }
    }
    else
    {
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

//####################################################################################################################
void
AmrCoreGP::function_F_04_function_A()
{
	BL_PROFILE("function_F_04_function_A()");
    readA();
    A = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    Aeq = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    for (int a=0; a < nump; a++) 
    {
        for(int i=0; i<numcom-1; i++) 
        {
            for(int j=0; j<numcom-1; j++) 
            {
                if (i==j) 
                {
                    A[a][i][j] = 0.5*findA(a,T,i,j);
                    Aeq[a][i][j] = 0.5*findA(a,Teq,i,j);
                } else 
                {
                    A[a][i][j] = findA(a,T,i,j);
                    Aeq[a][i][j] = findA(a,Teq,i,j);
                }
                //Print()<<"A["<<a<<","<<i<<","<<j<<"]="<<A[a][i][j]<<"\n";
            }
        }
    }
    //Print()<< "I am in function_F_04_function_A"<<'\n';
}

//####################################################################################################################
void
AmrCoreGP::function_F_04_function_B()
{
    getc();
    B = Vector<Vector<Real>>(nump,Vector<Real>(numcom-1,0.0));
    double sum_c = 0.0;
    for(int a=0; a< nump; a++)
    {
        if (a != (nump-1)) 
        {
            for(int i = 0; i< numcom-1; i++)
            {    
                for (int k=0; k < numcom-1; k++) 
                {
                    if (k!=i) 
                    {
                        sum_c += A[nump-1][k][i]*conc[nump-1][k] - A[a][k][i]*conc[a][k];
                    }
                }
                B[a][i] = (2.0*(A[nump-1][i][i]*conc[nump-1][i] - A[a][i][i]*conc[a][i]) + sum_c);
                sum_c=0.0;
                //Print()<<"B["<<a<<" , "<<i<<"] = "<<B[a][i]<<"\n";
            }
        }
    }
    //Print()<< "I am in function_F_04_function_B"<<'\n';

}

//####################################################################################################################
void
AmrCoreGP::function_F_04_function_C()
{
    C=Vector<Real>(nump,0.0);
    for(int a =0; a<nump;a++)
    {
        if (a != (nump-1)) 
        {
            for(int i=0; i<numcom-1; i++)
            {
                for(int j=0; j<numcom-1; j++)
                {
                    if (i <= j) 
                    {
                        C[a] += (A[a][i][j]*conc[a][i]*conc[a][j] - A[nump-1][i][j]*conc[nump-1][i]*conc[nump-1][j]);
                    }
                }
            }
        }
        //Print()<<"C["<<a<<"] = "<<C[a]<<"\n";
    }
    //Print()<< "I am in function_F_04_function_C"<<'\n';
}

//####################################################################################################################
void
AmrCoreGP::function_F_04_dc_dmu()
{
    dcdmu = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    for(int a=0; a<nump; a++)
    {
        Array2D <Real,0,compcount-2,0,compcount-2,Order::C> muc{};
        Array2D <Real,0,compcount-2,0,compcount-2,Order::C> inv_muc{};
        for (int l=0;l<numcom-1;l++) 
        {
            for (int m=0;m<numcom-1;m++) 
            {
                if (l==m) 
                {
                    muc(l,m)=2.0*A[a][l][m];
                } 
                else 
                {
                    muc(l,m)=A[a][l][m];
                       
                }
            }
        }

        if(numcom==2)
        {
            dcdmu[a][0][0] = 1/muc(0,0);
        }

        if(numcom==3)
        {
            double det = muc(0,0)*muc(1,1) - muc(0,1)*muc(1,0);
            dcdmu[a][0][0] = muc(1,1)/det;
            dcdmu[a][0][1] = -muc(0,1)/det;
            dcdmu[a][1][0] = -muc(1,0)/det;
            dcdmu[a][1][1] = muc(0,0)/det;
        }

        if(numcom>3)
        {
            mat_inv(muc,inv_muc,numcom);
            for(int m=0; m<numcom-1; m++)
            {
                for(int n=0; n< numcom-1; n++)
                {
                    dcdmu[a][m][n] = inv_muc(m,n);
                }
            }
        }       
     }   
     //Print()<< "I am in function_F_04_function_dc_dmu"<<'\n';
}

//####################################################################################################################
