#ifndef _READINPUT_H_
#define _READINPUT_H_

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <iostream>
#include <fstream>
#include <string>

#include "Variables.h"
using namespace std;

using namespace amrex;

void convert_string(vector<string> &x, vector<Real> &y) {
transform(x.begin(),x.end(), back_inserter(y), [](const string & astr){return stod(astr); });
}


void readinput(){

BL_PROFILE("readinput()");

//Reading from infile starts-----------------------------------------------------------------------------------------------

{	    Print()<<"Entering read input\n";
        ParmParse pp;
        pp.get("DIMENSION",dim);
        pp.get("MESH_X",ncellx);
        pp.get("MESH_Y",ncelly);
        pp.get("MESH_Z",ncellz);
        pp.get("DELTA_X",dx);
        pp.get("DELTA_Y",dy);
        pp.get("DELTA_Z",dz);
        pp.get("DELTA_t",dt);
        Print()<<"dt\n";
        pp.get("NUMPHASES",nump);
        pp.get("NUMCOMPONENTS",numcom);
        pp.get("NTIMESTEPS",nsteps);
        pp.query("NSMOOTH",nsmooth);
        pp.get("SAVET",savet);
        pp.get("STARTTIME",startt);
        pp.get("RESTART",restart);
        pp.get("numworkers",numworkers);
        pp.getarr("COMPONENTS",comp);
        pp.getarr("PHASES",phase);
        pp.getarr("GAMMA",gammaa);
        Print()<<"gammaa\n";
        int ss = pp.countname("DIFFUSIVITY");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("DIFFUSIVITY", i, n1, 0); 
    	diffu.push_back(n1);
    	}
        ss=0;
        pp.query("R",R);
        pp.get("V",Vm);
        pp.query("ELASTICITY",ELASTICITY);
        pp.query("rho",ro);
        pp.query("damping_factor",dampfac);
        pp.query("max_iterations",MAX_ITERATIONS);
        ss = pp.countname("EIGEN_STRAIN");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("EIGEN_STRAIN", i, n1, 0); 
    	egstr.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("VOIGT_ISOTROPIC");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("VOIGT_ISOTROPIC", i, n1, 0); 
    	voigiso.push_back(n1);
    	}
        ss=0;
        Print()<<"voig\n";

        ss = pp.countname("BOUNDARY");
    	for(int i=0;i<ss;i++){
        Vector<std::string> n2;
    	pp.getktharr ("BOUNDARY", i, n2, 0); 
    	bound.push_back(n2);
    	}
        ss=0;

        ss = pp.countname("BOUNDARY_VALUE");
    	for(int i=0;i<ss;i++){
        Vector<std::string> n2;
    	pp.queryktharr ("BOUNDARY_VALUE", i, n2, 0); 
    	boundval.push_back(n2);
    	}
        ss=0;

        pp.get("ISOTHERMAL",isothermal);
        pp.query("BINARY",binary);
        pp.query("TERNARY",ternary);
        pp.query("DILUTE",dilute);
        pp.get("T",T);
        pp.query("WRITEFORMAT",writeformat);
        pp.get("WRITEHDF5",writehdf5);
        pp.get("TRACK_PROGRESS",trackprog);
        pp.get("epsilon",eps);
        pp.get("tau",tau);
        pp.get("Tau",Tau);
        pp.get("Function_anisotropy",funcANI);
        pp.get("Anisotropy_type",ANItype);
        pp.getarr("dab",dab);
        //pp.get("Rotation_matrix",s_rotmat);
        Print()<<"dab\n";
        ss = pp.countname("Rotation_matrix");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("Rotation_matrix", i, n1, 0); 
    	rotmat.push_back(n1);
    	}
        ss=0;

        pp.get("Function_W",funcW);
      
        pp.queryarr("Gamma_abc",gammaa_abc);
        
        pp.query("Shift",shiftdom);
        pp.query("Shiftj",shiftj);
        pp.get("Writecomposition",writecomp);
        pp.query("Noise_phasefield",noise_pf);
        pp.query("Amp_Noise_Phase",amp_noise_phase);
        pp.get("Equilibrium_temperature",Teq);
        pp.get("Filling_temperature",Tfill);
        pp.queryarr("Tempgrady",tempgrady);
        pp.get("Function_F",funcf);

        ss = pp.countname("A");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("A", i, n1, 0); 
    	A1.push_back(n1);
    	}
        ss=0;
        Print()<<"A\n";

        ss = pp.countname("ceq");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("ceq", i, n1, 0); 
    	ceq.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("cfill");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("cfill", i, n1, 0); 
    	cfill.push_back(n1);
    	}
        ss=0;

        if (funcf != 4){
        ss = pp.countname("c_guess");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("c_guess", i, n1, 0); 
    	cguess.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("slopes");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.queryktharr ("slopes", i, n1, 0); 
    	slopes.push_back(n1);
    	}
        ss=0;
        }

        pp.get("num_thermo_phases",ntp);
        pp.get("tdbfname",tdbname);
        pp.getarr("tdb_phases",tdb_phase);
        pp.getarr("phase_map",phasemap);
        Print()<<"PM\n";

        ss = pp.countname("FILLCUBE");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("FILLCUBE", i, n1, 0); 
    	cube.push_back(n1);
    	}
        ss=0;
        
        ss = pp.countname("FILLCYLINDER");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("FILLCYLINDER", i, n1, 0); 
    	cylinder.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLSPHERE");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("FILLSPHERE", i, n1, 0); 
    	sphere.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLELLIPSE");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("FILLELLIPSE", i, n1, 0); 
    	ellipse.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLCYLINDERRANDOM");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("FILLCYLINDERRANDOM", i, n1, 0); 
    	cylrand.push_back(n1);
    	}
        ss=0;

        ss = pp.countname("FILLSPHERERANDOM");
    	for(int i=0;i<ss;i++){
        Vector<Real> n1;
    	pp.getktharr ("FILLSPHERERANDOM", i, n1, 0); 
    	sphrand.push_back(n1);
    	}
        ss=0;

        pp.queryarr("FILLCUBERANDOM",cuberand);

        pp.queryarr("FILLVORONOI2D",voronoi2D);

        pp.queryarr("FILLCUBEPATTERN",cubepat);

    }

    Print()<<"Done Parse\n";
    //Print()<<"Done read pp\n";

//Reading from infile finished-----------------------------------------------------------------------------------------------

//Extracting data from the strings----------------------------------------------------------------------------------------------- 

    //dim = stoi(s_dim);
    //ncellx = stoi(s_ncellx);
    //ncelly = stoi(s_ncelly);
    //ncellz = stoi(s_ncellz);
    //dx = stod(s_dx);
    //dy = stod(s_dy);
    //dz = stod(s_dz);
    //dt = stod(s_dt);
    //nump = stoi(s_nump);
    //numcom = stoi(s_numcom);
    //nsteps = stoi(s_nsteps);
    //nsmooth = stoi(s_nsmooth);
    //savet = stoi(s_savet);
    //startt = stoi(s_startt);
    //restart = stoi(s_restart);
    //numworkers = stoi(s_numworkers);

    //Print()<<"till numwork\n";

//Defining the characters that must be removed from the string-----------------------------------------------------------------------------------------------

    //string chars = ";( )";

//Storing name of the components-----------------------------------------------------------------------------------------------

    // for(char c:chars){
    //     s_comp.erase(remove(s_comp.begin(),s_comp.end(),c),s_comp.end());
    // }
     std::stringstream ss;
    // ss.str(s_comp);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     comp.push_back(substr);
    // }
    // ss.str("");
    // ss.clear(); 

//Storing name of the phases-----------------------------------------------------------------------------------------------

    // for(char c:chars){
    //     s_phase.erase(remove(s_phase.begin(),s_phase.end(),c),s_phase.end());
    // }
    // ss.str(s_phase);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     phase.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();

    unsigned int count{0};

//Storing the value for gamma-----------------------------------------------------------------------------------------------

    //Print()<<"gamma_start\n";

    // for(char c:chars){
    //         s_gammaa.erase(remove(s_gammaa.begin(),s_gammaa.end(),c),s_gammaa.end());
    //     }
    
    // ss.str("");
    // ss.clear();
    // ss.str(s_gammaa);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     val.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();
    // for (int i = 0; i < val.size(); i++)
    // {
    //     gammaa.push_back(stod(val[i]));
    // }
    // val.clear();

    gam = Vector<Vector<Real>>(nump,Vector<Real>(nump,0.0));
    int m=0;
    for (int a=0; a<nump; a++){
        for(int b=0; b<nump; b++){
            if(a!=b && a<b){
                gam[a][b] = gammaa[m];
                m++;
                //Print()<<"a: "<<a<<"b: "<<b<<" = "<<gam[a][b]<<"\n";
           }
        }
    }

    Print()<<"gamma_end\n";


//Storing Diffusivities of each phase-----------------------------------------------------------------------------------------------

    // while(count<s_diff.size()){
    //     for(char c:chars){
    //         s_diff[count].erase(std::remove(s_diff[count].begin(),s_diff[count].end(),c),s_diff[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // int dul = 0;
    // while(dul<s_diff.size()){
    //     ss.str(s_diff[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     diffu.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }

    diff = Vector<Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(numcom-1,Vector<Real>(numcom-1,0.0)));
    int cnt=0;
    for (int a=0; a<nump; a++){
        for(int k=0; k<numcom-1; k++){
            for(int l=0; l<numcom-1; l++){
                if(k==l){
                    diff[a][k][l] = diffu[a][2+cnt];
                    cnt++;
                    //Print()<<diffu[a][0]<<" , "<<diffu[a][1]<<" , "<<diffu[a][2]<<"\n";
                } 
            }
        }
        cnt=0;
        //diff[diffu[a][1]] = diffu[a][2];
    }
    
    count=0;

    Print()<<"Diff done\n";
    
//Storing R and Vm----------------------------------------------------------------------------------------------- 

    // R=stod(s_R);
    // Vm=stod(s_Vm);

//Elastic properties-----------------------------------------------------------------------------------------------

    // ELASTICITY = stoi(s_elast);
    // if(ELASTICITY==1){
    //     ro = stod(s_ro);
    //     dampfac = stod(s_dampfac);
    //     MAX_ITERATIONS = stoi(s_maxit);
    // }

    // while(count<s_egstr.size()){
    //     for(char c:chars){
    //         s_egstr[count].erase(std::remove(s_egstr[count].begin(),s_egstr[count].end(),c),s_egstr[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul<s_egstr.size()){
    //     ss.str(s_egstr[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     egstr.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    
    //count=0;
   

    // while(count<s_voigiso.size()){
    //     for(char c:chars){
    //         s_voigiso[count].erase(std::remove(s_voigiso[count].begin(),s_voigiso[count].end(),c),s_voigiso[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul<s_voigiso.size()){
    //     ss.str(s_voigiso[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     voigiso.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    
    //count=0;

//Storing boundary conditions-----------------------------------------------------------------------------------------------

    // while(count<s_bound.size()){
    //     for(char c:chars){
    //         s_bound[count].erase(std::remove(s_bound[count].begin(),s_bound[count].end(),c),s_bound[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul<s_bound.size()){
    //     ss.str(s_bound[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     bound.push_back(val);
    //     dul++;
    //     val.clear();
    // }
    
    // count=0;

//Storing boundary values---------------------------------------------------------------------------------------------------------

    // while(count<s_boundval.size()){
    //     for(char c:chars){
    //         s_boundval[count].erase(std::remove(s_boundval[count].begin(),s_boundval[count].end(),c),s_boundval[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul<s_boundval.size()){
    //     ss.str(s_boundval[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     boundval.push_back(val);
    //     dul++;
    //     val.clear();
    // }
    // count=0;

    //Print()<<"bval\n";

//Storing iso, bin, dil, T, writeformat, trackprog, eps, tau, Tau, funcANI, ANItype, dab----------------------------------------------------------------------------------------------- 

    // isothermal=stoi(s_isothermal);
    // if(s_binary!=""){
    //     binary=stoi(s_binary);
    // }
    // if(s_ternary!=""){
    //     ternary=stoi(s_ternary);
    // }
    // if(s_dilute!=""){
    //     dilute=stoi(s_dilute);
    // }
    // T=stoi(s_T);
    // for(char c:chars){
    //     s_writeformat.erase(remove(s_writeformat.begin(),s_writeformat.end(),c),s_writeformat.end());
    // }
    // writeformat = s_writeformat;
    // writehdf5=stoi(s_writehdf5);
    // trackprog=stoi(s_trackprog);
    // eps=stod(s_eps);
    // tau=stod(s_tau);
    // for(char c:chars){
    //     s_Tau.erase(remove(s_Tau.begin(),s_Tau.end(),c),s_Tau.end());
    // }
    // Tau=stod(s_Tau);

    // funcANI=stoi(s_funcANI);
    // ANItype=stoi(s_ANItype);
    // for(char c:chars){
    //     s_dab.erase(remove(s_dab.begin(),s_dab.end(),c),s_dab.end());
    // }
    
    // ss.str("");
    // ss.clear();
    // ss.str(s_dab);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     val.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();
    // for (int i = 0; i < val.size(); i++)
    // {
    //     dab.push_back(stod(val[i]));
    // }
    // val.clear();

    dabb = Vector<Vector<Real>>(nump,Vector<Real>(nump,0.0));
    m=0;
    for (int a=0; a<nump; a++){
        for(int b=0; b<nump; b++){
            if(a!=b && a<b){
                dabb[a][b] = dab[m];
                m++;
           }
        }
    }

    Print()<<"dab\n";

//Storing rotation matrix-----------------------------------------------------------------------------------------------

    count=0;

    // while(count<s_rotmat.size()){
    //     for(char c:chars){
    //         s_rotmat[count].erase(remove(s_rotmat[count].begin(),s_rotmat[count].end(),c),s_rotmat[count].end());
    //     }
    // count++;
    // }

    // count=0;
    // dul = 0;
    // while(dul<s_rotmat.size()){
    // ss.str(s_rotmat[dul]);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     val.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();
    // convert_string(val, dbl);
    // rotmat.push_back(dbl);
    // dul++;
    // val.clear();
    // dbl.clear();
    // }

    //Print()<<"rotmat\n";

//Storing funcW, gamma_abc, shift of domain, writecomp, noise, Teq, Tfill-----------------------------------------------------------------------------------------------
    
    //funcW=stoi(s_funcW);

    //Print()<<"gamma_abc_start\n";

    // if(nump==3){
    //     for(char c:chars){
    //         gamma_abc.erase(remove(gamma_abc.begin(),gamma_abc.end(),c),gamma_abc.end());
    //     }
    //     gammaa_abc.push_back(stod(gamma_abc));
    // }

    if(nump>=3){
    // for(char c:chars){
    //         gamma_abc.erase(remove(gamma_abc.begin(),gamma_abc.end(),c),gamma_abc.end());
    //     }

    gam_abc = Vector <Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(nump,0.0)));
    m=0;


    for(int i=0; i < nump; i++) {
        for (int j=i+1; j < nump; j++) {
            for (int k=j+1; k < nump; k++) {
                gam_abc[i][i][i] = 0.0;
                gam_abc[i][j][j] = 0.0;
                gam_abc[i][k][k] = 0.0;
                
                gam_abc[i][j][k] = gammaa_abc[m];
                gam_abc[i][k][j] = gam_abc[i][j][k];
                gam_abc[j][i][k] = gam_abc[i][j][k];
                gam_abc[j][k][i] = gam_abc[i][j][k];
                gam_abc[k][i][j] = gam_abc[i][j][k];
                gam_abc[k][j][i] = gam_abc[i][j][k];
                
                m++;
            }
        }
    }

    Print()<<"gamma_abc_end_all\n";
    }
   
    // shiftdom=stod(s_shiftdom);
    // shiftj=stod(s_shiftj);
    // writecomp=stoi(s_writecomp);
    // noise_pf=stoi(s_noise_pf);
    // amp_noise_phase=stod(s_amp_noise_phase);
    // Teq=stod(s_Teq);
    // Tfill=stod(s_Tfill);

//Storing tempgrad-----------------------------------------------------------------------------------------------

    // for(char c:chars){
    //     s_tempgrady.erase(remove(s_tempgrady.begin(),s_tempgrady.end(),c),s_tempgrady.end());
    // }
    // ss.str("");
    // ss.clear();
    // ss.str(s_tempgrady);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     val.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();
    // for (int i = 0; i < val.size(); i++)
    // {
    //     tempgrady.push_back(stod(val[i]));
    // }
    // val.clear();

//Storing value for functionF-----------------------------------------------------------------------------------------------

    //funcf=stoi(s_funcf);

//Storing A matrix-----------------------------------------------------------------------------------------------

    // while(count< s_A.size()){
    //     for(char c:chars){
    //         s_A[count].erase(std::remove(s_A[count].begin(),s_A[count].end(),c),s_A[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul< s_A.size()){
    //     ss.str(s_A[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     A1.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

//Storing ceq-----------------------------------------------------------------------------------------------

    // while(count< s_ceq.size()){
    //     for(char c:chars){
    //         s_ceq[count].erase(std::remove(s_ceq[count].begin(),s_ceq[count].end(),c),s_ceq[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul = 0; 
    // while(dul< s_ceq.size()){
    //     ss.str(s_ceq[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
        
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     ceq.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

    c_eq = Vector <Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(numcom,0)));
    for(auto k=2; k<ceq[0].size(); k++){
    	for(auto i=0; i<(int)ceq.size(); i++){
            auto v = (long unsigned int)ceq[i][0];
            auto w = (long unsigned int)ceq[i][1]; 
    		c_eq[v][w][count] = ceq[i][k];
    
    }
    count++;
    }	  
    count=0;

    Print()<<"ceq done\n";

//Storing cfill---------------------------------------------------------------------------------------------------------

    // while(count< s_cfill.size()){
    //     for(char c:chars){
    //         s_cfill[count].erase(std::remove(s_cfill[count].begin(),s_cfill[count].end(),c),s_cfill[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0; 
    // while(dul< s_cfill.size()){
    //     ss.str(s_cfill[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     cfill.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

    c_fill = Vector <Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(numcom,0)));
    for(auto k=2; k<cfill[0].size(); k++){
    	for(auto i=0; i<(int)cfill.size(); i++){
            auto v = (long unsigned int)cfill[i][0];
            auto w = (long unsigned int)cfill[i][1]; 
    		c_fill[v][w][count] = cfill[i][k];
    
    }
    count++;
    }	  
    count=0;

    Print()<<"cfill done\n";

//Storing cguess and slope if not using function F4-----------------------------------------------------------------------------------------------

     if (funcf != 4){
    // while(count< s_cguess.size()){
    //     for(char c:chars){
    //         s_cguess[count].erase(std::remove(s_cguess[count].begin(),s_cguess[count].end(),c),s_cguess[count].end());
    //     }
    //     count++;
    // }
    // count=0; 
    // dul=0;
    // while(dul< s_cguess.size()){
    //     ss.str(s_cguess[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str(""); 
    //     ss.clear();
    //     convert_string(val, dbl);
    //     cguess.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

    c_guess = Vector <Vector<Vector<Real>>>(nump,Vector<Vector<Real>>(nump,Vector<Real>(numcom,0)));
    for(auto k=2; k<cguess[0].size(); k++){
    	for(auto i=0; i<(int)cguess.size(); i++){
            auto v = (int)cguess[i][0];
            auto w = (int)cguess[i][1]; 
    		c_guess[v][w][count] = cguess[i][k];
    
    }
    count++;
    }	  
    count=0;

    Print()<<"cguess done\n";
    // while(count< s_slopes.size()){
    //     for(char c:chars){
    //         s_slopes[count].erase(std::remove(s_slopes[count].begin(),s_slopes[count].end(),c),s_slopes[count].end());
    //     }
    //     count++;
    // }
    // count=0; 
    // dul=0;
    // while(dul< s_slopes.size()){
    //     ss.str(s_slopes[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     slopes.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;
     }

//Storing number of thermodynamic phases, tdbfname, tdbphase and phase map-----------------------------------------------------------------------------------------------

    //ntp = stoi(s_ntp);

    // for(char c:chars){
    //     s_tdbname.erase(remove(s_tdbname.begin(),s_tdbname.end(),c),s_tdbname.end());
    // }
    // tdbname=s_tdbname;

    // for(char c:chars){
    //     s_tdbphase.erase(remove(s_tdbphase.begin(),s_tdbphase.end(),c),s_tdbphase.end());
    // }
    
    // ss.str(s_tdbphase);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     tdb_phase.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();

    // for(char c:chars){
    //     s_phasemap.erase(remove(s_phasemap.begin(),s_phasemap.end(),c),s_phasemap.end());
    // }
   
    // ss.str(s_phasemap);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     phasemap.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();

    // while(count< s_cube.size()){
    //     for(char c:chars){
    //         s_cube[count].erase(std::remove(s_cube[count].begin(),s_cube[count].end(),c),s_cube[count].end());
    //     }
    //     count++;
    // }
    // count=0;

//Storing filling data-----------------------------------------------------------------------------------------------

    // dul=0;
    // while(dul< s_cube.size()){
    //     ss.str(s_cube[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     cube.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    
    // count=0;

    // while(count< s_cylinder.size()){
    //     for(char c:chars){
    //         s_cylinder[count].erase(std::remove(s_cylinder[count].begin(),s_cylinder[count].end(),c),s_cylinder[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul< s_cylinder.size()){
    //     ss.str(s_cylinder[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     cylinder.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    
    // count=0;
    

    // while(count< s_ellipse.size()){
    //     for(char c:chars){
    //         s_ellipse[count].erase(std::remove(s_ellipse[count].begin(),s_ellipse[count].end(),c),s_ellipse[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul< s_ellipse.size()){
    //     ss.str(s_ellipse[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     ellipse.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // //Print()<<ellipse.size()<<"\n";
    // count=0;

    // while(count< s_sphere.size()){
    //     for(char c:chars){
    //         s_sphere[count].erase(std::remove(s_sphere[count].begin(),s_sphere[count].end(),c),s_sphere[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul< s_sphere.size()){
    //     ss.str(s_sphere[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     sphere.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

    // while(count< s_cylrand.size()){
    //     for(char c:chars){
    //         s_cylrand[count].erase(std::remove(s_cylrand[count].begin(),s_cylrand[count].end(),c),s_cylrand[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul< s_cylrand.size()){
    //     ss.str(s_cylrand[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     cylrand.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

    // while(count< s_sphrand.size()){
    //     for(char c:chars){
    //         s_sphrand[count].erase(std::remove(s_sphrand[count].begin(),s_sphrand[count].end(),c),s_sphrand[count].end());
    //     }
    //     count++;
    // }
    // count=0;
    // dul=0;
    // while(dul< s_sphrand.size()){
    //     ss.str(s_sphrand[dul]);
    //     while(ss.good()){
    //         std::string substr;
    //         getline(ss,substr,',');
    //         val.push_back(substr);
    //     }
    //     ss.str("");
    //     ss.clear();
    //     convert_string(val, dbl);
    //     sphrand.push_back(dbl);
    //     dul++;
    //     val.clear();
    //     dbl.clear();
    // }
    // count=0;

    // if(s_cuberand != ""){
    // for(char c:chars){
    //     s_cuberand.erase(remove(s_cuberand.begin(),s_cuberand.end(),c),s_cuberand.end());
    // }
    // ss.str("");
    // ss.clear();
    // ss.str(s_cuberand);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     val.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();
    // for (int i = 0; i < val.size(); i++)
    // {
    //     cuberand.push_back(stod(val[i]));
    // }
    // val.clear();

    // }

    // if(s_cubepat != ""){
    // for(char c:chars){
    //     s_cubepat.erase(remove(s_cubepat.begin(),s_cubepat.end(),c),s_cubepat.end());
    // }
    // ss.str("");
    // ss.clear();
    // ss.str(s_cubepat);
    // while(ss.good()){
    //     std::string substr;
    //     getline(ss,substr,',');
    //     val.push_back(substr);
    // }
    // ss.str("");
    // ss.clear();
    // for (int i = 0; i < val.size(); i++)
    // {
    //     cubepat.push_back(stod(val[i]));
    // }
    // val.clear();
    // }

}

#endif