#ifndef INITIALISE_FUNCTIONS_H_
#define INITIALISE_FUNCTIONS_H_

#include "Variables.H"
#include "Function_F4.H"

void init_functions(MultiFab& phi_new){

    if(funcf == 1){
    Print()<<"Function F1 will soon be added"<<"\n";
    }

    if(funcf == 2){
    Print()<<"Function F2 will soon be added"<<"\n";
    }

    if(funcf == 3){
    Print()<<"Function F3 will soon be added"<<"\n";
    }

    if(funcf == 4){
    dc_dmu = function_F_04_dc_dmu;
    c_mu = function_F_04_c_mu;
    Mu = function_F_04_Mu;
    function_A = function_F_04_function_A;
    function_B = function_F_04_function_B;
    function_D = function_F_04_function_D;
    free_energy = function_F_04_free_energy;
    }

    Vector<int> flag(4,0);
    //Print()<<ellipse.size()<<"\n";
    //Print()<<"ellipse : "<<ellipse[0][0]<<" "<<ellipse[0][1]<<" "<<ellipse[0][2]<<" "<<ellipse[0][4]<<" "<<ellipse[0][5]<<" "<<"\n";
    if(cylinder.size()>0){
        flag[0]=1;
    }
    if(sphere.size()>0){
        flag[1]=1;
    }
    if(cube.size()>0){
        flag[2]=1;
    }
    if(ellipse.size()>0){
        flag[3]=1;
        //Print()<<"Flag 3 : "<<flag[3]<<"\n";
    }
    
    if(flag[0]==1){
        Initialise_phi = init_phi_cyl;
        Initialise_phi(phi_new);
    }
    if(flag[1]==1){
        Initialise_phi = init_phi_sph;
        Initialise_phi(phi_new);
    }
    if(flag[2]==1){
        Initialise_phi = init_phi_cube;
        Initialise_phi(phi_new);
    }
    if(flag[3]==1){
        Initialise_phi = init_phi_ellip;
        //Print()<<"Inside ellipse func point"<<"\n";
        Initialise_phi(phi_new);
    }

    Init_liq(phi_new);
    
}
#endif