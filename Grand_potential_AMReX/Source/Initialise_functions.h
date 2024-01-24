#ifndef INITIALISE_FUNCTIONS_H_
#define INITIALISE_FUNCTIONS_H_

#include "Variables.h"
#include "Function_F4.h"
#include "Function_A.h"
#include "Function_W.h"
#include "Anisotropy_01.h"

void init_functions(MultiFab& phi_new){

    //Function pointers for Function F = 1----------------------------------------------------------------------------------
    if(funcf == 1){
        Print()<<"Function F1 will soon be added\n";
    }

    //Function pointers for Function F = 2----------------------------------------------------------------------------------
    if(funcf == 2){
        Print()<<"Function F2 will soon be added\n";
    }

    //Function pointers for Function F = 3----------------------------------------------------------------------------------
    if(funcf == 3){
        Print()<<"Function F3 will soon be added\n";
    }

    //Function pointers for Function F = 4----------------------------------------------------------------------------------
    if(funcf == 4){
        dc_dmu = function_F_04_dc_dmu;
        //c_mu = function_F_04_c_mu;
        Mu = function_F_04_Mu;
        function_A = function_F_04_function_A;
        function_B = function_F_04_function_B;
        function_C = function_F_04_function_C;
        //free_energy = function_F_04_free_energy;
    }

    //Function pointers for Function W = 1----------------------------------------------------------------------------------
    if(funcW == 1){
        Print()<<"This model currently supports only double well potential\n";
        abort();
    }

    //Function pointers for Function W = 2----------------------------------------------------------------------------------
    if(funcW == 2){
        dwdphi = function_W_02_dwdphi;
    }

    //Function pointers for Function A = 0----------------------------------------------------------------------------------
    if(funcANI == 0){
        if(dim == 2){
            aniso_term = function_A_00_iso_2D;
        }
        if(dim == 3){
            aniso_term = function_A_00_iso_3D;
        }
    }
    
    //Function pointers for Function A = 1----------------------------------------------------------------------------------
    if(funcANI == 1){
        if(dim == 2){
            aniso_term = function_A_01_ani_2D;
        }
        if(dim == 3){
            aniso_term = function_A_01_ani_3D;
        }
    }

    if(dim == 2){
        Chem_pot = dmudt_2D;
    }

    if(dim == 3){
        Chem_pot = dmudt_3D;
    }

    //Function pointers for initialising the domain----------------------------------------------------------------------------------
    
    if(cylinder.size()>0){
        Initialise_phi = init_phi_cyl;
        Initialise_phi(phi_new);
    }
    if(sphere.size()>0){
        Initialise_phi = init_phi_sph;
        Initialise_phi(phi_new);
    }
    if(cube.size()>0){
        Initialise_phi = init_phi_cube;
        Initialise_phi(phi_new);
    }
    if(ellipse.size()>0){
        Initialise_phi = init_phi_ellip;
        Initialise_phi(phi_new);
    }

    Init_liq(phi_new);
    
}
#endif