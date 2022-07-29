#ifndef INITIALIZE_FUNCTIONS_SOLVERLOOP_H_
#define INITIALIZE_FUNCTIONS_SOLVERLOOP_H_

void initialize_functions_solverloop(){

  /* if (FUNCTION_F == 1) {    
       free_energy               = function_F_01_free_energy;
       dc_dmu                    = function_F_01_dc_dmu;
       c_mu                      = function_F_01_c_mu;
       Mu                        = function_F_01_Mu;
       dpsi                      = function_F_01_dpsi;
//        function_A                = function_F_01_function_A;
       function_B                = function_F_01_function_B;
       function_C                = function_F_01_function_C;
       init_propertymatrices     = function_F_01_init_propertymatrices;
  } */
  if ((FUNCTION_F == 2)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a<NUMPHASES; a++) {
      for (b=0; b<NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_02_free_energy;
    dc_dmu                    = function_F_02_dc_dmu;
    c_mu                      = function_F_02_c_mu;
    Mu                        = function_F_02_Mu;
    dpsi                      = function_F_02_dpsi;
    CL_Solve_phi_com_Function = CL_Solve_phi_com_Function_F_2;
  }
  if ((FUNCTION_F == 3)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a<NUMPHASES; a++) {
      for (b=0; b<NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_03_free_energy;
    dc_dmu                    = function_F_03_dc_dmu;
    c_mu                      = function_F_03_c_mu;
    Mu                        = function_F_03_Mu;
    dpsi                      = function_F_03_dpsi;
    function_A                = function_F_03_function_A;
    function_B                = function_F_03_function_B;
    function_C                = function_F_03_function_C;
    init_propertymatrices     = function_F_03_init_propertymatrices;
    CL_Solve_phi_com_Function = CL_Solve_phi_com_Function_F_3;
  }
   if ((FUNCTION_F == 4)) {
    thermo_phase = (long*)malloc(NUMPHASES*sizeof(long));
    long a, b;
    for (a=0; a<NUMPHASES; a++) {
      for (b=0; b<NUM_THERMO_PHASES; b++) {
        if (strcmp(phase_map[a],Phases_tdb[b])==0) {
          thermo_phase[a] = b;
          break;
        }
      }
    }
    free_energy               = function_F_04_free_energy;
    dc_dmu                    = function_F_04_dc_dmu;
    c_mu                      = function_F_04_c_mu;
    Mu                        = function_F_04_Mu;
    dpsi                      = function_F_04_dpsi;
    function_A                = function_F_04_function_A;
    function_B                = function_F_04_function_B;
    function_C                = function_F_04_function_C;
    init_propertymatrices     = function_F_04_init_propertymatrices;
    CL_Solve_phi_com_Function = CL_Solve_phi_com_Function_F_4;
  } 
}
#endif
