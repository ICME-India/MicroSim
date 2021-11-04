#ifndef INITIALIZE_FUNCTIONS_SOLVERLOOP_H_
#define INITIALIZE_FUNCTIONS_SOLVERLOOP_H_

void initialize_functions_solverloop(){
  if (FUNCTION_F == 1) {    
       free_energy               = function_F_01_free_energy;
       dc_dmu                    = function_F_01_dc_dmu;
       c_mu                      = function_F_01_c_mu;
       Mu                        = function_F_01_Mu;
       dpsi                      = function_F_01_dpsi;
       function_A                = function_F_01_function_A;
       function_B                = function_F_01_function_B;
       function_C                = function_F_01_function_C;
       compute_chemicalpotential = function_F_01_compute_chemicalpotential;
  }
  if (FUNCTION_W == 1) {
    dwdphi        = function_W_01_dwdphi;
    dwdphi_smooth = function_W_01_dwdphi_smooth;
    OBSTACLE = 1;
    WELL = 0;
  }
  if (FUNCTION_W == 2) {
    dwdphi        = function_W_02_dwdphi;
    dwdphi_smooth = function_W_02_dwdphi_smooth;
    OBSTACLE = 0;
    WELL = 1;
  }
  if(FUNCTION_ANISOTROPY == 0) {
    dAdphi               = function_A_00_dAdphi;
    divdAdgradphi        = function_A_00_divdAdgradphi;
    dAdphi_smooth        = function_A_00_dAdphi;
    divdAdgradphi_smooth = function_A_00_divdAdgradphi;
    ANISOTROPY = 0;
  }
  if(FUNCTION_ANISOTROPY == 1) {
    dAdphi               = function_A_01_dAdphi;
    divdAdgradphi        = function_A_01_divdAdgradphi;
    dAdphi_smooth        = function_A_01_dAdphi_smooth;
    divdAdgradphi_smooth = function_A_01_divdAdgradphi_smooth;
    ANISOTROPY = 1;
  }
  if (FOLD == 4) {
    dAdq        = anisotropy_01_dAdq;
    function_ac = anisotropy_01_function_ac;
  }
  if (DIMENSION==2) {
    calculate_gradients                                     = calculate_gradients_2D;
    calculate_gradients_phasefield                          = calculate_gradients_phasefield_2D;
    calculate_gradients_concentration                       = calculate_gradients_concentration_2D;
    calculate_fluxes_concentration                          = calculate_fluxes_concentration_2D;
    calculate_divergence_concentration                      = calculate_divergence_concentration_2D;
    calculate_divergence_concentration_smooth               = calculate_divergence_concentration_smooth_2D;
    calculate_divergence_concentration_smooth_concentration = calculate_divergence_concentration_smooth_concentration_2D;
    calculate_divergence_phasefield                         = calculate_divergence_phasefield_2D;
    calculate_divergence_phasefield_smooth                  = calculate_divergence_phasefield_smooth_2D;
  }
}
#endif