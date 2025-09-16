#include <mpi.h>
#include <hdf5.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <endian.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include "tdbs/Thermo.h"
#include "tdbs/Thermo.c"
#include "functions/global_vars.h"
#include "functions/functions.h"
#include "functions/matrix.h"
#include "functions/utility_functions.h"
#include "functions/utility_io.h" // newly added
#include "functions/functionH.h"
#include "functions/functionF_01.h"
#include "functions/functionF_02.h"
#include "functions/functionF_03.h"
#include "functions/functionF_04.h"
#include "functions/functionF_05.h"
#include "functions/functionF_06.h"
#include "functions/functionF_elast.h"
#include "functions/functionQ.h"
#include "functions/anisotropy_01.h"
#include "functions/anisotropy_02.h"
#include "functions/anisotropy_gen.h" // newly added
#include "functions/functionW_01.h"
#include "functions/functionW_02.h"
#include "functions/function_A_00.h"
#include "functions/function_A_01.h"
#include "functions/functionTau.h"
#include "functions/functionD.h"
#include "functions/filling.h"
#include "functions/reading_input_parameters.h"
#include "functions/reading_input_parameters_new.h" // newly added
#include "functions/dynamic_input_parameters.h" // newly added
#include "functions/read_boundary_conditions.h"
#include "functions/print_input_parameters.h"
#include "functions/print_boundary_conditions.h"
#include "functions/initialize_variables.h"
#include "functions/free_variables.h"
#include "functions/free_variables_new.h" // newly added
#include "functions/fill_domain.h"
#include "functions/shift.h"
#include "functions/Temperature_gradient.h"
#include "solverloop/serialinfo_xy.h"
#include "solverloop/gradients.h"
#include "solverloop/simplex_projection.h"
#include "solverloop/calculate_gradients.h"
#include "solverloop/calculate_fluxes_concentration.h"
#include "solverloop/calculate_divergence_phasefield.h"
#include "solverloop/calculate_divergence_concentration.h"
#include "solverloop/calculate_divergence_stress.h"
#include "solverloop/calculate_LBM_distribution.h"
#include "solverloop/calculate_LBM_distribution3D.h" //  newly added
#include "solverloop/solverloop.h"
#include "solverloop/file_writer.h"
#include "solverloop/file_writer_3D.h"
#include "solverloop/file_writer_xmf.h" // newly added
#include "solverloop/mpiinfo_xyz.h"
#include "solverloop/mpiinfo_xyz_new.h"
#include "solverloop/boundary_mpi.h"
#include "solverloop/initialize_functions_solverloop.h"

int main(int argc, char * argv[]) 
{
  /// Start the clock
  clock_t clock_start = clock(), 
          clock_now;
  double  time_passed = 0;

  /// Initialize MPI Environment
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  /// assert the existence of Input and Filling Files
  ms_Assert_Terminate(!ms_FILE_IsRegular(argv[1]));
  ms_Assert_Terminate(!ms_FILE_IsRegular(argv[2]));

  // Read and Initialize Input Parameters
  if (USE_NEW_INPFILE) {
    ms_ReadInputParameters(argv[1]);
  }
  else {
    reading_input_parameters(argv);
    read_boundary_conditions(argv);
    msDyInp_Initialize(argv);
  }

  /// automatic mpi decomposition if decomposition isnt provided
  ms_mpi_Init_numworkers(argc, argv);

  /// Initialize the global tensors
  initialize_variables();
  
  /// Initialize the function pointers
  initialize_functions_solverloop();

  if (USE_GSL_MATINV)
  {
    matinv_gsl_Malloc(&dcdmu_gsl_data, NUMCOMPONENTS-1);
  }
  
  if (!((FUNCTION_F == 2) || (FUNCTION_F==5) || (GRAIN_GROWTH)))
  {
    init_propertymatrices(T);
  }

  if(taskid == MASTER)
  {
    for (a=0; a<NUMPHASES-1; a++) {
      for(k=0; k<NUMCOMPONENTS-1; k++) {
        printf("\n- [%d] slopes[%s][%s][%s]=%le", taskid, Phases[a], Phases[NUMPHASES-1], Components[k], slopes[a][NUMPHASES-1][k]);
        printf("\n- [%d] slopes[%s][%s][%s]=%le", taskid, Phases[a], Phases[a],           Components[k], slopes[a][a][k]);
      }
    }
  }

  // Build MPI derives types for send recieve calls
  Build_derived_type(gridinfo_instance, &MPI_gridinfo);

  if (ELASTICITY) {
    Build_derived_type_stress(iter_gridinfo_w_instance, &MPI_iter_gridinfo);
  }

  if (LBM) {
    Build_derived_type_lbm(&MPI_lbm_gridinfo);
  }

  if(taskid==MASTER) {
    // Create output directory
    mkdir("DATA",0777);
    if (!WRITEHDF5){
      for (n=0; n < numtasks; n++) {
        // create sub-directories for .vtk files
        sprintf(dirname,"DATA/Processor_%d",n);
        mkdir(dirname, 0777);
      }
    }

    print_input_parameters(argv);
    print_boundary_conditions(argv);

    serialinfo_xy();

    // If starting from 0, fill the domain using fill file
    if ((STARTTIME == 0) && (RESTART == 0)) {
      fill_domain(argv);
      // If fill file is not found it leads to BAD TERMINATION
    }
  }

  // names of all the relavent fields
  if (!USE_NEW_IO_BUFF)
  {
    populate_table_names();  
  }
  
  if(ELASTICITY) {
    for (int b=0; b<NUMPHASES; b++) {
      stiffness_phase_n[b].C11 = stiffness_phase[b].C11/stiffness_phase[NUMPHASES-1].C44;
      stiffness_phase_n[b].C12 = stiffness_phase[b].C12/stiffness_phase[NUMPHASES-1].C44;
      stiffness_phase_n[b].C44 = stiffness_phase[b].C44/stiffness_phase[NUMPHASES-1].C44;
    }
  }  

  // Initialize domain decomposition, allocate memory and do everything

  if (USE_NEW_MPIINFO)
  {
    ms_mpi_CreateUniverse();
  }
  else
  {
    Mpiinfo(taskid);
  }

  if (USE_NEW_IO_BUFF)
  {
    IO_Initialize();
    IO_MallocGlobal();
  }

  
  if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
    Calculate_Tau();
  } else {
    if ((FUNCTION_F ==5) && (!GRAIN_GROWTH)) {
      for (a=0; a < NUMPHASES-1; a++) {
        tau_ab[a][NUMPHASES-1] = (Lf*Lf)*epsilon*0.2222/(V*V*therm_cond*Teq);
        tau_ab[NUMPHASES-1][a] = tau_ab[a][NUMPHASES-1];
      }
      for (a=0; a < NUMPHASES-1; a++) {
        for (b=0; b < NUMPHASES-1; b++) {
          tau_ab[a][b] = (Lf*Lf)*epsilon*0.2222/(V*V*therm_cond*Teq);
        }
      }
    }
    if (taskid == MASTER){
      printf("\n - [%d] tau[0][NUMPHASES-1]=%le\n",taskid, tau_ab[0][NUMPHASES-1]);
    }
  }


  
  //Checking tdb functions
  
//   if (taskid == MASTER) {
//     double c_x;
//     double c[NUMCOMPONENTS-1];
//     double c_calc[NUMCOMPONENTS-1];
//     double mu[NUMCOMPONENTS-1];
//     double dpsi;    
//     char filename[1000];
//     double fe;
//     FILE *fp_check;
//     for (a=0; a<NUMPHASES; a++) {
//       sprintf(filename, "Thermodynamic_functions_%ld.dat", a);
//       fp_check = fopen(filename, "w");
//       for(c_x=0.01; c_x < 0.99;) {
//         c[0] = c_x;
//         Mu(c, T, a, mu);
//         dc_dmu(mu, c, T, a, dcdmu);
//         fe = free_energy(c, T, a);
//         dpsi = fe - mu[0]*c[0];
//         c_mu(mu, c, T, a, ceq[a][a]);
//         fprintf(fp_check, "%le %le %le %le %le %le\n", c_x, mu[0], dcdmu[0][0], fe,  dpsi, c_calc[0]);
//         c_x += 0.05;
//       }
//       fclose(fp_check);
//     }
//   }

  // restarting the simulation
  if ( RESTART && (STARTTIME > 0) ) {
    // read from data files
    if (WRITEHDF5) {
      readfromfile_mpi_hdf5(gridinfo_w, argv, numworkers, STARTTIME);
    } else {
      if (ASCII) {
        readfromfile_mpi(gridinfo_w, argv, STARTTIME);
      } else {
        readfromfile_mpi_binary(gridinfo_w, argv, STARTTIME);
      }
    }
    if (SHIFT) {
      shift_position = 0 ;
      FILE *fp;
      fp = fopen("DATA/shift.dat","r");
      char tempbuff[10000];
      while(fgets(tempbuff,10000,fp)) {
        sscanf(tempbuff ,"%ld %ld\n",&time_file, &position);
        if(time_file == STARTTIME){
          shift_position = position ;
          if(taskid ==  MASTER){
            printf("\n- shift_position = %ld  **FOUND**", position) ;
          }
          break ;
        }
      }
      // for(file_iter=0; file_iter <= STARTTIME/saveT; file_iter++) {
      //   fscanf(fp,"%ld %ld\n",&time_file, &position);
      // }
      fclose(fp);
      if(taskid ==  MASTER){
          printf("\n- shift_position = %ld  **FOUND**", position) ;
      }

    }
    //if (!ISOTHERMAL) {
      // temperature_gradientY.gradient_OFFSET = 0 ;
      // if (SHIFT) {
      //    temperature_gradientY.gradient_OFFSET = (temperature_gradientY.gradient_OFFSET) + (temperature_gradientY.velocity)*(STARTTIME*deltat) - shift_position*deltay;
      //    if(taskid ==  MASTER){
      //       printf("\n- [S] temperature_gradientY.gradient_OFFSET = %le", temperature_gradientY.gradient_OFFSET) ;
      //     }
      // } else {
      //   temperature_gradientY.gradient_OFFSET  = (temperature_gradientY.gradient_OFFSET) + (temperature_gradientY.velocity)*(STARTTIME*deltat);
      // }
    //}
  }

  // Apply the initial boundary conditions
  if(boundary_worker) {
   apply_boundary_conditions(taskid);
  }

  mpiexchange_left_right(taskid);
  mpiexchange_top_bottom(taskid);
  if (DIMENSION==3) {
    mpiexchange_front_back(taskid);
  }

  if (LBM) {
    // sync
    mpiexchange_top_bottom_lbm(taskid);
    mpiexchange_left_right_lbm(taskid);
    if (DIMENSION == 3) { 
      mpiexchange_front_back_lbm(taskid);
    }
    // apply boundary conditions
    if (boundary_worker) { 
      apply_boundary_conditions_lbm(taskid); 
    }
  }

  if (!RESTART && TEMPGRADY)
  {
      // printf("\n- restart = %ld, tgrady = %d", RESTART, TEMPGRADY);
      BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET;
      GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
      temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT + (workers_mpi.offset[Y]-workers_mpi.offset_y)*GRADIENT;
      
      fill_temperature_gradientY(gridinfo_w, shift_OFFSET, 0);
  }
  

  // if(RESTART && STARTTIME > 0 ){
  //   // write the intial files
  //   if (!WRITEHDF5) {
  //     if (ASCII == 0) {
  //       writetofile_mpi_binary(gridinfo_w, argv, 0 + STARTTIME +1);
  //     } else {
  //       writetofile_mpi(gridinfo_w, argv, 0 + STARTTIME +1);
  //     }
  //   } else {
  //     LBM_ALL_SAVE = LBM ;
  //     writetofile_mpi_hdf5(gridinfo_w, argv, 0 + STARTTIME + 1);
  //     LBM_ALL_SAVE = false ;
  //   }
  // }



  // Section 1 done!
  
   
//  exit(0);


/**
 * Section 2
 * Smoothing
 */

  if(!RESTART){
    if(taskid == MASTER){
      printf("\n- [%d] Smooth start", taskid) ;
    }
    

    for(t=1; t<nsmooth; t++) {

      // perform smoothing
      smooth(workers_mpi.start, workers_mpi.end);

      // apply boundary conditions
      if(boundary_worker) {
        apply_boundary_conditions(taskid);
      }

      // sync the workers data
      mpiexchange_left_right(taskid);
      mpiexchange_top_bottom(taskid);
      if (DIMENSION == 3) {
        mpiexchange_front_back(taskid);
      }

    } // end of smooth loop

    if (!WRITEHDF5) {
      if (ASCII == 0) {
        writetofile_mpi_binary(gridinfo_w, argv, 0 + STARTTIME);
      } else {
        writetofile_mpi(gridinfo_w, argv, 0 + STARTTIME);
      }
    } else {
      writetofile_mpi_hdf5(gridinfo_w, argv, 0 + STARTTIME);
      if (taskid == MASTER)
      {
        ms_write_xmf(argv[3], 0+STARTTIME);  
      }
    }

    if(taskid == MASTER){printf("\n- [%d] Smooth complete", taskid) ;}
  }else{
    if(taskid == MASTER){printf("\n- [%d] No Smoothing", taskid) ;}
  }

   


/**
 * Section 3
 * Solverloop
 */
  if(taskid == MASTER){printf("\n- [%d] Solverloop start", taskid) ;}

  for(t=1; t<=ntimesteps; t++){

    // sync the fields in workers
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    if (DIMENSION == 3) {
      mpiexchange_front_back(taskid);
    }
    
    // solve for phase field
    solverloop_phasefield(workers_mpi.start, workers_mpi.end);
    
    // apply boundary conditions
    if(boundary_worker) {
      apply_boundary_conditions(taskid);
    }

    // sync
    mpiexchange_left_right(taskid);
    mpiexchange_top_bottom(taskid);
    if (DIMENSION == 3) {
      mpiexchange_front_back(taskid);
    }
    
    // if Lattice boltzman is enabled
    // solve for fluid effects
    if (LBM) {
      // calculate the f distribution functions
      collision_lbm(workers_mpi.start, workers_mpi.end);
      
      // sync
      mpiexchange_top_bottom_lbm(taskid);
      mpiexchange_left_right_lbm(taskid);
      if (DIMENSION == 3) { 
        mpiexchange_front_back_lbm(taskid);
      }
      // apply boundary conditions
      if (boundary_worker) { 
        apply_boundary_conditions_lbm(taskid); 
      }
      
      // perform streaming operations
      streaming_row_forward(workers_mpi.start, workers_mpi.end);
      streaming_row_back(workers_mpi.start, workers_mpi.end);
      streaming_col_forward(workers_mpi.start, workers_mpi.end);
      streaming_col_back(workers_mpi.start, workers_mpi.end);
      if(DIMENSION ==3){
        streaming_Z_forward(workers_mpi.start, workers_mpi.end);
        streaming_Z_back(workers_mpi.start, workers_mpi.end);
      }
      
      // calculate velocities and rho
      field_variables_lbm(workers_mpi.start, workers_mpi.end);
    }
    
    if ((FUNCTION_F != 5) && (!GRAIN_GROWTH)) {
      solverloop_concentration(workers_mpi.start,workers_mpi.end);
    }
    
    if (TEMPGRADY) {
      BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
      GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
      temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT + (workers_mpi.offset[Y]-workers_mpi.offset_y)*GRADIENT;
      apply_temperature_gradientY(gridinfo_w, shift_OFFSET, t);
    }
    
    if (ELASTICITY) {
      for(iter=1; iter < MAX_ITERATIONS; iter++) {		//elasticity solver
		     mpiexchange_top_bottom_stress(taskid);
		     mpiexchange_left_right_stress(taskid);
         if (DIMENSION ==3) {
           mpiexchange_front_back_stress(taskid);
         }
		     iterative_stress_solver(workers_mpi.start, workers_mpi.end);
         calculate_avg_strain();
         
//          printf("avg_strain.xx=%le\n, avg_strain.yy=%le\n", avg_strain.xx, avg_strain.yy);
         
         if (boundary_worker) {
		        apply_boundary_conditions_stress(taskid);
         }
		     
// 		     if ((iter%100)==0) {
//            error = 0.0;
// 		       for(x=workers_mpi.start[X]; x<=workers_mpi.end[X]; x++) {
// 		         compute_error(x, &error);
// 		       }
// 		       printf("error=%le\n", error);
// 		       MPI_Reduce(&error,  &global_error,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
// 		       if (fabs(global_error) < tolerance) {
// 		        break;
// 		       }
//          }
      }
    }
    
    // Check for shift
    if (t%10 == 0)
    {
      if(SHIFT)
      {
        MPI_Iallreduce(&workers_max_min.INTERFACE_POS_MAX,  &INTERFACE_POS_GLOBAL,  1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD, &request);
        shift_ON = 0;
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        if(INTERFACE_POS_GLOBAL > shiftj)
        {
          shift_ON = 1;
        }
        if (shift_ON)
        {
          EIDT.interface_velocity = (INTERFACE_POS_GLOBAL - shiftj)*deltay / (EIDT.shift_interval*deltat);
          apply_shiftY(gridinfo_w, INTERFACE_POS_GLOBAL); 
          shift_OFFSET += (INTERFACE_POS_GLOBAL - shiftj);
          mpiexchange_top_bottom(taskid);

          if(taskid == rand() % numtasks)
          {
            printf("\n- step[%d] shift: y=%ld v=%le d=%ld int=%le", t, INTERFACE_POS_GLOBAL, EIDT.interface_velocity, (INTERFACE_POS_GLOBAL - shiftj), EIDT.shift_interval*deltat);
          }

          EIDT.shift_interval = 0 ;
        }
      }
    }

    EIDT.shift_interval += 1;

    if(boundary_worker) { 
      apply_boundary_conditions(taskid);
      if(LBM){   // safety net to accounnt for the possible non local interactions in streaming in lbm fields
        apply_boundary_conditions_lbm(taskid);
      }
    }
    

    if (t%saveT == 0) {

      if (USE_NEW_INPFILE) {
        ms_ReadInputParameters(argv[1]);
      }
      else {
        if(msDyInp_IsModified()){
          msDyInp_ReInitializeInpParameters() ;
        }
      }
      
      if (!WRITEHDF5) {
        if (ASCII == 0) {
          writetofile_mpi_binary(gridinfo_w, argv, t + STARTTIME);
        } else {
          writetofile_mpi(gridinfo_w, argv, t + STARTTIME);
        }
      } else {
        LBM_ALL_SAVE = LBM && (t % lbmSaveFreq == 0) ;

        writetofile_mpi_hdf5(gridinfo_w, argv, t + STARTTIME);

        if (taskid == MASTER)
        {
          ms_write_xmf(argv[3], t+STARTTIME);  
        }
        
        LBM_ALL_SAVE = false ;
      }

      if(SHIFT) {
        if (taskid == MASTER) {
          fp=fopen("DATA/shift.dat","a");
          fprintf(fp,"%ld %ld\n",t + STARTTIME, shift_OFFSET + shift_position);
          fclose(fp);
        }
      }
    }

    /// update the boiundary vaklues after shifting
    if(EIDT_SWITCH)
    {
      if (EIDT.mode == 1)
      {
        for(int i = 0 ; i < NUMCOMPONENTS-1; i++)
        {
          EIDT.comp_far_field[i] += EIDT.comp_rate[i] * deltat;
        }        
      }
      else if (EIDT.mode == 2)
      {
        EIDT_Update_FFComp(STARTTIME + t);
      }
    }


    //printf("Iteration=%d\n",t);
    if (t%time_output == 0 || t == 1)
    { 
      for (b=0; b<NUMPHASES; b++) {
        MPI_Reduce(&workers_max_min.phi_max[b],        &global_max_min.phi_max[b],          1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.phi_min[b],        &global_max_min.phi_min[b],          1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rel_change_phi[b], &global_max_min.rel_change_phi[b],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rel_change_phi[b],&global_max_min.rel_change_phi_max[b],1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      }
      for (k=0; k<(NUMCOMPONENTS-1); k++) {
        MPI_Reduce(&workers_max_min.mu_max[k],         &global_max_min.mu_max[k],          1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.mu_min[k],         &global_max_min.mu_min[k],          1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rel_change_mu[k],  &global_max_min.rel_change_mu[k],   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rel_change_mu[k],&global_max_min.rel_change_mu_max[k], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      }
      if(LBM) {
        MPI_Reduce(&workers_max_min.rho_max,         &global_max_min.rho_max,          1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&workers_max_min.rho_min,         &global_max_min.rho_min,          1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        workers_max_min.rho_max = -10 ; workers_max_min.rho_min = 10 ;
      }


      clock_now = clock();
      time_passed = ((double)(clock_now - clock_start)) / CLOCKS_PER_SEC;

      if (taskid == MASTER)
      {
        fprintf(stdout, "\n");
        for (b=0; b<NUMPHASES; b++)
        {
          fprintf(stdout, "\nstep[%ld] (%5.2f) %16s, Max = %le, Min = %le, Relative_Change=%le", 
            t+STARTTIME, time_passed, Phases[b], global_max_min.phi_max[b], global_max_min.phi_min[b], sqrt(global_max_min.rel_change_phi[b]));
        }
        for (k=0; k<NUMCOMPONENTS-1; k++)
        {
          fprintf(stdout, "\nstep[%ld] (%5.2f) %16s, Max = %le, Min = %le, Relative_Change=%le", 
            t+STARTTIME, time_passed, Components[k], global_max_min.mu_max[k], global_max_min.mu_min[k], sqrt(global_max_min.rel_change_mu[k]));
        }

        if(LBM)
        {
          fprintf(stdout, "\nstep[%ld] (%5.2f) %16s, Max = %le, Min = %le",
            t+STARTTIME, time_passed, "rho", global_max_min.rho_max, global_max_min.rho_min);
        }

        if(EIDT_SWITCH)
        {
          print_double_array("comp_ff", EIDT.comp_far_field, NUMCOMPONENTS-1, t+STARTTIME);
          if (EIDT.mode == 1)
          {
            print_double_array("comp_ff_rate", EIDT.comp_rate, NUMCOMPONENTS-1, t+STARTTIME);  
          }else if (EIDT.mode == 2)
          {
            print_double_array("comp_center", EIDT.comp_center, NUMCOMPONENTS-1, t+STARTTIME);  
          }
          
        }
      }
    }
  }// end of solverloop

  if(taskid == MASTER){printf("\n- [%d] Solverloop complete", taskid) ;}



/**
 * Section 3
 * Freeing allocated memory
 */

  if (USE_GSL_MATINV) {
    matinv_gsl_Free(&dcdmu_gsl_data);
  }

  if (USE_NEW_IO_BUFF) {
    IO_FreeGlobal();
  }
  
  // if (USE_NEW_MALLOC) {
  //   ms_FreeAllGlobalVariables(); 
  // }
  // else {
    free_variables();
  // }
  
  if (taskid == MASTER) {
    index_count = layer_size*rows_x;
    for (index_=0; index_ < index_count; index_++) {
      if ((&gridinfo[index_]) !=NULL) {
        free_memory_fields(&gridinfo[index_]);
      }
    }
    free(gridinfo);
    if (LBM) {
      // free(lbm_gridinfo);
      for (index_=0; index_ < index_count; index_++) {
        if ((&lbm_gridinfo[index_]) !=NULL) {
          free_memory_lbm_fields(&lbm_gridinfo[index_]);
        }
      }
      free(lbm_gridinfo);
    }
  }
  // for(i=0; i<size_fields; i++) {
  //   free(coordNames[i]);
  // }

  MPI_Type_free(&MPI_gridinfo_vector_b);
  MPI_Type_free(&MPI_gridinfo_vector_c);
  MPI_Type_free(&MPI_gridinfo);

  if (ELASTICITY) {
    MPI_Type_free(&MPI_gridinfo_vector_b_stress);
    MPI_Type_free(&MPI_gridinfo_vector_c_stress);
    MPI_Type_free(&MPI_iter_gridinfo);
  }
  if (LBM) {
    MPI_Type_free(&MPI_gridinfo_vector_b_lbm);
    MPI_Type_free(&MPI_gridinfo_vector_c_lbm);
    MPI_Type_free(&MPI_lbm_gridinfo);
  }

  // Finalize MPI
  MPI_Finalize();

  if(taskid == MASTER){printf("\n\n"); } ;

  return 0 ;
}



