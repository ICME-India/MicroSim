void CL_create_kernel_args() {

  printf("In CL_create_kernel_args\n");

  ker_SolverStress_iterative = clCreateKernel(program, "SolverStress_iterative", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverStress_iterative error %d in rank = %d \n", ret, rank);
    exit(1);
  }

  ker_SolverStress_iterative_2D = clCreateKernel(program, "SolverStress_iterative_2D", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverStress_iterative_2D error %d in rank = %d \n", ret, rank);
    exit(1);
  }
  
  ker_SolverCsClEq_F2 = clCreateKernel(program, "SolverCsClEq_F2", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCsClEq_F2 error %d in rank = %d \n", ret, rank);
    exit(1);
  }

  ker_SolverCsClEq_F3 = clCreateKernel(program, "SolverCsClEq_F3", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCsClEq_F3 error %d in rank = %d \n", ret, rank);
    exit(1);
  }

  ker_SolverCsClEq_F4 = clCreateKernel(program, "SolverCsClEq_F4", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCsClEq_F4 error %d in rank = %d \n", ret, rank);
    exit(1);
  }

  ker_SolverPhi_F2_smooth = clCreateKernel(program, "SolverPhi_F2_smooth", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F2_smooth error %d\n", ret);
  }

  ker_SolverPhi_F3_smooth = clCreateKernel(program, "SolverPhi_F3_smooth", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F3_smooth error %d\n", ret);
  }

  ker_SolverPhi_F4_smooth = clCreateKernel(program, "SolverPhi_F4_smooth", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F4_smooth error %d\n", ret);
  }
  
  ker_SolverPhi_F2 = clCreateKernel(program, "SolverPhi_F2", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F2 error %d\n", ret);
  }

  ker_SolverPhi_F3 = clCreateKernel(program, "SolverPhi_F3", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F3 error %d\n", ret);
  }

  ker_SolverPhi_F4 = clCreateKernel(program, "SolverPhi_F4", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverPhi_F4 error %d\n", ret);
  }

  ker_SolverCatr_F2_smooth = clCreateKernel(program, "SolverCatr_F2_smooth", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F2_smooth error %d\n", ret);
  }

  ker_SolverCatr_F3_smooth = clCreateKernel(program, "SolverCatr_F3_smooth", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F3_smooth error %d\n", ret);
  }

  ker_SolverCatr_F4_smooth = clCreateKernel(program, "SolverCatr_F4_smooth", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F4_smooth error %d\n", ret);
  }
  
  ker_SolverCatr_F2 = clCreateKernel(program, "SolverCatr_F2", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F2 error %d\n", ret);
  }
  
  ker_SolverCatr_F3 = clCreateKernel(program, "SolverCatr_F3", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F3 error %d\n", ret);
  }

  ker_SolverCatr_F4 = clCreateKernel(program, "SolverCatr_F4", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_SolverCatr_F4 error %d\n", ret);
  }
  
  ker_apply_BC_phi_y0_noflux = clCreateKernel(program, "apply_BC_phi_y0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_y0_noflux error %d\n", ret);
  }
  
  ker_apply_BC_phi_yn_noflux = clCreateKernel(program, "apply_BC_phi_yn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_yn_noflux error %d\n", ret);
  }
  
  ker_apply_BC_phi_y0_periodic = clCreateKernel(program, "apply_BC_phi_y0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_y0_periodic error %d\n", ret);
  }
  
  ker_apply_BC_phi_yn_periodic = clCreateKernel(program, "apply_BC_phi_yn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_yn_periodic error %d\n", ret);
  }
  
  ker_apply_BC_phi_z0_noflux = clCreateKernel(program, "apply_BC_phi_z0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_z0_noflux error %d\n", ret);
  }
  
  ker_apply_BC_phi_zn_noflux = clCreateKernel(program, "apply_BC_phi_zn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_zn_noflux error %d\n", ret);
  }
  
  ker_apply_BC_phi_z0_periodic = clCreateKernel(program, "apply_BC_phi_z0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_z0_periodic error %d\n", ret);
  }
  
  ker_apply_BC_phi_zn_periodic = clCreateKernel(program, "apply_BC_phi_zn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_zn_periodic error %d\n", ret);
  }
  
  ker_apply_BC_phi_x0_noflux = clCreateKernel(program, "apply_BC_phi_x0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_x0_noflux error %d\n", ret);
  }
  
  ker_apply_BC_phi_xn_noflux = clCreateKernel(program, "apply_BC_phi_xn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_xn_noflux error %d\n", ret);
  }
  
  ker_apply_BC_phi_x0_periodic = clCreateKernel(program, "apply_BC_phi_x0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_x0_periodic error %d\n", ret);
  }
  
  ker_apply_BC_phi_xn_periodic = clCreateKernel(program, "apply_BC_phi_xn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_phi_xn_periodic error %d\n", ret);
  }
  
  ker_apply_BC_com_y0_noflux = clCreateKernel(program, "apply_BC_com_y0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_y0_noflux error %d\n", ret);
  }
  
  ker_apply_BC_com_yn_noflux = clCreateKernel(program, "apply_BC_com_yn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_yn_noflux error %d\n", ret);
  }
  
  ker_apply_BC_com_y0_periodic = clCreateKernel(program, "apply_BC_com_y0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_y0_periodic error %d\n", ret);
  }
  
  ker_apply_BC_com_yn_periodic = clCreateKernel(program, "apply_BC_com_yn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_yn_periodic error %d\n", ret);
  }
  
  ker_apply_BC_com_z0_noflux = clCreateKernel(program, "apply_BC_com_z0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_z0_noflux error %d\n", ret);
  }
  
  ker_apply_BC_com_zn_noflux = clCreateKernel(program, "apply_BC_com_zn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_zn_noflux error %d\n", ret);
  }
  
  ker_apply_BC_com_z0_periodic = clCreateKernel(program, "apply_BC_com_z0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_z0_periodic error %d\n", ret);
  }
  
  ker_apply_BC_com_zn_periodic = clCreateKernel(program, "apply_BC_com_zn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_zn_periodic error %d\n", ret);
  }
  
  ker_apply_BC_com_x0_noflux = clCreateKernel(program, "apply_BC_com_x0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_x0_noflux error %d\n", ret);
  }
  
  ker_apply_BC_com_xn_noflux = clCreateKernel(program, "apply_BC_com_xn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_xn_noflux error %d\n", ret);
  }
  
  ker_apply_BC_com_x0_periodic = clCreateKernel(program, "apply_BC_com_x0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_x0_periodic error %d\n", ret);
  }
  
  ker_apply_BC_com_xn_periodic = clCreateKernel(program, "apply_BC_com_xn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_com_xn_periodic error %d\n", ret);
  }
  
  ker_addNoise = clCreateKernel(program, "addNoise", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_addNoise error %d\n", ret);
  }

  ker_copy_New_To_Old = clCreateKernel(program, "copy_New_To_Old", &ret);
  if (ret!=CL_SUCCESS) {
    printf("ker_copy_New_To_Old error %d\n", ret);
  }
  kernel10[0] = clCreateKernel(program, "update_temp_UC", &ret);
  if (ret!=CL_SUCCESS) {
    printf("kernel10 update_temp_UC error %d\n", ret);
  }
  kernel10[1] = clCreateKernel(program, "update_temp_DS", &ret);
  if (ret!=CL_SUCCESS) {
    printf("kernel10 update_temp_DS error %d\n", ret);
  }
  kernel10[2] = clCreateKernel(program, "update_temp_CR", &ret);
  if (ret!=CL_SUCCESS) {
    printf("kernel10 update_temp_CR error %d\n", ret);
  }
  kernel11[0] = clCreateKernel(program, "apply_BC_temp_it_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("kernel11 apply_BC_temp_it_noflux error %d\n", ret);
  }
  kernel11[1] = clCreateKernel(program, "apply_BC_temp_ib_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("kernel11 apply_BC_temp_ib_noflux error %d\n", ret);
  }



  ker_apply_BC_ela_y0_noflux = clCreateKernel(program, "apply_BC_ela_y0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_y0_noflux error %d\n", ret);
  }

  ker_apply_BC_ela_yn_noflux = clCreateKernel(program, "apply_BC_ela_yn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_yn_noflux error %d\n", ret);
  }

  ker_apply_BC_ela_y0_periodic = clCreateKernel(program, "apply_BC_ela_y0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_y0_periodic error %d\n", ret);
  }

  ker_apply_BC_ela_yn_periodic = clCreateKernel(program, "apply_BC_ela_yn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_yn_periodic error %d\n", ret);
  }

  ker_apply_BC_ela_z0_noflux = clCreateKernel(program, "apply_BC_ela_z0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_z0_noflux error %d\n", ret);
  }

  ker_apply_BC_ela_zn_noflux = clCreateKernel(program, "apply_BC_ela_zn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_zn_noflux error %d\n", ret);
  }

  ker_apply_BC_ela_z0_periodic = clCreateKernel(program, "apply_BC_ela_z0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_z0_periodic error %d\n", ret);
  }

  ker_apply_BC_ela_zn_periodic = clCreateKernel(program, "apply_BC_ela_zn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_zn_periodic error %d\n", ret);
  }

  ker_apply_BC_ela_x0_noflux = clCreateKernel(program, "apply_BC_ela_x0_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_x0_noflux error %d\n", ret);
  }

  ker_apply_BC_ela_xn_noflux = clCreateKernel(program, "apply_BC_ela_xn_noflux", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_xn_noflux error %d\n", ret);
  }

  ker_apply_BC_ela_x0_periodic = clCreateKernel(program, "apply_BC_ela_x0_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_x0_periodic error %d\n", ret);
  }

  ker_apply_BC_ela_xn_periodic = clCreateKernel(program, "apply_BC_ela_xn_periodic", &ret);
  if (ret!=CL_SUCCESS) {
    printf("apply_BC_ela_xn_periodic error %d\n", ret);
  }

  
  /* Set OpenCL kernel arguments */

  ret = clSetKernelArg(ker_SolverStress_iterative, 0,  sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverStress_iterative, 1,  sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_SolverStress_iterative, 2,  sizeof(cl_mem), &d_eigen_strain_phase);
  ret = clSetKernelArg(ker_SolverStress_iterative, 3,  sizeof(cl_mem), &d_stiffness_phase);
  ret = clSetKernelArg(ker_SolverStress_iterative, 4,  sizeof(cl_mem), &d_stiffness_phase_n);
  ret = clSetKernelArg(ker_SolverStress_iterative, 5,  sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverStress_iterative, 6,  sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverStress_iterative, 7,  sizeof(cl_mem), &d_tstep);

  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 0,  sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 1,  sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 2,  sizeof(cl_mem), &d_eigen_strain_phase);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 3,  sizeof(cl_mem), &d_stiffness_phase);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 4,  sizeof(cl_mem), &d_stiffness_phase_n);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 5,  sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 6,  sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverStress_iterative_2D, 7,  sizeof(cl_mem), &d_tstep);
  
  ret = clSetKernelArg(ker_SolverCsClEq_F2, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCsClEq_F2, 1, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCsClEq_F2, 2, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCsClEq_F2, 3, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCsClEq_F2, 4, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCsClEq_F2, 5, sizeof(cl_mem), &d_tstep);

  ret = clSetKernelArg(ker_SolverCsClEq_F3, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCsClEq_F3, 1, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCsClEq_F3, 2, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCsClEq_F3, 3, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCsClEq_F3, 4, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCsClEq_F3, 5, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverCsClEq_F3, 6, sizeof(cl_mem), &d_propf3);

  ret = clSetKernelArg(ker_SolverCsClEq_F4, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 1, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 2, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 3, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 4, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 5, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 6, sizeof(cl_mem), &d_propf4);
  ret = clSetKernelArg(ker_SolverCsClEq_F4, 7, sizeof(cl_mem), &d_propf4spline);

  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverPhi_F2_smooth, 6, sizeof(cl_mem), &d_tstep);

  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverPhi_F3_smooth, 7, sizeof(cl_mem), &d_propf3);

  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 7, sizeof(cl_mem), &d_propf4);
  ret = clSetKernelArg(ker_SolverPhi_F4_smooth, 8, sizeof(cl_mem), &d_propf4spline);
  

  ret = clSetKernelArg(ker_SolverPhi_F2, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverPhi_F2, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverPhi_F2, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverPhi_F2, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverPhi_F2, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverPhi_F2, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverPhi_F2, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverPhi_F2, 7, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_SolverPhi_F2, 8, sizeof(cl_mem), &d_eigen_strain_phase);
  ret = clSetKernelArg(ker_SolverPhi_F2, 9, sizeof(cl_mem), &d_stiffness_phase);
  ret = clSetKernelArg(ker_SolverPhi_F2, 10, sizeof(cl_mem), &d_stiffness_phase_n);

  ret = clSetKernelArg(ker_SolverPhi_F3, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverPhi_F3, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverPhi_F3, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverPhi_F3, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverPhi_F3, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverPhi_F3, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverPhi_F3, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverPhi_F3, 7, sizeof(cl_mem), &d_propf3);
  ret = clSetKernelArg(ker_SolverPhi_F3, 8, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_SolverPhi_F3, 9, sizeof(cl_mem), &d_eigen_strain_phase);
  ret = clSetKernelArg(ker_SolverPhi_F3, 10, sizeof(cl_mem), &d_stiffness_phase);
  ret = clSetKernelArg(ker_SolverPhi_F3, 11, sizeof(cl_mem), &d_stiffness_phase_n);

  ret = clSetKernelArg(ker_SolverPhi_F4, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverPhi_F4, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverPhi_F4, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverPhi_F4, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverPhi_F4, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverPhi_F4, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverPhi_F4, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverPhi_F4, 7, sizeof(cl_mem), &d_propf4);
  ret = clSetKernelArg(ker_SolverPhi_F4, 8, sizeof(cl_mem), &d_propf4spline);
  ret = clSetKernelArg(ker_SolverPhi_F4, 9, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_SolverPhi_F4, 10, sizeof(cl_mem), &d_eigen_strain_phase);
  ret = clSetKernelArg(ker_SolverPhi_F4, 11, sizeof(cl_mem), &d_stiffness_phase);
  ret = clSetKernelArg(ker_SolverPhi_F4, 12, sizeof(cl_mem), &d_stiffness_phase_n);

  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCatr_F2_smooth, 6, sizeof(cl_mem), &d_tstep);

  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverCatr_F3_smooth, 7, sizeof(cl_mem), &d_propf3);

  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 7, sizeof(cl_mem), &d_propf4);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 8, sizeof(cl_mem), &d_propf4spline);
  ret = clSetKernelArg(ker_SolverCatr_F4_smooth, 9, sizeof(cl_mem), &d_propf4spline1);

  ret = clSetKernelArg(ker_SolverCatr_F2, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCatr_F2, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverCatr_F2, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCatr_F2, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCatr_F2, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCatr_F2, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCatr_F2, 6, sizeof(cl_mem), &d_tstep);

  ret = clSetKernelArg(ker_SolverCatr_F3, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCatr_F3, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverCatr_F3, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCatr_F3, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCatr_F3, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCatr_F3, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCatr_F3, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverCatr_F3, 7, sizeof(cl_mem), &d_propf3);

  ret = clSetKernelArg(ker_SolverCatr_F4, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_SolverCatr_F4, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_SolverCatr_F4, 2, sizeof(cl_mem), &d_cscl);
  ret = clSetKernelArg(ker_SolverCatr_F4, 3, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_SolverCatr_F4, 4, sizeof(cl_mem), &d_pfmvar);
  ret = clSetKernelArg(ker_SolverCatr_F4, 5, sizeof(cl_mem), &d_temp);
  ret = clSetKernelArg(ker_SolverCatr_F4, 6, sizeof(cl_mem), &d_tstep);
  ret = clSetKernelArg(ker_SolverCatr_F4, 7, sizeof(cl_mem), &d_propf4);
  ret = clSetKernelArg(ker_SolverCatr_F4, 8, sizeof(cl_mem), &d_propf4spline);
  ret = clSetKernelArg(ker_SolverCatr_F4, 9, sizeof(cl_mem), &d_propf4spline1);
  
  ret = clSetKernelArg(ker_apply_BC_phi_y0_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_y0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_y0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_yn_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_yn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_yn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_y0_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_y0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_y0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_yn_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_yn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_yn_periodic, 2, sizeof(cl_mem), &d_cscl);
  
  ret = clSetKernelArg(ker_apply_BC_phi_z0_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_z0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_z0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_zn_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_zn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_zn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_z0_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_z0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_z0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_zn_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_zn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_zn_periodic, 2, sizeof(cl_mem), &d_cscl);
  
  ret = clSetKernelArg(ker_apply_BC_phi_x0_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_x0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_x0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_xn_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_xn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_xn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_x0_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_x0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_x0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_phi_xn_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_phi_xn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_phi_xn_periodic, 2, sizeof(cl_mem), &d_cscl);
//

  ret = clSetKernelArg(ker_apply_BC_com_y0_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_y0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_y0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_yn_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_yn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_yn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_y0_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_y0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_y0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_yn_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_yn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_yn_periodic, 2, sizeof(cl_mem), &d_cscl);
  
  ret = clSetKernelArg(ker_apply_BC_com_z0_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_z0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_z0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_zn_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_zn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_zn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_z0_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_z0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_z0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_zn_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_zn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_zn_periodic, 2, sizeof(cl_mem), &d_cscl);
  
  ret = clSetKernelArg(ker_apply_BC_com_x0_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_x0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_x0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_xn_noflux, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_xn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_xn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_x0_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_x0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_x0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_com_xn_periodic, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_apply_BC_com_xn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_com_xn_periodic, 2, sizeof(cl_mem), &d_cscl);
  

  ret = clSetKernelArg(ker_addNoise, 0, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_addNoise, 1, sizeof(cl_mem), &d_pfmdat);

  ret = clSetKernelArg(ker_copy_New_To_Old, 0, sizeof(cl_mem), &d_gridinfomO);
  ret = clSetKernelArg(ker_copy_New_To_Old, 1, sizeof(cl_mem), &d_gridinfomN);
  ret = clSetKernelArg(ker_copy_New_To_Old, 2, sizeof(cl_mem), &d_pfmdat);

  for ( i = 0; i < 3; i++ ) { 
    ret = clSetKernelArg(kernel10[i], 0, sizeof(cl_mem), &d_gridinfomN);
    ret = clSetKernelArg(kernel10[i], 1, sizeof(cl_mem), &d_tstep);
    ret = clSetKernelArg(kernel10[i], 2, sizeof(cl_mem), &d_pfmdat);
    ret = clSetKernelArg(kernel10[i], 3, sizeof(cl_mem), &d_pfmvar);
  }

  for ( i = 0; i < 2; i++ ) {
    ret = clSetKernelArg(kernel11[i], 0, sizeof(cl_mem), &d_gridinfomN);
    ret = clSetKernelArg(kernel11[i], 1, sizeof(cl_mem), &d_tstep);
    ret = clSetKernelArg(kernel11[i], 2, sizeof(cl_mem), &d_pfmdat);
    ret = clSetKernelArg(kernel11[i], 3, sizeof(cl_mem), &d_pfmvar);
  }



  ret = clSetKernelArg(ker_apply_BC_ela_y0_noflux, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_y0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_y0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_yn_noflux, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_yn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_yn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_y0_periodic, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_y0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_y0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_yn_periodic, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_yn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_yn_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_z0_noflux, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_z0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_z0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_zn_noflux, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_zn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_zn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_z0_periodic, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_z0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_z0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_zn_periodic, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_zn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_zn_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_x0_noflux, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_x0_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_x0_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_xn_noflux, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_xn_noflux, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_xn_noflux, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_x0_periodic, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_x0_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_x0_periodic, 2, sizeof(cl_mem), &d_cscl);

  ret = clSetKernelArg(ker_apply_BC_ela_xn_periodic, 0, sizeof(cl_mem), &d_iter_gridinfom);
  ret = clSetKernelArg(ker_apply_BC_ela_xn_periodic, 1, sizeof(cl_mem), &d_pfmdat);
  ret = clSetKernelArg(ker_apply_BC_ela_xn_periodic, 2, sizeof(cl_mem), &d_cscl);

  printf("Exit CL_create_kernel_args\n");



}
