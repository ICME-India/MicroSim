void free_variables(){

  ret = clFlush(cmdQ);
  ret = clFinish(cmdQ);
  ret = clReleaseKernel(ker_copy_New_To_Old);
  for ( i = 0; i < 3; i++ ) {
    ret = clReleaseKernel(kernel10[i]);
  }
  for ( i = 0; i < 2; i++ ) {
    ret = clReleaseKernel(kernel11[i]);
  }
  ret = clReleaseKernel(kernel12);
  ret = clReleaseProgram(program);    
  //ret = clReleaseMemObject(d_gridOld);
  //ret = clReleaseMemObject(d_gridNew);
  ret = clReleaseMemObject(d_gridinfomN);
  ret = clReleaseMemObject(d_gridinfomO);
  ret = clReleaseMemObject(d_cscl);
  ret = clReleaseMemObject(d_pfmdat);
  ret = clReleaseMemObject(d_pfmvar);
  ret = clReleaseMemObject(d_temp);
  ret = clReleaseMemObject(d_tstep);
  ret = clReleaseCommandQueue(cmdQ);
  ret = clReleaseContext(context);

  free(source_str); 

  //free(gridOld);
  //free(gridNew);
  free(cscl);
  free(temp);
  free(tstep);


  long index, gidy, i;  
  long index_count;
  long layer;
  
  index_count = layer_size*rows_x;
  
  Free3M(Diffusivity, NUMPHASES, NUMCOMPONENTS-1);
  Free3M(ceq,         NUMPHASES, NUMPHASES);
  Free3M(cfill,       NUMPHASES, NUMPHASES);
 
  
  for (index=0; index < index_count; index++) {
    free_memory_fields(&gridinfo[index]);
  }
  for (i=0; i<6; i++) {
    free(boundary[i]);
  }


  free(global_max_min.phi_max);
  free(global_max_min.phi_min);
  free(global_max_min.com_max);
  free(global_max_min.com_min);
  free(global_max_min.rel_change_phi);
  free(global_max_min.rel_change_com);

}
