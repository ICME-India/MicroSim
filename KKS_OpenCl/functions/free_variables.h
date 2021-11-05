void free_variables(){

  ret = clFlush(cmdQ);
  ret = clFinish(cmdQ);
  ret = clReleaseKernel(kernel1);
  ret = clReleaseKernel(kernel2);
  ret = clReleaseKernel(kernel3);
  ret = clReleaseKernel(kernel4);
  for ( i = 0; i < 9; i++ ) {
    ret = clReleaseKernel(kernel5[i]);
    ret = clReleaseKernel(kernel6[i]);
    ret = clReleaseKernel(kernel7[i]);
  }
  ret = clReleaseKernel(kernel8);
  ret = clReleaseKernel(kernel9);
  for ( i = 0; i < 3; i++ ) {
    ret = clReleaseKernel(kernel10[i]);
  }
  for ( i = 0; i < 2; i++ ) {
    ret = clReleaseKernel(kernel11[i]);
  }
  ret = clReleaseKernel(kernel12);
  ret = clReleaseProgram(program);    
  ret = clReleaseMemObject(d_gridOld);
  ret = clReleaseMemObject(d_gridNew);
  ret = clReleaseMemObject(d_cscl);
  ret = clReleaseMemObject(d_pfmdat);
  ret = clReleaseMemObject(d_pfmvar);
  ret = clReleaseMemObject(d_temp);
  ret = clReleaseMemObject(d_tstep);
  ret = clReleaseCommandQueue(cmdQ);
  ret = clReleaseContext(context);

  free(source_str); 

  free(gridOld);
  free(gridNew);
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


}
