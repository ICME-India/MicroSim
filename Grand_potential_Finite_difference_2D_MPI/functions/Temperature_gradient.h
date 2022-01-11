#ifndef TEMPERATURE_GRADIENT_H_
#define TEMPERATURE_GRADIENT_H_

void apply_temperature_gradientY(struct fields* gridinfo_w, long shift_OFFSET, long t) {
  long x, y, z, endx;
  long gidy;
//   double BASE_POS=0;
//   double GRADIENT;
//   double temp_bottom;
//    
//   BASE_POS    = (temperature_gradientY.gradient_OFFSET/deltay) - shift_OFFSET + ((temperature_gradientY.velocity/deltay)*(t*deltat));
//   GRADIENT    = (temperature_gradientY.DeltaT)*deltay/(temperature_gradientY.Distance);
//   temp_bottom = temperature_gradientY.base_temp - BASE_POS*GRADIENT + (workers_mpi.offset[Y]-workers_mpi.offset_y)*GRADIENT;
  
  for(x=0; x < workers_mpi.rows_x; x++) {
   for(z=0; z < workers_mpi.rows_z; z++) {
    for(y=0; y < workers_mpi.rows_y; y++) {
      gidy                          = x*workers_mpi.layer_size  + z*workers_mpi.rows_y + y;
      gridinfo_w[gidy].temperature  = temp_bottom               + y*GRADIENT;
//       gridinfo_w[gidy].temperature  = temperature_gradientY.base_temp;
//       if (gridinfo1[gidy].temperature > 1.0) {
//          gridinfo1[gidy].temperature = 1.0;
//       } else if (gridinfo1[gidy].temperature < temperature_gradientY.base_temp) {
//          gridinfo[gidy].temperature = temperature_gradientY.base_temp;
//       }
      }
    }
  }
}
#endif