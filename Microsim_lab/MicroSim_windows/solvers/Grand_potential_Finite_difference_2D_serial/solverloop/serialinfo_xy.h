void serialinfo_xy() {
  
  long a;
  long k;
  long index_count;
  long index;
  start[X]   = 3;
  start[Y]   = 3;
  start[Z]   = 3;
  
  rows_x     = MESH_X + 6;
  rows_y     = MESH_Y + 6;
  rows_z     = MESH_Z + 6;
  end[X]     = rows_x - 4;
  end[Y]     = rows_y - 4;
  end[Z]     = rows_z - 4;
  
  layer_size = rows_y*rows_z;
  
  if (DIMENSION == 2) {
    rows_z     = 1;
    start[Z]   = 0; 
    end[Z]     = 0;
    layer_size = rows_y;
  }
  
  index_count = layer_size*rows_x;  
  
//   grad  = (struct gradlayer*)malloc(sizeof(*grad));
//   allocate_memory_gradlayer(grad);
  
  int layer;
//   gradient1 = (struct gradlayer **)malloc((layer_size)*(sizeof(**gradient1)));
  for (layer=0; layer < 4; layer++) {
    gradient1[layer] = (struct gradlayer *)malloc((layer_size)*(sizeof(*gradient1[layer])));
    for (index=0; index < layer_size; index++) {
      allocate_memory_gradlayer(&gradient1[layer][index]);
    }
  }
  
  gradient = gradient1+1;
  
  gridinfo_instance = (struct fields *)malloc(sizeof(*gridinfo_instance));
  allocate_memory_fields(gridinfo_instance);
  
  gridinfo = (struct fields* )malloc((index_count)*sizeof(*gridinfo));
  
  for (index=0; index < index_count; index++) {
    allocate_memory_fields(&gridinfo[index]);
  }
  
  global_max_min.phi_max        = (double*)malloc(NUMPHASES*sizeof(double));
  global_max_min.phi_min        = (double*)malloc(NUMPHASES*sizeof(double));
  global_max_min.mu_max         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  global_max_min.mu_min         = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  global_max_min.rel_change_phi = (double*)malloc((NUMPHASES)*sizeof(double));
  global_max_min.rel_change_mu  = (double*)malloc((NUMCOMPONENTS-1)*sizeof(double));
  
  for (a=0; a<NUMPHASES; a++) {
    global_max_min.phi_max[a] = 1.0;
    global_max_min.phi_min[a] = 0.0;
  }
  for (k=0; k<NUMCOMPONENTS-1; k++) {
    global_max_min.mu_max[k] = 1.0;
    global_max_min.mu_min[k] = 1.0;
  }
}