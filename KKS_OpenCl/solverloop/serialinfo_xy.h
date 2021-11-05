void serialinfo_xy() {
  
  long a;
  long k;
  long index_count;
  long index;
  start[X]   = 1;
  start[Y]   = 1;
  start[Z]   = 1;
  
  rows_x     = MESH_X + 2;
  rows_y     = MESH_Y + 2;
  rows_z     = MESH_Z + 2;
  end[X]     = rows_x - 2;
  end[Y]     = rows_y - 2;
  end[Z]     = rows_z - 2;
  
  layer_size = rows_y*rows_z;
  
  if (DIMENSION == 2) {
    rows_z     = 1;
    start[Z]   = 0; 
    end[Z]     = 0;
    layer_size = rows_y;
  }
  
  index_count = layer_size*rows_x; 
  
  gridinfo = (struct fields* )malloc((index_count)*sizeof(*gridinfo));
  
  for (index=0; index < index_count; index++) {
    allocate_memory_fields(&gridinfo[index]);
  }
  
}