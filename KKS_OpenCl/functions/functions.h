#ifndef FUNCTIONS_DEFS_H_
#define FUNCTIONS_DEFS_H_

double* MallocV(long m);
double** MallocM(long m,long n);
double*** Malloc3M(long m, long n, long k);
double**** Malloc4M(long m, long n, long k, long l);
void FreeM(double **a, long m);
void Free3M(double ***a, long m, long n);
void Free4M(double ****Mat, long m, long n, long k);
void populate_matrix(double **Mat, char *tmpstr, long NUMPHASES);
void populate_matrix3M(double ***Mat, char *tmpstr, long NUMPHASES);
void populate_diffusivity_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_A_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_thermodynamic_matrix(double ***Mat, char *tmpstr, long NUMCOMPONENTS);
void populate_string_array(char**string, char *tmpstr, long size);
void get_Rotation_Matrix(double **R, double theta, int axis);
void populate_rotation_matrix(double ****Mat, double ****Inv_Rotation_matrix, char *tmpstr);
void initialize_boundary_conditions(char *tmpstr);
void initialize_boundary_points_values(char *tmpstr);
void assign_boundary_points_values(long j, double VALUE_LEFT, double VALUE_RIGHT, double VALUE_FRONT, double VALUE_BACK, double VALUE_TOP, double VALUE_BOTTOM);
void assign_buffer_points_conditions(long j, int BOUNDARY_LEFT, int BOUNDARY_RIGHT, int BOUNDARY_FRONT, int BOUNDARY_BACK, int BOUNDARY_TOP, int BOUNDARY_BOTTOM);
void PRINT_INT(char *key, int value, FILE *fp);
void PRINT_LONG(char *key, long value, FILE *fp);
void PRINT_DOUBLE(char *key, double value, FILE *fp);
void PRINT_MATRIX(char *key, double **Mat, long m, long n, FILE *fp);
void PRINT_VECTOR(char *key, double *Mat, long m, FILE *fp);
void PRINT_STRING_ARRAY(char *key, char **str, long m, FILE *fp);
void PRINT_STRING(char *key, char *str, FILE *fp);
void PRINT_BOUNDARY_CONDITIONS(FILE *fp);
void allocate_memory_fields(struct fields *ptr);
void free_memory_fields(struct fields *ptr);
double hphi(double *phi, long a);
double dhphi(double *phi, long b, long a);
void fill_phase_cube(struct fill_cube fill_cube_parameters, struct fields* gridinfo, long b);
void fill_phase_cylinder (struct fill_cylinder fill_cylinder_parameters, struct fields* gridinfo, long b);
void fill_phase_sphere(struct fill_sphere fill_sphere_parameters, struct fields* gridinfo, long b);
void fill_phase_ellipse (struct fill_ellipse fill_ellipse_parameters, struct fields* gridinfo, long b);
//void init_propertymatrices();
void fill_composition_cube(struct fields* gridinfo);
void reading_input_parameters(char *argv[]);
void free_variables();
void read_cells_vtk_2D(FILE *fp, struct fields *gridinfo);
void readfromfile_serial2D(struct fields* gridinfo, char *argv[], long t);
void read_cells_vtk_2D_binary(FILE *fp, struct fields *gridinfo);
void readfromfile_serial2D_binary(struct fields* gridinfo, char *argv[], long t);
void apply_temperature_gradientY(struct fields* gridinfo, long shift_OFFSET, long t);
void apply_shiftY_cscl(struct csle *cscl, long INTERFACE_POS_GLOBAL);
void mpi_xy(int rank, char *argv[]);

void read_cells_vtk_mpi(FILE *fp, struct fields *gridinfo);
void readfromfile_serialmpi(struct fields* gridinfo, char *argv[], long t);
void read_cells_vtk_mpi_binary(FILE *fp, struct fields *gridinfo);
void readfromfile_serialmpi_binary(struct fields* gridinfo, char *argv[], long t);

#endif






// // #include"constants1.h"
// double chempot(double *phi, double *c, int i);
// double freeenergy(double *phi, double *c, int i);
// double grandchempot(double *phi, double *c, int i);
// double D(struct variables*, double T, long index, long k, long l);
// double df_liquid(double *c, int i);
// double f_liquid (double *c, int i);
// 
// // double free_energy(double *c, double T, long a);
// // double function_A(double T, long i, long j, long a);
// // double function_B(double T, long i, long a);
// // double function_C(double T, long a);
// // double Mu(double *c,double T, long a, long i);
// // // double free_liquid(double *c);
// // double c_mu(double *mu, double T, long a, long i);
// // double dc_dmu(double *mu, double T, long a, long i, long j);
// 
// // double c_l(double *mu, int i);
// double df_phi(double *phi, double *c, int i);
// 
// double chempot_phi(double *phi, double *c, int i);
// double chempot_struct(struct variables* gridinfo, long index, double *mu);
// 
// // double dpsi(double *mu, double T, double *phi, long a);
// 
// double chempot_c(double *phi, double *c, int i);
// void writetofile(double *phi,double *c, int t);
// void writetofile_struct(struct variables* gridinfo,int t);
// void rdfrmfile(double *phi,double *c,int t);
// void rdfrmfile_struct(struct variables* gridinfo,int t);
// void initdomain_struct(struct variables *gridinfo);
// void initdomain(double *phi,double *c);
// void initialize_variables();
// 
// void calculate_d3gradients(long x, struct gradlayer **gradient);
// void calculate_divergence_concentration(long x, struct gradlayer **gradient);
// void calculate_fluxes_concentration(long x, struct gradlayer **gradient);
// void calculate_divergence_concentration_smooth(long x, struct gradlayer **gradient);
// void calculate_divergence_concentration_smooth_concentration(long x, struct gradlayer **gradient);
// void calculate_divergence_phasefield(long x, struct gradlayer **gradient);
// void calculate_divergence_phasefield_smooth(long x, struct gradlayer **gradient);
// void calculate_gradients(long x, struct gradlayer **gradient);
// void calculate_gradients_phasefield(long x, struct gradlayer **gradient);
// void calculate_gradients_concentration(long x, struct gradlayer **gradient);
// 
// // void compute_chemicalpotential(struct variables* gridinfo);
// 
// void swaplayers();
// void solverloop(long *start, long *end);
// void solverloop_phasefield(long *start, long *end);
// void solverloop_concentration(long *start, long *end);
// void  smooth(long *start, long *end);
// void  smooth_concentration(long *start, long *end);
// double FunctionTau(double *phi);
// 
// // void dAdq (double *qab, double* dadq, long a, long b);
// // double function_ac(double *qab, long a, long b);
// 
// 
// 
// // double dAdphi(double *phi, struct gradlayer **gradient, long gidy, long a);
// // double dAdphi_smooth(double *phi, struct gradlayer **gradient, long gidy, long a);
// // double divdAdgradphi(struct gradlayer **gradient, long index, long gidy, long a);
// // double divdAdgradphi_smooth(struct gradlayer **gradient, long index, long gidy, long a);
// 
// double d2gradphi(double *phi, struct gradlayer **gradient, long index, long gidy, long a);
// // void fill_phase_ellipse (double length, double e, double angle, long x_center, long y_center, struct variables* gridinfo, long b);
// // void rotate_vector(double *q, double *rot_q, long a, long b);
// // double dwdphi(double *phi, double *divphi, long a);
// // double dwdphi(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a);
// // double dwdphi_smooth(double *phi, double *divphi, struct gradlayer **gradient, long gidy, long a);
// double calculate_equilibrium_composition(double T, long i, long a);
// double calculate_der_equilibrium_composition(double T, long i, long a);
// // double func_dcbdT_phase(double *mu, double T, long a, long k);
// // double func_dBbdT(double T, long a, long k);
// 
// //Function to build derived MPI datatype
// // void Build_derived_type(struct variables* myNode, MPI_Datatype* mesg_mpi_t_ptr);
// void copyXZ(long copy_to, long copy_from, long x_start, long x_end, struct variables* gridinfo);
// void copyYZ(long copy_to, long copy_from, struct variables* gridinfo);
// // void Build_derived_type(struct variables* myNode, MPI_Datatype* MPI_gridinfo) {
// // #ifdef ISOTHERMAL
// // 	int block_lengths[3];
// // 	MPI_Aint displacements[3];
// // 	MPI_Datatype typelists[3];
// // 	MPI_Aint start_address;
// // 	MPI_Aint address;
// // 	//fill block lengths
// //         block_lengths[0] = NUMPHASES;
// //         block_lengths[1] = NUMCOMPONENTS-1;
// //         block_lengths[2] = NUMPHASES;
// //         //Fill Typelists
// //         typelists[0]     = MPI_DOUBLE;
// //         typelists[1]     = MPI_DOUBLE;
// //         typelists[2]     = MPI_DOUBLE;
// // 	//Calculate Displacements
// //         displacements[0] = 0;
// //         MPI_Get_address(&((*myNode).phia),&start_address);
// //         MPI_Get_address(&((*myNode).compi),&address);
// //         displacements[1] = address - start_address;
// //         MPI_Get_address(&((*myNode).deltaphi),&address);
// //         displacements[2] = address - start_address;
// //         //Create and Commit new mpi type
// // 	MPI_Type_create_struct(3, block_lengths, displacements, typelists, MPI_gridinfo);
// // 	MPI_Type_commit(MPI_gridinfo);
// // #else
// // 	int block_lengths[4];
// // 	MPI_Aint displacements[4];
// // 	MPI_Datatype typelists[4];
// // 	MPI_Aint start_address;
// // 	MPI_Aint address;
// // 	//fill block lengths
// // 	block_lengths[0] = NUMPHASES;
// // 	block_lengths[1] = NUMCOMPONENTS-1;
// // 	block_lengths[2] = NUMPHASES;
// //         block_lengths[3] = 1;
// //         //Fill Typelists
// // 	typelists[0] = MPI_DOUBLE;
// // 	typelists[1] = MPI_DOUBLE;
// // 	typelists[2] = MPI_DOUBLE;
// //         typelists[3] = MPI_DOUBLE;
// // 	//Calculate Displacements
// // 	displacements[0] = 0;
// // 	MPI_Get_address(&((*myNode).phia),&start_address);
// // 	MPI_Get_address(&((*myNode).compi),&address);
// // 	displacements[1] = address - start_address;
// // 	MPI_Get_address(&((*myNode).deltaphi),&address);
// // 	displacements[2] = address - start_address;
// //         MPI_Get_address(&((*myNode).temperature),&address);
// //         displacements[3] = address - start_address;
// //         //Create and Commit new mpi type
// // 	MPI_Type_create_struct(4,block_lengths,displacements,typelists,MPI_gridinfo);
// // 	MPI_Type_commit(MPI_gridinfo);
// // #endif
// // }
// void Mpiinfo(long taskid);
// void mpiexchange_left_right(long taskid);
// void mpiexchange_top_bottom(long taskid);
// void mpiboundary_left_right(long taskid);
// void mpiboundary_top_bottom(long taskid);
// void sendtomaster();
// void receivefrmworker();
// void init_filling ();
// // void fill_phase_cube(FILE *fp, struct fill_cube fill_cube_paramaters);
// void fill_random();
// long assign_phase (double rand_phase);
// int check_FLAG(int *FLAG_FILLED, long xo, long yo, long radius);
// void assign_FLAG(int *FLAG_FILLED, long xo, long yo, long radius);
// 
// 
// // void fill_composition_cube(FILE *fp, long j, struct fill_cube* fill_cube_paramaters);
// 
// void fill_cells_phase_field (FILE *fp, struct variables*, long b);
// void fill_cells_phase_field_liquid(FILE *fp, struct variables*, long b);
// void fill_cells_composition (FILE *fp, struct variables*, long k);
// void find_numcubes(char *pattern, struct filling_type* filling_type_phase);
// void fillpattern(char *pattern, struct filling_type* filling_type_phase, long start_cube, long numperiods, struct variables* gridinfo);
// void write_cells_phasefield (FILE *fp, struct variables *gridinfo, long b);
// void write_cells_phasefield_MATRIX(FILE *fp, struct variables *gridinfo, long b);
// void write_cells_composition (FILE *fp, struct variables *gridinfo, double T, long k);
// void write_cells_temperature (FILE *fp, struct variables *gridinfo);
// void write_cells_composition_file0(FILE *fp, struct variables *gridinfo, double T, long k);
// void write_cells_composition_MATRIX(FILE *fp, struct variables *gridinfo, double T, long k);
// // void apply_shiftY(struct variables* gridinfo1, long INTERFACE_POS_GLOBAL, long offset, long rows);
// void apply_shiftY(struct variables* gridinfo1, long INTERFACE_POS_GLOBAL, int *offset, int *rows);
// long check_SHIFT(long x);
// // void apply_temperature_gradientY(struct variables* gridinfo1, long offset, long rows, long shift_OFFSET, long t);
// void apply_temperature_gradientY(struct variables* gridinfo1, int *offset, int *rows, long shift_OFFSET, long t);
// void apply_boundary_conditions(long taskid);


// void Build_derived_type(struct variables* myNode, MPI_Datatype* MPI_gridinfo) {
//   if(ISOTHERMAL) {
//     int block_lengths[3];
//     MPI_Aint displacements[3];
//     MPI_Datatype typelists[3];
//     MPI_Aint start_address;
//     MPI_Aint address;
//     //fill block lengths
//     block_lengths[0] = NUMPHASES;
//     block_lengths[1] = NUMCOMPONENTS-1;
//     block_lengths[2] = NUMPHASES;
//     //Fill Typelists
//     typelists[0]     = MPI_DOUBLE;
//     typelists[1]     = MPI_DOUBLE;
//     typelists[2]     = MPI_DOUBLE;
//     //Calculate Displacements
//     displacements[0] = 0;
//     MPI_Get_address(&((*myNode).phia),&start_address);
//     MPI_Get_address(&((*myNode).compi),&address);
//     displacements[1] = address - start_address;
//     MPI_Get_address(&((*myNode).deltaphi),&address);
//     displacements[2] = address - start_address;
//     //Create and Commit new mpi type
//     MPI_Type_create_struct(3, block_lengths, displacements, typelists, MPI_gridinfo);
//     MPI_Type_commit(MPI_gridinfo);
//   } else {
//     int block_lengths[4];
//     MPI_Aint displacements[4];
//     MPI_Datatype typelists[4];
//     MPI_Aint start_address;
//     MPI_Aint address;
//     //fill block lengths
//     block_lengths[0] = NUMPHASES;
//     block_lengths[1] = NUMCOMPONENTS-1;
//     block_lengths[2] = NUMPHASES;
//     block_lengths[3] = 1;
//     //Fill Typelists
//     typelists[0] = MPI_DOUBLE;
//     typelists[1] = MPI_DOUBLE;
//     typelists[2] = MPI_DOUBLE;
//     typelists[3] = MPI_DOUBLE;
//     //Calculate Displacements
//     displacements[0] = 0;
//     MPI_Get_address(&((*myNode).phia),&start_address);
//     MPI_Get_address(&((*myNode).compi),&address);
//     displacements[1] = address - start_address;
//     MPI_Get_address(&((*myNode).deltaphi),&address);
//     displacements[2] = address - start_address;
//     MPI_Get_address(&((*myNode).temperature),&address);
//     displacements[3] = address - start_address;
//     //Create and Commit new mpi type
//     MPI_Type_create_struct(4,block_lengths,displacements,typelists,MPI_gridinfo);
//     MPI_Type_commit(MPI_gridinfo);
//   }
// }

// double* MallocV(double *a,long m) {
//   double *a;
//   a=(double *)malloc(sizeof(*a));
//   return a;
// }
// double** MallocM(long m,long n) {
//   int i;
//   double **a;
//   a=(double **)malloc(m*sizeof(**a));
//     for(i=0;i<m;i++) {
//       a[i]=(double *)malloc(n*sizeof(*a[i]));
//     }
//    return a;
// }
// double*** Malloc3M(long m, long n, long k) {
//   int i;
//   double*** a;
//   a=(double ***)malloc(m*sizeof(***a));
//     for (i=0;i<m;i++) {
//     *a[i]=(double **)malloc(n*sizeof(**a[i]));
//     for (i=0; i<n; i++) {
//       a[i][j]  = (double *)malloc(k*sizeof(*a[i][j]));
//     }
//   }
//   return a;
// }
// // double*** arr3dAlloc(const int ind1, const int ind2, const int ind3)
// // {
// //   int i;
// //   int j;
// //   double*** array = (double***) malloc( (ind1 * sizeof(double*)) + (ind1*ind2 * sizeof(double**)) + (ind1*ind2*ind3 * sizeof(double)) );
// //   for(i = 0; i < ind1; ++i) {
// //     array[i] = (double**)(array + ind1) + i * ind2;
// //     for(j = 0; j < ind2; ++j) {
// //       array[i][j] = (double*)(array + ind1 + ind1*ind2) + i*ind2*ind3 + j*ind3;
// //     }
// //   }
// //   return array;
// // }
// void FreeM(double **a, long m, long n) {
//   int i;
//   for(i=0;i<m;i++) {
//     free(a[i]);
//   }
// }
