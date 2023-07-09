#ifdef __APPLE__           
#include <OpenCL/opencl.h>           
#else          
#include <CL/cl.h>           
#endif

#define MAX_SOURCE_SIZE (0x100000)


void CL_device_kernel_build() {

  printf("In CL_device_kernel_build\n");
  
  // Get all platforms
  ret = clGetPlatformIDs(0, NULL, &ret_num_platforms);
  if (ret_num_platforms == 0) {
    printf("Found 0 platforms!\n");
  }

  cl_platform_id platform_id[ret_num_platforms];
  if ( rank == 0 ) { 
    printf("No. of platforms %d\n", ret_num_platforms);
  }

  if ( ret_num_platforms == 0 ) { 
    printf("\n-------------------------------------------\n");
    printf("ERROR:\n");
    printf("\t OpenCL installation is not available on this machine \t\n");
    printf("\n-------------------------------------------\n");
    exit(1);
  }

  if (platformnumber > ret_num_platforms-1) {
    printf("\n-------------------------------------------\n");
    printf("ERROR:\n");
    printf("Choosen platformID is greater than available platforms\n");
    printf("\n-------------------------------------------\n");
    exit(1);
  }

  ret = clGetPlatformIDs(ret_num_platforms, platform_id, NULL);
  if ( rank == 0 ) { 
    printf("\n# Choosen platformid = %d\n", platformnumber);
  }

  ret = clGetDeviceIDs(platform_id[platformnumber], CL_DEVICE_TYPE_ALL, 0, NULL, &ret_num_devices);
  if ( rank == 0 ) { 
    printf("\nNo. of devices = %d in platformId = %d\n", ret_num_devices, platformnumber);
  }

  cl_device_id device_id[ret_num_devices];

  ret = clGetDeviceIDs(platform_id[platformnumber], CL_DEVICE_TYPE_ALL, ret_num_devices, device_id, NULL);

  // Each device information
  for (j = 0; j < ret_num_devices; j++) {
    if ( rank == 0 ) { 
      printf("\n  DeviceId = %d \n", j);
    }
    // Get device name
    ret = clGetDeviceInfo(device_id[j], CL_DEVICE_NAME, sizeof(string), &string, NULL);
    if ( rank == 0 ) { 
      printf("\t\tName: %s\n", string);
    }
    // Get device OpenCL version
    ret = clGetDeviceInfo(device_id[j], CL_DEVICE_OPENCL_C_VERSION, sizeof(string), &string, NULL);

    // Get Max. Compute units
    cl_uint num;
    ret = clGetDeviceInfo(device_id[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &num, NULL);
    if ( rank == 0 ) { 
      printf("\t\tMax. Compute Units: %d\n", num);
    }

    ret = clGetDeviceInfo(device_id[j], CL_DEVICE_NAME, sizeof(device_name), &device_name, NULL);
    if (ret != CL_SUCCESS) {
      printf("Error: Failed to access device name!\n");
    }
    if ( rank == 0 ) { 
      printf(" \t\tDevice is  %s ",device_name);
    }
    ret = clGetDeviceInfo(device_id[j], CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL);
    if (ret != CL_SUCCESS) {
      printf("Error: Failed to access device type information!\n");
    }
    if ( rank == 0 ) { 
      if (device_type  == CL_DEVICE_TYPE_GPU)
        printf(" GPU from ");
      else if (device_type == CL_DEVICE_TYPE_CPU)
        printf(" CPU from ");
      else 
        printf(" non  CPU or GPU processor from ");
    }
    ret = clGetDeviceInfo(device_id[j], CL_DEVICE_VENDOR, sizeof(vendor_name), &vendor_name, NULL);
    if (ret != CL_SUCCESS) {
      printf("Error: Failed to access device vendor name!\n");
    }
    if ( rank == 0 ) { 
      printf(" %s \n",vendor_name);
    }

    clGetDeviceInfo(device_id[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(local_size), &local_size, &local_size_size);
    if ( rank == 0 ) { 
      printf("\t\tCL_DEVICE_LOCAL_MEM_SIZE = %d\n", (int)local_size);
    }
  }

  if (numtasks == 1) {
    device_run = devicenumber;
  }
  else {
    device_run = rank%ret_num_devices;
  }

  if (device_run > ret_num_devices-1) {
    printf("\n-------------------------------------------\n");
    printf("ERROR:\n");
    printf("Choosen deviceID is greater than no. of available devices\n");
    printf("\n-------------------------------------------\n");
    exit(1);
  }


  //device_run = rank%ret_num_devices;

  //device_run = devicenumber;


  /*if ( (rank%2 == 0) ) {
    device_run = 0;
  }
  else {
    device_run = 1;
  }*/


  printf("\n*******************************************************************\n");
  printf("*    DeviceId = %d is choosen         *\n", device_run);
  printf("*******************************************************************\n\n");


  /* Create OpenCL Context */
  context = clCreateContext(NULL, 1, &device_id[device_run], NULL, NULL, &ret);
  //printf("context create %d\n", ret);

  /* Create command queue */
  cmdQ = clCreateCommandQueue(context, device_id[device_run], 0, &ret);
  //printf("command0 q %d\n", ret);


  FILE *fp;
  const char fileName[] = "solverloop/CL_Kim_Kernel.cl";
  size_t source_size;
  //char *source_str;
  /* Load kernel source file */
  fp = fopen(fileName, "rb");
  if (!fp) {
    fprintf(stderr, "Failed to load kernel file.\n");
  }

  source_str = (char *)malloc(MAX_SOURCE_SIZE);
  source_size = fread(source_str, 1, MAX_SOURCE_SIZE, fp);
  fclose(fp);

  /* Create kernel from source */
  program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
  if (ret!=CL_SUCCESS) {
    printf("program build with source error %d\n", ret);
  }
  ret = clBuildProgram(program, 1, &device_id[device_run], NULL, NULL, NULL);
  if (ret!=CL_SUCCESS) {
    printf("program build error %d\n", ret);
  }
  else { 
    printf("Rank = %d: \t Program built \n", rank);
  }

  char buildString[1000000];
  ret = clGetProgramBuildInfo(program, device_id[device_run], CL_PROGRAM_BUILD_LOG, sizeof(char) * 1000000, buildString, NULL);
  if(ret != CL_SUCCESS) {
    printf("Could not get build log from rank = %d.\n", rank);
    exit(1);
  }
  printf("%s\n", buildString);


}
