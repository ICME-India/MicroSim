#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <string.h>
#include<time.h>
//#include<mpi.h>
#include<stddef.h>

#ifdef __APPLE__           
#include <OpenCL/opencl.h>           
#else          
#include <CL/cl.h>           
#endif

#define MAX_SOURCE_SIZE (0x100000)

#define LOCAL_SIZE 16

#define devicenumber 0
#define platformnumber 0


void main(int argc, char *argv[]) {
    
    int i;
    int j;
    int k;
    int index;
    int i1;
    int ig;
    int i0;
    
     int rank =0;
     
    int device_run;
    
    
    /*************************************/
    
    //cl_platform_id platform_id = NULL;
    //cl_device_id device_id = NULL;
    cl_context context = NULL;
    cl_command_queue cmdQ = NULL;
    cl_program program;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
    cl_ulong local_size;
    cl_int cl_local_size;
    size_t kernel_code_size, local_size_size;
    cl_event kerndone;
    cl_char device_name[1024] = {0}; 
    cl_char vendor_name[1024] = {0};
    cl_device_type device_type; 
    cl_char string[10240] = {0};
    cl_uint work_dim=2;
    
    /****************************************************************/
    // Get all platforms
    ret = clGetPlatformIDs(0, NULL, &ret_num_platforms);
    if (ret_num_platforms == 0) {
      printf("Found 0 platforms!\n");
    }
    
    
    cl_platform_id platform_id[ret_num_platforms];
    if ( rank == 0 ) { 
      printf("No. of platforms %d\n", ret_num_platforms);
    }
    
    ret = clGetPlatformIDs(ret_num_platforms, platform_id, NULL);
    if ( rank == 0 ) { 
      printf("Choosen platformid = %d\n", platformnumber);
    }
    
    ret = clGetDeviceIDs(platform_id[platformnumber], CL_DEVICE_TYPE_ALL, 0, NULL, &ret_num_devices);
    if ( rank == 0 ) { 
      printf("\nNo. of devices = %d in platformId = %d\n", ret_num_devices, platformnumber);
    }
    //printf("\nAll the devices in platform %d:\n", platformnumber);
    
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
	if(device_type  == CL_DEVICE_TYPE_GPU)
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
      // /* Get platform/device information */
      // ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
      // ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &ret_num_devices);
      // printf("no  of dev %d\n", ret_num_devices);
      //    clGetDeviceInfo(device_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(local_size), &local_size, &local_size_size);
      clGetDeviceInfo(device_id[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(local_size), &local_size, &local_size_size);
      if ( rank == 0 ) { 
	printf("\t\tCL_DEVICE_LOCAL_MEM_SIZE = %d\n", (int)local_size);
      }
    }
    
    
}
