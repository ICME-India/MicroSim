#pragma once

#include "utilities/Logger.hpp"
#include <sycl/sycl.hpp>
#include <mpi.h>

enum Device{CPU, GPU_A, FPGA};

class Device_Info {
    private:
        
        Device Accelarator;
        
        std::vector<sycl::device> Available_device;
        std::string name;
        std::string vendor;
        
        int max_compute_units;
        int workgroup_size;

        IO::Logger &logger;

        void Aspect_selector();
        sycl::device Select_Device();
        void getInfo();
        
    public:

        Device_Info(IO::Logger &logger);
        sycl::device selected_device;
        void LogInfo();
        ~Device_Info();
};