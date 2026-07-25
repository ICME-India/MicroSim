#include "comm/device_info.hpp"

#ifdef GPU
    Device_Info::Device_Info(IO::Logger &logg):
            Accelarator(GPU_A), 
            name(""), 
            vendor(""), 
            max_compute_units{0}, 
            workgroup_size{0},
            logger(logg){
        sycl::device selected_device = Device_Info::Select_Device();
        std::cout<<"Defined GPU device\n";
    }
#elif FPGA
    /// yet to be defined
#else
    Device_Info::Device_Info(IO::Logger &logg):
            Accelarator(CPU), 
            name(""), 
            vendor(""), 
            max_compute_units{0}, 
            workgroup_size{0},
            logger(logg) {
        sycl::device selected_device = Device_Info::Select_Device();
    }
#endif


void Device_Info::Aspect_selector() {

    int best_score = -1;

    for (const sycl::device &dev : Available_device) {
        int score = 0;
        std::string dev_name = dev.get_info<sycl::info::device::name>();
        std::string dev_vendor = dev.get_info<sycl::info::device::vendor>();

        if (!name.empty() && dev_name.find(name) != std::string::npos) {score++;}
        if (!vendor.empty() && dev_vendor.find(vendor) != std::string::npos) {score++;}

        if (score > best_score) {
            best_score = score;
            selected_device = dev;
        }
    }
    if (best_score <= 0) {
        logger.log("No matching aspect found", LogLevel::ERROR);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
}


void Device_Info::getInfo(){
    name = selected_device.get_info<sycl::info::device::name>();
    vendor = selected_device.get_info<sycl::info::device::vendor>();
    max_compute_units = selected_device.get_info<sycl::info::device::max_compute_units>();
    workgroup_size = selected_device.get_info<sycl::info::device::max_work_group_size>();
}


sycl::device Device_Info::Select_Device(){
    if(name.empty() && vendor.empty()){
        switch (Accelarator){
        case CPU:
            Available_device = sycl::device::get_devices(sycl::info::device_type::cpu);
            selected_device = Available_device[0];
            getInfo();
            return selected_device;

        case GPU_A:
            Available_device = sycl::device::get_devices(sycl::info::device_type::gpu);
            if(Available_device.empty()){
                logger.log("No GPU device available, switching to CPU",LogLevel::WARNING);
                Available_device = sycl::device::get_devices(sycl::info::device_type::cpu);
            }
            selected_device = Available_device[0];
            getInfo();
            return selected_device;

        default:
            logger.log("No matching device found", LogLevel::ERROR);
            MPI_Abort(MPI_COMM_WORLD, 0);
            return selected_device;

        }
    }else{
        Available_device = sycl::device::get_devices();
        Aspect_selector();
        getInfo();
        return selected_device;
    }
}


void Device_Info::LogInfo(){
    getInfo();
    logger.log("--------------------------------------------------------------", LogLevel::INFO);
    logger.log("                  Device Information                          ", LogLevel::INFO);
    logger.log("--------------------------------------------------------------", LogLevel::INFO);
    logger.log("selected device : "+Device_Info::name);
    logger.log("Vendor : "+Device_Info::vendor);
    logger.log("Maximum compute units : "+std::to_string(max_compute_units));
    logger.log("Maximum workgroup size : "+std::to_string(workgroup_size));

    logger.log("--------------------------------------------------------------", LogLevel::INFO);
}


Device_Info::~Device_Info(){};