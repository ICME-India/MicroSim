#include "utilities/Logger.hpp"
#include <iostream>
#include <sstream>
#include <ctime>
#include <mpi.h>

IO::Logger::Logger(const string& filename, LogLevel res_level, bool opt){
    print_level = res_level;
    to_screen = opt;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    std::ostringstream final_filename;
    final_filename << "DATA/rank_" << mpi_rank << "_" << filename;
    log_file.open(final_filename.str(), std::ios::out | std::ios::trunc);
}

IO::Logger::~Logger() {
    if (log_file.is_open())
        log_file.close();
}

void IO::Logger::check_status(LogLevel level){
    switch (level){
        case LogLevel::INFO:
            print = (print_level == LogLevel::INFO || print_level == LogLevel::ERROR);
            break;
        case LogLevel::WARNING:
            print = (print_level == LogLevel::DEBUG);
            break;
        case LogLevel::ERROR:
            print = true;
            break;
        case LogLevel::DEBUG:
            print = true;
            break;
        default:
            print = false;
            break;
    }
}

void IO::Logger::log(const string& message, LogLevel level) {
    check_status(level);
    if(print){
        std::lock_guard<std::mutex> lock(log_mutex);
        string level_str = levelToString(level);
        string timestamp = getCurrentTime();

        std::ostringstream log_line;
        //log_line << "[" << timestamp << "] [Rank " << mpi_rank << "] [" << level_str << "] " << message;

        log_line<<message;

        if(to_screen)std::cout << message << std::endl;
        if (log_file.is_open()){
            log_file << log_line.str() << std::endl;
        }
    }
}


string IO::Logger::getCurrentTime() {
    std::time_t now = std::time(nullptr);
    char buf[20];
    if (std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&now)))
        return buf;
    return "unknown-time";
}

string IO::Logger::levelToString(LogLevel level) {
    switch (level) {
        case LogLevel::INFO: return "INFO";
        case LogLevel::WARNING: return "WARNING";
        case LogLevel::ERROR: return "ERROR";
        case LogLevel::DEBUG: return "DEBUG";
        default: return "UNKNOWN";
    }
}