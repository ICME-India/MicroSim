#pragma once
#include <fstream>
#include <string>
#include <mutex>

using std::string;

enum class LogLevel{INFO, WARNING, ERROR, DEBUG};

namespace IO{

    class Logger {
        public:
            explicit Logger(const string& filename = "log.txt", LogLevel res_level = LogLevel::INFO, bool to_screen = false);
            ~Logger();

            void log(const string& message, LogLevel level = LogLevel::INFO);

        private:
            std::ofstream log_file;
            std::mutex log_mutex;
            int mpi_rank = 0;
            LogLevel print_level;
            bool print;
            bool to_screen;

            void check_status(LogLevel level);
            string getCurrentTime();
            string levelToString(LogLevel level);
    };
    
};