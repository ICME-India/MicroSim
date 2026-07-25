#include "thermo/Function_F4.hpp"

#include <gsl/gsl_spline.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

namespace microsim {

void ThermoDynamicsF4::extract_spline_coefficients(
    gsl_spline *spline,
    SplineCoefficients &coeffs,
    int n_points
) {
    coeffs.n_points = n_points;

    // Validate input
    if (!spline || n_points <= 0) {
        std::cerr << "ERROR: Invalid spline or n_points in extract_spline_coefficients" << std::endl;
        coeffs.c0 = coeffs.c1 = coeffs.c2 = coeffs.c3 = coeffs.x = nullptr;
        coeffs.n_points = 0;
        return;
    }

    // Allocate memory for coefficients
    coeffs.c0 = new double[n_points];
    coeffs.c1 = new double[n_points];
    coeffs.c2 = new double[n_points];
    coeffs.c3 = new double[n_points];
    coeffs.x = new double[n_points];  // Changed: only n_points, not n_points+1

    // Extract breakpoints
    // GSL spline has exactly n_points x-values (indices 0 to n_points-1)
    for(int i = 0; i < n_points; ++i) {  // Changed: < instead of <=
        coeffs.x[i] = spline->x[i];
        if(std::isnan(coeffs.x[i]) || std::isinf(coeffs.x[i])) {
            std::cerr << "ERROR: NaN/Inf in spline->x[" << i << "] = " << coeffs.x[i] << std::endl;
        }
    }

    // Extract cubic spline coefficients by evaluating derivatives
    // Instead of accessing GSL internals, use the GSL API to compute coefficients
    // For a cubic spline: f(x) = a + b*(x-xi) + c*(x-xi)^2 + d*(x-xi)^3
    // where a = f(xi), b = f'(xi), c = f''(xi)/2, d depends on boundary conditions

    gsl_interp_accel *acc = gsl_interp_accel_alloc();

    for(int i = 0; i < n_points; ++i) {
        double xi = coeffs.x[i];

        // c0 = f(xi)
        coeffs.c0[i] = gsl_spline_eval(spline, xi, acc);

        // c1 = f'(xi) - first derivative
        if (i < n_points - 1) {
            coeffs.c1[i] = gsl_spline_eval_deriv(spline, xi, acc);

            // c2 = f''(xi)/2 - second derivative / 2
            coeffs.c2[i] = gsl_spline_eval_deriv2(spline, xi, acc) / 2.0;

            // c3 needs to be computed from continuity with next interval
            // For simplicity, we'll compute it from the next point
            double xi1 = coeffs.x[i + 1];
            double dx = xi1 - xi;
            double fi = coeffs.c0[i];
            double fi1 = gsl_spline_eval(spline, xi1, acc);
            double dfi = coeffs.c1[i];
            double dfi1 = gsl_spline_eval_deriv(spline, xi1, acc);

            // Using Hermite interpolation formula to get c3
            coeffs.c3[i] = (2.0 * (fi - fi1) + (dfi + dfi1) * dx) / (dx * dx * dx);
        } else {
            // Last point - extrapolation not needed
            coeffs.c1[i] = 0.0;
            coeffs.c2[i] = 0.0;
            coeffs.c3[i] = 0.0;
        }

        // Validate extracted coefficients
        if(std::isnan(coeffs.c0[i]) || std::isinf(coeffs.c0[i]) ||
           std::isnan(coeffs.c1[i]) || std::isinf(coeffs.c1[i]) ||
           std::isnan(coeffs.c2[i]) || std::isinf(coeffs.c2[i]) ||
           std::isnan(coeffs.c3[i]) || std::isinf(coeffs.c3[i])) {
            std::cerr << "ERROR: NaN/Inf in extracted coefficients at index " << i << std::endl;
            std::cerr << "  x=" << xi << " c0=" << coeffs.c0[i] << " c1=" << coeffs.c1[i]
                      << " c2=" << coeffs.c2[i] << " c3=" << coeffs.c3[i] << std::endl;
        }
    }

    gsl_interp_accel_free(acc);
}

void ThermoDynamicsF4::Initialize_host(phases_components *ref) {
    std::cout << "Initializing ThermoDynamics on host..." << std::endl;
    if (!queue_ptr) {
        std::cerr << "Error: SYCL queue not set before Initialize_host" << std::endl;
        return;
    }
    sycl::queue &q = *queue_ptr;

    // Initialize thermo_phase mapping using SYCL shared memory
    thermo_phase = sycl::malloc_shared<int>(nPhase, q);
    for(int i = 0; i < nPhase; ++i) {
        // Find the index of phases_map[i] in phases_tdb
        std::string phase_name = ref->phases_map[i];
        for(int j = 0; j < ref->phases_tdb.size(); ++j) {
            if(ref->phases_tdb[j] == phase_name) {
                thermo_phase[i] = j;
                break;
            }
        }
        std::cout << "Phase mapping: phase[" << i << "] = " << ref->phases_map[i]
                  << " -> thermo_phase[" << i << "] = " << thermo_phase[i] << std::endl;
    }
    q.wait(); // Ensure initialization is complete

    char composition_file[1024];
    char hessian_file[1024];
    sprintf(composition_file, "tdbs_encrypted/Composition_%s.csv", ref->phases_tdb[0].c_str());
    std::ifstream file(composition_file);

    if (!file.is_open()) {
        std::cerr << "Error: file not found: " << composition_file << std::endl;
        return;
    }

    long num_lines{0};
    std::string line;

    /// @brief Count the number of lines in the file
    while (std::getline(file, line)){
        num_lines++;
    }

    std::cout << "Number of lines : " << num_lines << " \n";

    file.clear();
    file.seekg(0, std::ios::beg);

    /// @brief stores the Composition and hessian
    auto Malloc4M = [](int a, int b, int c, int d) {
        double ****array = new double***[a];
        for(int i = 0; i < a; ++i) {
            array[i] = new double**[b];
            for(int j = 0; j < b; ++j) {
                array[i][j] = new double*[c];
                for(int k = 0; k < c; ++k) {
                    array[i][j][k] = new double[d];
                }
            }
        }
        return array;
    };

    auto MallocM = [](int a, int b) {
        double **array = new double*[a];
        for(int i = 0; i < a; ++i) {
            array[i] = new double[b];
        }
        return array;
    };

    auto parsecsvList = [](std::vector<double>& result, const std::string& str) {
        std::stringstream ss(str);
        std::string item;
        while (std::getline(ss, item, ',')) {
            item.erase(0, item.find_first_not_of(" \t"));
            item.erase(item.find_last_not_of(" \t") + 1);
            result.push_back(std::stod(item));
        }
    };

    double ****Comp_ES{nullptr};
    double ****ThF{nullptr};
    double **T_ES{nullptr};
    double **T_Thf{nullptr};

    gsl_interp_accel**** acc_ES{nullptr};
    gsl_spline ****spline_ES{nullptr};
    gsl_interp_accel ****acc_Thf{nullptr};
    gsl_spline ****spline_Thf{nullptr};

    std::cout << "Allocating memory..." << std::endl;
    Comp_ES = Malloc4M(nPhase, nComp-1, 2, num_lines-1);
    ThF = Malloc4M(nPhase, nComp-1, nComp-1, num_lines-1);
    T_ES = MallocM(nPhase, num_lines);
    T_Thf = MallocM(nPhase, num_lines);

    // Calculate flat array sizes
    int ES_coeffs_size = nPhase * (nComp - 1) * 2;
    int Thf_coeffs_size = nPhase * (nComp - 1) * (nComp - 1);

    std::cout << "Allocating flat coefficient arrays..." << std::endl;
    std::cout << "  ES_coeffs_size: " << ES_coeffs_size << std::endl;
    std::cout << "  Thf_coeffs_size: " << Thf_coeffs_size << std::endl;

    // Allocate flat arrays on host
    ES_coeffs_flat = new SplineCoefficients[ES_coeffs_size];
    Thf_coeffs_flat = new SplineCoefficients[Thf_coeffs_size];

    // Initialize all structures
    for(int i = 0; i < ES_coeffs_size; ++i) {
        ES_coeffs_flat[i] = SplineCoefficients();
    }
    for(int i = 0; i < Thf_coeffs_size; ++i) {
        Thf_coeffs_flat[i] = SplineCoefficients();
    }

    // Read composition files
    for(int a = 0; a < nPhase - 1; ++a){
        sprintf(composition_file, "tdbs_encrypted/Composition_%s.csv", ref->phases_tdb[a].c_str());
        std::cout << "Reading " << composition_file << std::endl;
        std::ifstream file(composition_file);

        if (!file.is_open()) {
            std::cerr << "Error: file not found: " << composition_file << std::endl;
            continue;
        }

        std::string header;
        std::getline(file, header);

        for(int i = 0; i < num_lines - 1; ++i){
            std::string line;
            std::getline(file, line);
            std::vector<double> tempString;
            parsecsvList(tempString, line);

            T_ES[a][i] = tempString[0];

            // Check for NaN/Inf in CSV data
            if(std::isnan(T_ES[a][i]) || std::isinf(T_ES[a][i])) {
                std::cerr << "ERROR: NaN/Inf in T_ES[" << a << "][" << i << "] = " << T_ES[a][i] << std::endl;
            }

            for(int j = 0; j < nComp - 1; ++j){
                Comp_ES[a][j][0][i] = tempString[j + 1];
                Comp_ES[a][j][1][i] = tempString[nComp + j];

                if(std::isnan(Comp_ES[a][j][0][i]) || std::isinf(Comp_ES[a][j][0][i])) {
                    std::cerr << "ERROR: NaN/Inf in Comp_ES[" << a << "][" << j << "][0][" << i << "]" << std::endl;
                }
                if(std::isnan(Comp_ES[a][j][1][i]) || std::isinf(Comp_ES[a][j][1][i])) {
                    std::cerr << "ERROR: NaN/Inf in Comp_ES[" << a << "][" << j << "][1][" << i << "]" << std::endl;
                }
            }
        }
        file.close();
    }

    // Read hessian files
    for(int a = 0; a < nPhase; ++a){
        sprintf(hessian_file, "tdbs_encrypted/HSN_%s.csv", ref->phases_tdb[a].c_str());
        std::cout << "Reading " << hessian_file << std::endl;
        std::ifstream file(hessian_file);

        if (!file.is_open()) {
            std::cerr << "Error: file not found: " << hessian_file << std::endl;
            continue;
        }

        std::string header;
        std::getline(file, header);

        for(int i = 0; i < num_lines - 1; ++i){
            std::getline(file, line);
            std::vector<double> tempString;
            parsecsvList(tempString, line);

            T_Thf[a][i] = tempString[0];

            // Check for NaN/Inf in temperature data
            if(std::isnan(T_Thf[a][i]) || std::isinf(T_Thf[a][i])) {
                std::cerr << "ERROR: NaN/Inf in T_Thf[" << a << "][" << i << "] = " << T_Thf[a][i] << std::endl;
            }

            int count{nComp};

            for(int j = 0; j < nComp - 1; ++j){
                ThF[a][j][j][i] = tempString[j + 1];

                // Check for NaN/Inf in Hessian data
                if(std::isnan(ThF[a][j][j][i]) || std::isinf(ThF[a][j][j][i])) {
                    std::cerr << "ERROR: NaN/Inf in ThF[" << a << "][" << j << "][" << j << "][" << i
                              << "] = " << ThF[a][j][j][i] << std::endl;
                }

                for(int k = j + 1; k < nComp - 1; ++k){
                    ThF[a][j][k][i] = tempString[count];
                    ThF[a][k][j][i] = ThF[a][j][k][i];

                    if(std::isnan(ThF[a][j][k][i]) || std::isinf(ThF[a][j][k][i])) {
                        std::cerr << "ERROR: NaN/Inf in ThF[" << a << "][" << j << "][" << k << "][" << i
                                  << "] = " << ThF[a][j][k][i] << std::endl;
                    }

                    count++;
                }
            }
        }

        file.close();
    }

    std::cout << "Initializing temporary GSL splines..." << std::endl;
    acc_ES = new gsl_interp_accel***[nPhase];
    spline_ES = new gsl_spline***[nPhase];
    acc_Thf = new gsl_interp_accel***[nPhase];
    spline_Thf = new gsl_spline***[nPhase];

    // Create and extract Thf spline coefficients
    std::cout << "Creating Thf splines and extracting coefficients..." << std::endl;
    for(long a = 0; a < nPhase; ++a){
        acc_Thf[a] = new gsl_interp_accel**[nComp - 1];
        spline_Thf[a] = new gsl_spline**[nComp - 1];
        for(long k = 0; k < nComp - 1; ++k){
            acc_Thf[a][k] = new gsl_interp_accel*[nComp - 1];
            spline_Thf[a][k] = new gsl_spline*[nComp - 1];
            for(long l = 0; l < nComp - 1; ++l){
                acc_Thf[a][k][l] = gsl_interp_accel_alloc();
                spline_Thf[a][k][l] = gsl_spline_alloc(gsl_interp_cspline, num_lines - 1);
                gsl_spline_init(spline_Thf[a][k][l], T_Thf[a], ThF[a][k][l], num_lines - 1);
                
                // Extract coefficients and store in flat array
                int flat_idx = get_Thf_flat_index(a, k, l);
                extract_spline_coefficients(spline_Thf[a][k][l], Thf_coeffs_flat[flat_idx], num_lines - 1);
                std::cout << "  Thf[" << a << "][" << k << "][" << l << "] -> flat[" << flat_idx << "]" << std::endl;
            }
        }
    }

    // Create and extract ES spline coefficients
    std::cout << "Creating ES splines and extracting coefficients..." << std::endl;
    for(long a = 0; a < nPhase - 1; ++a){
        std::cout << "  Phase " << a << "/" << nPhase - 1 << std::endl;
        acc_ES[a] = new gsl_interp_accel**[nComp - 1];
        spline_ES[a] = new gsl_spline**[nComp - 1];
        for(int k = 0; k < nComp - 1; ++k){
            acc_ES[a][k] = new gsl_interp_accel*[2];
            spline_ES[a][k] = new gsl_spline*[2];
            for(int l = 0; l < 2; ++l){
                acc_ES[a][k][l] = gsl_interp_accel_alloc();
                spline_ES[a][k][l] = gsl_spline_alloc(gsl_interp_cspline, num_lines - 1);
                gsl_spline_init(spline_ES[a][k][l], T_ES[a], Comp_ES[a][k][l], num_lines - 1);
                
                // Extract coefficients and store in flat array
                int flat_idx = get_ES_flat_index(a, k, l);
                extract_spline_coefficients(spline_ES[a][k][l], ES_coeffs_flat[flat_idx], num_lines - 1);
                std::cout << "  ES[" << a << "][" << k << "][" << l << "] -> flat[" << flat_idx << "]" << std::endl;
            }
        }
    }
    
    // Handle the last phase - copy coefficients from first phase
    if (nPhase > 1) {
        std::cout << "Copying ES coefficients for last phase..." << std::endl;
        long last = nPhase - 1;
        for(int k = 0; k < nComp - 1; ++k){
            for(int l = 0; l < 2; ++l){
                int last_idx = get_ES_flat_index(last, k, l);
                int first_idx = get_ES_flat_index(0, k, l);
                
                // Share the same coefficient arrays (shallow copy)
                ES_coeffs_flat[last_idx] = ES_coeffs_flat[first_idx];
                std::cout << "  ES[" << last << "][" << k << "][" << l << "] = ES[0][" << k << "][" << l << "]" << std::endl;
            }
        }
    }

    std::cout << "Freeing GSL spline objects..." << std::endl;
    // Free GSL objects - we only need the coefficients
    for (int a = 0; a < nPhase; ++a) {
        for (int k = 0; k < nComp - 1; ++k) {
            for (int l = 0; l < nComp - 1; ++l) {
                gsl_interp_accel_free(acc_Thf[a][k][l]);
                gsl_spline_free(spline_Thf[a][k][l]);
            }
            delete[] acc_Thf[a][k];
            delete[] spline_Thf[a][k];
        }
        delete[] acc_Thf[a];
        delete[] spline_Thf[a];
    }
    delete[] acc_Thf;
    delete[] spline_Thf;

    for(int a = 0; a < nPhase - 1; ++a) {
        for(int k = 0; k < nComp - 1; ++k) {
            for(int l = 0; l < 2; ++l) {
                gsl_interp_accel_free(acc_ES[a][k][l]);
                gsl_spline_free(spline_ES[a][k][l]);
            }
            delete[] acc_ES[a][k];
            delete[] spline_ES[a][k];
        }
        delete[] acc_ES[a];
        delete[] spline_ES[a];
    }
    delete[] acc_ES;
    delete[] spline_ES;

    // ============================================================
    // SYCL: Copy coefficient data to device
    // ============================================================
    std::cout << "Copying spline coefficients to SYCL device..." << std::endl;
    
    // Calculate total memory needed for all coefficient arrays
    int total_ES_points = 0;
    int total_Thf_points = 0;
    
    for(int i = 0; i < ES_coeffs_size; ++i) {
        total_ES_points += ES_coeffs_flat[i].n_points;
    }
    for(int i = 0; i < Thf_coeffs_size; ++i) {
        total_Thf_points += Thf_coeffs_flat[i].n_points;
    }
    
    std::cout << "  Total ES points: " << total_ES_points << std::endl;
    std::cout << "  Total Thf points: " << total_Thf_points << std::endl;

    // Allocate shared memory using SYCL USM (accessible from both host and device)
    d_ES_coeffs_flat = sycl::malloc_shared<SplineCoefficients>(ES_coeffs_size, q);
    d_Thf_coeffs_flat = sycl::malloc_shared<SplineCoefficients>(Thf_coeffs_size, q);

    // Allocate contiguous shared memory for all coefficient arrays
    d_ES_c0 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_c1 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_c2 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_c3 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_x = sycl::malloc_shared<double>(total_ES_points, q);  // Fixed: just total_ES_points

    d_Thf_c0 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_c1 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_c2 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_c3 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_x = sycl::malloc_shared<double>(total_Thf_points, q);  // Fixed: just total_Thf_points

    // Pack coefficients into contiguous arrays and copy to device
    int ES_offset = 0;
    int ES_x_offset = 0;
    SplineCoefficients *temp_ES_coeffs = new SplineCoefficients[ES_coeffs_size];
    
    for(int i = 0; i < ES_coeffs_size; ++i) {
        int n = ES_coeffs_flat[i].n_points;
        if (n > 0) {
            // Copy coefficient data to contiguous arrays
            q.memcpy(d_ES_c0 + ES_offset, ES_coeffs_flat[i].c0, n * sizeof(double)).wait();
            q.memcpy(d_ES_c1 + ES_offset, ES_coeffs_flat[i].c1, n * sizeof(double)).wait();
            q.memcpy(d_ES_c2 + ES_offset, ES_coeffs_flat[i].c2, n * sizeof(double)).wait();
            q.memcpy(d_ES_c3 + ES_offset, ES_coeffs_flat[i].c3, n * sizeof(double)).wait();
            q.memcpy(d_ES_x + ES_x_offset, ES_coeffs_flat[i].x, n * sizeof(double)).wait();  // Fixed: n not n+1
            
            // Update structure to point to device memory
            temp_ES_coeffs[i].n_points = n;
            temp_ES_coeffs[i].c0 = d_ES_c0 + ES_offset;
            temp_ES_coeffs[i].c1 = d_ES_c1 + ES_offset;
            temp_ES_coeffs[i].c2 = d_ES_c2 + ES_offset;
            temp_ES_coeffs[i].c3 = d_ES_c3 + ES_offset;
            temp_ES_coeffs[i].x = d_ES_x + ES_x_offset;
            
            ES_offset += n;
            ES_x_offset += n;  // Fixed: n not n+1
        }
    }
    
    // Copy structures to device
    q.memcpy(d_ES_coeffs_flat, temp_ES_coeffs, ES_coeffs_size * sizeof(SplineCoefficients)).wait();
    delete[] temp_ES_coeffs;

    // Do the same for Thf coefficients
    int Thf_offset = 0;
    int Thf_x_offset = 0;
    SplineCoefficients *temp_Thf_coeffs = new SplineCoefficients[Thf_coeffs_size];
    
    for(int i = 0; i < Thf_coeffs_size; ++i) {
        int n = Thf_coeffs_flat[i].n_points;
        if (n > 0) {
            q.memcpy(d_Thf_c0 + Thf_offset, Thf_coeffs_flat[i].c0, n * sizeof(double)).wait();
            q.memcpy(d_Thf_c1 + Thf_offset, Thf_coeffs_flat[i].c1, n * sizeof(double)).wait();
            q.memcpy(d_Thf_c2 + Thf_offset, Thf_coeffs_flat[i].c2, n * sizeof(double)).wait();
            q.memcpy(d_Thf_c3 + Thf_offset, Thf_coeffs_flat[i].c3, n * sizeof(double)).wait();
            q.memcpy(d_Thf_x + Thf_x_offset, Thf_coeffs_flat[i].x, n * sizeof(double)).wait();  // Fixed: n not n+1
            
            temp_Thf_coeffs[i].n_points = n;
            temp_Thf_coeffs[i].c0 = d_Thf_c0 + Thf_offset;
            temp_Thf_coeffs[i].c1 = d_Thf_c1 + Thf_offset;
            temp_Thf_coeffs[i].c2 = d_Thf_c2 + Thf_offset;
            temp_Thf_coeffs[i].c3 = d_Thf_c3 + Thf_offset;
            temp_Thf_coeffs[i].x = d_Thf_x + Thf_x_offset;
            
            Thf_offset += n;
            Thf_x_offset += n;  // Fixed: n not n+1
        }
    }
    
    q.memcpy(d_Thf_coeffs_flat, temp_Thf_coeffs, Thf_coeffs_size * sizeof(SplineCoefficients)).wait();
    delete[] temp_Thf_coeffs;

    std::cout << "SYCL device memory transfer complete!" << std::endl;

    // Diagnostic: Check temperature ranges in spline data
    std::cout << "\n=== Spline Temperature Range Diagnostics ===" << std::endl;
    std::cout << "First ES spline detailed check:" << std::endl;
    if (ES_coeffs_flat[0].n_points > 0) {
        std::cout << "  n_points = " << ES_coeffs_flat[0].n_points << std::endl;
        std::cout << "  x ptr = " << (void*)ES_coeffs_flat[0].x << std::endl;
        std::cout << "  First 5 x values: ";
        for(int j = 0; j < std::min(5, ES_coeffs_flat[0].n_points); ++j) {
            std::cout << ES_coeffs_flat[0].x[j] << " ";
        }
        std::cout << std::endl;
        std::cout << "  Last 5 x values: ";
        int start = std::max(0, ES_coeffs_flat[0].n_points - 5);
        for(int j = start; j < ES_coeffs_flat[0].n_points; ++j) {
            std::cout << ES_coeffs_flat[0].x[j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nSummary of temperature ranges (showing first 10 only):" << std::endl;
    for(int i = 0; i < std::min(10, ES_coeffs_size); ++i) {
        if (ES_coeffs_flat[i].n_points > 0 && ES_coeffs_flat[i].x != nullptr) {
            std::cout << "ES_coeffs[" << i << "] n=" << ES_coeffs_flat[i].n_points
                      << " T range: [" << ES_coeffs_flat[i].x[0] << ", "
                      << ES_coeffs_flat[i].x[ES_coeffs_flat[i].n_points - 1] << "]" << std::endl;
        }
    }
    for(int i = 0; i < std::min(10, Thf_coeffs_size); ++i) {
        if (Thf_coeffs_flat[i].n_points > 0 && Thf_coeffs_flat[i].x != nullptr) {
            std::cout << "Thf_coeffs[" << i << "] n=" << Thf_coeffs_flat[i].n_points
                      << " T range: [" << Thf_coeffs_flat[i].x[0] << ", "
                      << Thf_coeffs_flat[i].x[Thf_coeffs_flat[i].n_points - 1] << "]" << std::endl;
        }
    }

    // Check if T is within range
    std::cout << "\nRequested temperature T = " << T << std::endl;

    // Diagnostic: Check for NaN/Inf in spline coefficients
    std::cout << "\n=== Checking Spline Coefficients for NaN/Inf ===" << std::endl;
    bool coeff_has_nan = false;
    for(int i = 0; i < Thf_coeffs_size; ++i) {
        if (Thf_coeffs_flat[i].n_points > 0) {
            for(int j = 0; j < Thf_coeffs_flat[i].n_points; ++j) {
                if(std::isnan(Thf_coeffs_flat[i].c0[j]) || std::isinf(Thf_coeffs_flat[i].c0[j]) ||
                   std::isnan(Thf_coeffs_flat[i].c1[j]) || std::isinf(Thf_coeffs_flat[i].c1[j]) ||
                   std::isnan(Thf_coeffs_flat[i].c2[j]) || std::isinf(Thf_coeffs_flat[i].c2[j]) ||
                   std::isnan(Thf_coeffs_flat[i].c3[j]) || std::isinf(Thf_coeffs_flat[i].c3[j])) {
                    std::cerr << "NaN/Inf in Thf_coeffs[" << i << "][" << j << "]" << std::endl;
                    coeff_has_nan = true;
                }
            }
        }
    }
    if (!coeff_has_nan) {
        std::cout << "No NaN/Inf found in Thf spline coefficients." << std::endl;
    }

    // Rest of the function continues...
    std::cout << "\nCalling function_A with T=" << T << std::endl;

    // After extracting all coefficients, call function_A to populate A matrix
    function_A(T);

    // Print A matrix for debugging
    std::cout << "A matrix after initialization:" << std::endl;
    for(int a = 0; a < nPhase; ++a) {
        std::cout << "Phase " << a << ":" << std::endl;
        for(int i = 0; i < nComp - 1; ++i) {
            for(int j = 0; j < nComp - 1; ++j) {
                std::cout << A[a_ind(a, i, j)] << " ";
            }
            std::cout << std::endl;
        }
    }

    //populate mu_c and c_mu with initial guesses
    for(int a = 0; a < nPhase; ++a) {
        for(int i = 0; i < nComp - 1; ++i) {
            for(int j = 0; j < nComp - 1; ++j){
                int muc_idx = muc_ind(a, i, j);
                int a_indx = a_ind(a, i, j);
                if(i==j){
                    muc[muc_idx] = 2.0*A[a_indx];
                }else{
                    muc[muc_idx] = A[a_indx];
                }
            }
        }
        MatrixInversion_device(&muc[a * (nComp-1) * (nComp-1)], &cmu[a * (nComp-1) * (nComp-1)], nComp - 1);
    }
    
    // Populate B and C
    for(int a = 0; a < nPhase; ++a) {
        for(int i = 0; i < nComp - 1; ++i) {
            int b_idx = b_ind(a, i);
            B[b_idx] = functionB(T, i, a);
        }
        C[c_ind(a)] = functionC(T, a);
    }

    // Debug: Print some A, B, C values to verify
    std::cout << "Sample A[0][0][0] = " << A[a_ind(0, 0, 0)] << std::endl;
    std::cout << "Sample B[0][0] = " << B[b_ind(0, 0)] << std::endl;
    std::cout << "Sample C[0] = " << C[c_ind(0)] << std::endl;

    // Check for NaN
    bool has_nan = false;
    for(int a = 0; a < nPhase; ++a) {
        for(int i = 0; i < nComp - 1; ++i) {
            for(int j = 0; j < nComp - 1; ++j) {
                if(std::isnan(A[a_ind(a, i, j)])) {
                    std::cerr << "NaN detected in A[" << a << "][" << i << "][" << j << "]" << std::endl;
                    has_nan = true;
                }
            }
            if(std::isnan(B[b_ind(a, i)])) {
                std::cerr << "NaN detected in B[" << a << "][" << i << "]" << std::endl;
                has_nan = true;
            }
        }
        if(std::isnan(C[c_ind(a)])) {
            std::cerr << "NaN detected in C[" << a << "]" << std::endl;
            has_nan = true;
        }
    }

    if(!has_nan) {
        std::cout << "No NaN values detected in A, B, C matrices." << std::endl;
    }

    std::cout << "Host initialization complete." << std::endl;
}

void ThermoDynamicsF4::copy_to_device(sycl::queue &q) {
    std::cout << "Copying spline coefficients to device..." << std::endl;
    
    queue_ptr = &q;
    
    // Calculate sizes
    int ES_coeffs_size = nPhase * (nComp - 1) * 2;
    int Thf_coeffs_size = nPhase * (nComp - 1) * (nComp - 1);
    
    // Allocate shared memory (accessible from both host and device)
    d_ES_coeffs_flat = sycl::malloc_shared<SplineCoefficients>(ES_coeffs_size, q);
    d_Thf_coeffs_flat = sycl::malloc_shared<SplineCoefficients>(Thf_coeffs_size, q);

    // Count total points
    int total_ES_points = 0;
    int total_Thf_points = 0;

    for(int i = 0; i < ES_coeffs_size; ++i) {
        total_ES_points += ES_coeffs_flat[i].n_points;
    }
    for(int i = 0; i < Thf_coeffs_size; ++i) {
        total_Thf_points += Thf_coeffs_flat[i].n_points;
    }

    // Allocate contiguous shared memory
    d_ES_c0 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_c1 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_c2 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_c3 = sycl::malloc_shared<double>(total_ES_points, q);
    d_ES_x = sycl::malloc_shared<double>(total_ES_points + ES_coeffs_size, q);

    d_Thf_c0 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_c1 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_c2 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_c3 = sycl::malloc_shared<double>(total_Thf_points, q);
    d_Thf_x = sycl::malloc_shared<double>(total_Thf_points, q);  // Fixed: just total_Thf_points
    
    // Pack coefficients into contiguous arrays and copy to device
    int ES_offset = 0;
    int ES_x_offset = 0;
    SplineCoefficients *temp_ES_coeffs = new SplineCoefficients[ES_coeffs_size];
    
    for(int i = 0; i < ES_coeffs_size; ++i) {
        int n = ES_coeffs_flat[i].n_points;
        if (n > 0) {
            // Copy coefficient data to contiguous arrays
            q.memcpy(d_ES_c0 + ES_offset, ES_coeffs_flat[i].c0, n * sizeof(double)).wait();
            q.memcpy(d_ES_c1 + ES_offset, ES_coeffs_flat[i].c1, n * sizeof(double)).wait();
            q.memcpy(d_ES_c2 + ES_offset, ES_coeffs_flat[i].c2, n * sizeof(double)).wait();
            q.memcpy(d_ES_c3 + ES_offset, ES_coeffs_flat[i].c3, n * sizeof(double)).wait();
            q.memcpy(d_ES_x + ES_x_offset, ES_coeffs_flat[i].x, n * sizeof(double)).wait();  // Fixed: n not n+1
            
            // Update structure to point to device memory
            temp_ES_coeffs[i].n_points = n;
            temp_ES_coeffs[i].c0 = d_ES_c0 + ES_offset;
            temp_ES_coeffs[i].c1 = d_ES_c1 + ES_offset;
            temp_ES_coeffs[i].c2 = d_ES_c2 + ES_offset;
            temp_ES_coeffs[i].c3 = d_ES_c3 + ES_offset;
            temp_ES_coeffs[i].x = d_ES_x + ES_x_offset;
            
            ES_offset += n;
            ES_x_offset += n;  // Fixed: n not n+1
        }
    }
    
    // Copy structures to device
    q.memcpy(d_ES_coeffs_flat, temp_ES_coeffs, ES_coeffs_size * sizeof(SplineCoefficients)).wait();
    delete[] temp_ES_coeffs;

    // Do the same for Thf coefficients
    int Thf_offset = 0;
    int Thf_x_offset = 0;
    SplineCoefficients *temp_Thf_coeffs = new SplineCoefficients[Thf_coeffs_size];
    
    for(int i = 0; i < Thf_coeffs_size; ++i) {
        int n = Thf_coeffs_flat[i].n_points;
        if (n > 0) {
            q.memcpy(d_Thf_c0 + Thf_offset, Thf_coeffs_flat[i].c0, n * sizeof(double)).wait();
            q.memcpy(d_Thf_c1 + Thf_offset, Thf_coeffs_flat[i].c1, n * sizeof(double)).wait();
            q.memcpy(d_Thf_c2 + Thf_offset, Thf_coeffs_flat[i].c2, n * sizeof(double)).wait();
            q.memcpy(d_Thf_c3 + Thf_offset, Thf_coeffs_flat[i].c3, n * sizeof(double)).wait();
            q.memcpy(d_Thf_x + Thf_x_offset, Thf_coeffs_flat[i].x, n * sizeof(double)).wait();  // Fixed: n not n+1
            
            temp_Thf_coeffs[i].n_points = n;
            temp_Thf_coeffs[i].c0 = d_Thf_c0 + Thf_offset;
            temp_Thf_coeffs[i].c1 = d_Thf_c1 + Thf_offset;
            temp_Thf_coeffs[i].c2 = d_Thf_c2 + Thf_offset;
            temp_Thf_coeffs[i].c3 = d_Thf_c3 + Thf_offset;
            temp_Thf_coeffs[i].x = d_Thf_x + Thf_x_offset;
            
            Thf_offset += n;
            Thf_x_offset += n;  // Fixed: n not n+1
        }
    }
    
    q.memcpy(d_Thf_coeffs_flat, temp_Thf_coeffs, Thf_coeffs_size * sizeof(SplineCoefficients)).wait();
    delete[] temp_Thf_coeffs;

    std::cout << "SYCL device memory transfer complete!" << std::endl;

}

void ThermoDynamicsF4::cleanup_spline_coefficients() {
    std::cout << "Cleaning up spline coefficients..." << std::endl;

    if(queue_ptr) {
        // Free device memory
        sycl::free(d_ES_coeffs_flat, *queue_ptr);
        sycl::free(d_Thf_coeffs_flat, *queue_ptr);
        sycl::free(d_ES_c0, *queue_ptr);
        sycl::free(d_ES_c1, *queue_ptr);
        sycl::free(d_ES_c2, *queue_ptr);
        sycl::free(d_ES_c3, *queue_ptr);
        sycl::free(d_ES_x, *queue_ptr);
        sycl::free(d_Thf_c0, *queue_ptr);
        sycl::free(d_Thf_c1, *queue_ptr);
        sycl::free(d_Thf_c2, *queue_ptr);
        sycl::free(d_Thf_c3, *queue_ptr);
        sycl::free(d_Thf_x, *queue_ptr);
    }

    // Free host memory
    std::set<void*> freed_pointers;

    int ES_coeffs_size = nPhase * (nComp - 1) * 2;
    int Thf_coeffs_size = nPhase * (nComp - 1) * (nComp - 1);

    for(int i = 0; i < ES_coeffs_size; ++i) {
        if (ES_coeffs_flat[i].c0 && freed_pointers.find(ES_coeffs_flat[i].c0) == freed_pointers.end()) {
            delete[] ES_coeffs_flat[i].c0;
            delete[] ES_coeffs_flat[i].c1;
            delete[] ES_coeffs_flat[i].c2;
            delete[] ES_coeffs_flat[i].c3;
            delete[] ES_coeffs_flat[i].x;
            freed_pointers.insert(ES_coeffs_flat[i].c0);
        }
    }

    for(int i = 0; i < Thf_coeffs_size; ++i) {
        if (Thf_coeffs_flat[i].c0 && freed_pointers.find(Thf_coeffs_flat[i].c0) == freed_pointers.end()) {
            delete[] Thf_coeffs_flat[i].c0;
            delete[] Thf_coeffs_flat[i].c1;
            delete[] Thf_coeffs_flat[i].c2;
            delete[] Thf_coeffs_flat[i].c3;
            delete[] Thf_coeffs_flat[i].x;
            freed_pointers.insert(Thf_coeffs_flat[i].c0);
        }
    }

    delete[] ES_coeffs_flat;
    delete[] Thf_coeffs_flat;

    // Free thermo_phase mapping
    if(thermo_phase && queue_ptr) {
        sycl::free(thermo_phase, *queue_ptr);
        thermo_phase = nullptr;
    }

    std::cout << "Cleanup complete." << std::endl;
}

} // namespace microsim