#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <memory>
#include <stdexcept>
#if __has_include("Thermo.h")
#include "Thermo.h"
#else
#include "../tdbs/Thermo.h"
#endif

namespace mcch {

class FreeEnergy {
public:
    virtual ~FreeEnergy() = default;
    virtual int num_components() const = 0;
    virtual double density(const std::vector<double>& c) const = 0;
    virtual void chemical_derivatives(const std::vector<double>& c, std::vector<double>& df_dc) const = 0;
};

// 1. Regular Solution Model:
// f0 = R*T * sum_i(c_i * ln(c_i + eps)) + sum_{i < j} (Omega_ij * c_i * c_j)
class RegularSolution : public FreeEnergy {
private:
    int n_comp;
    double RT;
    std::vector<std::vector<double>> omega; // interaction parameters Omega_ij
    double eps; // small regularization to avoid log(0)

public:
    RegularSolution(int num_c, double rt_param, const std::vector<std::vector<double>>& omega_matrix, double epsilon = 1e-12)
        : n_comp(num_c), RT(rt_param), omega(omega_matrix), eps(epsilon) {
        if ((int)omega.size() != n_comp) {
            throw std::invalid_argument("Omega matrix dimension does not match number of components.");
        }
    }

    int num_components() const override {
        return n_comp;
    }

    double get_RT() const { return RT; }
    double get_eps() const { return eps; }
    const std::vector<std::vector<double>>& get_omega() const { return omega; }

    double density(const std::vector<double>& c) const override {
        double f = 0.0;
        // Ideal mixing entropy term
        for (int i = 0; i < n_comp; ++i) {
            double ci = std::max(c[i], eps);
            f += RT * ci * std::log(ci);
        }
        // Excess enthalpy (interaction) term
        for (int i = 0; i < n_comp; ++i) {
            for (int j = i + 1; j < n_comp; ++j) {
                f += omega[i][j] * c[i] * c[j];
            }
        }
        return f;
    }

    void chemical_derivatives(const std::vector<double>& c, std::vector<double>& df_dc) const override {
        df_dc.resize(n_comp);
        for (int i = 0; i < n_comp; ++i) {
            double ci = std::max(c[i], eps);
            double d_entropy = RT * (std::log(ci) + 1.0);
            double d_enthalpy = 0.0;
            for (int j = 0; j < n_comp; ++j) {
                if (i != j) {
                    d_enthalpy += omega[i][j] * c[j];
                }
            }
            df_dc[i] = d_entropy + d_enthalpy;
        }
    }
};

// 2. Polynomial Multi-Well Model:
// f0 = sum_{i < j} W_ij * c_i^2 * c_j^2 + sum_i A_i * c_i^2 * (1 - c_i)^2
class PolynomialMultiWell : public FreeEnergy {
private:
    int n_comp;
    std::vector<std::vector<double>> W; // pairwise barrier heights
    std::vector<double> A; // self double-well barrier heights

public:
    PolynomialMultiWell(int num_c, const std::vector<std::vector<double>>& W_matrix, const std::vector<double>& A_vec = {})
        : n_comp(num_c), W(W_matrix), A(A_vec) {
        if ((int)W.size() != n_comp) {
            throw std::invalid_argument("W matrix dimension does not match number of components.");
        }
        if (A.empty()) {
            A.assign(n_comp, 0.0);
        }
    }

    int num_components() const override {
        return n_comp;
    }

    const std::vector<std::vector<double>>& get_W() const { return W; }
    const std::vector<double>& get_A() const { return A; }

    double density(const std::vector<double>& c) const override {
        double f = 0.0;
        for (int i = 0; i < n_comp; ++i) {
            for (int j = i + 1; j < n_comp; ++j) {
                f += W[i][j] * c[i] * c[i] * c[j] * c[j];
            }
        }
        for (int i = 0; i < n_comp; ++i) {
            if (A[i] != 0.0) {
                double term = c[i] * (1.0 - c[i]);
                f += A[i] * term * term;
            }
        }
        return f;
    }

    void chemical_derivatives(const std::vector<double>& c, std::vector<double>& df_dc) const override {
        df_dc.resize(n_comp);
        for (int i = 0; i < n_comp; ++i) {
            double d = 0.0;
            for (int j = 0; j < n_comp; ++j) {
                if (i != j) {
                    d += 2.0 * W[i][j] * c[i] * c[j] * c[j];
                }
            }
            if (A[i] != 0.0) {
                // Derivative of A_i * c_i^2 * (1 - c_i)^2:
                // = 2 * A_i * c_i * (1 - c_i) * (1 - 2*c_i)
                d += 2.0 * A[i] * c[i] * (1.0 - c[i]) * (1.0 - 2.0 * c[i]);
            }
            df_dc[i] = d;
        }
    }
};

// 3. Binary Double Well (for classic Cahn-Hilliard verification):
// f0(c) = W * c_0^2 * (1 - c_0)^2
class BinaryDoubleWell : public FreeEnergy {
private:
    double W;

public:
    BinaryDoubleWell(double barrier = 1.0) : W(barrier) {}

    int num_components() const override {
        return 2;
    }

    double get_barrier() const { return W; }
    double get_W() const { return W; }

    double density(const std::vector<double>& c) const override {
        double c0 = c[0];
        return W * c0 * c0 * (1.0 - c0) * (1.0 - c0);
    }

    void chemical_derivatives(const std::vector<double>& c, std::vector<double>& df_dc) const override {
        df_dc.resize(2);
        double c0 = c[0];
        df_dc[0] = 2.0 * W * c0 * (1.0 - c0) * (1.0 - 2.0 * c0);
        df_dc[1] = 0.0;
    }
};

// 4. Redlich-Kister CALPHAD Model:
// Evaluates molar Gibbs energy and chemical derivatives via CALPHAD database functions
class Redlich_Kister : public FreeEnergy {
private:
    int n_comp;
    double T;
public:
    Redlich_Kister(int num_c = 2, double Temp = 800.0)
        : n_comp(num_c), T(Temp) {
        if (n_comp != 2) {
            throw std::invalid_argument("Redlich_Kister thermodynamic model currently only supports binary systems (num_components = 2).");
        }
    }

    int num_components() const override {
        return n_comp;
    }

    double density(const std::vector<double>& c) const override {
        if (c.size() < 2) {
            throw std::invalid_argument("Redlich_Kister::density requires at least 2 components.");
        }
        double f = 0.0;
        GES(T, c.data(), &f);
        return f;
    }

    void chemical_derivatives(const std::vector<double>& c, std::vector<double>& df_dc) const override {
        if (c.size() < 2) {
            throw std::invalid_argument("Redlich_Kister::chemical_derivatives requires at least 2 components.");
        }
        df_dc.resize(n_comp);
        dGES(T, c.data(), df_dc.data());
    }

    double temperature() const {
        return T;
    }

    void set_temperature(double Temp) {
        T = Temp;
    }
};


} // namespace mcch
