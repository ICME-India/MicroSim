#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <memory>

namespace mcch {

enum class MobilityType {
    CONSTANT,
    DEGENERATE,
    MATRIX
};

class MobilityModel {
public:
    virtual ~MobilityModel() = default;
    virtual MobilityType type() const = 0;
    virtual int num_components() const = 0;
    virtual bool is_constant() const = 0;

    // Evaluate mobility matrix M_ij for given concentration vector c
    // Output: flattened (n_comp x n_comp) matrix
    virtual void evaluate_matrix(const std::vector<double>& c, std::vector<double>& M_out) const = 0;
};

// 1. Constant scalar or matrix mobility: M_ij = M0 * delta_ij or full matrix M_ij
class ConstantMobility : public MobilityModel {
protected:
    int n_comp;
    double M0;
    std::vector<std::vector<double>> M_matrix;
    std::vector<double> M_flat;
    bool is_diag;

public:
    // Constructor with scalar mobility (e.g. for binary systems)
    ConstantMobility(int num_c, double mobility = 1.0)
        : n_comp(num_c), M0(mobility), is_diag(true) {
        M_matrix.assign(n_comp, std::vector<double>(n_comp, 0.0));
        M_flat.assign(n_comp * n_comp, 0.0);
        for (int i = 0; i < n_comp; ++i) {
            M_matrix[i][i] = M0;
            M_flat[i * n_comp + i] = M0;
        }
    }

    // Constructor with full mobility matrix (for multicomponent systems)
    ConstantMobility(int num_c, const std::vector<std::vector<double>>& matrix)
        : n_comp(num_c), is_diag(true) {
        M_matrix.assign(n_comp, std::vector<double>(n_comp, 0.0));
        M_flat.assign(n_comp * n_comp, 0.0);
        M0 = 1.0;
        for (int i = 0; i < n_comp; ++i) {
            for (int j = 0; j < n_comp; ++j) {
                double val = (i < (int)matrix.size() && j < (int)matrix[i].size()) ? matrix[i][j] : 0.0;
                M_matrix[i][j] = val;
                M_flat[i * n_comp + j] = val;
                if (i != j && std::abs(val) > 1e-14) {
                    is_diag = false;
                }
            }
        }
        if (n_comp > 0) {
            M0 = M_matrix[0][0];
        }
    }

    MobilityType type() const override { return is_diag ? MobilityType::CONSTANT : MobilityType::MATRIX; }
    int num_components() const override { return n_comp; }
    bool is_constant() const override { return true; }
    double value() const { return M0; }
    double value(int i, int j) const { return M_matrix[i][j]; }
    const std::vector<std::vector<double>>& matrix() const { return M_matrix; }
    const std::vector<double>& flat_matrix() const { return M_flat; }
    bool is_diagonal() const { return is_diag; }

    void evaluate_matrix(const std::vector<double>& c, std::vector<double>& M_out) const override {
        (void)c;
        M_out = M_flat;
    }
};

// 2. Degenerate Onsager mobility: M_ij = M0 * (c_i * delta_ij - c_i * c_j)
// Guarantees zero flux at c_i = 0 and sum_i M_ij = 0
class DegenerateMobility : public MobilityModel {
private:
    int n_comp;
    double M0;

public:
    DegenerateMobility(int num_c, double mobility = 1.0)
        : n_comp(num_c), M0(mobility) {}

    MobilityType type() const override { return MobilityType::DEGENERATE; }
    int num_components() const override { return n_comp; }
    bool is_constant() const override { return false; }
    double base_value() const { return M0; }

    void evaluate_matrix(const std::vector<double>& c, std::vector<double>& M_out) const override {
        M_out.resize(n_comp * n_comp);
        for (int i = 0; i < n_comp; ++i) {
            double ci = std::max(0.0, std::min(1.0, c[i]));
            for (int j = 0; j < n_comp; ++j) {
                double cj = std::max(0.0, std::min(1.0, c[j]));
                if (i == j) {
                    M_out[i * n_comp + j] = M0 * ci * (1.0 - ci);
                } else {
                    M_out[i * n_comp + j] = -M0 * ci * cj;
                }
            }
        }
    }
};

// 3. User-defined constant mobility matrix (n_comp x n_comp)
class MatrixMobility : public ConstantMobility {
public:
    MatrixMobility(int num_c, const std::vector<std::vector<double>>& matrix)
        : ConstantMobility(num_c, matrix) {}

    MobilityType type() const override { return MobilityType::MATRIX; }
};

} // namespace mcch
