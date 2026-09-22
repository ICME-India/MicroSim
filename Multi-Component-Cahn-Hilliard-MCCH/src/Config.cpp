#include "Config.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <algorithm>

namespace mcch {

// Simple lightweight JSON tokenizer/parser
namespace json_simple {

enum TokenType { TOK_LBRACE, TOK_RBRACE, TOK_LBRACKET, TOK_RBRACKET, TOK_COLON, TOK_COMMA, TOK_STRING, TOK_NUMBER, TOK_BOOL, TOK_NULL, TOK_EOF };

struct Token {
    TokenType type;
    std::string str_val;
    double num_val;
};

class Lexer {
    std::string src;
    size_t pos = 0;
public:
    Lexer(const std::string& input) : src(input) {}

    Token next() {
        while (pos < src.size() && (std::isspace(src[pos]) || src[pos] == '\r' || src[pos] == '\n')) {
            pos++;
        }
        if (pos >= src.size()) return {TOK_EOF, "", 0.0};

        char c = src[pos];
        if (c == '{') { pos++; return {TOK_LBRACE, "{", 0.0}; }
        if (c == '}') { pos++; return {TOK_RBRACE, "}", 0.0}; }
        if (c == '[') { pos++; return {TOK_LBRACKET, "[", 0.0}; }
        if (c == ']') { pos++; return {TOK_RBRACKET, "]", 0.0}; }
        if (c == ':') { pos++; return {TOK_COLON, ":", 0.0}; }
        if (c == ',') { pos++; return {TOK_COMMA, ",", 0.0}; }

        if (c == '"') {
            pos++;
            std::string s;
            while (pos < src.size() && src[pos] != '"') {
                if (src[pos] == '\\' && pos + 1 < src.size()) {
                    pos++;
                }
                s += src[pos++];
            }
            if (pos < src.size()) pos++; // skip closing quote
            return {TOK_STRING, s, 0.0};
        }

        if (std::isdigit(c) || c == '-' || c == '+') {
            size_t start = pos++;
            while (pos < src.size() && (std::isdigit(src[pos]) || src[pos] == '.' || src[pos] == 'e' || src[pos] == 'E' || src[pos] == '-' || src[pos] == '+')) {
                pos++;
            }
            std::string num_str = src.substr(start, pos - start);
            double val = std::stod(num_str);
            return {TOK_NUMBER, num_str, val};
        }

        // Keywords (true, false, null)
        size_t start = pos;
        while (pos < src.size() && std::isalpha(src[pos])) pos++;
        std::string kw = src.substr(start, pos - start);
        if (kw == "true" || kw == "false") return {TOK_BOOL, kw, kw == "true" ? 1.0 : 0.0};
        return {TOK_NULL, kw, 0.0};
    }
};

} // namespace json_simple

void SimulationConfig::load_from_json(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open JSON config file: " + filepath);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string text = buffer.str();

    json_simple::Lexer lexer(text);
    json_simple::Token tok = lexer.next();
    if (tok.type != json_simple::TOK_LBRACE) {
        throw std::runtime_error("Invalid JSON: root must be an object");
    }

    while (true) {
        tok = lexer.next();
        if (tok.type == json_simple::TOK_RBRACE || tok.type == json_simple::TOK_EOF) break;
        if (tok.type == json_simple::TOK_COMMA) continue;

        if (tok.type != json_simple::TOK_STRING) {
            throw std::runtime_error("Expected string key in JSON object");
        }
        std::string key = tok.str_val;

        tok = lexer.next();
        if (tok.type != json_simple::TOK_COLON) {
            throw std::runtime_error("Expected ':' after key " + key);
        }

        tok = lexer.next();

        if (key == "dim") dim = (int)tok.num_val;
        else if (key == "nx") nx = (int)tok.num_val;
        else if (key == "ny") ny = (int)tok.num_val;
        else if (key == "nz") nz = (int)tok.num_val;
        else if (key == "dx") dx = tok.num_val;
        else if (key == "dy") dy = tok.num_val;
        else if (key == "dz") dz = tok.num_val;
        else if (key == "bc_x") bc_x = tok.str_val;
        else if (key == "bc_y") bc_y = tok.str_val;
        else if (key == "bc_z") bc_z = tok.str_val;
        else if (key == "num_components") num_components = (int)tok.num_val;
        else if (key == "free_energy_type") free_energy_type = tok.str_val;
        else if (key == "RT") RT = tok.num_val;
        else if (key == "mobility_type") mobility_type = tok.str_val;
        else if (key == "mobility_val" || key == "mobility" || key == "mobility_matrix" || key == "M") {
            if (tok.type == json_simple::TOK_NUMBER) {
                mobility_val = tok.num_val;
                mobility_matrix.assign(num_components, std::vector<double>(num_components, 0.0));
                for (int i = 0; i < num_components; ++i) mobility_matrix[i][i] = tok.num_val;
            } else if (tok.type == json_simple::TOK_LBRACKET) {
                mobility_matrix.clear();
                while (true) {
                    tok = lexer.next();
                    if (tok.type == json_simple::TOK_RBRACKET) break;
                    if (tok.type == json_simple::TOK_COMMA) continue;
                    if (tok.type == json_simple::TOK_LBRACKET) {
                        std::vector<double> row;
                        while (true) {
                            tok = lexer.next();
                            if (tok.type == json_simple::TOK_RBRACKET) break;
                            if (tok.type == json_simple::TOK_COMMA) continue;
                            row.push_back(tok.num_val);
                        }
                        mobility_matrix.push_back(row);
                    }
                }
                if (!mobility_matrix.empty() && !mobility_matrix[0].empty()) {
                    mobility_val = mobility_matrix[0][0];
                }
            }
        }
        else if (key == "integrator") integrator = tok.str_val;
        else if (key == "device") device = tok.str_val;
        else if (key == "use_gpu") device = (tok.num_val > 0.5 || tok.str_val == "true") ? "gpu" : "cpu";
        else if (key == "dt") dt = tok.num_val;
        else if (key == "stabilization") stabilization = tok.num_val;
        else if (key == "total_steps") total_steps = (int)tok.num_val;
        else if (key == "output_interval") output_interval = (int)tok.num_val;
        else if (key == "diag_interval") diag_interval = (int)tok.num_val;
        else if (key == "initial_condition") initial_condition = tok.str_val;
        else if (key == "noise_amp") noise_amp = tok.num_val;
        else if (key == "seed") seed = (unsigned int)tok.num_val;
        else if (key == "droplet_comp") droplet_comp = (int)tok.num_val;
        else if (key == "droplet_radius") droplet_radius = tok.num_val;
        else if (key == "droplet_diffuse_width") droplet_diffuse_width = tok.num_val;
        else if (key == "Temperature" || key == "temperature" || key == "T") T = tok.num_val;
        else if (key == "output_dir") output_dir = tok.str_val;
        else if (key == "c_mean" && tok.type == json_simple::TOK_LBRACKET) {
            c_mean.clear();
            while (true) {
                tok = lexer.next();
                if (tok.type == json_simple::TOK_RBRACKET) break;
                if (tok.type == json_simple::TOK_COMMA) continue;
                c_mean.push_back(tok.num_val);
            }
        }
        else if (key == "A" && tok.type == json_simple::TOK_LBRACKET) {
            A.clear();
            while (true) {
                tok = lexer.next();
                if (tok.type == json_simple::TOK_RBRACKET) break;
                if (tok.type == json_simple::TOK_COMMA) continue;
                A.push_back(tok.num_val);
            }
        }
        else if (key == "W" && tok.type == json_simple::TOK_LBRACKET) {
            W.clear();
            while (true) {
                tok = lexer.next();
                if (tok.type == json_simple::TOK_RBRACKET) break;
                if (tok.type == json_simple::TOK_COMMA) continue;
                if (tok.type == json_simple::TOK_LBRACKET) {
                    std::vector<double> row;
                    while (true) {
                        tok = lexer.next();
                        if (tok.type == json_simple::TOK_RBRACKET) break;
                        if (tok.type == json_simple::TOK_COMMA) continue;
                        row.push_back(tok.num_val);
                    }
                    W.push_back(row);
                }
            }
        }
        else if (key == "omega" && tok.type == json_simple::TOK_LBRACKET) {
            omega.clear();
            while (true) {
                tok = lexer.next();
                if (tok.type == json_simple::TOK_RBRACKET) break;
                if (tok.type == json_simple::TOK_COMMA) continue;
                if (tok.type == json_simple::TOK_LBRACKET) {
                    std::vector<double> row;
                    while (true) {
                        tok = lexer.next();
                        if (tok.type == json_simple::TOK_RBRACKET) break;
                        if (tok.type == json_simple::TOK_COMMA) continue;
                        row.push_back(tok.num_val);
                    }
                    omega.push_back(row);
                }
            }
        }
        else if (key == "kappa") {
            if (tok.type == json_simple::TOK_NUMBER) {
                // Scalar kappa: will be expanded to diagonal matrix in finalize_defaults
                kappa.assign(num_components, std::vector<double>(num_components, 0.0));
                for (int i = 0; i < num_components; ++i) kappa[i][i] = tok.num_val;
            } else if (tok.type == json_simple::TOK_LBRACKET) {
                kappa.clear();
                while (true) {
                    tok = lexer.next();
                    if (tok.type == json_simple::TOK_RBRACKET) break;
                    if (tok.type == json_simple::TOK_COMMA) continue;
                    if (tok.type == json_simple::TOK_LBRACKET) {
                        std::vector<double> row;
                        while (true) {
                            tok = lexer.next();
                            if (tok.type == json_simple::TOK_RBRACKET) break;
                            if (tok.type == json_simple::TOK_COMMA) continue;
                            row.push_back(tok.num_val);
                        }
                        kappa.push_back(row);
                    }
                }
            }
        }
        else {
            // Skip unrecognized value / object
            if (tok.type == json_simple::TOK_LBRACKET) {
                int depth = 1;
                while (depth > 0) {
                    tok = lexer.next();
                    if (tok.type == json_simple::TOK_LBRACKET) depth++;
                    if (tok.type == json_simple::TOK_RBRACKET) depth--;
                    if (tok.type == json_simple::TOK_EOF) break;
                }
            } else if (tok.type == json_simple::TOK_LBRACE) {
                int depth = 1;
                while (depth > 0) {
                    tok = lexer.next();
                    if (tok.type == json_simple::TOK_LBRACE) depth++;
                    if (tok.type == json_simple::TOK_RBRACE) depth--;
                    if (tok.type == json_simple::TOK_EOF) break;
                }
            }
        }
    }

    finalize_defaults();
}

void SimulationConfig::finalize_defaults() {
    if (dim == 2) {
        nz = 1;
        dz = 1.0;
    }

    if (c_mean.empty() || c_mean.size() != (size_t)num_components) {
        c_mean.assign(num_components, 1.0 / num_components);
    }

    if (kappa.empty() || kappa.size() != (size_t)num_components) {
        double val = (!kappa.empty() && !kappa[0].empty()) ? kappa[0][0] : 1.0;
        kappa.assign(num_components, std::vector<double>(num_components, 0.0));
        for (int i = 0; i < num_components; ++i) {
            kappa[i][i] = val;
        }
    }

    if (W.empty() || W.size() != (size_t)num_components) {
        W.assign(num_components, std::vector<double>(num_components, 0.0));
        for (int i = 0; i < num_components; ++i) {
            for (int j = 0; j < num_components; ++j) {
                if (i != j) W[i][j] = 1.0;
            }
        }
    }

    if (omega.empty() || omega.size() != (size_t)num_components) {
        omega.assign(num_components, std::vector<double>(num_components, 0.0));
        for (int i = 0; i < num_components; ++i) {
            for (int j = 0; j < num_components; ++j) {
                if (i != j) omega[i][j] = 3.0;
            }
        }
    }

    if (A.empty() || A.size() != (size_t)num_components) {
        A.assign(num_components, 0.0);
    }

    if (mobility_matrix.empty() || mobility_matrix.size() != (size_t)num_components) {
        double val = (!mobility_matrix.empty() && !mobility_matrix[0].empty()) ? mobility_matrix[0][0] : mobility_val;
        mobility_matrix.assign(num_components, std::vector<double>(num_components, 0.0));
        for (int i = 0; i < num_components; ++i) {
            mobility_matrix[i][i] = val;
        }
        mobility_val = val;
    }
}

} // namespace mcch
