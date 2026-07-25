import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


# Define the symbolic variables representing physical quantities
c = sp.symbols('c')           # Composition (independent variable, 0 to 1)
T = sp.symbols('T')           # Current Temperature
T_eq = sp.symbols('T_eq')     # Reference Equilibrium Temperature

# Define symbolic parameters for a Generic Phase
# A: Curvature
# Beq: Equilibrium linear coefficient
# dBdT: Slope of B with temperature
# C_off: Vertical energy offset
A_sym, Beq_sym, dBdT_sym, C_off_sym = sp.symbols('A B_eq dBdT C_offset')



# The temperature dependent coefficient B(T)
B_T_expr = Beq_sym + dBdT_sym * (T - T_eq)

# The Parabolic Gibbs Free Energy Density f(c, T)
f_sym_expr = A_sym * (c**2) + B_T_expr * c + C_off_sym

# The Chemical Potential mu(c, T)
# Calculated dynamically by taking the derivative of f with respect to c
mu_sym_expr = sp.diff(f_sym_expr, c)


print("--- Symbolic Thermodynamic Equations ---")
print(f"Free Energy f(c,T): {f_sym_expr}")
print(f"Chem. Potential mu(c,T) = df/dc: {mu_sym_expr}")
print("----------------------------------------\n")



# "Solid" (alpha) and "Liquid" (beta) with different properties.
Ref_Temp = 1000.0  # Reference temperature (e.g., Melting Point)

# Parameters for Solid Phase (Alpha)
# Narrow parabola (high A), stable at lower T
params_solid = {
    A_sym: 80000,        # High curvature (narrow well)
    Beq_sym: 80000,     # Positions minimum near c=0.5 at T_eq
    dBdT_sym: 0,        # Positive entropy effect
    C_off_sym: 0,    # Vertical shift
    T_eq: Ref_Temp
}

# Parameters for Liquid Phase (Beta)
# Wider parabola (lower A), stable at higher T
params_liquid = {
    A_sym: 40000,        # Lower curvature (wider well)
    Beq_sym: 40000,     # Positions minimum near c=0.5 at T_eq
    dBdT_sym: -0,       # Negative entropy effect relative to solid
    C_off_sym: 0,    # Vertical shift
    T_eq: Ref_Temp
}

# Substitute the specific phase parameters into the general symbolic expressions
f_solid_sym = f_sym_expr.subs(params_solid)
mu_solid_sym = mu_sym_expr.subs(params_solid)

f_liquid_sym = f_sym_expr.subs(params_liquid)
mu_liquid_sym = mu_sym_expr.subs(params_liquid)

# Convert symbolic expressions to highly efficient numpy-callable functions
# These functions now take (c, T) as numerical inputs.
f_solid_func = sp.lambdify((c, T), f_solid_sym, 'numpy')
mu_solid_func = sp.lambdify((c, T), mu_solid_sym, 'numpy')

f_liquid_func = sp.lambdify((c, T), f_liquid_sym, 'numpy')
mu_liquid_func = sp.lambdify((c, T), mu_liquid_sym, 'numpy')



# Define composition range for plotting
c_range = np.linspace(0.01, 0.99, 500)

# Define temperatures to visualize: Below, At, and Above Reference Temperature
temperatures = [Ref_Temp - 150, Ref_Temp, Ref_Temp + 150]
temp_labels = ["Solid Stable", "Equilibrium", "Liquid Stable"]

fig, axes = plt.subplots(3, 2, figsize=(14, 15))
fig.suptitle("Parabolic Free Energy Model", fontsize=16, y=0.99)

for i, current_T in enumerate(temperatures):
    # Calculate numerical data for this temperature
    f_s_data = f_solid_func(c_range, current_T)
    f_l_data = f_liquid_func(c_range, current_T)
    mu_s_data = mu_solid_func(c_range, current_T)
    mu_l_data = mu_liquid_func(c_range, current_T)

    # --- Left Column: Free Energy Curves (f vs c) ---
    ax_f = axes[i, 0]
    ax_f.plot(c_range, f_s_data, 'b-', linewidth=2, label='Solid Phase ($f_S$)')
    ax_f.plot(c_range, f_l_data, 'r--', linewidth=2, label='Liquid Phase ($f_L$)')
    
    # Visual aid: Highlight minimums
    min_s_idx = np.argmin(f_s_data)
    min_l_idx = np.argmin(f_l_data)
    ax_f.scatter(c_range[min_s_idx], f_s_data[min_s_idx], color='blue', marker='o')
    ax_f.scatter(c_range[min_l_idx], f_l_data[min_l_idx], color='red', marker='o')

    ax_f.set_title(f"Free Energy Density at T = {current_T:.0f} K ({temp_labels[i]})")
    ax_f.set_ylabel("Free Energy ($J/m^3$)")
    ax_f.set_xlabel("Composition ($c$)")
    ax_f.legend()
    ax_f.grid(True, alpha=0.3)

    # --- Right Column: Chemical Potential Curves (mu vs c) ---
    # Note: For a parabolic f, mu is exactly linear with c.
    ax_mu = axes[i, 1]
    ax_mu.plot(c_range, mu_s_data, 'b-', linewidth=2, label='Solid Chem Pot ($\mu_S$)')
    ax_mu.plot(c_range, mu_l_data, 'r--', linewidth=2, label='Liquid Chem Pot ($\mu_L$)')

    # driving forces for diffusion balance *if* an interface existed.
    # (KKS model condition: mu_S = mu_L)
    idx = np.argwhere(np.diff(np.sign(mu_s_data - mu_l_data))).flatten()
    if len(idx) > 0:
         ax_mu.plot(c_range[idx], mu_s_data[idx], 'ko', markersize=8, label='$\mu_S = \mu_L$')

    ax_mu.set_title(f"Chemical Potential ($\partial f / \partial c$) at T = {current_T:.0f} K")
    ax_mu.set_ylabel("Chemical Potential ($\mu$)")
    ax_mu.set_xlabel("Composition ($c$)")
    ax_mu.legend()
    ax_mu.grid(True, alpha=0.3)

plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.show()