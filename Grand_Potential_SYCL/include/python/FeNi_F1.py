import numpy as np
import matplotlib.pyplot as plt


class ThermoDynamics:
    def __init__(self, n_comp, n_phase, T_eq_ref):
        self.nComp = n_comp
        self.nPhase = n_phase
        self.T_eq = T_eq_ref
        
        # Dimensions for matrices
        # A: [phase][comp][comp] (Interaction/Curvature)
        self.A = np.zeros((n_phase, n_comp-1, n_comp-1))
        
        # Equilibrium data containers
        # ceq[a][b][k]: Conc of comp k in phase b, when in equilibrium with phase a
        self.ceq = np.zeros((n_phase, n_phase, n_comp-1))
        
        # slopes[a][b][k]: Slope dT/dc of phase b boundary with phase a
        self.slopes = np.zeros((n_phase, n_phase, n_comp-1))
        
        # Computed Parameters
        self.Beq = np.zeros((n_phase, n_comp-1))
        self.B = np.zeros((n_phase, n_comp-1))
        self.dBbdT = np.zeros((n_phase, n_comp-1))
        self.C_const = np.zeros(n_phase)
        self.DELTA_C = np.zeros((n_phase, n_comp-1))
        self.DELTA_T = np.zeros((n_phase, n_phase)) # [a][b]
        self.Mu = np.zeros((n_phase, n_comp-1))
        
        # Linear Algebra cache (Inverse of A matrix)
        self.cmu = np.zeros((n_phase, n_comp-1, n_comp-1))

    # --- Loading Data ---
    def set_A(self, phase_idx, val):
        # For binary (1 independent comp), A is scalar.
        self.A[phase_idx, 0, 0] = val

    def set_ceq(self, phase_a, phase_b, val):
        self.ceq[phase_a, phase_b, 0] = val

    def set_slope(self, phase_a, phase_b, val):
        self.slopes[phase_a, phase_b, 0] = val


    def initialize_host(self):
        # NOTE: ISOTHERMAL CASE
        # Invert A matrices (muc -> cmu)
        # NOTE: calculates 'muc' then inverts it. 
        # Free energy F = A*ci*cj + B*cj + C. 
        # chemical potential \frac{\partial F}{\parial Ci} mu = 2*A*cj + B  for(diagonal elements) else A*cj + B.
        # dmu/dc = 2*A (for diagonal elements) else A.
        
        for a in range(self.nPhase):
            muc = np.zeros((self.nComp-1, self.nComp-1))
            for i in range(self.nComp-1):
                for j in range(self.nComp-1):
                    if i == j:
                        muc[i,j] = 2.0 * self.A[a, i, j] # Hessian is 2*A
                    else:
                        muc[i,j] = self.A[a, i, j]
            
            # Invert to get cmu (susceptibility)
            try:
                self.cmu[a] = np.linalg.inv(muc)
            except:
                self.cmu[a] = 0.0

        # 2. Calculate DELTA_T, DELTA_C, dBbdT
        ref_phase = self.nPhase - 1 # Liquid
        
        for a in range(self.nPhase):
            # Reset
            self.DELTA_T[a, ref_phase] = 0.0
            self.DELTA_T[a, a] = 0.0
            
            # Calculate DELTA parameters based on slopes
            for k in range(self.nComp-1):
                c_liq_eq = self.ceq[a, ref_phase, k]
                c_sol_eq = self.ceq[a, a, k]
                slope_liq = self.slopes[a, ref_phase, k]
                slope_sol = self.slopes[a, a, k]
                
                delta_c = c_liq_eq - c_sol_eq
                
                self.DELTA_C[a, k] = delta_c
                self.DELTA_T[a, a]         += slope_sol * (-delta_c)
                self.DELTA_T[a, ref_phase] += slope_liq * (-delta_c) 
                

            # Calculate dB/dT (Entropy difference)
            dcbdT_sol = np.zeros(self.nComp-1)
            dcbdT_liq = np.zeros(self.nComp-1)
            
            for k in range(self.nComp-1):
                if abs(self.DELTA_T[a, a]) > 1e-12:
                    dcbdT_sol[k] = self.DELTA_C[a, k] / self.DELTA_T[a, a]
                if abs(self.DELTA_T[a, ref_phase]) > 1e-12:
                    dcbdT_liq[k] = self.DELTA_C[a, k] / self.DELTA_T[a, ref_phase]

                # Formula for dBbdT
                # dBbdT = 2 * (A_liq * dcbdT_liq - A_sol * dcbdT_sol)
                term_liq = self.A[ref_phase, k, k] * dcbdT_liq[k]
                term_sol = self.A[a, k, k] * dcbdT_sol[k]
                self.dBbdT[a, k] = 2.0 * (term_liq - term_sol)

                for i in range(self.nComp-1):
                    if i != k:
                        term_liq = self.A[ref_phase, k, i] * dcbdT_liq[i]
                        term_sol = self.A[a, k, i] * dcbdT_sol[i]
                        self.dBbdT[a, k] += (term_liq - term_sol)

        # 3. Initialize Beq (at T_eq)
        for a in range(self.nPhase):
            for i in range(self.nComp-1):
                self.Beq[a, i] = self.functionB(self.T_eq, i, a, use_Beq_only=True)
            self.C_const[a] = self.functionC(self.T_eq, a)

    
    def functionB(self, T, i, a, use_Beq_only=False):
        # In runtime, B = Beq + dBdT * dT
        if not use_Beq_only:
            # Runtime update (linear approximation)
             return self.Beq[a, i] + self.dBbdT[a, i] * (T - self.T_eq)
        
        # Calculation from Phase Diagram (The "Construction" step)
        if a == (self.nPhase - 1): # Liquid
            return 0.0 # Reference
        
        ref_phase = self.nPhase - 1
        c_liq = np.zeros(self.nComp-1)
        c_sol = np.zeros(self.nComp-1)

        sum_c = 0.0
        
        # Interpolate equilibrium C based on T (Linear Phase Diagram approx)
        for k in range(self.nComp-1):
            if abs(self.DELTA_T[a, ref_phase]) > 1e-9:
                c_liq[k] = self.ceq[a, ref_phase, k] - (self.DELTA_C[a, k] * (self.T_eq - T) / self.DELTA_T[a, ref_phase])
            else: 
                c_liq[k] = self.ceq[a, ref_phase, k]
                
            if abs(self.DELTA_T[a, a]) > 1e-9:
                c_sol[k] = self.ceq[a, a, k] - (self.DELTA_C[a, k] * (self.T_eq - T) / self.DELTA_T[a, a])
            else:
                c_sol[k] = self.ceq[a, a, k]

            if k != i:
                sum_c += self.A[ref_phase, i, i] * c_liq[i] - self.A[a, i, i] * c_sol[i]
            
        
        # Calculate B_ai
        # B = 2*(A_liq*c_liq - A_sol*c_sol) for binary
        term_liq = self.A[ref_phase, i, i] * c_liq[i]
        term_sol = self.A[a, i, i] * c_sol[i]
        
        return 2.0 * (term_liq - term_sol) + sum_c


    def functionC(self, T, a):
        # Shifts the curve vertically to match Grand Potentials
        if a == (self.nPhase - 1): return 0.0
        
        ref_phase = self.nPhase - 1
        c_liq = np.zeros(self.nComp-1)
        c_sol = np.zeros(self.nComp-1)

        sum_c = 0.0
        
        # Same interpolation
        for k in range(self.nComp-1):
            if abs(self.DELTA_T[a, ref_phase]) > 1e-9:
                c_liq[k] = self.ceq[a, ref_phase, k] - (self.DELTA_C[a, k] * (self.T_eq - T) / self.DELTA_T[a, ref_phase])
            else: c_liq[k] = self.ceq[a, ref_phase, k]
            
            if abs(self.DELTA_T[a, a]) > 1e-9:
                c_sol[k] = self.ceq[a, a, k] - (self.DELTA_C[a, k] * (self.T_eq - T) / self.DELTA_T[a, a])
            else: c_sol[k] = self.ceq[a, a, k]

            for i in range(self.nComp-1):
                for j in range(self.nComp-1):
                    if i <= j:
                        sum_c += (self.A[a, i, j] * c_sol[i] * c_sol[j] - self.A[ref_phase, i, j] * c_liq[i] * c_liq[j])

        return sum_c

    # --- Free Energy Calculation ---
    def get_free_energy(self, c_val, T, phase_idx):
        # Update parameters for this T
        sum = 0.0
        for i in range(self.nComp-1):
            for j in range(self.nComp-1):
                if i <= j:
                    sum += self.A[phase_idx, i, j] * c_val[i] * c_val[j]
        
            self.B[phase_idx, i] = self.Beq[phase_idx, i] + self.dBbdT[phase_idx, i] * (T - self.T_eq)
            sum += self.B[phase_idx, i] * c_val[i]
        
        sum += self.functionC(T, phase_idx)
        return sum


    def get_chemical_potential(self, c_val, T, phase_idx):
        sum = 0.0
        for k in range(self.nComp-1):
            sum += 2.0 * self.A[phase_idx, k, k] * c_val[k] + self.Beq[phase_idx, k] + self.dBbdT[phase_idx, k] * (T - self.T_eq)
            for i in range(self.nComp-1):
                if i != k:
                    sum += self.A[phase_idx, k, i] * c_val[i]
        self.Mu[phase_idx] = sum
        return sum

    def get_grand_potential(self, c_val, T, phase_idx):
        f = self.get_free_energy(c_val, T, phase_idx)
        mu_term = 0.0
        for i in range(self.nComp-1):
            mu_term += self.get_chemical_potential(c_val, T, phase_idx) * c_val[i]

        return f - mu_term

# ==========================================
#               Ni-Fe SYSTEM
# ==========================================
# System Constants
T_EQ_REF = 1768.0 # Reference Temp (e.g. Pure Fe melting or Alloy ref)
N_COMP = 2
N_PHASE = 3 # 0:FCC, 1:BCC, 2:LIQ

td = ThermoDynamics(N_COMP, N_PHASE, T_EQ_REF)

# --- 1. Load A (Curvature) ---
td.set_A(0, 196602.716)   # FCC
td.set_A(1, 213864.84785) # BCC
td.set_A(2, 155449.2845)  # LIQ

# --- 2. Load ceq (Equilibrium Concentrations) ---
# Pair 0-2 (FCC-Liq)
td.set_ceq(0, 0, 0.043456)    # C_sol_fcc
td.set_ceq(0, 2, 0.04905849)  # C_liq (w.r.t fcc)

# Pair 1-2 (BCC-Liq)
td.set_ceq(1, 1, 0.03825473)  # C_sol_bcc
td.set_ceq(1, 2, 0.04905849)  # C_liq (w.r.t bcc)

# Pair 2 (Liquid reference)
td.set_ceq(2, 0, 0.04905849)
td.set_ceq(2, 1, 0.04905849)
td.set_ceq(2, 2, 0.04905849)

# --- 3. Load Slopes ---
# Pair 0-2 (FCC-Liq)
td.set_slope(0, 0, -244.62189) # Slope Solidus
td.set_slope(0, 2, -217.7516)  # Slope Liquidus

# Pair 1-2 (BCC-Liq)
td.set_slope(1, 1, -618.8046)  # Slope Solidus
td.set_slope(1, 2, -466.03965) # Slope Liquidus

# Redundant liquid entries
td.set_slope(2, 0, -217.7516)
td.set_slope(2, 1, -466.03965)

# --- 4. Initialize ---
td.initialize_host()

# ==========================================
#               VISUALIZATION
# ==========================================
# Define Plotting Range
c_range = np.linspace(0.00, 0.08, 500) # Zoom in on the relevant low C range (0 to 10%)
plot_temps = [T_EQ_REF]
labels = ["FCC (0)", "BCC (1)", "LIQUID (2)"]
colors = ['b', 'g', 'r']

fig, axes = plt.subplots(3, 3, figsize=(18, 15))
fig.suptitle(f"Thermodynamics of Ni-Fe System  (T_eq = {T_EQ_REF} K)", fontsize=14)

cols = ["Free Energy Density ($f$)", "Chemical Potential ($\mu$)", "Grand Potential ($\Psi$)"]

for i, T in enumerate(plot_temps):
    ax_row = axes[i]
    
    # Storage for finding intersections
    mus = []
    psis = []

    for p in range(N_PHASE):
        f_vals = []
        mu_vals = []
        psi_vals = []
        
        for c in c_range:
            c_vec = np.array([c])
            f_vals.append(td.get_free_energy(c_vec, T, p))
            mu_vals.append(td.get_chemical_potential(c_vec, T, p))
            psi_vals.append(td.get_grand_potential(c_vec, T, p))
            
        f_vals = np.array(f_vals)
        mu_vals = np.array(mu_vals)
        psi_vals = np.array(psi_vals)
        
        mus.append(mu_vals)
        psis.append(psi_vals)
        
        # Plot 1: Free Energy
        ax_row[0].plot(c_range, f_vals, color=colors[p], label=labels[p], linewidth=2)
        
        # Plot 2: Chemical Potential
        ax_row[1].plot(c_range, mu_vals, color=colors[p], label=labels[p], linestyle='--')
        
        # Plot 3: Legendre (Psi vs Mu)
        ax_row[2].plot(mu_vals, psi_vals, color=colors[p], label=labels[p])

    # Decoration
    ax_row[0].set_ylabel(f"T = {T:.0f} K\nEnergy ($J/m^3$)")
    ax_row[0].legend()
    ax_row[0].set_ylim(bottom=min(td.get_free_energy(np.array([0.05]), T, 2) - 1000, 0)) # Auto scale helper
    
    ax_row[1].set_ylabel("Chem Pot ($\mu$)")
    
    ax_row[2].set_ylabel("Grand Pot ($\Psi$)")
    ax_row[2].set_xlabel("Chem Pot ($\mu$)")
    
    # Mark Equilibrium (Intersection in Legendre plot)
    # Check intersection between FCC(0) and LIQ(2)
    # 1. Find overlapping range of mu
    mu_min = max(np.min(mus[0]), np.min(mus[2]))
    mu_max = min(np.max(mus[0]), np.max(mus[2]))
    
    if mu_max > mu_min:
        common_mu = np.linspace(mu_min, mu_max, 1000)
        
        # Ensure sorted for interp
        sort_idx_0 = np.argsort(mus[0])
        sort_idx_2 = np.argsort(mus[2])
        
        psi_0_int = np.interp(common_mu, mus[0][sort_idx_0], psis[0][sort_idx_0])
        psi_2_int = np.interp(common_mu, mus[2][sort_idx_2], psis[2][sort_idx_2])
        
        idx = np.argwhere(np.diff(np.sign(psi_0_int - psi_2_int))).flatten()
        if len(idx) > 0:
            mu_eq = common_mu[idx[0]]
            psi_eq = psi_0_int[idx[0]]
            ax_row[2].plot(mu_eq, psi_eq, 'k*', markersize=12, label='FCC/LIQ Eq')
            
            # Optional: Project back to Chemical Potential Plot
            c_eq_fcc = np.interp(mu_eq, mus[0][sort_idx_0], c_range[sort_idx_0])
            c_eq_liq = np.interp(mu_eq, mus[2][sort_idx_2], c_range[sort_idx_2])
            
            # Draw tie-line on Chem Pot plot
            ax_row[1].axhline(mu_eq, color='k', linestyle=':', alpha=0.4)
            ax_row[1].plot([c_eq_fcc, c_eq_liq], [mu_eq, mu_eq], 'k-o', markersize=5, alpha=0.6)
            ax_row[1].text(c_eq_fcc, mu_eq, f" {c_eq_fcc:.4f}", fontsize=8, verticalalignment='bottom')
            ax_row[1].text(c_eq_liq, mu_eq, f" {c_eq_liq:.4f}", fontsize=8, verticalalignment='bottom')

            # Draw tangent on Free Energy plot
            f_eq_fcc = td.get_free_energy(np.array([c_eq_fcc]), T, 0)
            f_eq_liq = td.get_free_energy(np.array([c_eq_liq]), T, 2)
            #ax_row[0].plot([c_eq_fcc, c_eq_liq], [f_eq_fcc, f_eq_liq], 'k--', linewidth=1.5, alpha=0.7)
            #ax_row[0].plot([c_eq_fcc, c_eq_liq], [f_eq_fcc, f_eq_liq], 'ko', markersize=5)
    
    # Check intersection between BCC(1) and LIQ(2)
    mu_min = max(np.min(mus[1]), np.min(mus[2]))
    mu_max = min(np.max(mus[1]), np.max(mus[2]))
    
    if mu_max > mu_min:
        common_mu = np.linspace(mu_min, mu_max, 1000)
        
        sort_idx_1 = np.argsort(mus[1])
        sort_idx_2 = np.argsort(mus[2])

        psi_1_int = np.interp(common_mu, mus[1][sort_idx_1], psis[1][sort_idx_1])
        psi_2_int = np.interp(common_mu, mus[2][sort_idx_2], psis[2][sort_idx_2])
        
        idx = np.argwhere(np.diff(np.sign(psi_1_int - psi_2_int))).flatten()
        if len(idx) > 0:
            mu_eq = common_mu[idx[0]]
            psi_eq = psi_1_int[idx[0]]
            ax_row[2].plot(mu_eq, psi_eq, 'm*', markersize=12, label='BCC/LIQ Eq')

            # Optional: Project back to Chemical Potential Plot
            c_eq_bcc = np.interp(mu_eq, mus[1][sort_idx_1], c_range[sort_idx_1])
            c_eq_liq = np.interp(mu_eq, mus[2][sort_idx_2], c_range[sort_idx_2])
            
            # Draw tie-line on Chem Pot plot
            ax_row[1].axhline(mu_eq, color='m', linestyle=':', alpha=0.4)
            ax_row[1].plot([c_eq_bcc, c_eq_liq], [mu_eq, mu_eq], 'm-o', markersize=5, alpha=0.6)
            ax_row[1].text(c_eq_bcc, mu_eq, f" {c_eq_bcc:.4f}", color='m', fontsize=8, verticalalignment='top')
            ax_row[1].text(c_eq_liq, mu_eq, f" {c_eq_liq:.4f}", color='m', fontsize=8, verticalalignment='top')

            # Draw tangent on Free Energy plot
            f_eq_bcc = td.get_free_energy(np.array([c_eq_bcc]), T, 1)
            f_eq_liq = td.get_free_energy(np.array([c_eq_liq]), T, 2)
            ax_row[0].plot([c_eq_bcc, c_eq_liq], [f_eq_bcc, f_eq_liq], 'm--', linewidth=1.5, alpha=0.7)
            ax_row[0].plot([c_eq_bcc, c_eq_liq], [f_eq_bcc, f_eq_liq], 'mo', markersize=5)

    ax_row[2].legend()

    for ax in ax_row:
        ax.grid(True, alpha=0.3)
        
axes[0, 0].set_title("Free Energy ($f$) vs Comp ($c$)")
axes[0, 1].set_title("Chem Pot ($\mu$) vs Comp ($c$)")
axes[0, 2].set_title("Legendre Transform ($\Psi$) vs $\mu$")

plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.show()