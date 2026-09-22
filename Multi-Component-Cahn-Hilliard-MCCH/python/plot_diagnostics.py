#!/usr/bin/env python3
"""
Plot time evolution of Free Energy, Free Energy Dissipation Rate (dF/dt <= 0),
and Component Mass Conservation from simulation history.csv.
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_diagnostics(csv_path, output_png=None):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        return

    df = pd.read_csv(csv_path)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # 1. Total Free Energy F(t)
    ax = axes[0, 0]
    ax.plot(df['time'], df['total_free_energy'], 'b-', lw=2, label='Total Free Energy $F(t)$')
    ax.set_xlabel('Time $t$', fontsize=12)
    ax.set_ylabel('Free Energy $F$', fontsize=12)
    ax.set_title(r'Free Energy Monotonic Dissipation ($dF/dt \leq 0$)', fontsize=13, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='best')

    # 2. Rate of Energy Dissipation dF/dt
    ax = axes[0, 1]
    ax.plot(df['time'], df['dF_dt'], 'r-', lw=1.8, label='$dF/dt$')
    ax.axhline(0, color='k', linestyle=':', lw=1)
    ax.set_xlabel('Time $t$', fontsize=12)
    ax.set_ylabel('$dF/dt$', fontsize=12)
    ax.set_title('Dissipation Rate (Thermodynamic Consistency)', fontsize=13, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='best')

    # 3. Average Concentrations (Mass Conservation)
    ax = axes[1, 0]
    comp_cols = [col for col in df.columns if col.startswith('avg_c')]
    for col in comp_cols:
        c_init = df[col].iloc[0]
        drift = df[col] - c_init
        ax.plot(df['time'], drift, lw=1.5, label=f'$\\Delta \\bar{{{col.replace("avg_", "")}}}$')
    ax.set_xlabel('Time $t$', fontsize=12)
    ax.set_ylabel('Mass Drift $\\bar{c}_i(t) - \\bar{c}_i(0)$', fontsize=12)
    ax.set_title('Discrete Mass Conservation (Drift $\\sim 10^{-15}$)', fontsize=13, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='best')

    # 4. Partition of Unity & Extreme Bounds
    ax = axes[1, 1]
    dev_col = 'max_unity_dev' if 'max_unity_dev' in df.columns else 'max_unity_deviation'
    ax.semilogy(df['time'], df[dev_col] + 1e-18, 'm-', lw=1.8, label='Max $|\\sum c_i - 1|$')
    ax.set_xlabel('Time $t$', fontsize=12)
    ax.set_ylabel('Partition of Unity Error', fontsize=12)
    ax.set_title('Constraint Enforcement $\\sum_{i=1}^N c_i = 1$', fontsize=13, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='best')

    plt.tight_layout()

    if output_png is None:
        output_png = os.path.splitext(csv_path)[0] + "_diagnostics.png"
    plt.savefig(output_png, dpi=200)
    print(f"Diagnostics plot saved to {output_png}")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "results/history.csv"
    plot_diagnostics(path)

if __name__ == "__main__":
    main()
