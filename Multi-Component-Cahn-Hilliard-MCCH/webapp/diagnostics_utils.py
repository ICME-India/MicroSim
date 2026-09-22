"""
Parser and statistical processor for Cahn-Hilliard simulation history.csv files.
Extracts energy dissipation, mass conservation, constraint enforcement, and performance metrics.
"""

import os
import pandas as pd
import numpy as np

def parse_history_csv(csv_path, max_chart_points=1000):
    """
    Parses history.csv and returns structured data for plotting and thermodynamic analysis.
    """
    if not os.path.exists(csv_path):
        return None

    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        return {"error": f"Failed to read CSV: {str(e)}"}

    if df.empty or 'time' not in df.columns:
        return {"error": "history.csv is empty or invalid format"}

    # Standardize column names
    unity_col = 'max_unity_dev' if 'max_unity_dev' in df.columns else ('max_unity_deviation' if 'max_unity_deviation' in df.columns else None)
    comp_cols = [c for c in df.columns if c.startswith('avg_c')]

    # Compute mass drifts
    mass_drifts = {}
    for c in comp_cols:
        c_init = float(df[c].iloc[0])
        mass_drifts[c] = (df[c] - c_init).tolist()

    # Downsampling for responsive chart rendering if necessary
    total_rows = len(df)
    stride = max(1, total_rows // max_chart_points)
    df_plot = df.iloc[::stride].copy()

    # Ensure last row is included
    if df_plot.index[-1] != df.index[-1]:
        df_plot = pd.concat([df_plot, df.iloc[[-1]]])

    # Summary thermodynamics
    f_init = float(df['total_free_energy'].iloc[0])
    f_final = float(df['total_free_energy'].iloc[-1])
    delta_f = f_final - f_init
    pct_f_drop = (delta_f / abs(f_init) * 100.0) if abs(f_init) > 1e-12 else 0.0

    max_drift_val = 0.0
    for c in comp_cols:
        c_init = float(df[c].iloc[0])
        c_max_d = float(np.max(np.abs(df[c] - c_init)))
        if c_max_d > max_drift_val:
            max_drift_val = c_max_d

    max_unity_val = float(df[unity_col].max()) if unity_col else 0.0
    max_df_dt = float(df['dF_dt'].max()) if 'dF_dt' in df.columns else 0.0
    min_df_dt = float(df['dF_dt'].min()) if 'dF_dt' in df.columns else 0.0

    result = {
        "num_records": total_rows,
        "initial_free_energy": f_init,
        "final_free_energy": f_final,
        "delta_free_energy": delta_f,
        "pct_free_energy_decrease": -pct_f_drop,
        "max_mass_drift": max_drift_val,
        "max_unity_deviation": max_unity_val,
        "max_df_dt": max_df_dt,
        "min_df_dt": min_df_dt,
        "is_monotone_dissipative": max_df_dt <= 1e-8,
        "components": [c.replace('avg_', '') for c in comp_cols],
        "series": {
            "steps": df_plot['step'].tolist() if 'step' in df_plot.columns else list(range(len(df_plot))),
            "times": [round(t, 6) for t in df_plot['time'].tolist()],
            "total_free_energy": df_plot['total_free_energy'].tolist(),
            "dF_dt": df_plot['dF_dt'].tolist() if 'dF_dt' in df_plot.columns else [],
            "unity_deviation": df_plot[unity_col].tolist() if unity_col else [],
            "min_c": df_plot['min_c'].tolist() if 'min_c' in df_plot.columns else [],
            "max_c": df_plot['max_c'].tolist() if 'max_c' in df_plot.columns else [],
            "avg_concentrations": {c: df_plot[c].tolist() for c in comp_cols},
            "mass_drifts": {c: [drift[i] for i in df_plot.index] for c, drift in mass_drifts.items()}
        }
    }

    return result
