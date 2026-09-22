"""
Flask Web Application for the Multi-Component Cahn-Hilliard Phase-Field Solver Studio.
Provides RESTful APIs for simulation control, real-time logging, VTK slicing, diagnostics, and visual analytics.
"""

import os
import sys
import glob
import json
import re
from typing import Dict, Any, List
from flask import Flask, render_template, request, jsonify, send_file, Response, send_from_directory
import numpy as np
import io

# Add project root to path
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from webapp.solver_runner import SolverRunner
from webapp.vtk_utils import (
    parse_vtk_file,
    parse_vtk_header,
    extract_slice_2d,
    render_slice_png,
    render_rgb_composite_png,
    extract_line_cut
)
from webapp.diagnostics_utils import parse_history_csv

app = Flask(__name__,
            template_folder=os.path.join(os.path.dirname(__file__), "templates"),
            static_folder=os.path.join(os.path.dirname(__file__), "static"))

runner = SolverRunner(WORKSPACE_ROOT)

# Friendly preset names
PRESET_NAMES = {
    "binary_spinodal.json": {
        "title": "Classic Binary Spinodal (N=2)",
        "description": "Two-phase separation with double-well free energy functional, periodic boundaries.",
        "icon": "bi-circle-half"
    },
    "ternary_spinodal.json": {
        "title": "Ternary Spinodal (N=3)",
        "description": "Three-phase separation forming diffuse interface networks and tri-junction points.",
        "icon": "bi-triangle"
    },
    "ternary_droplet.json": {
        "title": "Ternary Droplet / Ripening (N=3)",
        "description": "Isolated circular droplet with Neumann (no-flux) boundary conditions.",
        "icon": "bi-droplet"
    },
    "regular_solution_3comp.json": {
        "title": "Regular Solution / Degenerate Mobility (N=3)",
        "description": "Logarithmic Flory-Huggins entropy with degenerate Onsager mobility tensor.",
        "icon": "bi-diagram-3"
    },
    "3d_spinodal.json": {
        "title": "Full 3D Spinodal Decomposition (48x48x48)",
        "description": "Three-dimensional bicontinuous interconnected phase structure.",
        "icon": "bi-box"
    }
}


def find_available_runs() -> List[Dict[str, Any]]:
    """Scan workspace for simulation result directories."""
    runs = []
    # Check top-level directories starting with results_ or containing history.csv
    candidates = glob.glob(os.path.join(WORKSPACE_ROOT, "results*"))
    for d in sorted(candidates):
        if os.path.isdir(d):
            base_name = os.path.basename(d)
            history_path = os.path.join(d, "history.csv")
            vtk_files = glob.glob(os.path.join(d, "solution_*.vtk"))

            if os.path.exists(history_path) or len(vtk_files) > 0:
                # Pretty title
                title = base_name.replace("results_", "").replace("_", " ").title()
                if "3D" in title or "3d" in title:
                    title = title.replace("3D", "3D").replace("3d", "3D")

                runs.append({
                    "id": base_name,
                    "title": title,
                    "path": d,
                    "has_history": os.path.exists(history_path),
                    "vtk_count": len(vtk_files),
                    "is_active": (runner.status == "running" and runner.run_id == base_name)
                })
    return runs


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/runs', methods=['GET'])
def get_runs():
    runs = find_available_runs()
    return jsonify({"runs": runs, "active_run": runner.run_id if runner.status == "running" else None})


@app.route('/api/presets', methods=['GET'])
def get_presets():
    examples_dir = os.path.join(WORKSPACE_ROOT, "examples")
    presets = []
    if os.path.exists(examples_dir):
        for f in sorted(os.listdir(examples_dir)):
            if f.endswith('.json'):
                fpath = os.path.join(examples_dir, f)
                try:
                    with open(fpath, 'r') as jf:
                        cfg = json.load(jf)
                    info = PRESET_NAMES.get(f, {
                        "title": f.replace('.json', '').replace('_', ' ').title(),
                        "description": "Simulation preset configuration.",
                        "icon": "bi-gear"
                    })
                    presets.append({
                        "file": f,
                        "title": info["title"],
                        "description": info["description"],
                        "icon": info["icon"],
                        "config": cfg
                    })
                except Exception as e:
                    pass
    return jsonify({"presets": presets})


@app.route('/api/preset/<filename>', methods=['GET'])
def get_preset(filename):
    fpath = os.path.join(WORKSPACE_ROOT, "examples", filename)
    if not os.path.exists(fpath):
        return jsonify({"error": f"Preset {filename} not found"}), 404
    with open(fpath, 'r') as f:
        cfg = json.load(f)
    return jsonify(cfg)


@app.route('/api/run/<run_id>/metadata', methods=['GET'])
def get_run_metadata(run_id):
    run_dir = os.path.join(WORKSPACE_ROOT, run_id)
    if not os.path.exists(run_dir):
        return jsonify({"error": f"Run {run_id} does not exist"}), 404

    # Look for config
    config = {}
    config_paths = [
        os.path.join(run_dir, "config.json"),
        os.path.join(WORKSPACE_ROOT, f"config_{run_id}.json"),
        os.path.join(WORKSPACE_ROOT, "examples", f"{run_id.replace('results_', '')}.json")
    ]
    for cp in config_paths:
        if os.path.exists(cp):
            try:
                with open(cp, 'r') as f:
                    config = json.load(f)
                break
            except Exception:
                pass

    # Find all VTK files and extract step numbers
    vtk_files = sorted(glob.glob(os.path.join(run_dir, "solution_*.vtk")))
    timesteps = []
    step_pattern = re.compile(r'solution_(\d+)\.vtk')

    first_vtk_info = None
    for vf in vtk_files:
        fname = os.path.basename(vf)
        m = step_pattern.search(fname)
        if m:
            step_num = int(m.group(1))
            timesteps.append({
                "step": step_num,
                "file": fname
            })

    # Sort timesteps by step
    timesteps.sort(key=lambda x: x["step"])

    # Parse header of first VTK file to get grid dimensions & component fields
    nx, ny, nz, dim = 128, 128, 1, 2
    fields = ['c0', 'c1']
    if vtk_files:
        hdr = parse_vtk_header(vtk_files[0])
        if hdr:
            nx, ny, nz = hdr["nx"], hdr["ny"], hdr["nz"]
            dim = 3 if nz > 1 else 2
            fields = hdr["fields"]

    comp_fields = [f for f in fields if f.startswith('c')]
    mu_fields = [f for f in fields if f.startswith('mu')]

    has_history = os.path.exists(os.path.join(run_dir, "history.csv"))

    return jsonify({
        "run_id": run_id,
        "nx": nx,
        "ny": ny,
        "nz": nz,
        "dim": dim,
        "components": comp_fields,
        "chemical_potentials": mu_fields,
        "num_components": len(comp_fields),
        "total_frames": len(timesteps),
        "timesteps": timesteps,
        "has_history": has_history,
        "config": config
    })


@app.route('/api/run/<run_id>/slice_img', methods=['GET'])
def get_slice_image(run_id):
    run_dir = os.path.join(WORKSPACE_ROOT, run_id)
    step = request.args.get('step', type=int, default=0)
    field = request.args.get('field', type=str, default='c0')
    plane = request.args.get('plane', type=str, default='xy')
    slice_idx = request.args.get('slice_idx', type=int, default=None)
    colormap = request.args.get('colormap', type=str, default='viridis')
    slabs = request.args.get('overlay_slabs', type=int, default=0)
    target_size = request.args.get('target_size', type=int, default=512)

    # Locate VTK file
    vtk_path = os.path.join(run_dir, f"solution_{step:05d}.vtk")
    if not os.path.exists(vtk_path):
        # Fallback to closest available file
        files = sorted(glob.glob(os.path.join(run_dir, "solution_*.vtk")))
        if not files:
            return "No VTK files in run", 404
        vtk_path = files[0]

    parsed = parse_vtk_file(vtk_path)
    if not parsed:
        return "Failed to parse VTK file", 500

    headers = {
        "X-Step": str(parsed["step"]),
        "X-Time": str(parsed["time"]),
        "X-Nx": str(parsed["nx"]),
        "X-Ny": str(parsed["ny"]),
        "X-Nz": str(parsed["nz"])
    }

    # Handle RGB composite request
    if field == 'rgb':
        comp_keys = [k for k in sorted(parsed["fields"].keys()) if k.startswith('c')]
        if len(comp_keys) >= 3:
            s0, _ = extract_slice_2d(parsed, 'c0', plane, slice_idx)
            s1, _ = extract_slice_2d(parsed, 'c1', plane, slice_idx)
            s2, _ = extract_slice_2d(parsed, 'c2', plane, slice_idx)
            png_bytes = render_rgb_composite_png(s0, s1, s2, overlay_slabs=slabs, orig_nx=parsed["nx"], target_size=target_size)
            resp = Response(png_bytes, mimetype='image/png')
            for k, v in headers.items():
                resp.headers[k] = v
            resp.headers["X-Is-RGB"] = "1"
            return resp
        else:
            # Fallback to c0 with coolwarm
            field = 'c0'
            colormap = 'coolwarm'

    # Single field extraction
    if field not in parsed["fields"]:
        field = 'c0' if 'c0' in parsed["fields"] else list(parsed["fields"].keys())[0]

    slice_2d, extent = extract_slice_2d(parsed, field, plane, slice_idx)
    if slice_2d is None:
        return "Slice extraction error", 500

    png_bytes, vmin, vmax = render_slice_png(
        slice_2d,
        colormap_name=colormap,
        overlay_slabs=slabs,
        orig_nx=parsed["nx"],
        target_size=target_size
    )

    resp = Response(png_bytes, mimetype='image/png')
    for k, v in headers.items():
        resp.headers[k] = v
    resp.headers["X-Field-Min"] = f"{vmin:.4f}"
    resp.headers["X-Field-Max"] = f"{vmax:.4f}"
    resp.headers["X-Slice-Plane"] = plane
    resp.headers["X-Slice-Pos"] = str(extent.get("slice_pos", 0))

    return resp


@app.route('/api/run/<run_id>/hover_info', methods=['GET'])
def get_hover_info(run_id):
    """Returns exact local composition and chemical potential values at (x, y, z)."""
    run_dir = os.path.join(WORKSPACE_ROOT, run_id)
    step = request.args.get('step', type=int, default=0)
    x = request.args.get('x', type=int, default=0)
    y = request.args.get('y', type=int, default=0)
    z = request.args.get('z', type=int, default=0)

    vtk_path = os.path.join(run_dir, f"solution_{step:05d}.vtk")
    if not os.path.exists(vtk_path):
        return jsonify({"error": "File not found"}), 404

    parsed = parse_vtk_file(vtk_path)
    if not parsed:
        return jsonify({"error": "Parse error"}), 500

    nx, ny, nz = parsed["nx"], parsed["ny"], parsed["nz"]
    x = max(0, min(nx - 1, x))
    y = max(0, min(ny - 1, y))
    z = max(0, min(nz - 1, z))

    values = {}
    sum_c = 0.0
    for k, arr in parsed["fields"].items():
        val = float(arr[z, y, x])
        values[k] = round(val, 5)
        if k.startswith('c'):
            sum_c += val

    return jsonify({
        "step": step,
        "time": parsed["time"],
        "x": x,
        "y": y,
        "z": z,
        "values": values,
        "sum_c": round(sum_c, 6)
    })


@app.route('/api/run/<run_id>/linecut', methods=['GET'])
def get_linecut(run_id):
    run_dir = os.path.join(WORKSPACE_ROOT, run_id)
    step = request.args.get('step', type=int, default=0)
    axis = request.args.get('axis', type=str, default='x')
    coord = request.args.get('coord', type=int, default=None)
    z_slice = request.args.get('z_slice', type=int, default=0)

    vtk_path = os.path.join(run_dir, f"solution_{step:05d}.vtk")
    if not os.path.exists(vtk_path):
        return jsonify({"error": "File not found"}), 404

    parsed = parse_vtk_file(vtk_path)
    if not parsed:
        return jsonify({"error": "Parse error"}), 500

    cut_data = extract_line_cut(parsed, axis=axis, coord=coord, z_slice=z_slice)
    return jsonify(cut_data)


@app.route('/api/run/<run_id>/diagnostics', methods=['GET'])
def get_diagnostics(run_id):
    csv_path = os.path.join(WORKSPACE_ROOT, run_id, "history.csv")
    diag = parse_history_csv(csv_path)
    if not diag:
        return jsonify({"error": "No history.csv available for this run"}), 404
    return jsonify(diag)


@app.route('/api/simulation/start', methods=['POST'])
def start_simulation():
    data = request.get_json() or {}
    config = data.get('config', {})
    backend = data.get('backend', 'cpp')
    device = data.get('device', 'cpu')
    mpi_ranks = data.get('mpi_ranks', 2)

    result = runner.start_simulation(config, backend=backend, mpi_ranks=mpi_ranks, device=device)
    return jsonify(result)


@app.route('/api/simulation/status', methods=['GET'])
def simulation_status():
    return jsonify(runner.get_status())


@app.route('/api/simulation/logs', methods=['GET'])
def simulation_logs():
    since = request.args.get('since', type=int, default=0)
    lines, next_idx = runner.get_logs(since)
    return jsonify({"lines": lines, "next_index": next_idx})


@app.route('/api/simulation/stop', methods=['POST'])
def stop_simulation():
    return jsonify(runner.stop_simulation())


@app.route('/api/download/vtk/<run_id>/<step>', methods=['GET'])
def download_vtk(run_id, step):
    step_num = int(step)
    fpath = os.path.join(WORKSPACE_ROOT, run_id, f"solution_{step_num:05d}.vtk")
    if not os.path.exists(fpath):
        return "File not found", 404
    return send_file(fpath, as_attachment=True)


@app.route('/api/download/csv/<run_id>', methods=['GET'])
def download_csv(run_id):
    fpath = os.path.join(WORKSPACE_ROOT, run_id, "history.csv")
    if not os.path.exists(fpath):
        return "File not found", 404
    return send_file(fpath, as_attachment=True, download_name=f"{run_id}_history.csv")


@app.route('/api/download/config/<run_id>', methods=['GET'])
def download_config(run_id):
    fpath = os.path.join(WORKSPACE_ROOT, run_id, "config.json")
    if not os.path.exists(fpath):
        return "File not found", 404
    return send_file(fpath, as_attachment=True, download_name=f"{run_id}_config.json")


if __name__ == '__main__':
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    print(f"Starting Cahn-Hilliard Phase-Field Studio on http://localhost:{port}")
    app.run(host='0.0.0.0', port=port, debug=False, threaded=True)
