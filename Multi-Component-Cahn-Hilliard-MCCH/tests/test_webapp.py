"""
Automated End-to-End Test Suite for Cahn-Hilliard Web Application.
Verifies all REST API endpoints, HTML rendering, VTK slicing, diagnostics parsing, and simulation execution.
"""

import sys
import os
import json
import time

# Ensure project root is in path
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from webapp.app import app

def run_tests():
    print("==================================================================")
    print("  Running Cahn-Hilliard WebApp Test Suite")
    print("==================================================================")

    client = app.test_client()

    # 1. Test Home Page
    res = client.get('/')
    assert res.status_code == 200, f"Expected 200 for /, got {res.status_code}"
    assert b"Cahn-Hilliard Phase-Field Studio" in res.data
    print(" [PASS] GET / -> HTTP 200 (HTML rendered with brand title)")

    # 2. Test Static Assets
    for asset in ['/static/vendor/bootstrap.min.css', '/static/vendor/chart.umd.min.js', '/static/css/style.css', '/static/js/app.js', '/static/js/viewer.js', '/static/js/charts.js']:
        res = client.get(asset)
        assert res.status_code == 200, f"Expected 200 for {asset}, got {res.status_code}"
        assert len(res.data) > 0
    print(" [PASS] Static assets verified (Bootstrap, Chart.js, styles, JS modules)")

    # 3. Test /api/runs
    res = client.get('/api/runs')
    assert res.status_code == 200
    runs = res.json.get('runs', [])
    assert len(runs) >= 5, f"Expected at least 5 preset runs, got {len(runs)}"
    run_ids = [r['id'] for r in runs]
    print(f" [PASS] GET /api/runs -> Found runs: {run_ids}")

    # 4. Test /api/presets
    res = client.get('/api/presets')
    assert res.status_code == 200
    presets = res.json.get('presets', [])
    assert len(presets) >= 5, f"Expected at least 5 presets, got {len(presets)}"
    print(f" [PASS] GET /api/presets -> Found {len(presets)} simulation presets")

    # 5. Test /api/run/<run_id>/metadata for 2D run
    res = client.get('/api/run/results_ternary_spinodal/metadata')
    assert res.status_code == 200
    meta = res.json
    assert meta['dim'] == 2
    assert meta['nx'] == 128
    assert 'c0' in meta['components']
    assert len(meta['timesteps']) > 0
    print(f" [PASS] GET /api/run/results_ternary_spinodal/metadata -> 2D, nx={meta['nx']}, frames={meta['total_frames']}")

    # 6. Test /api/run/<run_id>/metadata for 3D run
    res = client.get('/api/run/results_3d_spinodal/metadata')
    assert res.status_code == 200
    meta_3d = res.json
    assert meta_3d['dim'] == 3
    assert meta_3d['nz'] > 1
    print(f" [PASS] GET /api/run/results_3d_spinodal/metadata -> 3D, nx={meta_3d['nx']}, nz={meta_3d['nz']}")

    # 7. Test /api/run/<run_id>/slice_img
    res = client.get('/api/run/results_ternary_spinodal/slice_img?step=20000&field=c0&colormap=viridis&overlay_slabs=4')
    assert res.status_code == 200
    assert res.content_type == 'image/png'
    assert len(res.data) > 1000
    assert 'X-Field-Min' in res.headers
    print(f" [PASS] GET /api/run/.../slice_img -> Generated 2D PNG ({len(res.data)} bytes, min={res.headers['X-Field-Min']}, max={res.headers['X-Field-Max']})")

    # 8. Test RGB composite
    res = client.get('/api/run/results_ternary_spinodal/slice_img?step=20000&field=rgb&overlay_slabs=4')
    assert res.status_code == 200
    assert res.content_type == 'image/png'
    assert res.headers.get('X-Is-RGB') == '1'
    print(f" [PASS] GET /api/run/.../slice_img?field=rgb -> Generated RGB Phase Composite ({len(res.data)} bytes)")

    # 9. Test Hover Coordinates Info
    res = client.get('/api/run/results_ternary_spinodal/hover_info?step=20000&x=50&y=60')
    assert res.status_code == 200
    hover = res.json
    assert 'values' in hover
    assert 'c0' in hover['values']
    assert abs(hover['sum_c'] - 1.0) < 1e-4
    print(f" [PASS] GET /api/run/.../hover_info -> x=50, y=60: {hover['values']}, sum_c={hover['sum_c']}")

    # 10. Test Line Cut Profile
    res = client.get('/api/run/results_ternary_spinodal/linecut?step=20000&axis=x&coord=64')
    assert res.status_code == 200
    cut = res.json
    assert cut['axis'] == 'x'
    assert 'c0' in cut['profiles']
    assert len(cut['positions']) == 128
    print(f" [PASS] GET /api/run/.../linecut -> 128 data points along X at fixed Y=64")

    # 11. Test Diagnostics CSV Parser
    res = client.get('/api/run/results_ternary_spinodal/diagnostics')
    assert res.status_code == 200
    diag = res.json
    assert diag['is_monotone_dissipative'] == True
    assert diag['max_mass_drift'] < 1e-10
    print(f" [PASS] GET /api/run/.../diagnostics -> Initial F={diag['initial_free_energy']:.2f}, Final F={diag['final_free_energy']:.2f}, Monotonic: {diag['is_monotone_dissipative']}")

    # 12. Test Live Simulation Execution via API (Fast test run with 50 steps using C++ MPI)
    test_sim_config = {
        "dim": 2,
        "nx": 64,
        "ny": 64,
        "nz": 1,
        "dx": 1.0,
        "dy": 1.0,
        "bc_x": "periodic",
        "bc_y": "periodic",
        "num_components": 2,
        "free_energy_type": "binary_double_well",
        "mobility_type": "constant",
        "mobility_val": 1.0,
        "kappa": 1.0,
        "integrator": "rk2",
        "dt": 0.01,
        "total_steps": 50,
        "output_interval": 25,
        "diag_interval": 25,
        "initial_condition": "random",
        "c_mean": [0.5, 0.5],
        "noise_amp": 0.05,
        "seed": 42,
        "output_dir": "results_webapp_test"
    }

    print("\n [TEST] Launching background simulation (C++ MPI 2 ranks)...")
    res = client.post('/api/simulation/start', json={
        "config": test_sim_config,
        "backend": "cpp",
        "mpi_ranks": 2
    })
    assert res.status_code == 200
    assert res.json['success'] == True

    # Poll status until completion (up to 15 seconds)
    completed = False
    for _ in range(30):
        time.sleep(0.5)
        status_res = client.get('/api/simulation/status')
        st = status_res.json
        print(f"   Status: {st['status']}, Step: {st['current_step']}/{st['total_steps']}, MLUPS: {st['current_mlups']}")
        if st['status'] == 'completed':
            completed = True
            break
        elif st['status'] in ('failed', 'stopped'):
            break

    assert completed, f"Simulation did not complete successfully. Status: {st}"
    print(" [PASS] Simulation completed successfully via webapp process manager!")

    # Verify generated output files
    assert os.path.exists('results_webapp_test/history.csv')
    assert os.path.exists('results_webapp_test/solution_00050.vtk')

    # Clean up test output folder
    os.system('rm -rf results_webapp_test config_results_webapp_test.json')
    print(" [PASS] Cleanup of test results completed.")

    print("\n==================================================================")
    print("  ALL 12 TESTS PASSED! Cahn-Hilliard Web App is 100% Functional.")
    print("==================================================================")

if __name__ == "__main__":
    run_tests()
