#!/usr/bin/env python3
"""
Launcher script for the Cahn-Hilliard Multi-Component Phase-Field Studio Web App.
Usage:
    python3 run_app.py [--port 5000] [--host 0.0.0.0] [--open-browser]
"""

import sys
import os
import argparse
import webbrowser
import threading
import time

# Ensure project root is in sys.path
WORKSPACE_ROOT = os.path.dirname(os.path.abspath(__file__))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

from webapp.app import app

def print_banner(host, port):
    print("=" * 74)
    print("  ⚛️  Generalized Multi-Component Cahn-Hilliard Phase-Field Studio")
    print("  High-Performance MPI Slab (Pencil) Parallelized C++17 & Python Solvers")
    print("=" * 74)
    print(f"  • Local Web Server URL:    http://localhost:{port}")
    if host == "0.0.0.0":
        print(f"  • Network Accessible URL:  http://0.0.0.0:{port}")
    print("  • Solvers Supported:       C++17 MPI Slab HPC Solver & Python Reference")
    print("  • Features Available:      Interactive 2D/3D Slicing, RGB Phase Composite,")
    print("                             MPI Slab Boundaries Overlay, 1D Cross-Sections,")
    print("                             Thermodynamic Diagnostics (dF/dt <= 0),")
    print("                             Real-Time Simulation Monitoring & Live Console.")
    print("=" * 74)
    print("  Press Ctrl+C to terminate the application server.\n")


def open_browser_delayed(url, delay=1.0):
    time.sleep(delay)
    try:
        webbrowser.open(url)
    except Exception:
        pass


def main():
    parser = argparse.ArgumentParser(description="Launch Cahn-Hilliard Studio Web App")
    parser.add_argument("-p", "--port", type=int, default=5000, help="Port to bind (default: 5000)")
    parser.add_argument("-H", "--host", type=str, default="0.0.0.0", help="Host address (default: 0.0.0.0)")
    parser.add_argument("--open-browser", action="store_true", help="Automatically open browser")
    parser.add_argument("--debug", action="store_true", help="Enable Flask debug mode")
    args = parser.parse_args()

    print_banner(args.host, args.port)

    if args.open_browser:
        t = threading.Thread(target=open_browser_delayed, args=(f"http://localhost:{args.port}", 1.2), daemon=True)
        t.start()

    try:
        app.run(host=args.host, port=args.port, debug=args.debug, threaded=True)
    except KeyboardInterrupt:
        print("\n[Studio] Web server stopped.")


if __name__ == "__main__":
    main()
