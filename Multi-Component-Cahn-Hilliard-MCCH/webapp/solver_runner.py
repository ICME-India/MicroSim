"""
Asynchronous Process Manager for Cahn-Hilliard Simulations.
Spawns, monitors, and terminates MPI C++ and Python solver jobs while capturing real-time metrics.
"""

import os
import sys
import json
import time
import signal
import subprocess
import threading
import re
import shutil
from typing import Dict, Any, Optional

class SolverRunner:
    def __init__(self, workspace_root: str):
        self.workspace_root = os.path.abspath(workspace_root)
        self.process: Optional[subprocess.Popen] = None
        self.thread: Optional[threading.Thread] = None
        self.lock = threading.Lock()

        # State attributes
        self.status = "idle"  # idle, running, completed, failed, stopped
        self.backend = "cpp"
        self.mpi_ranks = 2
        self.run_id = ""
        self.output_dir = ""
        self.config_path = ""
        self.start_time = 0.0
        self.end_time = 0.0

        # Live metrics
        self.total_steps = 1000
        self.current_step = 0
        self.current_time = 0.0
        self.current_energy = 0.0
        self.current_df_dt = 0.0
        self.current_mlups = 0.0
        self.error_message = ""
        self.logs = []
        self.max_logs = 2000

    def get_status(self) -> Dict[str, Any]:
        with self.lock:
            elapsed = 0.0
            if self.status == "running":
                elapsed = time.time() - self.start_time
            elif self.end_time > self.start_time:
                elapsed = self.end_time - self.start_time

            pct = 0.0
            if self.total_steps > 0:
                pct = min(100.0, round((self.current_step / self.total_steps) * 100.0, 1))

            return {
                "status": self.status,
                "backend": self.backend,
                "mpi_ranks": self.mpi_ranks,
                "run_id": self.run_id,
                "output_dir": self.output_dir,
                "current_step": self.current_step,
                "total_steps": self.total_steps,
                "progress_percent": pct,
                "current_time": self.current_time,
                "current_energy": self.current_energy,
                "current_df_dt": self.current_df_dt,
                "current_mlups": self.current_mlups,
                "elapsed_seconds": round(elapsed, 2),
                "error_message": self.error_message,
                "log_count": len(self.logs)
            }

    def get_logs(self, since_index: int = 0):
        with self.lock:
            if since_index < len(self.logs):
                return self.logs[since_index:], len(self.logs)
            return [], len(self.logs)

    def _find_cpp_solver(self) -> Optional[str]:
        candidates = []
        if "MCCH_SOLVER_BIN" in os.environ:
            candidates.append(os.environ["MCCH_SOLVER_BIN"])
        candidates.append(os.path.join(self.workspace_root, "build", "mcch_solver"))
        candidates.append(os.path.join(self.workspace_root, "bin", "mcch_solver"))
        which_path = shutil.which("mcch_solver")
        if which_path:
            candidates.append(which_path)
        candidates.append(os.path.expanduser("~/.local/bin/mcch_solver"))
        candidates.append("/usr/local/bin/mcch_solver")

        for p in candidates:
            if p and os.path.isfile(p) and os.access(p, os.X_OK):
                return os.path.abspath(p)
        return None

    def _find_mpi_launcher(self) -> Optional[str]:
        if "MPIEXEC" in os.environ:
            return os.environ["MPIEXEC"]
        for launcher in ["mpirun", "mpiexec", "srun"]:
            path = shutil.which(launcher)
            if path:
                return path
        return None

    def start_simulation(self, config: Dict[str, Any], backend: str = "cpp", mpi_ranks: int = 2, device: str = "cpu") -> Dict[str, Any]:
        with self.lock:
            if self.status == "running":
                return {"success": False, "error": "A simulation is already in progress"}

            self.backend = backend
            self.device = device
            self.mpi_ranks = max(1, int(mpi_ranks))

            # Determine output directory
            out_dir = config.get("output_dir", "")
            if not out_dir:
                timestamp = time.strftime("%Y%m%d_%H%M%S")
                out_dir = f"results_run_{timestamp}"
                config["output_dir"] = out_dir

            config["device"] = device
            self.run_id = os.path.basename(out_dir)
            self.output_dir = os.path.join(self.workspace_root, out_dir)
            os.makedirs(self.output_dir, exist_ok=True)

            # Save config JSON inside output_dir and project root
            config_filename = f"config_{self.run_id}.json"
            self.config_path = os.path.join(self.workspace_root, config_filename)
            with open(self.config_path, 'w') as f:
                json.dump(config, f, indent=2)

            # Also save a copy inside the run folder
            with open(os.path.join(self.output_dir, "config.json"), 'w') as f:
                json.dump(config, f, indent=2)

            self.total_steps = int(config.get("total_steps", 1000))
            self.current_step = 0
            self.current_time = 0.0
            self.current_energy = 0.0
            self.current_df_dt = 0.0
            self.current_mlups = 0.0
            self.error_message = ""
            self.logs = []
            self.status = "running"
            self.start_time = time.time()
            self.end_time = 0.0

        # Determine command line
        mpi_launcher = self._find_mpi_launcher()

        if backend == "cpp":
            executable = self._find_cpp_solver()
            if not executable:
                with self.lock:
                    self.status = "failed"
                    self.error_message = (
                        "C++ binary 'mcch_solver' not found. Please compile or install the solver "
                        "using './install.sh' or 'cmake -B build && cmake --build build'."
                    )
                return {"success": False, "error": self.error_message}

            if self.mpi_ranks > 1 or mpi_launcher:
                launcher = mpi_launcher if mpi_launcher else "mpirun"
                cmd = [launcher, "-np", str(self.mpi_ranks), executable, "-c", self.config_path]
            else:
                cmd = [executable, "-c", self.config_path]

            if device in ("gpu", "cuda"):
                cmd.append("--gpu")
        else:
            py_solver = os.path.join(self.workspace_root, "python", "mcch_solver.py")
            if not os.path.exists(py_solver):
                which_py = shutil.which("mcch_solver.py")
                if which_py:
                    py_solver = which_py
            if self.mpi_ranks > 1 and mpi_launcher:
                cmd = [mpi_launcher, "-np", str(self.mpi_ranks), sys.executable, py_solver, "-c", self.config_path]
            else:
                cmd = [sys.executable, py_solver, "-c", self.config_path]
            if device in ("gpu", "cuda"):
                cmd.append("--gpu")

        self.thread = threading.Thread(target=self._run_process_thread, args=(cmd,), daemon=True)
        self.thread.start()

        return {"success": True, "run_id": self.run_id, "output_dir": self.output_dir}

    def _run_process_thread(self, cmd):
        log_msg = f"[RUNNER] Launching: {' '.join(cmd)}"
        with self.lock:
            self.logs.append(log_msg)

        # Regex to parse simulation step rows
        # Format: Step Time Total Energy dF/dt Max |sum(c)-1| c_min c_max MLUPS
        step_pattern = re.compile(r'^\s*(\d+)\s+([\d\.\+eE\-]+)\s+([\d\.\+eE\-]+)\s+([\d\.\+eE\-]+)\s+([\d\.\+eE\-]+)\s+([\d\.\+eE\-]+)\s+([\d\.\+eE\-]+)\s+([\d\.\+eE\-]+)')

        try:
            self.process = subprocess.Popen(
                cmd,
                cwd=self.workspace_root,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                bufsize=1,
                preexec_fn=os.setsid
            )

            for line in self.process.stdout:
                clean_line = line.rstrip()
                with self.lock:
                    if len(self.logs) < self.max_logs:
                        self.logs.append(clean_line)
                    else:
                        self.logs.pop(0)
                        self.logs.append(clean_line)

                    # Match step progress
                    m = step_pattern.match(clean_line)
                    if m:
                        try:
                            self.current_step = int(m.group(1))
                            self.current_time = float(m.group(2))
                            self.current_energy = float(m.group(3))
                            self.current_df_dt = float(m.group(4))
                            self.current_mlups = float(m.group(8))
                        except (ValueError, IndexError):
                            pass

            self.process.wait()
            ret_code = self.process.returncode

            with self.lock:
                self.end_time = time.time()
                if self.status == "stopped":
                    self.logs.append("[RUNNER] Simulation was stopped by user.")
                elif ret_code == 0:
                    self.status = "completed"
                    self.current_step = self.total_steps
                    self.logs.append(f"[RUNNER] Simulation completed successfully in {round(self.end_time - self.start_time, 2)}s.")
                else:
                    self.status = "failed"
                    self.error_message = f"Process exited with non-zero code {ret_code}"
                    self.logs.append(f"[RUNNER] Simulation failed with code {ret_code}")

        except Exception as e:
            with self.lock:
                self.status = "failed"
                self.error_message = str(e)
                self.logs.append(f"[RUNNER EXCEPTION] {str(e)}")
        finally:
            self.process = None

    def stop_simulation(self) -> Dict[str, Any]:
        with self.lock:
            if not self.process or self.status != "running":
                return {"success": False, "message": "No active simulation to stop"}

            self.status = "stopped"
            try:
                # Send SIGTERM to process group
                os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            except Exception as e:
                try:
                    self.process.terminate()
                except Exception:
                    pass

        # Brief wait then kill if needed
        time.sleep(0.5)
        if self.process and self.process.poll() is None:
            try:
                os.killpg(os.getpgid(self.process.pid), signal.SIGKILL)
            except Exception:
                pass

        return {"success": True, "message": "Simulation stopped"}
