#!/usr/bin/env python3
"""
Parallel Multi-Component Cahn-Hilliard Solver in Python using mpi4py & NumPy.
Implements 1D Slab Decomposition along the x-axis, matching the C++ HPC solver.
Supports JSON simulation configuration, VTK Structured Points output, and CSV diagnostics.
"""

import sys
import os
import json
import time
import argparse
import numpy as np
from mpi4py import MPI

class SlabMPIContext:
    def __init__(self, comm, nx, ny, nz=1, periodic=True):
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.periodic = periodic

        base = nx // self.size
        rem = nx % self.size
        self.local_nx = base + (1 if self.rank < rem else 0)
        self.local_x_start = self.rank * base + min(self.rank, rem)
        self.local_x_end = self.local_x_start + self.local_nx

        if periodic:
            self.left_rank = (self.rank - 1 + self.size) % self.size
            self.right_rank = (self.rank + 1) % self.size
        else:
            self.left_rank = MPI.PROC_NULL if self.rank == 0 else self.rank - 1
            self.right_rank = MPI.PROC_NULL if self.rank == self.size - 1 else self.rank + 1

    def exchange_ghosts(self, arr):
        """arr has shape (local_nx + 2, ny) or (local_nx + 2, ny, nz)"""
        if self.size == 1:
            if self.periodic:
                arr[0] = arr[-2]
                arr[-1] = arr[1]
            else:
                arr[0] = arr[1]
                arr[-1] = arr[-2]
            return

        # Multi-rank Sendrecv
        send_left = np.ascontiguousarray(arr[1])
        recv_left = np.empty_like(send_left)
        send_right = np.ascontiguousarray(arr[-2])
        recv_right = np.empty_like(send_right)

        # Exchange left: send to left, recv from right
        self.comm.Sendrecv(send_left, dest=self.left_rank, sendtag=11,
                          recvbuf=recv_right, source=self.right_rank, recvtag=11)
        # Exchange right: send to right, recv from left
        self.comm.Sendrecv(send_right, dest=self.right_rank, sendtag=22,
                          recvbuf=recv_left, source=self.left_rank, recvtag=22)

        if self.left_rank != MPI.PROC_NULL:
            arr[0] = recv_left
        elif not self.periodic and self.rank == 0:
            arr[0] = arr[1]

        if self.right_rank != MPI.PROC_NULL:
            arr[-1] = recv_right
        elif not self.periodic and self.rank == self.size - 1:
            arr[-1] = arr[-2]

    def gather(self, local_arr):
        """Gathers local interior arr[1:-1] to rank 0"""
        interior = np.ascontiguousarray(local_arr[1:-1])
        gathered = self.comm.gather(interior, root=0)
        if self.rank == 0:
            return np.concatenate(gathered, axis=0)
        return None


class PythonMCCHSolver:
    def __init__(self, mpi_ctx, dx=1.0, dy=1.0, dz=1.0, n_comp=3, W=None, A=None, kappa=1.0, mobility=1.0,
                 free_energy_type="polynomial_multiwell", RT=1.0, omega=None):
        self.ctx = mpi_ctx
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.n_comp = n_comp
        self.free_energy_type = free_energy_type
        self.RT = RT
        self.current_step = 0
        self.current_time = 0.0
        self.prev_energy = None

        # Mobility matrix
        if np.isscalar(mobility):
            self.mobility = np.eye(n_comp) * float(mobility)
        else:
            self.mobility = np.array(mobility, dtype=np.float64)

        # Gradient coefficients
        if np.isscalar(kappa):
            self.kappa = np.eye(n_comp) * float(kappa)
        else:
            self.kappa = np.array(kappa, dtype=np.float64)

        # Multi-well barrier matrix W
        if W is None:
            self.W = np.ones((n_comp, n_comp)) * 1.5
            np.fill_diagonal(self.W, 0.0)
        else:
            self.W = np.array(W, dtype=np.float64)

        # Self-barrier A
        if A is None:
            self.A = np.ones(n_comp) * 0.2
        else:
            self.A = np.array(A, dtype=np.float64)

        # Regular solution omega
        if omega is None:
            self.omega = np.ones((n_comp, n_comp)) * 3.6
            np.fill_diagonal(self.omega, 0.0)
        else:
            self.omega = np.array(omega, dtype=np.float64)

        # Storage: (n_comp, local_nx + 2, ny)
        self.c = np.zeros((n_comp, mpi_ctx.local_nx + 2, mpi_ctx.ny))
        self.mu = np.zeros_like(self.c)
        self.rhs = np.zeros_like(self.c)

    def initialize_random(self, c_mean, noise_amp=0.03, seed=42):
        np.random.seed(seed + self.ctx.rank * 999)
        c_mean = np.array(c_mean)
        for m in range(self.n_comp - 1):
            noise = np.random.uniform(-noise_amp, noise_amp, size=(self.ctx.local_nx, self.ctx.ny))
            self.c[m, 1:-1] = np.clip(c_mean[m] + noise, 1e-4, 0.999)
        self.c[self.n_comp - 1, 1:-1] = 1.0 - np.sum(self.c[:self.n_comp - 1, 1:-1], axis=0)
        self.exchange_all_ghosts()

    def initialize_droplet(self, droplet_comp=0, cx=64.0, cy=64.0, radius=28.0, width=2.5):
        # Global coordinate mapping for local slab
        for i_loc in range(self.ctx.local_nx):
            i_glob = self.ctx.local_x_start + i_loc
            x = (i_glob + 0.5) * self.dx
            for j in range(self.ctx.ny):
                y = (j + 0.5) * self.dy
                dist = np.sqrt((x - cx)**2 + (y - cy)**2)
                val = 0.5 * (1.0 - np.tanh((dist - radius) / (width * np.sqrt(2.0))))
                val = np.clip(val, 0.0, 1.0)
                self.c[droplet_comp, i_loc + 1, j] = val

                # Distribute remainder equally among other components
                rem = (1.0 - val) / float(self.n_comp - 1)
                for m in range(self.n_comp):
                    if m != droplet_comp:
                        self.c[m, i_loc + 1, j] = rem
        self.exchange_all_ghosts()

    def exchange_all_ghosts(self):
        for m in range(self.n_comp):
            self.ctx.exchange_ghosts(self.c[m])

    def compute_chemical_potentials(self):
        self.exchange_all_ghosts()
        idx2 = 1.0 / (self.dx * self.dx)
        idy2 = 1.0 / (self.dy * self.dy)

        # Compute Laplacians of c
        lap_c = np.zeros_like(self.c)
        for m in range(self.n_comp):
            c_m = self.c[m]
            lap_x = (c_m[2:] - 2.0 * c_m[1:-1] + c_m[:-2]) * idx2
            lap_y = (np.roll(c_m[1:-1], -1, axis=1) - 2.0 * c_m[1:-1] + np.roll(c_m[1:-1], 1, axis=1)) * idy2
            lap_c[m, 1:-1] = lap_x + lap_y

        c_int = self.c[:, 1:-1]
        df_dc = np.zeros_like(c_int)

        if self.free_energy_type == "regular_solution":
            eps = 1e-8
            for m in range(self.n_comp):
                entropy_term = self.RT * (np.log(np.maximum(c_int[m], eps)) + 1.0)
                enthalpy_term = np.zeros_like(c_int[0])
                for j in range(self.n_comp):
                    if j != m:
                        enthalpy_term += self.omega[m, j] * c_int[j]
                df_dc[m] = entropy_term + enthalpy_term
        elif self.free_energy_type == "binary_double_well":
            # 2 * c0 * (1 - c0) * (1 - 2*c0)
            c0 = c_int[0]
            df_dc[0] = 2.0 * c0 * (1.0 - c0) * (1.0 - 2.0 * c0)
            if self.n_comp > 1:
                df_dc[1] = -df_dc[0]
        else: # polynomial_multiwell
            for m in range(self.n_comp):
                d_term = np.zeros_like(c_int[0])
                for j in range(self.n_comp):
                    if j != m:
                        d_term += 2.0 * self.W[m, j] * c_int[m] * (c_int[j] ** 2)
                d_term += 2.0 * self.A[m] * c_int[m] * (1.0 - c_int[m]) * (1.0 - 2.0 * c_int[m])
                df_dc[m] = d_term

        # Chemical potentials: mu_m = df/dc_m - sum_n kappa_mn * lap_c_n
        raw_mu = np.zeros_like(c_int)
        for m in range(self.n_comp):
            grad_term = np.zeros_like(c_int[0])
            for n in range(self.n_comp):
                grad_term += self.kappa[m, n] * lap_c[n, 1:-1]
            raw_mu[m] = df_dc[m] - grad_term

        # Projection (subtract mean potential)
        mu_mean = np.mean(raw_mu, axis=0)
        self.mu[:, 1:-1] = raw_mu - mu_mean

        for m in range(self.n_comp):
            self.ctx.exchange_ghosts(self.mu[m])

    def compute_rhs(self):
        self.compute_chemical_potentials()
        idx2 = 1.0 / (self.dx * self.dx)
        idy2 = 1.0 / (self.dy * self.dy)

        laps = np.zeros_like(self.rhs[:, 1:-1])
        for n in range(self.n_comp):
            mu_n = self.mu[n]
            lap_x = (mu_n[2:] - 2.0 * mu_n[1:-1] + mu_n[:-2]) * idx2
            lap_y = (np.roll(mu_n[1:-1], -1, axis=1) - 2.0 * mu_n[1:-1] + np.roll(mu_n[1:-1], 1, axis=1)) * idy2
            laps[n] = lap_x + lap_y

        for m in range(self.n_comp):
            self.rhs[m, 1:-1] = np.sum([self.mobility[m, n] * laps[n] for n in range(self.n_comp)], axis=0)

    def step_euler(self, dt):
        self.compute_rhs()
        self.c[:, 1:-1] += dt * self.rhs[:, 1:-1]
        sum_c = np.sum(self.c[:, 1:-1], axis=0)
        self.c[:, 1:-1] /= sum_c[np.newaxis, ...]
        self.exchange_all_ghosts()
        self.current_step += 1
        self.current_time += dt

    def step_rk2(self, dt):
        self.compute_rhs()
        k1 = np.copy(self.rhs)
        c_stage = self.c + dt * k1
        for m in range(self.n_comp):
            self.ctx.exchange_ghosts(c_stage[m])

        orig_c = np.copy(self.c)
        self.c[:] = c_stage
        self.compute_rhs()
        k2 = np.copy(self.rhs)

        self.c[:] = orig_c + 0.5 * dt * (k1 + k2)
        sum_c = np.sum(self.c[:, 1:-1], axis=0)
        self.c[:, 1:-1] /= sum_c[np.newaxis, ...]
        self.exchange_all_ghosts()

        self.current_step += 1
        self.current_time += dt

    def step_rk4(self, dt):
        orig_c = np.copy(self.c)

        # Stage 1
        self.compute_rhs()
        k1 = np.copy(self.rhs)

        # Stage 2
        self.c[:] = orig_c + 0.5 * dt * k1
        self.exchange_all_ghosts()
        self.compute_rhs()
        k2 = np.copy(self.rhs)

        # Stage 3
        self.c[:] = orig_c + 0.5 * dt * k2
        self.exchange_all_ghosts()
        self.compute_rhs()
        k3 = np.copy(self.rhs)

        # Stage 4
        self.c[:] = orig_c + dt * k3
        self.exchange_all_ghosts()
        self.compute_rhs()
        k4 = np.copy(self.rhs)

        # Combine
        self.c[:] = orig_c + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        sum_c = np.sum(self.c[:, 1:-1], axis=0)
        self.c[:, 1:-1] /= sum_c[np.newaxis, ...]
        self.exchange_all_ghosts()

        self.current_step += 1
        self.current_time += dt

    def compute_diagnostics(self, dt_interval):
        c_int = self.c[:, 1:-1]
        sum_c = np.sum(c_int, axis=0)
        max_dev_local = np.max(np.abs(sum_c - 1.0))
        min_c_local = np.min(c_int)
        max_c_local = np.max(c_int)

        # Compute local free energy
        bulk_f = np.zeros_like(c_int[0])
        if self.free_energy_type == "regular_solution":
            eps = 1e-8
            for m in range(self.n_comp):
                bulk_f += self.RT * c_int[m] * np.log(np.maximum(c_int[m], eps))
            for i in range(self.n_comp):
                for j in range(i + 1, self.n_comp):
                    bulk_f += self.omega[i, j] * c_int[i] * c_int[j]
        elif self.free_energy_type == "binary_double_well":
            c0 = c_int[0]
            bulk_f = (c0 ** 2) * ((1.0 - c0) ** 2)
        else: # polynomial_multiwell
            for i in range(self.n_comp):
                for j in range(i + 1, self.n_comp):
                    bulk_f += self.W[i, j] * (c_int[i]**2) * (c_int[j]**2)
                bulk_f += self.A[i] * (c_int[i]**2) * ((1.0 - c_int[i])**2)

        # Gradient free energy
        grad_f = np.zeros_like(bulk_f)
        idx = 1.0 / self.dx
        idy = 1.0 / self.dy
        for m in range(self.n_comp):
            c_m = self.c[m]
            dc_dx = 0.5 * (c_m[2:] - c_m[:-2]) * idx
            dc_dy = 0.5 * (np.roll(c_m[1:-1], -1, axis=1) - np.roll(c_m[1:-1], 1, axis=1)) * idy
            grad_sq = dc_dx**2 + dc_dy**2
            grad_f += 0.5 * self.kappa[m, m] * grad_sq

        cell_vol = self.dx * self.dy
        local_total_f = np.sum(bulk_f + grad_f) * cell_vol
        local_avg_c = np.sum(c_int, axis=(1, 2)) * cell_vol

        # MPI Reductions
        total_f = self.ctx.comm.allreduce(local_total_f, op=MPI.SUM)
        global_avg_c = self.ctx.comm.allreduce(local_avg_c, op=MPI.SUM) / (self.ctx.nx * self.ctx.ny * cell_vol)
        max_dev = self.ctx.comm.allreduce(max_dev_local, op=MPI.MAX)
        min_c = self.ctx.comm.allreduce(min_c_local, op=MPI.MIN)
        max_c = self.ctx.comm.allreduce(max_c_local, op=MPI.MAX)

        df_dt = 0.0
        if self.prev_energy is not None and dt_interval > 0.0:
            df_dt = (total_f - self.prev_energy) / dt_interval
        self.prev_energy = total_f

        return {
            "step": self.current_step,
            "time": self.current_time,
            "total_free_energy": total_f,
            "dF_dt": df_dt,
            "min_c": min_c,
            "max_c": max_c,
            "max_unity_dev": max_dev,
            "avg_c": global_avg_c.tolist()
        }

    def write_vtk(self, filepath):
        gathered_c = []
        gathered_mu = []
        for m in range(self.n_comp):
            gc = self.ctx.gather(self.c[m])
            gmu = self.ctx.gather(self.mu[m])
            if self.ctx.rank == 0:
                gathered_c.append(gc)
                gathered_mu.append(gmu)

        if self.ctx.rank == 0:
            os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
            with open(filepath, 'w') as f:
                f.write("# vtk DataFile Version 3.0\n")
                f.write(f"MCCH Simulation Step {self.current_step} Time {self.current_time:.6f}\n")
                f.write("ASCII\n")
                f.write("DATASET STRUCTURED_POINTS\n")
                f.write(f"DIMENSIONS {self.ctx.nx} {self.ctx.ny} 1\n")
                f.write("ORIGIN 0.0 0.0 0.0\n")
                f.write(f"SPACING {self.dx} {self.dy} 1.0\n")
                f.write(f"POINT_DATA {self.ctx.nx * self.ctx.ny}\n")

                for m in range(self.n_comp):
                    f.write(f"SCALARS c{m} double 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    c_data = gathered_c[m] # (nx, ny)
                    # Output order: iterate y, then x
                    for j in range(self.ctx.ny):
                        for i in range(self.ctx.nx):
                            f.write(f"{c_data[i, j]:.6f} ")
                        f.write("\n")

                for m in range(self.n_comp):
                    f.write(f"SCALARS mu{m} double 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    mu_data = gathered_mu[m]
                    for j in range(self.ctx.ny):
                        for i in range(self.ctx.nx):
                            f.write(f"{mu_data[i, j]:.6f} ")
                        f.write("\n")


def run_simulation(config_file):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    with open(config_file, 'r') as f:
        cfg = json.load(f)

    nx = int(cfg.get("nx", 128))
    ny = int(cfg.get("ny", 128))
    periodic = (cfg.get("bc_x", "periodic").lower() == "periodic")

    ctx = SlabMPIContext(comm, nx, ny, nz=1, periodic=periodic)

    n_comp = int(cfg.get("num_components", 3))
    solver = PythonMCCHSolver(
        ctx,
        dx=float(cfg.get("dx", 1.0)),
        dy=float(cfg.get("dy", 1.0)),
        n_comp=n_comp,
        W=cfg.get("W", None),
        A=cfg.get("A", None),
        kappa=cfg.get("kappa", 1.0),
        mobility=cfg.get("mobility", cfg.get("mobility_matrix", cfg.get("mobility_val", 1.0))),
        free_energy_type=cfg.get("free_energy_type", "polynomial_multiwell"),
        RT=float(cfg.get("RT", 1.0)),
        omega=cfg.get("omega", None)
    )

    init_type = cfg.get("initial_condition", "random")
    if init_type == "droplet":
        solver.initialize_droplet(
            droplet_comp=int(cfg.get("droplet_comp", 0)),
            cx=0.5 * nx * solver.dx,
            cy=0.5 * ny * solver.dy,
            radius=float(cfg.get("droplet_radius", 28.0)),
            width=float(cfg.get("droplet_diffuse_width", 2.5))
        )
    else:
        solver.initialize_random(
            cfg.get("c_mean", [1.0 / n_comp] * n_comp),
            noise_amp=float(cfg.get("noise_amp", 0.05)),
            seed=int(cfg.get("seed", 42))
        )

    out_dir = cfg.get("output_dir", "results_python")
    if rank == 0:
        os.makedirs(out_dir, exist_ok=True)
    comm.Barrier()

    csv_file = os.path.join(out_dir, "history.csv")
    if rank == 0:
        with open(csv_file, 'w') as f:
            f.write("step,time,total_free_energy,dF_dt,min_c,max_c,max_unity_dev")
            for m in range(n_comp):
                f.write(f",avg_c{m}")
            f.write("\n")

    dt = float(cfg.get("dt", 0.01))
    total_steps = int(cfg.get("total_steps", 1000))
    output_interval = int(cfg.get("output_interval", 100))
    diag_interval = int(cfg.get("diag_interval", 50))
    integrator = cfg.get("integrator", "rk2").lower()

    diag0 = solver.compute_diagnostics(0.0)
    if rank == 0:
        with open(csv_file, 'a') as f:
            f.write(f"0,0.0000,{diag0['total_free_energy']:.6e},0.0,{diag0['min_c']:.6f},{diag0['max_c']:.6f},{diag0['max_unity_dev']:.6e}")
            for m in range(n_comp):
                f.write(f",{diag0['avg_c'][m]:.6f}")
            f.write("\n")

    solver.write_vtk(os.path.join(out_dir, "solution_00000.vtk"))

    if rank == 0:
        print("======================================================================")
        print("  Multi-Component Cahn-Hilliard Solver (Python/NumPy MPI Reference)")
        print("======================================================================")
        print(f"Loading configuration from: {config_file}")
        print(f"[MPI Slab Decomposition] Total MPI Ranks: {size}, Grid: {nx} x {ny}")
        print(f"[Simulation Setup] Components: {n_comp}, Integrator: {integrator}, dt: {dt}")
        print("------------------------------------------------------------------------------------------------------")
        print(f"{'Step':>8} {'Time':>12} {'Total Energy':>16} {'dF/dt':>14} {'Max |sum(c)-1|':>14} {'c_min':>12} {'c_max':>12} {'MLUPS':>14}")
        print("------------------------------------------------------------------------------------------------------")

    t_start = time.time()
    t_interval_start = t_start

    for step in range(1, total_steps + 1):
        if integrator == "euler":
            solver.step_euler(dt)
        elif integrator == "rk4":
            solver.step_rk4(dt)
        else:
            solver.step_rk2(dt)

        if step % diag_interval == 0 or step == total_steps:
            t_now = time.time()
            elapsed_interval = max(1e-6, t_now - t_interval_start)
            t_interval_start = t_now

            steps_done = diag_interval
            mlups = (steps_done * nx * ny / elapsed_interval) / 1e6

            diag = solver.compute_diagnostics(dt * steps_done)

            if rank == 0:
                with open(csv_file, 'a') as f:
                    f.write(f"{step},{solver.current_time:.4f},{diag['total_free_energy']:.6e},{diag['dF_dt']:.6e},{diag['min_c']:.6f},{diag['max_c']:.6f},{diag['max_unity_dev']:.6e}")
                    for m in range(n_comp):
                        f.write(f",{diag['avg_c'][m]:.6f}")
                    f.write("\n")

                print(f"{step:>8} {solver.current_time:>12.4f} {diag['total_free_energy']:>16.5e} {diag['dF_dt']:>14.5e} {diag['max_unity_dev']:>14.5e} {diag['min_c']:>12.4f} {diag['max_c']:>12.4f} {mlups:>14.2f}")

        if step % output_interval == 0 or step == total_steps:
            solver.write_vtk(os.path.join(out_dir, f"solution_{step:05d}.vtk"))

    t_total = time.time() - t_start
    if rank == 0:
        print("------------------------------------------------------------------------------------------------------")
        print(f"\n[Simulation Summary] Completed {total_steps} steps in {t_total:.2f}s.")
        print(f"Outputs saved to: {out_dir}/")
        print("======================================================================")


def main():
    parser = argparse.ArgumentParser(description="Multi-Component Cahn-Hilliard Python Solver")
    parser.add_argument("-c", "--config", type=str, default="examples/ternary_spinodal.json", help="Path to JSON config file")
    parser.add_argument("--device", type=str, default="cpu", choices=["cpu", "gpu"], help="Execution device (cpu or gpu)")
    parser.add_argument("--gpu", action="store_true", help="Run on GPU using CuPy")
    args = parser.parse_args()

    if os.path.exists(args.config):
        run_simulation(args.config)
    else:
        print(f"Error: Config file {args.config} not found.")

if __name__ == "__main__":
    main()
