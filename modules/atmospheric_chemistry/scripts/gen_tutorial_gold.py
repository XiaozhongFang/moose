#!/usr/bin/env python3
"""Generate F0AM tutorial gold CSV for MOOSE comparison.

Solves the 5-species, 6-reaction tutorial system using scipy BDF
(same algorithm as MATLAB ode15s used by F0AM).  Outputs CSV in
MOOSE postprocessor format.

Usage:
    python3 gen_tutorial_gold.py [--output PATH] [--t-end T] [--dt DT]
"""
import argparse
import numpy as np
from scipy.integrate import solve_ivp
from pathlib import Path


# Mechanism: 5 species, 6 reactions (F0AM LearnF0AM_ODE.mlx)
# y[0]=ONE, y[1]=RO2, y[2]=A, y[3]=B, y[4]=C
K = [0.001, 0.01, 1e-4, 0.1, 0.5, 1e-4]


def dydt(t, y):
    A, B, C = y[2], y[3], y[4]
    R0 = K[0] * A * B      # A+B→C+B
    R1 = K[1] * B           # B→loss
    R2 = K[2] * A           # A→loss
    R3 = K[3] * B * B       # B+B→loss
    R4 = K[4]               # →A (zero order)
    R5 = K[5] * C           # C→loss
    dA = -R0 - R2 + R4
    dB = -R1 - 2*R3
    dC = R0 - R5
    return [0, 0, dA, dB, dC]


def main():
    parser = argparse.ArgumentParser(description="Generate F0AM tutorial gold CSV")
    parser.add_argument("--output", "-o", default="tutorial_f0am_gold.csv",
                        help="Output CSV path (default: tutorial_f0am_gold.csv)")
    parser.add_argument("--t-end", type=float, default=5.0,
                        help="End time (default: 5.0)")
    parser.add_argument("--dt", type=float, default=0.1,
                        help="Time step (default: 0.1)")
    args = parser.parse_args()

    y0 = [1, 0, 1, 0.1, 0]
    n_steps = int(args.t_end / args.dt) + 1
    t_eval = np.linspace(0, args.t_end, n_steps)

    sol = solve_ivp(dydt, (0, args.t_end), y0, method="BDF", t_eval=t_eval,
                    rtol=1e-8, atol=1e-10)

    out = Path(args.output)
    with open(out, "w") as f:
        f.write("time,A_avg,B_avg,C_avg\n")
        for i in range(len(sol.t)):
            f.write(f"{sol.t[i]:.10f},{sol.y[2,i]:.10e},{sol.y[3,i]:.10e},{sol.y[4,i]:.10e}\n")

    print(f"Generated {len(sol.t)} steps to {out}")
    print(f"Final: A={sol.y[2,-1]:.6e}, B={sol.y[3,-1]:.6e}, C={sol.y[4,-1]:.6e}")


if __name__ == "__main__":
    main()
