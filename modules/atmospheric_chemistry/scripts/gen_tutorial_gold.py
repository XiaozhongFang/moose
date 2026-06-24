#!/usr/bin/env python3
"""Generate F0AM tutorial gold CSV for MOOSE comparison.
Solves 5-species, 6-reaction system using scipy.odeint (same as MATLAB ode15s)."""
import numpy as np
from scipy.integrate import solve_ivp

# Mechanism: 5 species, 6 reactions (F0AM LearnF0AM_ODE.mlx)
# y[0]=ONE, y[1]=RO2, y[2]=A, y[3]=B, y[4]=C
k = [0.001, 0.01, 1e-4, 0.1, 0.5, 1e-4]

def dydt(t, y):
    A, B, C = y[2], y[3], y[4]
    R = np.zeros(6)
    R[0] = k[0] * A * B      # A+B→C+B
    R[1] = k[1] * B           # B→loss
    R[2] = k[2] * A           # A→loss
    R[3] = k[3] * B * B       # B+B→loss
    R[4] = k[4]               # →A (zero order)
    R[5] = k[5] * C           # C→loss
    dA = -R[0] - R[2] + R[4]
    dB = -R[1] - 2*R[3]        # B is catalyst in Rx1 (net 0), consumed in Rx2,Rx4
    dC = R[0] - R[5]
    return [0, 0, dA, dB, dC]

y0 = [1, 0, 1, 0.1, 0]
t_span = (0, 5)
t_eval = np.linspace(0, 5, 51)

sol = solve_ivp(dydt, t_span, y0, method='BDF', t_eval=t_eval, rtol=1e-8, atol=1e-10)

with open('/home/fangxiaozhong/git_repo/moose/modules/chemical_reactions/test/tests/actions/gold/tutorial_f0am_gold.csv', 'w') as f:
    f.write('time,A_avg,B_avg,C_avg\n')
    for i in range(len(sol.t)):
        f.write(f'{sol.t[i]:.10f},{sol.y[2,i]:.10e},{sol.y[3,i]:.10e},{sol.y[4,i]:.10e}\n')

print(f'Generated {len(sol.t)} time steps to tutorial_f0am_gold.csv')
print(f'Final: A={sol.y[2,-1]:.6e}, B={sol.y[3,-1]:.6e}, C={sol.y[4,-1]:.6e}')
