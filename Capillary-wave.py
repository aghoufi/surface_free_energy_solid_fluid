###Script Python Complet – `capillary_wave.py`**

import numpy as np
from scipy.special import sph_harm

# --- Paramètres globaux ---
N               # Number of atoms
l = 6           # degree of spherical harmonic 
r_cut = 3.0     # cutoff

# --- Positions ---
positions = np.zeros((N, 3))
with open('positions.xyz', 'r') as f:
    for i in range(N):
        x, y, z = map(float, f.readline().split())
        positions[i] = [x, y, z]

# --- find neighbors  ---
def find_neighbors(positions, r_cut):
    N = len(positions)
    neighbors = [[] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if i != j:
                rij = positions[j] - positions[i]
                dist = np.linalg.norm(rij)
                if dist < r_cut:
                    neighbors[i].append(j)
    return neighbors

neighbors = find_neighbors(positions, r_cut)
num_neighbors = [len(nlist) for nlist in neighbors]

# --- Calculation of q6m ---
def compute_q6(positions, neighbors, l=6):
    N = len(positions)
    q6 = np.zeros((N, 2 * l + 1), dtype=complex)

    for i in range(N):
        for j in neighbors[i]:
            rij = positions[j] - positions[i]
            r = np.linalg.norm(rij)
            if r == 0:
                continue
            theta = np.arccos(rij[2] / r)
            phi = np.arctan2(rij[1], rij[0])
            for m in range(-l, l + 1):
                q6[i, m + l] += sph_harm(m, l, phi, theta)
        if len(neighbors[i]) > 0:
            q6[i] /= len(neighbors[i])
    return q6

q6 = compute_q6(positions, neighbors, l)

# --- Calculation de Sij ---
def compute_Sij(q6, neighbors):
    N = len(q6)
    Sij = np.zeros(N)
    for i in range(N):
        norm_i = np.linalg.norm(q6[i])
        if norm_i == 0:
            continue
        for j in neighbors[i]:
            norm_j = np.linalg.norm(q6[j])
            if norm_j == 0:
                continue
            dot = np.vdot(q6[i], q6[j])
            Sij[i] += np.real(dot / (norm_i * norm_j))
        if len(neighbors[i]) > 0:
            Sij[i] /= len(neighbors[i])
    return Sij

Sij = compute_Sij(q6, neighbors)

# --- Calculation of order parameter phi ---
def compute_phi(Sij, neighbors, threshold=0.7):
    N = len(Sij)
    phi_vals = np.zeros(N)
    for i in range(N):
        count = 0
        for j in neighbors[i]:
            if Sij[i] > threshold and Sij[j] > threshold:
                count += 1
        if len(neighbors[i]) > 0:
            phi_vals[i] = count / len(neighbors[i])
    return phi_vals

phi_vals = compute_phi(Sij, neighbors)

# --- Save the results ---
with open('order_parameter.dat', 'w') as f:
    for i in range(N):
        f.write(f"{positions[i,2]:.6f} {phi_vals[i]:.6f}\n")

print("✅ end of calculation and results are stored in 'order_parameter.dat'")
```

#commentaries 
-The positions.xyz file must exist in the same folder and contain exactly N lines (one per atom).
-The script uses scipy for spherical harmonics.
