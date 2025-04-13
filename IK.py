import numpy as np
import matplotlib.pyplot as plt

# --- System parameters ---
N                                  # number of particles
box[1], box[2], box[3]             # Box lengths                 
box = np.array([box[1], box[2], box[3] ]) # Box dimensions
bins = 100                         # number of Bins according to z (normal to the interface)
dz = box[2] / bins
z_edges = np.linspace(0, box[2], bins + 1)
z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
rcut = 12.0                         # Cut-off pour LJ (Angstrom)

# --- Lennard-Jones parameters --- (data file)
#epsilon 
#sigma 

# --- positions and velocities ---  (data file)
#positions 
#velocities 
#masses 

# --- periodic boundary conditions ---
def minimum_image(rij, box):
    return rij - np.round(rij / box) * box

# --- Lennard Jones forces ---
def lj_force(rij, epsilon, sigma):
    r2 = np.dot(rij, rij)
    if r2 == 0 or np.sqrt(r2) > rcut:
        return np.zeros(3)
    sr2 = (sigma ** 2) / r2
    sr6 = sr2 ** 3
    sr12 = sr6 ** 2
    force_mag = 48 * epsilon * (sr12 - 0.5 * sr6) / r2
    return force_mag * rij

# --- Forces calculations ---
forces = np.zeros((N, 3))
for i in range(N - 1):
    for j in range(i + 1, N):
        rij = minimum_image(positions[i] - positions[j], box)
        fij = lj_force(rij, epsilon, sigma)
        forces[i] += fij
        forces[j] -= fij  # Newton

# --- Kinetic pressure ---
P_kin = np.zeros(bins)
for i in range(N):
    z = positions[i, 2]
    bin_index = int(z // dz)
    if bin_index < bins:
        P_kin[bin_index] += masses[i] * velocities[i, 2] ** 2

P_kin /= dz * box[0] * box[1]

# --- configurational pressure (Irving-Kirkwood) ---
P_vir = np.zeros(bins)
for i in range(N - 1):
    for j in range(i + 1, N):
        rij = minimum_image(positions[i] - positions[j], box)
        r = np.linalg.norm(rij)
        if r == 0 or r > rcut:
            continue
        fij = lj_force(rij, epsilon, sigma)

        # Integration on i → j
        n_steps = 20
        for s in np.linspace(0, 1, n_steps):
            r_s = positions[i] - s * rij
            z = r_s[2] % box[2]
            bin_index = int(z // dz)
            if bin_index < bins:
                contrib = rij[2] * fij[2] / n_steps
                P_vir[bin_index] += contrib

P_vir *= 0.5 / (dz * box[0] * box[1])

# --- Pression totale ---
P_total = P_kin + P_vir


