import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from proton_path import proton_path
from material import atom
from material import material
import seaborn
import sys

"""Initilise starting configoration"""
x0 = [0, 0, 0]
omega0 = [np.pi / 2, 0]
E0 = 100
E0_sd = 0.0001
dt = 0.05

N = 10000  # Number of replicate protons
# Length of air before beam enters water
air_gap = 0.0

"""Read in atoms and materials"""
atom_data = np.loadtxt("../Materials/atoms.txt", dtype=str)
atoms = []
for i in range(atom_data.shape[0]):
    if atom_data[i, 0] == "hydrogen":
        tmp = atom(
            int(atom_data[i, 2]),
            int(atom_data[i, 1]),
            "../Splines/" + atom_data[i, 0] + "_el_rate.txt",
            "",
            "",
            "../Splines/" + atom_data[i, 0] + "_el_angle_cdf.txt",
            "",
            "",
        )
    else:
        tmp = atom(
            int(atom_data[i, 2]),
            int(atom_data[i, 1]),
            "../Splines/" + atom_data[i, 0] + "_el_rate.txt",
            "../Splines/" + atom_data[i, 0] + "_ne_rate.txt",
            "../Splines/" + atom_data[i, 0] + "_ne_yield.txt",
            "../Splines/" + atom_data[i, 0] + "_el_angle_cdf.txt",
            "../Splines/" + atom_data[i, 0] + "_ne_angle_cdf.txt",
            "../Splines/" + atom_data[i, 0] + "_ne_energy_cdf.txt",
        )
    atoms.append(tmp)

material_names = np.loadtxt("../Materials/materials.txt", dtype=str)
materials = []
for name in material_names:
    tmp = material("../Materials/" + name + ".txt", atoms)
    materials.append(tmp)

"""Simulate paths"""
# Radius of the circle
nozzle_radius = 0.05
x_sd = 0.05

# Initial energies and points in the nozzle
x = np.zeros(N)
y = x
z = x
for i in range(N):
    y[i] = np.random.normal(0, x_sd)
    z[i] = np.random.normal(0, x_sd)
    while y[i] * y[i] + z[i] * z[i] > nozzle_radius * nozzle_radius:
        y[i] = np.random.normal(0, x_sd)
        z[i] = np.random.normal(0, x_sd)
E0 = np.random.normal(E0, E0_sd, N)

# Combine the coordinates
x = np.column_stack((x, y, z))

# Loop over points
result = proton_path(x[0], omega0, E0[0], dt, materials, air_gap)
x_res = [result.points[0][0][0]]
y_res = [result.points[0][0][1]]
s_res = [0]
for i in range(1, len(result.points)):
    x_res.append(result.points[i][0][0])
    y_res.append(result.points[i][0][1])
    s_res.append(result.points[i - 1][2] - result.points[i][2])
# Plot limits need to be set manually for sensible results
x_range = (0, 9)
y_range = (-5, 5)
heatmap, xedges, yedges = np.histogram2d(
    x_res, y_res, bins=90, range=[x_range, y_range], weights=s_res, density=False
)

for i in range(1, N):
    result = proton_path(x[i], omega0, E0[i], dt, materials, air_gap)
    x_res = [result.points[0][0][0]]
    y_res = [result.points[0][0][1]]
    s_res = [0]
    for i in range(1, len(result.points)):
        x_res.append(result.points[i][0][0])
        y_res.append(result.points[i][0][1])
        s_res.append(result.points[i - 1][2] - result.points[i][2])
    heatmap_tmp, xedges, yedges = np.histogram2d(
        x_res, y_res, bins=90, range=[x_range, y_range], weights=s_res, density=False
    )
    heatmap = heatmap + heatmap_tmp

np.divide(heatmap, N)

fig, ax = plt.subplots(figsize=(8, 8))  # Adjust the size as needed
seaborn.heatmap(
    np.transpose(heatmap),
    cmap="viridis",
    norm=mcolors.LogNorm(),
    yticklabels=[round(x, 1) for x in np.linspace(-5, 5, 90)],
    xticklabels=[round(x, 1) for x in np.linspace(0, 9, 90)],
)

plt.savefig("heatmap.pdf")

fig, ax = plt.subplots()
plt.plot(np.divide(np.sum(heatmap, axis=1), N))

plt.savefig("bragg-peak.pdf")
