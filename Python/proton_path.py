import math
import numpy as np
import random
from wright_fisher import wright_fisher


class proton_path:

    def __init__(self, x0, omega0, E0, dt, materials, air_gap):
        self.dt = dt  # time step
        self.points = [
            [x0, omega0, E0]
        ]  # initialise paths with start position, angle, energy
        self.scatter(materials, air_gap)

    def scatter(self, materials, air_gap):
        absorption_e = 0.037
        material_index = 0
        crossed = 0
        while self.points[-1][2] > absorption_e:
            self.spherical_BM(
                self.points[-1][0],
                self.points[-1][1],
                self.points[-1][2],
                self.dt,
                materials[material_index],
            )
            if crossed == 0 and self.points[-1][0][0] >= air_gap:
                material_index = material_index + 1
                crossed = 1
            if self.points[-1][2] > absorption_e:
                non_e_jump = materials[material_index].nonelastic_rate(
                    self.points[-1][2]
                )
                e_jump = materials[material_index].elastic_rate(self.points[-1][2])
                rutherford_lb = 2.5 * materials[material_index].multiple_scattering_sd(
                    self.points[-1][2], self.dt
                )
                rutherford_jump = 0
                if rutherford_lb < math.pi:
                    rutherford_jump = materials[material_index].rutherford_rate(
                        self.points[-1][2], rutherford_lb
                    )
                alpha = non_e_jump + e_jump + rutherford_jump
                if random.random() < 1 - math.exp(-self.dt * alpha):
                    u = random.random()
                    if u < rutherford_jump / alpha:
                        self.points[-1][1] = materials[
                            material_index
                        ].rutherford_scatter(self.points[-1][1], rutherford_lb)
                    elif u < (rutherford_jump + e_jump) / alpha:
                        self.points[-1][1] = materials[material_index].elastic_scatter(
                            self.points[-1][1], self.points[-1][2]
                        )
                    else:
                        self.points[-1][1], self.points[-1][2] = materials[
                            material_index
                        ].nonelastic_scatter(self.points[-1][1], self.points[-1][2])
        return self.points

    def spherical_BM(self, x, omega, e0, dt, material):
        wf = wright_fisher(material.multiple_scattering_sd(e0, self.dt) ** 2)
        theta = 2 * math.pi * random.random()
        Y = np.array((np.cos(theta), np.sin(theta)))

        # Set u = (e_d - z) / |e_d - z|
        z = [
            [
                np.sin(omega[0]) * np.cos(omega[1]),
                np.sin(omega[0]) * np.sin(omega[1]),
                np.cos(omega[0]),
            ]
        ]
        e_d = np.zeros(3)
        e_d[2] = 1.0
        dist = np.linalg.norm(e_d - z)
        if dist > 1e-9:
            u = (e_d - z) / dist
        else:
            u = np.ones(3)
            if z[0] > 0:
                u[0] = -1
            if z[1] > 0:
                u[1] = -1

        # Compute O(z) = I - 2uu_T
        O_z = np.identity(3) - 2 * np.outer(u, u)
        c = np.concatenate(
            (np.array([2 * np.sqrt(wf * (1 - wf)) * Y]), np.array([1 - 2 * wf])),
            axis=None,
        )
        w = np.dot(O_z, c)

        angle = [np.arccos(w[2]), np.arctan2(w[1], w[0])]
        pos = [
            x[0] + np.sin(angle[0]) * np.cos(angle[1]) * dt,
            x[1] + np.sin(angle[0]) * np.sin(angle[1]) * dt,
            x[2] + np.cos(angle[0]) * dt,
        ]
        delta = np.linalg.norm(np.array(self.points[-1][0]) - np.array(pos))
        energy_mu = material.bethe_bloch(e0) * delta
        energy_sd = math.sqrt(delta) * material.energy_straggling_sd(e0)
        energy = max(0, e0 - max(random.gauss(energy_mu, energy_sd), 0))

        self.points.append([pos, angle, energy])
