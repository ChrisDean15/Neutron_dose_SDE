import math
import numpy as np
import random
from cross_sections import cross_section_1d
from cross_sections import cross_section_2d
from cross_sections import cross_section_3d
import sys


class atom:

    def __init__(self, a0, z0, el_r, ne_r, ne_y, el_a, ne_a, ne_e):
        self.a = a0
        self.z = z0
        self.el_rate = cross_section_1d(el_r)
        self.ne_rate = cross_section_1d(ne_r)
        self.ne_yield = cross_section_1d(ne_y)
        self.el_angle_cdf = cross_section_2d(el_a)
        self.ne_angle_cdf = cross_section_2d(ne_a)
        self.ne_energy_cdf = cross_section_3d(ne_e)

    def cm_to_lab_frame(self, ang, e0, e_delta):
        mp = 938.346  # mass of proton * c^2, MeV
        mn = mp * self.a  # mass of colliding nucleus * c^2, MeV
        u = math.sqrt(e0 * (e0 + mp)) / (e0 + mp + mn)
        g = 1 / math.sqrt(1 - u * u)
        p = math.sqrt((e0 - e_delta) * (e0 - e_delta + 2 * mp))
        e = e0 - e_delta + mp
        v_ratio = u * (e - u * p) / (p - u * e)
        return math.atan(math.sin(ang) / (g * (math.cos(ang) + v_ratio)))


class material:

    def __init__(self, path, atoms):
        data = np.loadtxt(path, max_rows=2)
        self.density = data[0]
        self.I = data[1] / 1e6
        data = np.loadtxt(path, skiprows=2)
        self.x = np.zeros(data.shape[0])
        self.at = []
        for i in range(data.shape[0]):
            self.at.append(atoms[int(data[i, 0])])
            self.x[i] = data[i, 1]

    def bethe_bloch(self, e0):
        mecsq = 0.511  # mass of electron * speed of light squared, MeV
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        betasq = (2 * mpcsq + e0) * e0 / (mpcsq + e0) ** 2
        num = (
            0.3072
            * self.density
            * (math.log(2 * mecsq * betasq / self.I * (1 - betasq)) - betasq)
            * sum([self.x[i] * self.at[i].z for i in range(len(self.at))])
            / betasq
        )
        denom = sum([self.at[i].a * self.x[i] for i in range(len(self.at))])
        return num / denom

    def multiple_scattering_sd(self, e0, dt):
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        pv = (2 * mpcsq + e0) * e0 / (mpcsq + e0)
        betasq = (2 * mpcsq + e0) * e0 / (mpcsq + e0) ** 2
        c = 29979245800  # speed of light, cm / s
        vel = math.sqrt(betasq) * c
        p = pv / math.sqrt(betasq)  # momentum in MeV / c.
        # effective chi_c_sq is just the sum of individual elements
        chi_c_sq = (
            sum(
                [
                    self.x[i] * self.at[i].z * (self.at[i].z + 1) / self.at[i].a
                    for i in range(len(self.at))
                ]
            )
            * 0.157
            * dt
            * self.density
            / (pv * pv)
        )
        # effective chi_a_sq is a weighted average on the log-scale
        chi_a_sq_vec = [
            2.007e-5
            * self.at[i].z ** (2 / 3)
            * (1 + 3.34 * (self.at[i].z / (137 * vel)) ** 2)
            / (p * p)
            for i in range(len(self.at))
        ]
        chi_a_sq = sum(
            [
                self.x[i]
                * self.at[i].z
                * (self.at[i].z + 1)
                * math.log(chi_a_sq_vec[i])
                / self.at[i].a
                for i in range(len(self.at))
            ]
        )
        denom = sum(
            [
                self.x[i] * self.at[i].z * (self.at[i].z + 1.0) / self.at[i].a
                for i in range(len(self.at))
            ]
        )
        chi_a_sq = math.exp(chi_a_sq / denom)
        omega = chi_c_sq / chi_a_sq
        F = 0.98
        v = omega / (2 * (1 - F))
        ret = math.sqrt(chi_c_sq * ((1 + v) * math.log(1 + v) / v - 1)) / (1 + F * F)
        ret = ret * dt / math.sqrt(ret * ret + dt * dt)
        return ret

    def energy_straggling_sd(self, e0):
        alpha = 1 / 137.0
        log_hbar = (
            -21 * math.log(10) + math.log(4.136) - math.log(2 * math.pi)
        )  # MeV * s
        log_c = math.log(29979245800)  # cm / s
        log_avogadro = math.log(6) + 23 * math.log(10)
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        betasq = (2 * mpcsq + e0) * e0 / (mpcsq + e0) ** 2
        a = sum(
            [self.x[i] * self.at[i].a for i in range(len(self.at))]
        )  # average molar mass
        z = sum(
            [self.x[i] * self.at[i].z for i in range(len(self.at))]
        )  # average electrons per molecule
        log_molecule_density = (
            math.log(self.density) + log_avogadro - math.log(a)
        )  # molecules / cm^3
        ret = (
            4
            * math.pi
            * z
            * (1 - betasq / 2)
            / math.sqrt(1 - betasq)
            * math.exp(2 * (math.log(alpha) + log_hbar + log_c) + log_molecule_density)
        )
        return math.sqrt(ret)

    def nonelastic_rate(self, e0):
        log_avogadro = math.log(6) + 23 * math.log(10)
        log_barns_to_cmsq = -24 * math.log(10)
        a = sum(
            [self.x[i] * self.at[i].a for i in range(len(self.at))]
        )  # average molar mass
        ret = sum(
            [
                self.at[i].a * self.x[i] * self.at[i].ne_rate.evaluate(e0)
                for i in range(len(self.at))
            ]
        )
        log_molecule_density = (
            math.log(self.density) + log_avogadro - math.log(a)
        )  # molecules / cm^3
        ret *= math.exp(log_barns_to_cmsq + log_molecule_density) / a
        return ret  # rate per cm

    def elastic_rate(self, e0):
        log_avogadro = math.log(6) + 23 * math.log(10)
        log_barns_to_cmsq = -24 * math.log(10)
        a = sum(
            [self.x[i] * self.at[i].a for i in range(len(self.at))]
        )  # average molar mass
        ret = sum(
            [
                self.at[i].a * self.x[i] * self.at[i].el_rate.evaluate(e0)
                for i in range(len(self.at))
            ]
        )
        log_molecule_density = (
            math.log(self.density) + log_avogadro - math.log(a)
        )  # molecules / cm^3
        ret *= math.exp(log_barns_to_cmsq + log_molecule_density) / a
        return ret  # rate per cm

    def rutherford_rate(self, e0, lb):
        log_ahbarc = math.log(197.3 / 137) - 13 * math.log(10)
        log_avogadro = math.log(6) + 23 * math.log(10)
        mpcsq = 938.346  # mass of proton * speed of light squared, MeV
        pv = (2 * mpcsq + e0) * e0 / (mpcsq + e0)
        a = sum(
            [self.x[i] * self.at[i].a for i in range(len(self.at))]
        )  # average molar mass
        az = sum(
            [
                self.x[i] * self.at[i].a * self.at[i].z * self.at[i].z
                for i in range(len(self.at))
            ]
        )
        log_molecule_density = (
            math.log(self.density) + log_avogadro - math.log(a)
        )  # molecules / cm^3
        sig = (
            math.exp(
                2
                * (
                    log_ahbarc
                    + math.log(math.cos(lb / 2))
                    - math.log(pv)
                    - math.log(math.sin(lb / 2))
                )
                + log_molecule_density
            )
            * math.pi
        )
        ret = sig * az / a
        return ret

    def nonelastic_scatter(self, ang0, e0):
        beta = 2 * math.pi * np.random.uniform()
        rate = sum(
            [
                self.at[i].a * self.x[i] * self.at[i].ne_rate.evaluate(e0)
                for i in range(len(self.at))
            ]
        )
        u = np.random.uniform()
        ind = 0
        tmp = self.at[ind].a * self.x[ind] * self.at[ind].ne_rate.evaluate(e0) / rate
        while tmp < u:
            ind = ind + 1
            tmp = (
                tmp
                + self.at[ind].a
                * self.x[ind]
                * self.at[ind].ne_rate.evaluate(e0)
                / rate
            )
        if np.random.uniform() > self.at[ind].ne_yield.evaluate(e0):
            # proton absorbed & track ends
            ang = 0
            energy = 0
        else:
            alpha = self.at[ind].ne_angle_cdf.sample(e0)
            energy = self.at[ind].ne_energy_cdf.sample(e0, alpha)
            alpha = self.at[ind].cm_to_lab_frame(alpha, e0, e0 - energy)
            ang = ang0
            if 0 <= beta and beta < math.pi / 2:
                ang[0] -= math.atan(math.sin(beta) * math.tan(alpha))
                ang[1] -= math.atan(math.cos(beta) * math.tan(alpha))
            elif math.pi / 2 <= beta and beta < math.pi:
                ang[0] -= math.atan(math.sin(math.pi - beta) * math.tan(alpha))
                ang[1] += math.atan(math.cos(math.pi - beta) * math.tan(alpha))
            elif math.pi <= beta and beta < 3 * math.pi / 2:
                ang[0] += math.atan(math.sin(3 * math.pi / 2 - beta) * math.tan(alpha))
                ang[1] += math.atan(math.cos(3 * math.pi / 2 - beta) * math.tan(alpha))
            else:
                ang[0] += math.atan(math.sin(2 * math.pi - beta) * math.tan(alpha))
                ang[1] -= math.atan(math.cos(2 * math.pi - beta) * math.tan(alpha))
        return ang, energy

    def elastic_scatter(self, ang0, e0):
        beta = 2 * math.pi * np.random.uniform()
        rate = sum(
            [
                self.at[i].a * self.x[i] * self.at[i].el_rate.evaluate(e0)
                for i in range(len(self.at))
            ]
        )
        u = np.random.uniform()
        ind = 0
        tmp = self.at[ind].a * self.x[ind] * self.at[ind].el_rate.evaluate(e0) / rate
        while tmp < u:
            ind = ind + 1
            tmp = (
                tmp
                + self.at[ind].a
                * self.x[ind]
                * self.at[ind].el_rate.evaluate(e0)
                / rate
            )
        alpha = self.at[ind].el_angle_cdf.sample(e0)
        alpha = self.at[ind].cm_to_lab_frame(alpha, e0, 0)
        ang = ang0
        if 0 <= beta and beta < math.pi / 2:
            ang[0] -= math.atan(math.sin(beta) * math.tan(alpha))
            ang[1] -= math.atan(math.cos(beta) * math.tan(alpha))
        elif math.pi / 2 <= beta and beta < math.pi:
            ang[0] -= math.atan(math.sin(math.pi - beta) * math.tan(alpha))
            ang[1] += math.atan(math.cos(math.pi - beta) * math.tan(alpha))
        elif math.pi <= beta and beta < 3 * math.pi / 2:
            ang[0] += math.atan(math.sin(3 * math.pi / 2 - beta) * math.tan(alpha))
            ang[1] += math.atan(math.cos(3 * math.pi / 2 - beta) * math.tan(alpha))
        else:
            ang[0] += math.atan(math.sin(2 * math.pi - beta) * math.tan(alpha))
            ang[1] -= math.atan(math.cos(2 * math.pi - beta) * math.tan(alpha))
        return ang

    def rutherford_scatter(self, ang0, lb):
        beta = 2 * math.pi * np.random.uniform()
        u = np.random.uniform()
        alpha = math.acos(
            (math.cos(lb / 2) ** 2 * (1 - u) - math.sin(lb / 2) ** 2)
            / (1 - u * math.cos(lb / 2) ** 2))
        )
        ang = ang0
        if 0 <= beta and beta < math.pi / 2:
            ang[0] -= math.atan(math.sin(beta) * math.tan(alpha))
            ang[1] -= math.atan(math.cos(beta) * math.tan(alpha))
        elif math.pi / 2 <= beta and beta < math.pi:
            ang[0] -= math.atan(math.sin(math.pi - beta) * math.tan(alpha))
            ang[1] += math.atan(math.cos(math.pi - beta) * math.tan(alpha))
        elif math.pi <= beta and beta < 3 * math.pi / 2:
            ang[0] += math.atan(math.sin(3 * math.pi / 2 - beta) * math.tan(alpha))
            ang[1] += math.atan(math.cos(3 * math.pi / 2 - beta) * math.tan(alpha))
        else:
            ang[0] += math.atan(math.sin(2 * math.pi - beta) * math.tan(alpha))
            ang[1] -= math.atan(math.cos(2 * math.pi - beta) * math.tan(alpha))
        return ang
