import math
import numpy as np
import random


class cross_section_1d:

    def __init__(self, path):
        self.energy = []
        self.rate = []
        if path != "":
            data = np.loadtxt(path)
            self.energy = data[0]
            self.rate = data[1]

    def evaluate(self, e0):
        ret = 0
        if len(self.energy) > 0:
            if e0 < self.energy[0]:
                ret = self.rate[0]
            elif e0 > self.energy[-1]:
                ret = self.rate[-1]
            else:
                i = 0
                while self.energy[i] < e0:
                    i = i + 1
                ret = (
                    (self.energy[i] - e0) * self.rate[i - 1]
                    + (e0 - self.energy[i - 1]) * self.rate[i]
                ) / (self.energy[i] - self.energy[i - 1])
        return ret


class cross_section_2d:

    def __init__(self, path):
        self.energy = []
        self.angle = []
        self.cdf = []
        if path != "":
            self.energy = np.loadtxt(path, max_rows=1)
            self.angle = np.loadtxt(path, skiprows=1, max_rows=1)
            self.cdf = np.loadtxt(path, skiprows=2)

    def sample(self, e0):
        u = random.random()
        c = 1
        if e0 <= self.energy[0]:
            while u > self.cdf[0][c]:
                c = c + 1

            ret = (
                (self.cdf[0][c] - u) * self.angle[c - 1]
                + (u - self.cdf[0][c - 1]) * self.angle[c]
            ) / (self.cdf[0][c] - self.cdf[0][c - 1])
        elif e0 >= self.energy[-1]:
            while u > self.cdf[-1][c]:
                c = c + 1
            ret = (
                (self.cdf[-1][c] - u) * self.angle[c - 1]
                + (u - self.cdf[-1][c - 1]) * self.angle[c]
            ) / (self.cdf[-1][c] - self.cdf[-1][c - 1])
        else:
            r = 1
            while e0 > self.energy[r]:
                r = r + 1
            while u > self.cdf[r][c]:
                c = c + 1
            ret = (
                (self.cdf[r][c] - u) * self.angle[c - 1]
                + (u - self.cdf[r][c - 1]) * self.angle[c]
            ) / (self.cdf[r][c] - self.cdf[r][c - 1])
            c = 1
            while u > self.cdf[r - 1][c]:
                c = c + 1
            ret2 = (
                (self.cdf[r - 1][c] - u) * self.angle[c - 1]
                + (u - self.cdf[r - 1][c - 1]) * self.angle[c]
            ) / (self.cdf[r - 1][c] - self.cdf[r - 1][c - 1])
            ret = ((self.energy[r] - e0) * ret2 + (e0 - self.energy[r - 1]) * ret) / (
                self.energy[r] - self.energy[r - 1]
            )
        return ret


class cross_section_3d:

    def __init__(self, path):
        self.energy = []
        self.angle = []
        self.cdf = []
        self.exit_energy = []
        if path != "":
            self.energy = np.loadtxt(path, max_rows=1)
            self.angle = np.loadtxt(path, skiprows=1, max_rows=1)
            data = np.loadtxt(path, skiprows=2)
            self.cdf = data[0::2, :]
            self.exit_energy = data[1::2, :]

    def find_energy_index(self, e0):
        ret = 0
        while (ret < len(self.energy)) and (e0 > self.energy[ret]):
            ret = ret + 1
        return ret

    def find_angle_index(self, ang):
        ret = 0
        while (ret < len(self.angle)) and (ang > self.angle[ret]):
            ret = ret + 1
        return ret

    def sample(self, e0, ang):
        u = random.random()
        e_row = self.find_energy_index(e0)
        ang_row = self.find_angle_index(ang)
        c = 1
        if e_row == 0:
            if ang_row == 0:
                while u > self.cdf[0][c]:
                    c = c + 1
                ret = (
                    (self.cdf[0][c] - u) * min(e0, self.exit_energy[0][c - 1])
                    + (u - self.cdf[0][c - 1]) * fmin(e0, self.exit_energy[0][c])
                ) / (self.cdf[0][c] - self.cdf[0][c - 1])
            elif ang_row == len(self.angle):
                r = ang_row - 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
            else:
                r = ang_row
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                c = 1
                while u > self.cdf[r - 1][c]:
                    c = c + 1
                ret2 = (
                    (self.cdf[r - 1][c] - u) * min(e0, self.exit_energy[r - 1][c - 1])
                    + (u - self.cdf[r - 1][c - 1]) * min(e0, self.exit_energy[r - 1][c])
                ) / (self.cdf[r - 1][c] - self.cdf[r - 1][c - 1])
                ret = (
                    (self.angle[r] - ang) * ret2 + (ang - self.angle[r - 1]) * ret
                ) / (self.angle[r] - self.angle[r - 1])
        elif e_row == len(self.energy):
            if ang_row == 0:
                r = (e_row - 1) * len(self.angle)
                while u > cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
            elif ang_row == len(self.angle):
                r = (e_row - 1) * len(self.angle) + ang_row - 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
            else:
                r = (e_row - 1) * len(self.angle) + ang_row
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                c = 1
                while u > self.cdf[r - 1][c]:
                    c = c + 1
                ret2 = (
                    (self.cdf[r - 1][c] - u) * min(e0, self.exit_energy[r - 1][c - 1])
                    + (u - self.cdf[r - 1][c - 1]) * min(e0, self.exit_energy[r - 1][c])
                ) / (self.cdf[r - 1][c] - self.cdf[r - 1][c - 1])
                ret = (
                    (self.angle[ang_row] - ang) * ret2
                    + (ang - self.angle[ang_row - 1]) * ret
                ) / (self.angle[ang_row] - self.angle[ang_row - 1])
        else:
            if ang_row == 0:
                r = e_row * len(self.angle)
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                c = 1
                r = (e_row - 1) * len(self.angle)
                while u > self.cdf[r][c]:
                    c = c + 1
                ret2 = (
                    (self.cdf[r][c] - u) * min(e, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                ret = (
                    (self.energy[e_row] - e0) * ret2
                    + (e0 - self.energy[e_row - 1]) * ret
                ) / (self.energy[e_row] - self.energy[e_row - 1])
            elif ang_row == len(self.angle):
                r = e_row * len(self.angle) + ang_row - 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                c = 1
                r = (e_row - 1) * len(self.angle) + ang_row - 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret2 = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                ret = (
                    (self.energy[e_row] - e0) * ret2
                    + (e0 - self.energy[e_row - 1]) * ret
                ) / (self.energy[e_row] - self.energy[e_row - 1])
            else:
                r = (e_row - 1) * len(self.angle) + ang_row - 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                r = (e_row - 1) * len(self.angle) + ang_row
                c = 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret2 = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                ret = (
                    (self.angle[ang_row] - ang) * ret
                    + (ang - self.angle[ang_row - 1]) * ret2
                ) / (self.angle[ang_row] - self.angle[ang_row - 1])
                r = e_row * len(self.angle) + ang_row - 1
                c = 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret3 = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                r = e_row * len(self.angle) + ang_row
                c = 1
                while u > self.cdf[r][c]:
                    c = c + 1
                ret4 = (
                    (self.cdf[r][c] - u) * min(e0, self.exit_energy[r][c - 1])
                    + (u - self.cdf[r][c - 1]) * min(e0, self.exit_energy[r][c])
                ) / (self.cdf[r][c] - self.cdf[r][c - 1])
                ret3 = (
                    (self.angle[ang_row] - ang) * ret3
                    + (ang - self.angle[ang_row - 1]) * ret4
                ) / (self.angle[ang_row] - self.angle[ang_row - 1])
                ret = (
                    (self.energy[e_row] - e0) * ret
                    + (e0 - self.energy[e_row - 1]) * ret3
                ) / (self.energy[e_row] - self.energy[e_row - 1])
        return ret
