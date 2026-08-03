# !/usr/bin/env python

# Copied from https://arxiv.org/pdf/2103.02794.pdf

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from sys import argv

# kB = 1.38064852e-23; u = 1.660538921e-27; M = 3.90e34; L = 6.17e17; T = 3e14; RHO = M/(L**3)

kB = 1.38064852e-23  # Boltzmann constant in SI units
u = 1.660538921e-27  # atomic mass unit in kg
M = 3.90e34  # code mass unit in kg
L = 6.17e17  # code length unit in m
T = 3e14  # code time unit in s
RHO = M / (L ** 3)  # code density unit in kg/m^3
# convert from kg/m^3 to amu/m^3
# RHO = RHO * 6.022e26 * 1e-6
# L = L * 3.240756e-17   # convert from m to pc

def sod(t, temperature):
    # --------o - - Input Section - -o - - - - - - -
    # Sod Shock Tube Problem
    # PrL = 0.04; PrR = 0.005; rhoL = 0.48; rhoR = 0.06; uL = 0.0; uR = 0.0; t_end = t; mu = 0.35; # mu = dt / dx
    # PrL = 0.0192; PrR = 0.0024; rhoL = 0.48; rhoR = 0.06; uL = 0.0; uR = 0.0; t_end = t; mu = 0.35; # mu = dt / dx
    # PrL = 6.62e-15; PrR = 8.28e-16; rhoL = 7.97e-20; rhoR = 9.96e-21; uL = 0.0; uR = 0.0; t_end = t; mu = 0.35; # mu = dt / dx
    # PrL = 2.45e-14; PrR = 3.07e-15; rhoL = 2.95e-19; rhoR = 3.69e-20; uL = 0.0; uR = 0.0; t_end = t; mu = 0.35; # mu = dt / dx
    uL = 0.0
    uR = 0.0
    t_end = t
    mu = 0.35  # mu = dt / dx
    rhoL = 16.0 / 9.0 * RHO
    rhoR = 2.0 / 9.0 * RHO
    PrL = temperature * kB / u * rhoL * 10.0 / 9.0
    PrR = temperature * kB / u * rhoR * 8.0 / 9.0
    # print('rhoL =', rhoL, 'rhoR =', rhoR, 'PrL =', PrL, 'PrR =', PrR)
    # --------o - - END Input Section - -o - - - - - - - - - - - - - - - - -
    gamma = 1.4
    gammaA = gamma - 1.0
    gammaB = 1 / gammaA
    gammaC = gamma + 1.0
    PRL = PrR / PrL
    cR = np.sqrt(gamma * PrR / rhoR)
    cL = np.sqrt(gamma * PrL / rhoL)
    CRL = cR / cL
    machL = (uL - uR) / cL

    def func(p34):
        wortel = np.sqrt(2 * gamma * (gammaA + gammaC * p34))
        yy = (gammaA * CRL * (p34 - 1)) / wortel
        yy = (1 + machL * gammaA / 2 - yy) ** (2 * gamma / gammaA)
        y = yy / p34 - PRL
        return y

    p34 = fsolve(func, 3.0)  # p34 = p3 / p4
    # print(p34)
    p3 = p34 * PrR
    alpha = gammaC / gammaA
    rho3 = rhoR * (1 + alpha * p34) / (alpha + p34)
    rho2 = rhoL * (p34 * PrR / PrL) ** (1 / gamma)
    u2 = uL - uR + (2 / gammaA) * cL * (1 - (p34 * PrR / PrL) ** (gammaA / (2 * gamma)))
    c2 = np.sqrt(gamma * p3 / rho2)
    spos = (0.25 * L + t_end * cR * np.sqrt(
        gammaA / (2 * gamma) + gammaC / (2 * gamma) * p34) + t_end * uR)  # Shock position
    conpos = 0.25 * L + u2 * t_end + t_end * uR  # Position of contact discontinuity
    pos1 = 0.25 * L + (uL - cL) * t_end  # Start of expansion fan
    pos2 = 0.25 * L + (u2 + uR - c2) * t_end  # End of expansion fan
    xgrid = np.linspace(0, 0.5 * L, 500)
    PrE = np.zeros((1, len(xgrid)))
    uE = np.zeros((1, len(xgrid)))
    rhoE = np.zeros((1, len(xgrid)))
    machE = np.zeros((1, len(xgrid)))
    cE = np.zeros((1, len(xgrid)))
    xgrid = np.matrix(xgrid)
    for i in range(0, xgrid.size):
        if xgrid[0, i] <= pos1:
            PrE[0, i] = PrL
            rhoE[0, i] = rhoL
            uE[0, i] = uL
            cE[0, i] = np.sqrt(gamma * PrE[0, i] / rhoE[0, i])
            machE[0, i] = uE[0, i] / cE[0, i]
        elif xgrid[0, i] <= pos2:
            PrE[0, i] = (PrL * (1 + (pos1 - xgrid[0, i]) / (cL * alpha * t_end)) ** (2 * gamma / gammaA))
            rhoE[0, i] = (rhoL * (1 + (pos1 - xgrid[0, i]) / (cL * alpha * t_end)) ** (2 / gammaA))
            uE[0, i] = uL + (2 / gammaC) * (xgrid[0, i] - pos1) / t_end
            cE[0, i] = np.sqrt(gamma * PrE[0, i] / rhoE[0, i])
            machE[0, i] = uE[0, i] / cE[0, i]
        elif xgrid[0, i] <= conpos:
            PrE[0, i] = p3
            rhoE[0, i] = rho2
            uE[0, i] = u2 + uR
            cE[0, i] = np.sqrt(gamma * PrE[0, i] / rhoE[0, i])
            machE[0, i] = uE[0, i] / cE[0, i]
        elif xgrid[0, i] <= spos:
            PrE[0, i] = p3
            rhoE[0, i] = rho3
            uE[0, i] = u2 + uR
            cE[0, i] = np.sqrt(gamma * PrE[0, i] / rhoE[0, i])
            machE[0, i] = uE[0, i] / cE[0, i]
        else:
            PrE[0, i] = PrR
            rhoE[0, i] = rhoR
            uE[0, i] = uR
            cE[0, i] = np.sqrt(gamma * PrE[0, i] / rhoE[0, i])
            machE[0, i] = uE[0, i] / cE[0, i]
    entropy_E = np.log(PrE / rhoE ** gamma)
    # calculate temperature from ideal gas law
    # P = rho * kB * T / u
    # -> T = P * u / (rho * kB)
    tempE = PrE * u / (rhoE * kB)
    return np.array(xgrid - 0.25 * L) * 3.240756e-17, np.array(PrE) * 6.022e26 * 1e-6, np.array(uE), np.array(
        rhoE) * 6.022e26 * 1e-6, np.array(machE), np.array(entropy_E), np.array(tempE)


""" ---o--- File Write Section ---o--- """


def file_write():
    (xgrid, PrE, uE, rhoE, machE, entropy_E, tempE) = sod(0.25)

    file = open("Riemann.txt", "w")
    for i in range(0, xgrid.size):
        file.write(" %.6f " % (xgrid[0, i]) + ' ' + str(" %.6f " % rhoE[0, i]) + ' ' \
                   + (" %.6f " % PrE[0, i]) + ' ' + str(" %.6f " % uE[0, i]) + ' ' \
                   + (" %.6f " % machE[0, i]) + ' ' + str(" %.6f " % entropy_E[0, i]))
        file.write(" \n ")
    file.close()


# RHO *= 6.022e26 * 1e-6
# L *= 3.240756e-17

def plot():
    (xgrid, PrE, uE, rhoE, machE, entropy_E, tempE) = sod(0.01 * T)
    # plt.plot(xgrid[0, :], rhoE[0, :], 'r', label='Density')
    # plt.plot(xgrid[0,:], PrE[0,:], 'b', label='Pressure')
    # plt.plot(xgrid[0,:], uE[0,:], 'g', label='Velocity')
    # plt.plot(xgrid[0,:], machE[0,:], 'y', label='Mach')
    # plt.plot(xgrid[0,:], entropy_E[0,:], 'k', label='Entropy')
    plt.plot(xgrid[0,:], tempE[0,:], 'k', label='Temperature')
    plt.legend(loc='upper right')
    plt.show()


if __name__ == "__main__":
    plot()
