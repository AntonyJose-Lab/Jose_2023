# -*- coding: utf-8 -*-
"""
This is a program for simulating the plots that result when a part of a regulatory loop is inhibited by beginning with a system at steady state, providing a pulse of inhibition to one of the entities/sensors and then following the system over time. This program illustrates the effect of duration and strength of the inhibition. For all real biological systems, all parameter values are >0. The simulation is based on the differential equations for the heritable regulatory architecture A. Giving tp a value of 90 and p a value of 0 below simulates a genetic change
Author: Antony Jose
"""
## Import required libraries.
from scipy.integrate import solve_ivp
import numpy as np
import math
import matplotlib.pyplot as plt
import os
plt.close('all')

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures/2023_4_19_Perturbations_of_HRA_A', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables/2023_4_19_Perturbations_of_HRA_A', exist_ok=True)
os.chdir('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/')
curr_dir = os.getcwd()

## time step for the simulation and duration.
t = (0, 100)  # total time

ss = 10  # duration of initial steady state.
tp = 1  # duration of reduction.
d = 0.5  # perturbation factor needed to observe a defect (e.g., 2 = 2x upregulation, 0.5 = 0.5 downregulation)
p = 0.7  # multiplier for describing the extent of experimental perturbation induced (i.e., p*d gives the extent of perturbation)

## Architecture A
HRA = 'A'
x0 = 10
Tx = 0.05
Ty = 0.1

kxy = 0.07
y0 = Tx * x0 / kxy  # 7.1428...
kyx = Tx * Ty / kxy  # 0.071428...

param = [x0, y0, -Tx, -Ty, kxy, kyx]
print(HRA, param)


## Epigenetic perturbation in the amount of each entity in A (x and y)
def HRA_A(t, u):
    du = [0, 0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if t < ss:
        du[0] = param[2] * u[0] + param[4] * u[1]
        du[1] = param[5] * u[0] + param[3] * u[1]
    elif ss < t < (ss + tp):
        u[0] = x0 * p * d
        du[0] = 0
        du[1] = param[5] * u[0] + param[3] * u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[2] * u[0] + param[4] * u[1]
        du[1] = param[5] * u[0] + param[3] * u[1]
    return du


A0 = [param[0], param[1]]
sol = solve_ivp(HRA_A, t, A0, max_step=0.01)
time = sol.t

x = sol.y[0, :]
y = sol.y[1, :]

x_defect = x0 * d
y_defect = y0 * d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color='brown', linewidth=1, label='x')
plt.plot(time, y, color='green', linewidth=1, label='y')
plt.plot(time, x_d, color='brown', linewidth=1, label='x', linestyle='--')
plt.plot(time, y_d, color='green', linewidth=1, label='y', linestyle='--')
plt.xlabel('time', fontsize=8)
plt.ylabel('entities', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))  # do not show grid
ymax = max(max(x), max(y)) + 1
plt.axis([0, 100, 0, ymax])  # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2023_4_19_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2023_4_19_Perturbations_of_HRA_A/', file_name + '_A_tp_1.svg')
fig.savefig(path, dpi=300)

## time step for the simulation and duration.
t = (0,100) # total time

ss = 10 # duration of initial steady state.
tp = 1 # duration of reduction.
d = 0.2 # perturbation factor needed to observe a defect (e.g., 2 = 2x upregulation, 0.5 = 0.5 downregulation)
p = 0.7 # multiplier for describing the extent of experimental perturbation induced (i.e., p*d gives the extent of perturbation)

## Architecture A
HRA = 'A'
x0 = 10
Tx = 0.05
Ty = 0.1

kxy = 0.07
y0 = Tx*x0/kxy # 7.1428...
kyx = Tx*Ty/kxy # 0.071428...

param = [x0, y0, -Tx, -Ty, kxy, kyx]
print(HRA, param)
## Epigenetic perturbation in the amount of each entity in A (x and y)
def HRA_A(t, u):
    du = [0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if t < ss:
        du[0] = param[2]*u[0] + param[4]*u[1]
        du[1] = param[5]*u[0] + param[3]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[5]*u[0] + param[3]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[2]*u[0] + param[4]*u[1]
        du[1] = param[5]*u[0] + param[3]*u[1]
    return du
A0 = [param[0], param[1]]
sol = solve_ivp(HRA_A, t, A0, max_step = 0.01)
time = sol.t

x = sol.y[0, :]
y = sol.y[1, :]

x_defect = x0*0.5
y_defect = y0*0.5
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2023_4_19_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2023_4_19_Perturbations_of_HRA_A/', file_name + '_A_d_0.2.svg')
fig.savefig(path, dpi=300)

## time step for the simulation and duration.
t = (0, 100)  # total time

ss = 10  # duration of initial steady state.
tp = 10  # duration of reduction.
d = 0.5  # perturbation factor needed to observe a defect (e.g., 2 = 2x upregulation, 0.5 = 0.5 downregulation)
p = 0.7  # multiplier for describing the extent of experimental perturbation induced (i.e., p*d gives the extent of perturbation)

## Architecture A
HRA = 'A'
x0 = 10
Tx = 0.05
Ty = 0.1

kxy = 0.07
y0 = Tx * x0 / kxy  # 7.1428...
kyx = Tx * Ty / kxy  # 0.071428...

param = [x0, y0, -Tx, -Ty, kxy, kyx]
print(HRA, param)


## Epigenetic perturbation in the amount of each entity in A (x and y)
def HRA_A(t, u):
    du = [0, 0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if t < ss:
        du[0] = param[2] * u[0] + param[4] * u[1]
        du[1] = param[5] * u[0] + param[3] * u[1]
    elif ss < t < (ss + tp):
        u[0] = x0 * p * d
        du[0] = 0
        du[1] = param[5] * u[0] + param[3] * u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[2] * u[0] + param[4] * u[1]
        du[1] = param[5] * u[0] + param[3] * u[1]
    return du


A0 = [param[0], param[1]]
sol = solve_ivp(HRA_A, t, A0, max_step=0.01)
time = sol.t

x = sol.y[0, :]
y = sol.y[1, :]

x_defect = x0 * d
y_defect = y0 * d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color='brown', linewidth=1, label='x')
plt.plot(time, y, color='green', linewidth=1, label='y')
plt.plot(time, x_d, color='brown', linewidth=1, label='x', linestyle='--')
plt.plot(time, y_d, color='green', linewidth=1, label='y', linestyle='--')
plt.xlabel('time', fontsize=8)
plt.ylabel('entities', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))  # do not show grid
ymax = max(max(x), max(y)) + 1
plt.axis([0, 100, 0, ymax])  # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2023_4_19_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2023_4_19_Perturbations_of_HRA_A/', file_name + '_A_tp_10.svg')
fig.savefig(path, dpi=300)
