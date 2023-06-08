# -*- coding: utf-8 -*-
"""
This is a program for simulating the plots that result after each genetic perturbation by begining with a system at steady state and then following the system over time. For all real biological systems, all parameter values are >0.
The simulation is based on the differental equations for each heritable regulatory architecture. Giving tp a value of 90 and p a value of 0 below simulates a genetic change
Author: Antony Jose
"""
## Import required libraries.
from scipy.integrate import solve_ivp
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import csv
plt.close('all')

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures', exist_ok=True)
os.chdir('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/')
curr_dir = os.getcwd()

## Writing out the parameters from the simulation as a .xls file with different sheets.
f = open('Tables/2022_8_17_SteadyStates_of_Heritable_Architectures/2022_8_26_Illustrative_SteadyStates_for_26_HRAs_params_and_genetic_perturb.csv', 'w', encoding='UTF8', newline='')
writer = csv.writer(f)

## time step for the simulation and duration.
t = (0,300) # total time

ss = 10 # duration of initial steady state.
tp = 90 # duration of epigenetic perturbation. (with p=0, this makes it a genetic mutation)
d = 0.5 # perturbation factor needed to observe a defect (e.g., 2 = 2x upregulation, 0.5 = 0.5 downregulation)
p = 0 # multiplier fordescribing the extent of experimental perturbation induced (i.e., p*d gives the extent of perturbation)


## Reference exponential decay
def exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz):
    turnover_x = np.zeros(30000)
    turnover_y = np.zeros(30000)
    turnover_z = np.zeros(30000)
    time_decay = np.zeros(30000)
    result = [time_decay, turnover_x, turnover_y, turnover_z]
    if HRA == 'A':
        i = 0
        for t in np.arange(0,300,0.01):
            if t < ss:
                turnover_x[i] = x0
                turnover_y[i] = y0
                time_decay[i] = t
                i = i+1
            elif t < 101:
                #turnover_x[i] = x0*math.exp(-Tx*(t-ss))
                turnover_y[i] = y0*math.exp(-Ty*(t-ss))
                time_decay[i] = t
                i = i+1
            elif 100 < t < (ss + 101):
                turnover_x[i] = x0
                turnover_y[i] = y0
                time_decay[i] = t
                i = i+1
            elif t > (ss + 100):
                turnover_x[i] = x0*math.exp(-Tx*(t-ss-100))
                #turnover_y[i] = y0*math.exp(-Ty*(t-ss-100))
                time_decay[i] = t
                i = i+1
            result = [time_decay, turnover_x, turnover_y]

    else:
        i = 0
        for t in np.arange(0,300,0.01):
            if t < ss:
                turnover_x[i] = x0
                turnover_y[i] = y0
                turnover_z[i] = z0
                time_decay[i] = t
                i = i+1
            elif t < 101:
                #turnover_x[i] = x0*math.exp(-Tx*(t-ss))
                turnover_y[i] = y0*math.exp(-Ty*(t-ss))
                turnover_z[i] = z0*math.exp(-Tz*(t-ss))
                time_decay[i] = t
                i = i+1
            elif 100 < t < (ss + 101):
                turnover_x[i] = x0
                turnover_y[i] = y0
                turnover_z[i] = z0
                time_decay[i] = t
                i = i+1
            elif (ss + 100) < t < 201:
                turnover_x[i] = x0*math.exp(-Tx*(t-ss-100))
                #turnover_y[i] = y0*math.exp(-Ty*(t-ss-100))
                turnover_z[i] = z0*math.exp(-Tz*(t-ss-100))
                time_decay[i] = t
                i = i+1
            elif 200 < t < (ss + 201):
                turnover_x[i] = x0
                turnover_y[i] = y0
                turnover_z[i] = z0
                time_decay[i] = t
                i = i+1
            elif t > (ss + 200):
                turnover_x[i] = x0*math.exp(-Tx*(t-ss-200))
                turnover_y[i] = y0*math.exp(-Ty*(t-ss-200))
                #turnover_z[i] = z0*math.exp(-Tz*(t-ss-200))
                time_decay[i] = t
                i = i+1
            result = [time_decay, turnover_x, turnover_y, turnover_z]
    return result

## Architecture A
HRA = 'A'
x0 = 10
Tx = 0.05
Ty = 0.1

kxy = 0.07
y0 = Tx*x0/kxy # 7.1428...
kyx = Tx*Ty/kxy # 0.071428...

param = [x0, y0, -Tx, -Ty, kxy, kyx]
param_name = 'HRA', 'x0, y0, -Tx, -Ty, kxy, kyx'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
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
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[2]*u[0] + param[4]*u[1]
        du[1] = 0
    elif (100 + ss + tp) < t :
        du[0] = param[2]*u[0] + param[4]*u[1]
        du[1] = param[5]*u[0] + param[3]*u[1]
    return du
A0 = [param[0], param[1]]
sol = solve_ivp(HRA_A, t, A0, max_step = 0.01)
time = sol.t

x = sol.y[0, :]
y = sol.y[1, :]

x_defect = x0*d
y_defect = y0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))

A_exp = exp_decay(HRA, x0, y0, 0, Tx, Ty, 0)
time_decay = A_exp[0]
x_exp = A_exp[1]
y_exp = A_exp[2]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Axy.svg')
fig.savefig(path, dpi=300)

# ## Architecture B
HRA = 'B'
x0 = 10
Tx = 0.05
Ty = 0.1
Tz = 0.15

kxy = 0.07
kzy = 0.09

y0 = Tx*x0/kxy # 7.143...
kyx = Tx*Ty/kxy # 0.07143...
z0 = x0*Tx*kzy/(kxy*Tz) # 4.285...

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in B (x, y, z)
def HRA_B(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_B, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

B_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = B_exp[0]
x_exp = B_exp[1]
y_exp = B_exp[2]
z_exp = B_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Bxyz.svg')
fig.savefig(path, dpi=300)

## Architecture C
HRA = 'C'
x0 = 10
Tx = 0.05
Ty = 0.1
Tz = 0.15

kxy = 0.08
kzx = 0.05

y0 = Tx*x0/kxy # 6.25
kyz = Tx*Ty*Tz/(kxy*kzx) # 0.1875
z0 = x0*kzx/Tz # 3.33...

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyz, kzx]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyz, kzx'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in C (x, y, z)
def HRA_C(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_C, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

C_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = C_exp[0]
x_exp = C_exp[1]
y_exp = C_exp[2]
z_exp = C_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Cxyz.svg')
fig.savefig(path, dpi=300)

# ## Architecture D
HRA = 'D'
x0 = 10
Tx = 0.05
Ty = 0.1
Tz = 0.15

kxy = 0.08
kzx = 0.05
kzy = 0.03

y0 = Tx*x0/kxy # 6.25
kyx = Tx*Ty/kxy # 0.0625
z0 = (x0/Tz)*(kzx + kzy*kyx/Ty) # 4.5833...

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in D (x, y, z)
def HRA_D(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_D, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

D_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = D_exp[0]
x_exp = D_exp[1]
y_exp = D_exp[2]
z_exp = D_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Dxyz.svg')
fig.savefig(path, dpi=300)

# ## Architecture E
HRA = 'E'
x0 = 10
Tx = 0.05
Ty = 0.1
Tz = 0.15

kxy = 0.08
kzx = 0.01
kzy = 0.04

y0 = Tx*x0/kxy # 6.25
kyx = Tx*Ty/kxy # 0.0625
z0 = (x0/Tz)*(-kzx + kzy*kyx/Ty) # 1.33...

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kzx]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kzx'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in E (x, y, z)
def HRA_E(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_E, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

E_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = E_exp[0]
x_exp = E_exp[1]
y_exp = E_exp[2]
z_exp = E_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Exyz.svg')
fig.savefig(path, dpi=300)


# ## Architecture F
HRA = 'F'
x0 = 10
Tx = 0.05
Ty = 0.11
Tz = 0.14

kxy = 0.08
kyx = 0.05
kzy = 0.08

y0 = Tx*x0/kxy # 50
kyz = (Tx*Ty*Tz - Tz*kyx*kxy)/(Tx*kzy) # 0.005
z0 = x0*Tx*kzy/(kxy*Tz) # 10

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in F (x, y, z)
def HRA_F(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_F, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

F_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = F_exp[0]
x_exp = F_exp[1]
y_exp = F_exp[2]
z_exp = F_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Fxyz.svg')
fig.savefig(path, dpi=300)


## Architecture G
HRA = 'G'
x0 = 10
Tx = 0.08
Ty = 0.11
Tz = 0.05

kxy = 0.14
kyx = 0.08
kzy = 0.05

y0 = Tx*x0/kxy #
kyz = (Tz*kyx*kxy - Tx*Ty*Tz)/(Tx*kzy) #
z0 = x0*Tx*kzy/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in H (x, y, z)
def HRA_G(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_G, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

G_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = G_exp[0]
x_exp = G_exp[1]
y_exp = G_exp[2]
z_exp = G_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Gxyz.svg')
fig.savefig(path, dpi=300)

## Architecture H
HRA = 'H'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.02

y0 = x0*kyx/Ty #
kxz = (Tz*Ty*Tx - Tz*kxy*kyx)/(kyx*kzy) #
z0 = x0*kyx*kzy/(Ty*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in H (x, y, z)
def HRA_H(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_H, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

H_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = H_exp[0]
x_exp = H_exp[1]
y_exp = H_exp[2]
z_exp = H_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Hxyz.svg')
fig.savefig(path, dpi=300)

## Architecture I
HRA = 'I'
x0 = 10
Tx = 0.11
Ty = 0.05
Tz = 0.11

kxy = 0.05
kyx = 0.14
kzy = 0.14

y0 = x0*kyx/Ty # 15
kxz = (Tz*kxy*kyx - Tz*Ty*Tx)/(kyx*kzy) # 0.0033...
z0 = x0*kyx*kzy/(Ty*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in H (x, y, z)
def HRA_I(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_I, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

I_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = I_exp[0]
x_exp = I_exp[1]
y_exp = I_exp[2]
z_exp = I_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Ixyz.svg')
fig.savefig(path, dpi=300)

## Architecture J
HRA = 'J'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.01

kxy = 0.01
kyx = 0.08
kzy = 0.02

y0 = x0*kyx/Ty #
kxz = (Tz*kxy*kyx + Tz*Ty*Tx)/(kyx*kzy) #
z0 = x0*kyx*kzy/(Ty*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, -kxy, kyx, kzy, kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, -kxy, kyx, kzy, kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in H (x, y, z)
def HRA_J(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[9]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0]
        du[2] = param[5]*u[2] + param[8]*u[1]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_J, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

J_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = J_exp[0]
x_exp = J_exp[1]
y_exp = J_exp[2]
z_exp = J_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Jxyz.svg')
fig.savefig(path, dpi=300)

## Architecture K
HRA = 'K'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.001
kzx = 0.04

y0 = x0*Tx/kxy #
kyz = (Tz*Ty*Tx - Tz*kxy*kyx)/(kzx*kxy + Tx*kzy) #
z0 = x0*(kzx*kxy + kzy*Tx)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in K (x, y, z)
def HRA_K(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_K, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

K_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = K_exp[0]
x_exp = K_exp[1]
y_exp = K_exp[2]
z_exp = K_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Kxyz.svg')
fig.savefig(path, dpi=300)

## Architecture L
HRA = 'L'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.001
kzx = 0.04

y0 = x0*Tx/kxy #
kyz = (Tz*Ty*Tx - Tz*kxy*kyx)/(kzx*kxy - Tx*kzy) #
z0 = x0*(kzx*kxy - kzy*Tx)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, -kzy, kzx, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, -kzy, kzx, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in L (x, y, z)
def HRA_L(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_L, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

L_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = L_exp[0]
x_exp = L_exp[1]
y_exp = L_exp[2]
z_exp = L_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Lxyz.svg')
fig.savefig(path, dpi=300)

## Architecture M
HRA = 'M'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.02
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*Ty*Tx - Tz*kxy*kyx)/(Tx*kzy - kzx*kxy) #
z0 = x0*(kzy*Tx - kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kzx, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kzx, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in M (x, y, z)
def HRA_M(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_M, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

M_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = M_exp[0]
x_exp = M_exp[1]
y_exp = M_exp[2]
z_exp = M_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Mxyz.svg')
fig.savefig(path, dpi=300)

## Architecture N
HRA = 'N'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.02
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*kxy*kyx - Tz*Ty*Tx)/(Tx*kzy + kzx*kxy) #
z0 = x0*(kzy*Tx + kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, -kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, -kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in N (x, y, z)
def HRA_N(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_N, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

N_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = N_exp[0]
x_exp = N_exp[1]
y_exp = N_exp[2]
z_exp = N_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Nxyz.svg')
fig.savefig(path, dpi=300)

## Architecture O
HRA = 'O'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.02
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*kxy*kyx + Tz*Ty*Tx)/(Tx*kzy + kzx*kxy) #
z0 = x0*(kzy*Tx + kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, kzy, kzx, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, kzy, kzx, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in O (x, y, z)
def HRA_O(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_O, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

O_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = O_exp[0]
x_exp = O_exp[1]
y_exp = O_exp[2]
z_exp = O_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Oxyz.svg')
fig.savefig(path, dpi=300)

## Architecture P
HRA = 'P'
x0 = 10
Tx = 0.05
Ty = 0.02
Tz = 0.04

kxy = 0.02
kyx = 0.03
kzy = 0.01
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*kxy*kyx + Tz*Ty*Tx)/(-Tx*kzy + kzx*kxy) #
z0 = x0*(-kzy*Tx + kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, -kzy, kzx, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, -kzy, kzx, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in P (x, y, z)
def HRA_P(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_P, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

P_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = P_exp[0]
x_exp = P_exp[1]
y_exp = P_exp[2]
z_exp = P_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Pxyz.svg')
fig.savefig(path, dpi=300)

## Architecture Q
HRA = 'Q'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.05

kxy = 0.08
kyx = 0.05
kzy = 0.05
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*kxy*kyx - Tz*Ty*Tx)/(-Tx*kzy + kzx*kxy) #
z0 = x0*(-kzy*Tx + kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, -kzy, kzx, -kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, -kzy, kzx, -kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in Q (x, y, z)
def HRA_Q(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_Q, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

Q_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = Q_exp[0]
x_exp = Q_exp[1]
y_exp = Q_exp[2]
z_exp = Q_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Qxyz.svg')
fig.savefig(path, dpi=300)

## Architecture R
HRA = 'R'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.05

kxy = 0.08
kyx = 0.05
kzy = 0.11
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*kxy*kyx - Tz*Ty*Tx)/(Tx*kzy - kzx*kxy) #
z0 = x0*(kzy*Tx - kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kzx, -kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, -kzx, -kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in R (x, y, z)
def HRA_R(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_R, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

R_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = R_exp[0]
x_exp = R_exp[1]
y_exp = R_exp[2]
z_exp = R_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Rxyz.svg')
fig.savefig(path, dpi=300)

## Architecture S
HRA = 'S'
x0 = 10
Tx = 0.05
Ty = 0.002
Tz = 0.04

kxy = 0.01
kyx = 0.03
kzy = 0.02
kzx = 0.05

y0 = x0*Tx/kxy #
kyz = (Tz*kxy*kyx + Tz*Ty*Tx)/(Tx*kzy - kzx*kxy) #
z0 = x0*(kzy*Tx - kzx*kxy)/(kxy*Tz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, kzy, -kzx, kyz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, kzy, -kzx, kyz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in S (x, y, z)
def HRA_S(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_S, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

S_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = S_exp[0]
x_exp = S_exp[1]
y_exp = S_exp[2]
z_exp = S_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Sxyz.svg')
fig.savefig(path, dpi=300)

## Architecture T
HRA = 'T'
x0 = 10
Tx = 0.14
Ty = 0.14
Tz = 0.14

kxy = 0.14
kyx = 0.05
kyz = 0.05
kzx = 0.05
kxz = 0.05

kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx - Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz + kxy*kxz*kyx) #
y0 = x0*(Tx*Tz - kzx*kxz)/(kxy*Tz + kzy*kxz) #
z0 = x0*(Tx*Ty - kxy*kyx)/(kxy*kyz + Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, kyz, kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, kyz, kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in T (x, y, z)
def HRA_T(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_T, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

T_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = T_exp[0]
x_exp = T_exp[1]
y_exp = T_exp[2]
z_exp = T_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Txyz.svg')
fig.savefig(path, dpi=300)

## Architecture U
HRA = 'U'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.08

kxy = 0.08
kyx = 0.05
kyz = 0.08
kzx = 0.05
kxz = 0.14

kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz - kzy*kxz) #
z0 = x0*(Tx*Ty - kxy*kyx)/(kxy*kyz - Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, kyz, -kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, kyz, -kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in U (x, y, z)
def HRA_U(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_U, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

U_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = U_exp[0]
x_exp = U_exp[1]
y_exp = U_exp[2]
z_exp = U_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Uxyz.svg')
fig.savefig(path, dpi=300)

## Architecture V
HRA = 'V'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.11

kxy = 0.04
kyx = 0.10
kxz = 0.05
kzx = 0.05
kyz = 0.05

kzy =  kxy*(Tx*Ty*Tz - Tz*kxy*kyx + Ty*kxz*kzx + kxy*kyz*kzx)/(-Tx*kxy*kyz - kxy*kxz*kyx)#
y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz - kzy*kxz) #
z0 = x0*(kxy*kyx - Tx*Ty)/(kxy*kyz + Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, -kyz, -kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, kzy, kzx, -kyz, -kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in V (x, y, z)
def HRA_V(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_V, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

V_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = V_exp[0]
x_exp = V_exp[1]
y_exp = V_exp[2]
z_exp = V_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Vxyz.svg')
fig.savefig(path, dpi=300)

## Architecture W
HRA = 'W'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.05

kxy = 0.08
kyx = 0.05
kxz = 0.11
kzx = 0.14
kyz = 0.05

kzy =  -kxy*(Tx*Ty*Tz - Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz + kzy*kxz) #
z0 = x0*(Tx*Ty - kxy*kyx)/(kxy*kyz - Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, -kzy, kzx, kyz, -kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, kyx, -kzy, kzx, kyz, -kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in W (x, y, z)
def HRA_W(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_W, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

W_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = W_exp[0]
x_exp = W_exp[1]
y_exp = W_exp[2]
z_exp = W_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Wxyz.svg')
fig.savefig(path, dpi=300)

## Architecture X
HRA = 'X'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.05

kxy = 0.05
kyx = 0.11
kxz = 0.05
kzx = 0.08
kyz = 0.08

kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx - Ty*kxz*kzx + kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
y0 = x0*(Tx*Tz - kzx*kxz)/(-kxy*Tz + kzy*kxz) #
z0 = x0*(Tx*Ty - kxy*kyx)/(-kxy*kyz + Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, -kxy, -kyx, kzy, kzx, kyz, kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, -kxy, -kyx, kzy, kzx, kyz, kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in W (x, y, z)
def HRA_X(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_X, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

X_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = X_exp[0]
x_exp = X_exp[1]
y_exp = X_exp[2]
z_exp = X_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Xxyz.svg')
fig.savefig(path, dpi=300)

## Architecture Y
HRA = 'Y'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.05

kxy = 0.11
kyx = 0.08
kxz = 0.08
kzx = 0.05
kyz = 0.05

kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
y0 = x0*(Tx*Tz + kzx*kxz)/(-kxy*Tz + kzy*kxz) #
z0 = x0*(Tx*Ty - kxy*kyx)/(-kxy*kyz + Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, -kxy, -kyx, kzy, -kzx, kyz, kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, -kxy, -kyx, kzy, -kzx, kyz, kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in W (x, y, z)
def HRA_Y(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_Y, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

Y_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = Y_exp[0]
x_exp = Y_exp[1]
y_exp = Y_exp[2]
z_exp = Y_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Yxyz.svg')
fig.savefig(path, dpi=300)

## Architecture Z
HRA = 'Z'
x0 = 10
Tx = 0.05
Ty = 0.05
Tz = 0.05

kxy = 0.14
kyx = 0.05
kxz = 0.08
kzx = 0.05
kyz = 0.11

kzy = -kxy*(Tx*Ty*Tz + Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz + kxy*kxz*kyx)#
y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz + kzy*kxz) #
z0 = x0*(Tx*Ty + kxy*kyx)/(kxy*kyz - Ty*kxz) #

param = [x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, -kzy, kzx, kyz, -kxz]
param_name = 'HRA', 'x0, y0, z0, -Tx, -Ty, -Tz, kxy, -kyx, -kzy, kzx, kyz, -kxz'
writer.writerow(param_name)
param_val = [HRA, param]
writer.writerow(param_val)
## Epigenetic perturbation in the amount of each entity in W (x, y, z)
def HRA_Z(t, u):
    du = [0,0,0]
    if u[0] < 0:
        u[0] = 0
    if u[1] < 0:
        u[1] = 0
    if u[2] < 0:
        u[2] = 0
    if t < ss:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif ss < t < (ss + tp) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (ss + tp) < t < 101:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (100 + ss) < t < (100 + ss + tp) :
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif (100 + ss + tp) < t < 201:
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    elif 200 < t < (200 + ss):
        u[0] = x0
        u[1] = y0
        u[2] = z0
    elif (200 + ss) < t < (200 + ss + tp) :
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = 0
    elif (200 + ss + tp) < t < 301 :
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2]
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0]
    return du
A0 = [param[0], param[1], param[2]]
sol = solve_ivp(HRA_Z, t, A0, max_step = 0.01)

time = sol.t
x = sol.y[0, :]
y = sol.y[1, :]
z = sol.y[2, :]

x_defect = x0*d
y_defect = y0*d
z_defect = z0*d
x_d = np.repeat(x_defect, len(time))
y_d = np.repeat(y_defect, len(time))
z_d = np.repeat(z_defect, len(time))

Z_exp = exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz)
time_decay = Z_exp[0]
x_exp = Z_exp[1]
y_exp = Z_exp[2]
z_exp = Z_exp[3]

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time_decay, z_exp, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_18_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_18_Genetic_Perturbations_of_Heritable_Architectures/', file_name + '_Zxyz.svg')
fig.savefig(path, dpi=300)

f.close()
