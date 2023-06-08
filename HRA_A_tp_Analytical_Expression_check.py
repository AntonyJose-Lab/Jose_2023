# -*- coding: utf-8 -*-
"""
This is a program for checking the analytical solution with the numerically computed solution for epigenetic perturbations (loss-of-function) of the simplest HRA 'A'. For all real biological systems, all parameter values are >0. 
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
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures/2023_5_21_RegulatoryNetworks', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables/2023_5_21_RegulatoryNetworks', exist_ok=True)
os.chdir('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/')
curr_dir = os.getcwd()

## Writing out the parameters from the simulation as a .xls file with different sheets.
f = open('Tables/2023_5_21_RegulatoryNetworks/2023_5_21_HRA_A_reduction_for_tp_Comparisons_of_Analytic_Expression_and_Numerical_Solution.csv', 'w', encoding='UTF8', newline='')
writer = csv.writer(f)

## time step for the simulation and duration.
t = (0,300) # total time
ss = 10 # duration of initial steady state.

## Reference exponential decay
def exp_decay(HRA, x0, y0, z0, Tx, Ty, Tz):
    turnover_x = np.zeros(30000)
    turnover_y = np.zeros(30000)
    time_decay = np.zeros(30000)
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
    return result

## Architecture A
HRA = 'A'
x0 = 10
Tx = 0.05
Ty = 0.1

kxy = 0.07
y0 = Tx*x0/kxy # 7.1428...
kyx = Tx*Ty/kxy # 0.071428...

## calculated steady states with tp when x is perturbed to x0*d*p or y is perturbed to y0*d*p.
d = 0.5 # perturbation factor needed to observe a defect (e.g., 2 = 2x upregulation, 0.5 = 0.5 downregulation)
p = 0.8 # multiplier fordescribing the extent of experimental perturbation induced (i.e., p*d gives the extent of perturbation)

t_crit_x = (1/Ty)*math.log((1/(d*p) - 1)/((1/p - 1)*(1 + Ty/Tx)))
t_crit_y = (1/Tx)*math.log((1/(d*p) - 1)/((1/p - 1)*(1 + Tx/Ty)))

xp = x0*d*p
yp = kyx*xp/Ty + ((y0*Ty - xp*kyx)/Ty)*math.exp(-Ty*t_crit_x)
xps_x = (xp*Ty + yp*kxy)/(Tx + Ty)
yps_x = (yp*Tx + xp*kyx)/(Tx + Ty)

param = [x0, y0, -Tx, -Ty, kxy, kyx, t_crit_x, t_crit_y]
param_name = 'x0', 'y0', '-Tx', '-Ty', 'kxy', 'kyx', 't_crit_x', 't_crit_y'
writer.writerow(param_name)
writer.writerow(param)

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
    elif ss < t < (ss + t_crit_x) :
        u[0] = x0*p*d
        du[0] = 0
        du[1] = param[5]*u[0] + param[3]*u[1]
    elif (ss + t_crit_x) < t < 101:
        du[0] = param[2]*u[0] + param[4]*u[1]
        du[1] = param[5]*u[0] + param[3]*u[1]
    elif 100 < t < (100 + ss):
        u[0] = x0
        u[1] = y0
    elif (100 + ss) < t < (100 + ss + t_crit_y) :
        u[1] = y0*p*d
        du[0] = param[2]*u[0] + param[4]*u[1]
        du[1] = 0
    elif (100 + ss + t_crit_y) < t :
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
plt.plot(time, x_d, color = 'brown', linewidth = 2, label = 'x_def', linestyle='dotted')
plt.plot(time, y_d, color = 'green', linewidth = 2, label = 'y_def', linestyle='dotted')
plt.plot(time_decay, x_exp, color = 'brown', linewidth = 1, label = 'x_decay', linestyle='dashed')
plt.plot(time_decay, y_exp, color = 'green', linewidth = 1, label = 'y_decay', linestyle='dashed')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False)
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y)) + 1
plt.axis([0, 300, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# plt.show()
file_name = '2023_5_21_HRA_A_reduction_for_tp_Comparisons_of_Analytic_Expression_and_Numerical_Solution'
path = os.path.join(curr_dir, 'Figures/2023_5_21_RegulatoryNetworks/', file_name + '_Axy.svg')
fig.savefig(path, dpi=300)

f.close()
