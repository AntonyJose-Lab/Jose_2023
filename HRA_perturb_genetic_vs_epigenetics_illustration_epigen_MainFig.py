# -*- coding: utf-8 -*-
"""
This is a program for simulating the plots that result after each epigenetic perturbation by begining with a system at steady state, providing a pulse of epigenetic perturbation 
and then following the system over time. For all real biological systems, all parameter values are >0.
The simulation is based on the differental equations for each heritable regulatory architecture. Giving tp a value of 90 and p a value of 0 below simulates a genetic change
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
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration', exist_ok=True)
os.chdir('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/')
curr_dir = os.getcwd()

## time step for the simulation and duration.
t = (0,100) # total time

ss = 10 # duration of initial steady state.
tp = 5 # duration of epigenetic perturbation.
d = 0.5 # perturbation factor needed to observe a defect (e.g., 2 = 2x upregulation, 0.5 = 0.5 downregulation)
p = 0.5 # multiplier fordescribing the extent of experimental perturbation induced (i.e., p*d gives the extent of perturbation)

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

x_defect = x0*d
y_defect = y0*d
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
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Axy.svg')
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
print(HRA,param)
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
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[9]*u[2] 
        du[2] = 0    
    elif (ss + tp) < t < 101 :
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Gxyz.svg')
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
print(HRA,param)
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Pxyz.svg')
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
print(HRA,param)
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
        u[1] = y0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] 
        du[1] = 0
        du[2] = param[5]*u[2] + param[8]*u[1] + param[9]*u[0] 
    elif (ss + tp) < t < 101:
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Rxyz.svg')
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
print(HRA,param)
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Wxyz.svg')
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
print(HRA,param)
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
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2] 
        du[2] = 0    
    elif (ss + tp) < t < 101 :
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Xxyz.svg')
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
print(HRA,param)
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
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2] 
        du[2] = 0    
    elif (ss + tp) < t < 101 :
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Yxyz.svg')
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
print(HRA,param)
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
        u[2] = z0*p*d
        du[0] = param[3]*u[0] + param[6]*u[1] + param[11]*u[2]
        du[1] = param[4]*u[1] + param[7]*u[0] + param[10]*u[2] 
        du[2] = 0    
    elif (ss + tp) < t < 101 :
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

# plotting and saving the resulting figures.
fig, ax = plt.subplots()
plt.plot(time, x, color = 'brown', linewidth = 1, label = 'x')
plt.plot(time, y, color = 'green', linewidth = 1, label = 'y')
plt.plot(time, z, color = 'blue', linewidth = 1, label = 'z')
plt.plot(time, x_d, color = 'brown', linewidth = 1, label = 'x', linestyle='--')
plt.plot(time, y_d, color = 'green', linewidth = 1, label = 'y', linestyle='--')
plt.plot(time, z_d, color = 'blue', linewidth = 1, label = 'z', linestyle='--')
plt.xlabel('time', fontsize = 8)
plt.ylabel('entities', fontsize = 8)
plt.xticks(fontsize = 8)
plt.yticks(fontsize = 8)
plt.grid(False) 
plt.title('heritable regulatory architecture ' + str(HRA))                 # do not show grid
ymax = max(max(x), max(y), max(z)) + 1
plt.axis([0, 100, 0, ymax])     # show axes measures
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.show()
file_name = '2022_8_26_Perturbations_HRA'
path = os.path.join(curr_dir, 'Figures/2022_8_26_Epigenetic_Perturbations_of_Heritable_Architectures_Epi_Illustration/', file_name + '_Zxyz.svg')
fig.savefig(path, dpi=300)