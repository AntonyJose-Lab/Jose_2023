# -*- coding: utf-8 -*-
"""
This is a program for determining the parameter values that permit steady states for each of the 26 heritable regulatory architectures.
For all real biological systems, all parameter values are > 0. Each architecture allows a maximum of 3 constraints based on the steady state differential equations (dx/dt = 0, dy/dt = 0, dz/dt = 0).
This constrains the value of 3 variables in terms of the other free variables.
By arbitrarily setting the value of x0 at 10, and computing the values for y0, z0 and one of the regulatory interactions (e.g., kxy),
we can search for all the other parameters that should allow reaching steady state.
Author: Antony Jose
"""
## Import required libraries.

import numpy as np
import math
import os
import csv

# Directory for saving files, changing into, and setting working directory.
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Figures/2022_8_17_SteadyStates_of_Heritable_Architectures', exist_ok=True)
os.makedirs('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/Tables/2022_8_17_SteadyStates_of_Heritable_Architectures', exist_ok=True)
os.chdir('/Users/antonyjose/Desktop/JoseLab/Modeling/2022_2_9_onwards/Regulatory_networks/analyses/Python/')
curr_dir = os.getcwd()

## Writing out the parameters from the simulation as a .xls file with different sheets.
f = open('Tables/2022_8_17_SteadyStates_of_Heritable_Architectures/2022_8_17_SteadyStates_for_26_HRAs_params.csv', 'w', encoding='UTF8', newline='')
writer = csv.writer(f)
header = ['x0','y0','z0','Tx','Ty','Tz','kxy','kyx','kxz','kzx','kyz','kzy', 'HRA']
writer.writerow(header)

## Chose random numbers for the values of all other parameters and deter,ine cases values that.

## Architecture A
HRA = 'A'
x0 = 10
z0 = 0
Tz = 0
kzx = 0
kxz = 0
kyz = 0
kzy = 0

for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for kxy in np.arange(0.05,0.15,0.03):
            y0 = Tx*x0/kxy
            kyx = Tx*Ty/kxy
            if y0 > 0 and kyx > 0:
                param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                writer.writerow(param_val)

# ## Architecture B
HRA = 'B'
x0 = 10
kxz = 0
kzx = 0
kyz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kzy in np.arange(0.05,0.15,0.03):
                    y0 = Tx*x0/kxy
                    kyx = Tx*Ty/kxy
                    z0 = x0*Tx*kzy/(kxy*Tz)
                    if y0 > 0 and kyx > 0 and z0 > 0:
                        param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                        writer.writerow(param_val)

## Architecture C
HRA = 'C'
x0 = 10
kxz = 0
kyx = 0
kzy = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kzx in np.arange(0.05,0.15,0.03):
                    y0 = Tx*x0/kxy
                    kyz = Tx*Ty*Tz/(kxy*kzx)
                    z0 = x0*kzx/Tz
                    if y0 > 0 and kyz > 0 and z0 > 0:
                        param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                        writer.writerow(param_val)

# ## Architecture D
HRA = 'D'
x0 = 10
kxz = 0
kyz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kzx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        y0 = Tx*x0/kxy
                        kyx = Tx*Ty/kxy
                        z0 = (x0/Tz)*(kzx + kzy*kyx/Ty)
                        if y0 > 0 and kyx > 0 and z0 > 0:
                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                            writer.writerow(param_val)

# ## Architecture E
HRA = 'E'
x0 = 10
kxz = 0
kyz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kzx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        if kzx != kzy*kyx/Ty:
                            y0 = Tx*x0/kxy
                            kyx = Tx*Ty/kxy
                            z0 = (x0/Tz)*(-kzx + kzy*kyx/Ty)
                            if y0 > 0 and kyx > 0 and z0 > 0:
                                param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                writer.writerow(param_val)

# ## Architecture F
HRA = 'F'
x0 = 10
kzx = 0
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        y0 = Tx*x0/kxy
                        kyz = (Tx*Ty*Tz - Tz*kyx*kxy)/(Tx*kzy)
                        z0 = x0*Tx*kzy/(kxy*Tz)
                        if y0 > 0 and kyz > 0 and z0 > 0:
                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                            writer.writerow(param_val)

## Architecture G
HRA = 'G'
x0 = 10
kzx = 0
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        y0 = Tx*x0/kxy #
                        kyz = (Tz*kyx*kxy - Tx*Ty*Tz)/(Tx*kzy) #
                        z0 = x0*Tx*kzy/(kxy*Tz) #
                        if y0 > 0 and kyz > 0 and z0 > 0:
                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                            writer.writerow(param_val)

## Architecture H
HRA = 'H'
x0 = 10
kzx = 0
kyz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        y0 = x0*kyx/Ty #
                        kxz = (Tz*Ty*Tx - Tz*kxy*kyx)/(kyx*kzy) #
                        z0 = x0*kyx*kzy/(Ty*Tz) #
                        if y0 > 0 and kxz > 0 and z0 > 0:
                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                            writer.writerow(param_val)

## Architecture I
HRA = 'I'
x0 = 10
kzx = 0
kyz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        y0 = x0*kyx/Ty # 15
                        kxz = (Tz*kxy*kyx - Tz*Ty*Tx)/(kyx*kzy) # 0.0033...
                        z0 = x0*kyx*kzy/(Ty*Tz) #
                        if y0 > 0 and kxz > 0 and z0 > 0:
                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                            writer.writerow(param_val)

## Architecture J
HRA = 'J'
x0 = 10
kzx = 0
kyz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        y0 = x0*kyx/Ty #
                        kxz = (Tz*kxy*kyx + Tz*Ty*Tx)/(kyx*kzy) #
                        z0 = x0*kyx*kzy/(Ty*Tz) #
                        if y0 > 0 and kxz > 0 and z0 > 0:
                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                            writer.writerow(param_val)

## Architecture K
HRA = 'K'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            y0 = x0*Tx/kxy #
                            kyz = (Tz*Ty*Tx - Tz*kxy*kyx)/(kzx*kxy + Tx*kzy) #
                            z0 = x0*(kzx*kxy + kzy*Tx)/(kxy*Tz) #
                            if y0 > 0 and kyz > 0 and z0 > 0:
                                param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                writer.writerow(param_val)

## Architecture L
HRA = 'L'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            if kzx*kxy != Tx*kzy:
                                y0 = x0*Tx/kxy #
                                kyz = (Tz*Ty*Tx - Tz*kxy*kyx)/(kzx*kxy - Tx*kzy) #
                                z0 = x0*(kzx*kxy - kzy*Tx)/(kxy*Tz) #
                                if y0 > 0 and kyz > 0 and z0 > 0:
                                    param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                    writer.writerow(param_val)

## Architecture M
HRA = 'M'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            if Tx*kzy != kzx*kxy:
                                y0 = x0*Tx/kxy #
                                kyz = (Tz*Ty*Tx - Tz*kxy*kyx)/(Tx*kzy - kzx*kxy) #
                                z0 = x0*(kzy*Tx - kzx*kxy)/(kxy*Tz) #
                                if y0 > 0 and kyz > 0 and z0 > 0:
                                    param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                    writer.writerow(param_val)

## Architecture N
HRA = 'N'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            y0 = x0*Tx/kxy #
                            kyz = (Tz*kxy*kyx - Tz*Ty*Tx)/(Tx*kzy + kzx*kxy) #
                            z0 = x0*(kzy*Tx + kzx*kxy)/(kxy*Tz) #
                            if y0 > 0 and kyz > 0 and z0 > 0:
                                param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                writer.writerow(param_val)

## Architecture O
HRA = 'O'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            y0 = x0*Tx/kxy #
                            kyz = (Tz*kxy*kyx + Tz*Ty*Tx)/(Tx*kzy + kzx*kxy) #
                            z0 = x0*(kzy*Tx + kzx*kxy)/(kxy*Tz) #
                            if y0 > 0 and kyz > 0 and z0 > 0:
                                param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                writer.writerow(param_val)

## Architecture P
HRA = 'P'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            if Tx*kzy != kzx*kxy:
                                y0 = x0*Tx/kxy #
                                kyz = (Tz*kxy*kyx + Tz*Ty*Tx)/(-Tx*kzy + kzx*kxy) #
                                z0 = x0*(-kzy*Tx + kzx*kxy)/(kxy*Tz) #
                                if y0 > 0 and kyz > 0 and z0 > 0:
                                    param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                    writer.writerow(param_val)

## Architecture Q
HRA = 'Q'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            if kzx*kxy != Tx*kzy:
                                y0 = x0*Tx/kxy #
                                kyz = (Tz*kxy*kyx - Tz*Ty*Tx)/(-Tx*kzy + kzx*kxy) #
                                z0 = x0*(-kzy*Tx + kzx*kxy)/(kxy*Tz) #
                                if y0 > 0 and kyz > 0 and z0 > 0:
                                    param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                    writer.writerow(param_val)

## Architecture R
HRA = 'R'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            if Tx*kzy != kzx*kxy:
                                y0 = x0*Tx/kxy #
                                kyz = (Tz*kxy*kyx - Tz*Ty*Tx)/(Tx*kzy - kzx*kxy) #
                                z0 = x0*(kzy*Tx - kzx*kxy)/(kxy*Tz) #
                                if y0 > 0 and kyz > 0 and z0 > 0:
                                    param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                    writer.writerow(param_val)

## Architecture S
HRA = 'S'
x0 = 10
kxz = 0
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kzy in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            if Tx*kzy != kzx*kxy:
                                y0 = x0*Tx/kxy #
                                kyz = (Tz*kxy*kyx + Tz*Ty*Tx)/(Tx*kzy - kzx*kxy) #
                                z0 = x0*(kzy*Tx - kzx*kxy)/(kxy*Tz) #
                                if y0 > 0 and kyz > 0 and z0 > 0:
                                    param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                    writer.writerow(param_val)

## Architecture T
HRA = 'T'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx - Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz + kxy*kxz*kyx)#
                                if kzy > 0:
                                    y0 = x0*(Tx*Tz - kzx*kxz)/(kxy*Tz + kzy*kxz) #
                                    z0 = x0*(Tx*Ty - kxy*kyx)/(kxy*kyz + Ty*kxz) #
                                    if y0 > 0 and kzy > 0 and z0 > 0:
                                        param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                        writer.writerow(param_val)

## Architecture U
HRA = 'U'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                if Tx*kxy*kyz != kxy*kxz*kyx:
                                    kzy =  kxy*(Tx*Ty*Tz - Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
                                    if kzy > 0 and kxy*Tz != kzy*kxz and kxy*kyz != Ty*kxz:
                                        y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz - kzy*kxz) #
                                        z0 = x0*(Tx*Ty - kxy*kyx)/(kxy*kyz - Ty*kxz) #
                                        if y0 > 0 and z0 > 0:
                                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                            writer.writerow(param_val)

## Architecture V
HRA = 'V'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                kzy =  kxy*(Tx*Ty*Tz - Tz*kxy*kyx + Ty*kxz*kzx + kxy*kyz*kzx)/(-Tx*kxy*kyz - kxy*kxz*kyx)#
                                if kzy > 0 and kxy*Tz != kzy*kxz:
                                    y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz - kzy*kxz) #
                                    z0 = x0*(kxy*kyx - Tx*Ty)/(kxy*kyz + Ty*kxz) #
                                    if y0 > 0 and z0 > 0:
                                        param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                        writer.writerow(param_val)

## Architecture W
HRA = 'W'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                if Tx*kxy*kyz != kxy*kxz*kyx:
                                    kzy =  -kxy*(Tx*Ty*Tz - Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
                                    if kzy > 0  and kxy*kyz != Ty*kxz:
                                       y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz + kzy*kxz) #
                                       z0 = x0*(Tx*Ty - kxy*kyx)/(kxy*kyz - Ty*kxz) #
                                       if y0 > 0 and z0 > 0:
                                           param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                           writer.writerow(param_val)

## Architecture X
HRA = 'X'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                if Tx*kxy*kyz != kxy*kxz*kyx:
                                    kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx - Ty*kxz*kzx + kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
                                    if kzy > 0 and kxy*Tz != kzy*kxz and kxy*kyz != Ty*kxz:
                                        y0 = x0*(Tx*Tz - kzx*kxz)/(-kxy*Tz + kzy*kxz) #
                                        z0 = x0*(Tx*Ty - kxy*kyx)/(-kxy*kyz + Ty*kxz) #
                                        if y0 > 0 and kzy > 0 and z0 > 0:
                                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                            writer.writerow(param_val)

## Architecture Y
HRA = 'Y'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                if Tx*kxy*kyz != kxy*kxz*kyx:
                                    kzy =  kxy*(Tx*Ty*Tz -Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz - kxy*kxz*kyx)#
                                    if kzy > 0 and kxy*Tz != kzy*kxz and kxy*kyz != Ty*kxz:
                                        y0 = x0*(Tx*Tz + kzx*kxz)/(-kxy*Tz + kzy*kxz) #
                                        z0 = x0*(Tx*Ty - kxy*kyx)/(-kxy*kyz + Ty*kxz) #
                                        if y0 > 0 and kzy > 0 and z0 > 0:
                                            param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                            writer.writerow(param_val)

## Architecture Z
HRA = 'Z'
for Tx in np.arange(0.05,0.15,0.03):
    for Ty in  np.arange(0.05,0.15,0.03):
        for Tz in np.arange(0.05,0.15,0.03):
            for  kxy in np.arange(0.05,0.15,0.03):
                for  kyx in np.arange(0.05,0.15,0.03):
                    for  kyz in np.arange(0.05,0.15,0.03):
                        for kzx in np.arange(0.05,0.15,0.03):
                            for kxz in np.arange(0.05,0.15,0.03):
                                kzy =  -kxy*(Tx*Ty*Tz + Tz*kxy*kyx + Ty*kxz*kzx - kxy*kyz*kzx)/(Tx*kxy*kyz + kxy*kxz*kyx)#
                                if kzy > 0  and kxy*kyz != Ty*kxz:
                                    y0 = x0*(Tx*Tz + kzx*kxz)/(kxy*Tz + kzy*kxz) #
                                    z0 = x0*(Tx*Ty + kxy*kyx)/(kxy*kyz - Ty*kxz) #
                                    if y0 > 0 and kzy > 0 and z0 > 0:
                                        param_val = [x0,y0,z0,Tx,Ty,Tz,kxy,kyx,kxz,kzx,kyz,kzy,HRA]
                                        writer.writerow(param_val)

f.close()
