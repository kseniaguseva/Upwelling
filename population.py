#!/usr/bin/env python

# Copyright (C) 2019 Ksenia Guseva <ksenia@skewed.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from math import *
from numpy import *
from pylab import *
from scipy.integrate import ode
from space import *

#########################################################
#  Summary                                              #
#  -------                                              #
#  -------                                              #
#                                                       #
#  Classes:                                             #
#           1. Param:                                   #
#                    1. init                            #
#                    2. pop_model                       #
#                    3. initialize_population           # 
#                    4. pop_model_init                  #
#                                                       #
#           2. Pop_state:                               #
#                    1. init                            #
#                    2. interpolate                     #
#                    3. copy_and_boundary               #
#                    4. diffusion                       #
#                                                       #
#           3.  State_mean:                             #
#                    1. init                            #
#                    2. space_average                   #
#                    3. null time                       #
#                    4. make array                      #
#                    5. time average                    #
#                                                       #
#########################################################

class Param():
    # The default parameters are taken from Table 1 of  Chakraborty et al (2014)
    def __init__(self, R, t0, r1 = 0.15, r2 = 0.15, e1 = 0.03, e2 = 0.03, a= 0.2,
                 b = 0.2, c = 0.4, llambda = 0.6, mue = 0.035, beta = 0.33,
                 gamma =0.5, d = 0.075, theta1= 1., theta2 = 0.8, s = 0.04,
                 alpha1 =0.25, alpha2 = 0.2, phi = 0.4, k = 0.015, N0 = 2.):
        self.r1 = r1*t0
        self.r2 = r2*t0
        self.e1 = e1/R
        self.e2 = e2/R
        self.a = (a/b)*t0
        self.c = (c/b)*R
        self.llambda = llambda*t0
        self.mue = mue/R
        self.beta = beta
        self.gamma = gamma
        self.d = d*t0
        self.theta1 = theta1
        self.theta2 = theta2        
        self.s = s*t0
        self.alpha1 = alpha1
        self.alpha2 = alpha2        
        self.phi = phi
        self.k_re = k*t0
        self.N0_re = N0/R

    def pop_model(self, t, start, k, grid):
        # system of differential equations for the population dynamics
        # computed for a spatial grid

        start = reshape(start, (4, grid.Nx, grid.Ny))
        N = start[0]
        Np = start[1]
        Tp = start[2]
        Z = start[3]


        dNdt = zeros((grid.Nx, grid.Ny))
        dNpdt = zeros((grid.Nx, grid.Ny))
        dTpdt = zeros((grid.Nx, grid.Ny))
        dZdt = zeros((grid.Nx, grid.Ny))


        f1 = N[1:,1:-1]/(self.e1 + N[1:,1:-1])
        f2 = N[1:,1:-1]/(self.e2 + N[1:,1:-1])

        h1 = self.llambda*(1 - self.phi)*(Np[1:,1:-1]*Np[1:,1:-1])/(self.mue**2 + (Np[1:,1:-1]*Np[1:,1:-1]) + (Tp[1:,1:-1]*Tp[1:,1:-1]))
        h2 = self.llambda*self.phi*(Tp[1:,1:-1]*Tp[1:,1:-1])/(self.mue**2 + (Tp[1:,1:-1]*Tp[1:,1:-1]) + (Np[1:,1:-1]*Np[1:,1:-1]))
        g = self.a/(1 + self.c*Np[1:,1:-1] + self.c*Tp[1:,1:-1])

        dNdt[1:,1:-1] = k[1:,1:-1]*(self.N0_re - N[1:,1:-1]) - g*(f1*Np[1:,1:-1] +f2*Tp[1:,1:-1]) +self.r1*Np[1:,1:-1]+self.r2*Tp[1:,1:-1] + self.beta*(h1+h2)*Z[1:,1:-1] + self.gamma*self.d*Z[1:,1:-1]
        dNpdt[1:,1:-1] = self.theta1*f1*g*Np[1:,1:-1] - self.r1*Np[1:,1:-1] - h1*Z[1:,1:-1] - (self.s)*Np[1:,1:-1]
        dTpdt[1:,1:-1] = self.theta2*f2*g*Tp[1:,1:-1] - self.r2*Tp[1:,1:-1] - h2*Z[1:,1:-1] -(self.s)*Tp[1:,1:-1]
        dZdt[1:,1:-1] = (self.alpha1*h1 + self.alpha2*h2)*Z[1:,1:-1] - self.d*Z[1:,1:-1]


        dNdt = dNdt.flatten()
        dNpdt = dNpdt.flatten()
        dTpdt = dTpdt.flatten()
        dZdt = dZdt.flatten()

        return array([dNdt, dNpdt, dTpdt, dZdt]).flatten()

    # Sets up the influx condition
    def initialize_population(self, R):
        # computes the steady state solution of the system without space
        # this value is used for the influx concentration in the model
        
        pop_init = array([0.01, 0.01, 0.01, 0.01])/R
        solver = ode(self.pop_model_init).set_integrator('dopri5')
        solver.set_initial_value(pop_init, 0).set_f_params()
        solver.integrate(20)
        pop_init = solver.y 

        return pop_init

    def pop_model_init(self, t, start):
        # system of differential equations of the population model
        # used to compute the steady state for a homogeneous system (no spatial dependence)
        # Used as an influx
        
        N = start[0]
        Np = start[1]
        Tp = start[2]
        Z = start[3]

        f1 = N/(self.e1 + N)
        f2 = N/(self.e2 + N)
        h1 = self.llambda*(1 -  self.phi)*(Np**2)/( self.mue**2 + Np**2 + Tp**2)
        h2 =  self.llambda* self.phi*(Tp**2)/( self.mue**2 + Tp**2 + Np**2)
        g =  self.a/(1 +  self.c*Np +  self.c*Tp)

        dNdt = self.k_re*(self.N0_re - N) - g*(f1*Np +f2*Tp) +  self.r1*Np + self.r2*Tp +  self.beta*(h1+h2)*Z +  self.gamma* self.d*Z
        dNpdt =  self.theta1*f1*g*Np - self.r1*Np - h1*Z - (self.s)*Np
        dTpdt =  self.theta2*f2*g*Tp -  self.r2*Tp - h2*Z - ( self.s)*Tp
        dZdt = (self.alpha1*h1 +  self.alpha2*h2)*Z -  self.d*Z

        return [dNdt, dNpdt, dTpdt, dZdt]

    
# Sets up and updates  the species on the spatial grid
# Uses functions from space.py

class Pop_state():
    def __init__(self, Nx, Ny,  Nut_init, Np_init, Tp_init, Z_init, x, y, fr = 0.2):
        self.Nut = init_species(Nx, Ny, fr*Nut_init, x, y)
        self.Np = init_species(Nx, Ny, fr*Np_init, x, y)
        self.Tp = init_species(Nx, Ny, fr*Tp_init, x, y)
        self.Z = init_species(Nx, Ny, fr*Z_init, x, y)
        
    def interpolate(self, grid, xf, yf):
        Nut_t = interpolate(self.Nut, grid, xf, yf)
        Np_t = interpolate(self.Np, grid, xf, yf)
        Tp_t = interpolate(self.Tp, grid, xf, yf)
        Z_t = interpolate(self.Z, grid, xf, yf)
        return Nut_t, Np_t, Tp_t, Z_t
    
    def copy_and_boundary(self, Nut_new, Np_new, Tp_new, Z_new, grid):
        self.Nut = Nut_new
        self.Nut = island(self.Nut, grid.Nx, grid.Ny, grid.x, grid.y)
        self.Nut[:, 0] = self.Nut[:, 1].copy()
        self.Nut[:, -1] = self.Nut[:, -2].copy() 
        self.Np = Np_new
        self.Np[:, 0] = self.Np[:, 1].copy() 
        self.Np[:, -1] = self.Np[:, -2].copy() 
        self.Tp = Tp_new
        self.Tp[:, 0] = self.Tp[:, 1].copy() 
        self.Tp[:, -1] = self.Tp[:, -2].copy() 
        self.Z = Z_new
        self.Z[:, 0] = self.Z[:, 1].copy() 
        self.Z[:, -1] = self.Z[:, -2].copy()
 
    def diffusion(self, grid, xf, yf, D, dt):
        self.Nut = diffusion(self.Nut, grid.Nx, grid.Ny, grid.dx, grid.dy, D, dt, grid.x, grid.y)
        self.Np = diffusion(self.Np, grid.Nx, grid.Ny,  grid.dx, grid.dy, D, dt, grid.x, grid.y)
        self.Tp = diffusion(self.Tp, grid.Nx, grid.Ny,  grid.dx, grid.dy, D, dt, grid.x, grid.y)
        self.Z = diffusion(self.Z, grid.Nx, grid.Ny, grid.dx, grid.dy, D, dt, grid.x, grid.y)
       


# Computes spatial and time averages
# These averages are used in the output
    
class State_mean(object):
    def __init__(self):
        self.Nut_fr_av = 0
        self.Np_fr_av = 0
        self.Tp_fr_av = 0
        self.Z_fr_av = 0

        self.Nut_fr_min = 0
        self.Np_fr_min = 0
        self.Tp_fr_min = 0
        self.Z_fr_min = 0

        self.Nut_fr_max = 0
        self.Np_fr_max = 0
        self.Tp_fr_max = 0
        self.Z_fr_max = 0

        self.N_time = []
        self.Np_time = []
        self.Tp_time = []
        self.Z_time = []

    # Compute the spatial average, min and max for the concentrations
    # (it ignores the island region)
    def space_average(self, Nut, Np, Tp, Z, Nx, Ny, x, y):
        self.N_time.append(space_mean(Nut, Nx, Ny, x, y))
        self.Np_time.append(space_mean(Np, Nx, Ny, x, y))
        self.Tp_time.append(space_mean(Tp, Nx, Ny, x, y))
        self.Z_time.append(space_mean(Z, Nx, Ny, x, y))

    
    def null_time(self):
        self.N_time = []
        self.Np_time = []
        self.Tp_time = []
        self.Z_time = []

    def make_array(self):
        self.N_time =  array(self.N_time)
        self.Np_time =  array(self.Np_time)
        self.Tp_time =  array(self.Tp_time)
        self.Z_time =  array(self.Z_time)

    # Computes the time average of the space average saved previously
    def time_average(self):      
                
        self.Nut_fr_av = mean(self.N_time[:, 0]) 
        self.Np_fr_av = mean(self.Np_time[:, 0])
        self.Tp_fr_av = mean(self.Tp_time[:, 0])
        self.Z_fr_av = mean(self.Z_time[:, 0])
        
        self.Nut_fr_min = mean(self.N_time[:, 1])
        self.Np_fr_min = mean(self.Np_time[:, 1])
        self.Tp_fr_min = mean(self.Tp_time[:, 1])
        self.Z_fr_min = mean(self.Z_time[:, 1])

        self.Nut_fr_max = mean(self.N_time[:, 2])
        self.Np_fr_max = mean(self.Np_time[:, 2])
        self.Tp_fr_max = mean(self.Tp_time[:, 2])
        self.Z_fr_max = mean(self.Z_time[:, 2])


