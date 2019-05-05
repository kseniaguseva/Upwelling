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
from scipy.integrate import odeint
from random import randint
from scipy.integrate import ode

from flow import *
from population import * 
from set_output import *
from intermit import *
from upwelling import *
from space import *

#### Initializations

from sys2 import *

dt = 0.01 
months = 12.
Dl = 1.                # length of the upwelling pulse in flow periods 
k_upt = 500            # strength of the upwelling pulse  k_up = k_upt*k_d 

grid = Grid(500, 300)  
flow = Flow()
state = State_mean()

upw =  UpwellingRegion(param.k_re, param.N0_re, Dl, k_upt, grid, int(months/dt), dt)

## Compute inflow with the 0.2 (default value) of the steady state
# Initializes the concentration values on the grid

pop_init = param.initialize_population(R)

pop_state = Pop_state(grid.Nx, grid.Ny,  pop_init[0], pop_init[1],
                      pop_init[2], pop_init[3], grid.x, grid.y)

solver = ode(param.pop_model).set_integrator('dopri5')

#### Semi-Lagrangian algorithm

for ti in range(0, int(months/dt)):
    t = dt*ti
    ## Sets the upwelling region
    upw.k_re_range[upw.region[0]:upw.region[1], upw.region[2]:upw.region[3]] = upw.k_upwelling[ti]

    ## Integrate backwards 
    Xout = odeint(flow.advect, grid.Grid, array([t, t-dt]),
                  args = (grid.Nx, grid.Ny, grid.init_E))

    x_back, y_back = structure_grid(Xout, grid.Nx, grid.Ny,
                                    int(len(grid.Grid)/2))

    ### Interpolate to find the concentration value
    Nut_t, Np_t, Tp_t, Z_t = pop_state.interpolate(grid, x_back, y_back)

    ## Evolve population dynamics
    pop = array([Nut_t.flatten(), Np_t.flatten(), Tp_t.flatten(), Z_t.flatten()]).flatten()
    solver.set_initial_value(pop, t-dt).set_f_params(upw.k_re_range, grid)
    solver.integrate(t)
    pop = solver.y 
    pop = reshape(pop, (4, grid.Nx, grid.Ny))

    ## Boundary conditions + update concentration grid
    pop_state.copy_and_boundary(pop[0], pop[1], pop[2], pop[3], grid)
    
    for aux_dt in linspace(0, 1, 10):
        pop_state.diffusion(grid, grid.x, grid.y, flow.D, dt/10)

    ## Prints into a file
    if (ti%25 == 0):
        out_space(pop_state, R, upw.Dl, upw.k_upt, "Data_space", t)
    
    state.space_average(R*pop_state.Nut, R*pop_state.Np, R*pop_state.Tp,
                        R*pop_state.Z, grid.Nx, grid.Ny, grid.x, grid.y)
    out_timeseries(state, upw.Dl, upw.k_upt, upw.k_upwelling/t0, "Data_timeseries")


    
