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
import upw_loops 



# Defining a grid for the advection dynamics

class Grid():
    def __init__(self, Nx, Ny, init_x = -2, init_y = -3, final_x = 8, final_y = 3):
        self.Nx = Nx
        self.Ny = Ny
        self.init_x = init_x
        self.init_y = init_y
        self.final_x = final_x
        self.final_y = final_y
        self.X = linspace(init_x, final_x, Nx)
        self.Y = linspace(init_y, final_y, Ny)
        self.dx = self.X[1]- self.X[0]
        self.dy = self.Y[1]- self.Y[0]
        self.init_E = int(1./self.dx)
        x, y = meshgrid(self.X, self.Y,  indexing = "ij")
        self.x = x
        self.y = y
        self.Grid = concatenate((x.flatten(), y.flatten()), axis = 0)
        


        
def structure_grid(Xout, Nx, Ny, size):
    Xout = array(Xout)

    xf = reshape(Xout[-1, :size], (Nx, Ny))
    yf = reshape(Xout[-1, size:], (Nx, Ny))

    return xf, yf 

def space_mean(Mat, Nx, Ny, x, y):
    av = upw_loops.get_gridmean(Nx, Ny, Mat, x, y)
    return av
    
def island(Mat, Nx, Ny, x, y):
    upw_loops.get_island(Nx, Ny, Mat, x, y)
    return Mat

def check_circle(x, y):
    if (sqrt(x**2 + y**2) < 1.0):
        return 1
    else:
        return 0


def init_species(Nx, Ny, eqv, x, y):
    #used for the Pop_state()._init_
    Mat = eqv*ones((Nx, Ny))
    
    for i in range(0, Nx):
        for j in range(0, Ny):
            if (sqrt(x[i, j]**2 + y[i, j]**2) < 1.0):
                Mat[i,j] = 0
    return Mat
    


def diffusion(Mat, Nx, Ny, dx, dy, D, dt, x, y):
    #used for the Pop_state().diffusion()
    Mat_new = Mat.copy()
    upw_loops.diffusion(Nx, Ny, dx, dy, Mat, Mat_new, D, dt, x, y)

    return Mat_new

def interpolate(Mat, grid, xf, yf):
    #used for the Pop_state().interpolate()
    
    Mat_new = Mat.copy()
    upw_loops.get_inerpolation_c(grid.Nx, grid.Ny, xf, yf, grid.init_x,
                                 grid.init_y, grid.X, grid.Y, grid.dx,
                                 grid.dy, Mat, Mat_new)
    
    return Mat_new
