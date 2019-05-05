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
import random
from upw_loops import * 



class UpwellingRegion():
    def __init__(self, k_re, N0_re, Dl, k_upt, grid, t_end, dt):
        self.Dl = Dl
        self.k_upt = k_upt
        iumin, iumax, jumin, jumax = up_area(grid)
        self.region = array([iumin, iumax, jumin, jumax])
        self.k_re_range = k_re*ones((grid.Nx, grid.Ny))
        self.k_upwelling =  get_constant_upwelling(k_re, k_upt*k_re, Dl, t_end, dt)

class UpwellingRegion_fromseires():
    def __init__(self, k_re, N0_re, Dl, k_upt, grid, k_int):
        self.Dl = Dl
        self.k_upt = k_upt
        iumin, iumax, jumin, jumax = up_area(grid)
        self.region = array([iumin, iumax, jumin, jumax])
        self.k_re_range = k_re*ones((grid.Nx, grid.Ny))
        self.k_upwelling = zeros(len(k_int))
        for i in range(0, len(k_int)):
            if (k_int[i] == 0):
                self.k_upwelling[i] =  k_re
            else:
                self.k_upwelling[i] =  k_upt*(k_int[i])*k_re
        
class UpwellingRegion_intermittent():
    def __init__(self, k_re, N0_re, Dl, k_upt, grid, k_int):
        self.Dl = Dl
        self.k_upt = k_upt
        iumin, iumax, jumin, jumax = up_area(grid)
        self.region = array([iumin, iumax, jumin, jumax])
        self.k_re_range = k_re*ones((grid.Nx, grid.Ny))
        self.k_upwelling =  k_upt*abs(k_int)*k_re + k_re
    

        

def get_constant_upwelling(kmin, kmax, Dl, tend, dt):
    k_range = kmin*ones((tend))
    n = int(tend*dt)
    init_up = 0
    for i in range(0, n):
        if (i > 0):
            if (i%4 == 0):
                init_up = int(random.random()/dt)
                low = int((i/dt) + init_up)
                up = int(low + (Dl/dt))
                if (low > len(k_range)):
                    break
                if(up > len(k_range)):
                    up = -1
                    Dl = (len(k_range) - low -1)*dt
                k_range[low:up] = kmax*ones(int(Dl/dt))
    return k_range

def get_t0_upwelling(kmin, kmax, Dl, tend, dt, tl):
    k_range = kmin*ones((tend))
    low = int((2/dt) + (tl/dt))
    up = int(low + (Dl/dt))
    k_range[low:up] = kmax*ones(int(Dl/dt))
    return k_range


def up_area(grid):
    iumin = int((-1 - grid.init_x)/grid.dx)
    iumax = int((1 - grid.init_x)/grid.dx)

    jumin = int((2.0 - grid.init_y)/grid.dy)
    jumax = int((2.5 - grid.init_y)/grid.dy)

    return iumin, iumax, jumin, jumax


