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


class Param_bursting():
    def __init__(self, mu = 1.815, nu1 = 1., nu2= 1.815, nu3 = 0.44, nu4 = 2.86, nu5 =2.86):
        self.mu = mu
        self.nu1 = nu1
        self.nu2 = nu2
        self.nu3 = nu3
        self.nu4 = nu4
        self.nu5 = nu5



def bursting(x, t, mu, nu1, nu2, nu3 , nu4, nu5):
    dx1 = x[1]
    dx2 = -x[0]**3 - 2*x[0]*x[2] + x[0]*x[4] - mu*x[1]
    dx3 = x[3]
    dx4 = -(x[2]**3) - nu1*(x[0]**2) + x[4]*x[2] - nu2*x[3]
    dx5 = - nu3*x[4] - nu4*(x[0]**2) - nu5*(x[2]**2 - 1)

    return [dx1, dx2, dx3, dx4, dx5]




