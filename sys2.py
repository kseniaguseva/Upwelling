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


from population import *


R = 0.005088 #g C m^-3
t0 = 32. # d

param = Param(R, t0, e1 = 0.02, e2 = 0.1, theta1 = .4, theta2 = .4, phi =
              0.05, r1 = 0.1, r2 = 0.1, alpha1 = 0.5, alpha2 =0.2, a = 0.2, d
              = 0.065, b = 0.1, mue = 0.02, s = 0.08, llambda =1.3,
              k = 0.0045, N0 = 2.)
