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


# R = m0/V0
# m0 = 10^{12} mmolN = 7,9 x 10^{10} gC
#
# V0 = r^3 = (25000)^{3} = 1,5625 x 10^{13} m^3
# R = 0.005088

from population import *

              # Units
R = 0.005088  # g C m^-3
t0 = 32.      # d

param = Param(R, t0, e1 = 0.02, e2 = 0.1, theta1 = .4, theta2 = .8, phi = 0.5,
              r1 = 0.1, r2 = 0.05, alpha1 = 0.25, alpha2 =0.2, a = 0.2, d =
              0.065, b = 0.1, mue = 0.02, s = 0.08, llambda =0.65,
              k = 0.015, N0 = 2.)
