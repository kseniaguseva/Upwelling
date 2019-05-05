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


def out_tint(state, Dl, tl, k_upwelling, name):
    np.savetxt("%s/N_timeseries_%.4f_%.3f.txt" %(name, Dl, tl), array(state.N_time))
    np.savetxt("%s/N_p_timeseries_%.4f_%.3f.txt" %(name, Dl, tl), array(state.Np_time))
    np.savetxt("%s/T_p_timeseries_%.4f_%.3f.txt" %(name, Dl, tl), array(state.Tp_time))
    np.savetxt("%s/Z_timeseries_%.4f_%.3f.txt" %(name, Dl, tl), array(state.Z_time))
    np.savetxt("%s/k_timeseries_%.4f_%.3f.txt" %(name, Dl, tl), array(k_upwelling))

def out_timeseries(state, Dl, k_upt, k_upwelling, name):
    np.savetxt("%s/N_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.N_time))
    np.savetxt("%s/N_p_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.Np_time))
    np.savetxt("%s/T_p_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.Tp_time))
    np.savetxt("%s/Z_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.Z_time))
    np.savetxt("%s/k_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(k_upwelling))

def out_timeseries_grid(state, Dl, k_upt, k_upwelling, name):
    #np.savetxt("%s/N_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.N_time))                   
    np.savetxt("%s/N_p_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.Np_time))
    np.savetxt("%s/T_p_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.Tp_time))
    np.savetxt("%s/Z_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(state.Z_time))
    #np.savetxt("%s/k_timeseries_%.4f_%d.txt" %(name, Dl, k_upt), array(k_upwelling))       

    
def out_space(pop_state, R, Dl, k_upt, name, t):
      np.savetxt("%s/Np_t%.2f_Dl%.4f_k%d.txt" %(name, t, Dl, k_upt), R*pop_state.Np)
      np.savetxt("%s/Tp_t%.2f_Dl%.4f_k%d.txt" %(name, t, Dl, k_upt), R*pop_state.Tp)
      np.savetxt("%s/Z_t%.2f_Dl%.4f_k%d.txt" %(name, t, Dl, k_upt), R*pop_state.Z)