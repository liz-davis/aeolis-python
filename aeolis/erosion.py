'''This file is part of AeoLiS.

AeoLiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AeoLiS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AeoLiS.  If not, see <http://www.gnu.org/licenses/>.

AeoLiS  Copyright (C) 2015 Bas Hoonhout

bas.hoonhout@deltares.nl         b.m.hoonhout@tudelft.nl
Deltares                         Delft University of Technology
Unit of Hydraulic Engineering    Faculty of Civil Engineering and Geosciences
Boussinesqweg 1                  Stevinweg 1
2629 HVDelft                     2628CN Delft
The Netherlands                  The Netherlands

'''
from __future__ import absolute_import, division

import logging
from scipy import ndimage, misc
from scipy.stats import norm, mode
import numpy as np
import math
#import matplotlib.pyplot as plt

# package modules
import aeolis.wind

from aeolis.utils import *


# from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)



def run_ph12(s, p, t):
    ''' Calculates bed level change due to dune erosion
    
    Calculates bed level change due to dune erosion according to Palmsten and Holman (2012).
    
    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    t : float
        Model time

    Returns
    -------
    dict
        Spatial grids
        
    '''
        
    # Define variables
    Ho = np.interp(t, p['wave_file'][:, 0], p['wave_file'][:, 1])
    Tp = np.interp(t, p['wave_file'][:, 0], p['wave_file'][:, 2])
    wl = np.interp(t, p['tide_file'][:, 0], p['tide_file'][:, 1])

    zToe = p['dune_toe_elevation']
    beta = p['beach_slope']
    dt = p['dt_opt']
    ny = p['ny']
    Kd = 1.0 # Coefficient to account for hugher runup on dune
    Cs = p['Cs']

    # wave runup calcs
    for iy in range(ny + 1):
        twl = s['R'][iy][0] * Kd + wl


        if twl > zToe:
            # parameter set up
            x = s['x'][iy, :]
            zb = s['zb'][iy, :]
            eta = s['eta'][iy][0]
            R = s['R'][iy][0]
            sigma_s = s['sigma_s'][iy][0]
            dx = np.abs(x[1] - x[0])
            Bt = beta * 0.54  # dune toe slope trajectory
            dVResidual_prev = 0  # add this parameter in to be consistent with other codes  

            # find dune base location
            st0 = np.nanargmin(np.abs(zb - zToe))
            xToe = x[st0]

            # dune toe trajectory
            zbT = np.ones(len(x)) * np.nan
            zbT[st0:] = Bt * (x[st0:] - x[st0]) + zToe

            # correct toe trajectory that exceeds actual bed elevation
            ix = zbT > zb
            zbT[ix] = zb[ix]
            
            # Initial volume calculations
            Vc = np.cumsum(dx * (zb[st0:] - zbT[st0:]))
            Vc = Vc - Vc[0]
            
            # Collision calculations
            p_collision = 1 - norm.cdf(zToe, eta + wl, sigma_s)
            Nc = p_collision * dt / Tp
            
            # Preliminary dV calc using unadjusted Cs to find erosion windo
            dV = 4 * Cs * (np.max(twl - zToe, 0)) ** 2 * Nc
            dVT = dV - dVResidual_prev
            ds = 0 if dVT < 0 else np.nanargmin(np.abs(Vc - dVT))
            st = st0 + ds

            # Compute vegetation scaling ONLY over eroded zone
            bgCE = 0.8 - ((s['bgbiomass'][iy, :] - 42.380) / (2296 - 42.380)) * (0.8 - 0.2) # min-max normalization formula
            valid_bgCE = bgCE[st0:st][zb[st0:st] < twl]
            avg_bgCE = np.mean(valid_bgCE) if np.any(valid_bgCE) else 0.8

            # Adjust Cs
            Cs_adjusted = Cs * avg_bgCE

            # Recompute dV using adjusted Cs
            dV = 4 * Cs_adjusted * (np.max(twl - zToe, 0)) ** 2 * Nc
            dVT = dV - dVResidual_prev
            ds = 0 if dVT < 0 else np.nanargmin(np.abs(Vc - dVT))
            st = st0 + ds
            dVResidual = Vc[ds] - dVT
            dz = -dVResidual / dx # Add residual back to the dune toe so as to have mass conservation
        
        
            # Update Morphology
            zb_new = zb.copy()
            zb_new[st0:st] = zbT[st0:st]
            dz_erosion = zb - zb_new
            
            # Update spatial grid with new bed levels
            s['zb'][iy, :] = zb_new
            s['dzbveg_waves'][iy, :] = dz_erosion

            
            # Vegetation removal from eroded region (> 0.3 m erosion)
            ix_erosion = dz_erosion >= 0.3  # Erosion exceeds 30 cm
            
            if np.any(ix_erosion):  # Only act if erosion condition is met
                # Remove vegetation from spatial grid if erosion condition is met
                i_start = np.where(ix_erosion)[0][0]    # Most landward eroded cell (highest index)
                s['rhoveg'][iy, :i_start+1] = 0.        # Remove vegetation seaward (lower indices)
                s['hveg'][iy, :i_start+1] = 0.
                s['vegetated'][iy, :i_start+1] = False
                s['lateral'][iy, :i_start+1] = False


    return s
