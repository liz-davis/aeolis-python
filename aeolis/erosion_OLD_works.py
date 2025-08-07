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
        
        
    Ho = np.interp(t, p['wave_file'][:, 0], p['wave_file'][:, 1])
    Tp = np.interp(t, p['wave_file'][:, 0], p['wave_file'][:, 2])
    wl = np.interp(t, p['tide_file'][:, 0], p['tide_file'][:, 1])

    zToe = p['dune_toe_elevation']
    beta = p['beach_slope']
    dt = p['dt_opt']

    # wave runup calcs
    Kd = 1.0  # Coefficient to account for higher runup on dune
    ny = p['ny']
    wl = interp_circular(t, p['tide_file'][:, 0], p['tide_file'][:, 1])
    Tp = interp_circular(t, p['wave_file'][:, 0], p['wave_file'][:, 2])

    for iy in range(ny + 1):
        twl = s['R'][iy][0] * Kd + wl
        #print("Debug: Entering the loop")
        #print(f"Debug: Value of twl: {twl}, zToe: {zToe}")
        
        # # calculate scaling coefficient for bg biomass *** CAN MOVE INSIDE IF STATEMENT to reduce computation
        # print("Entering the new bgCE calculation block")
        # bgCE = 0.8 - ((s['bgbiomass'][iy, :] - 42.380) / (2296 - 42.380)) * (0.8 - 0.2) # calculated with maximum-minimum normalization formula
        #bgbiomass_calc = s['bgbiomass'][iy, :]
        #print(f"test calculated bgbiomass scaling coefficient: {bgCE}", flush = True) 
        #print(f"Debug: s['bgbiomass'][iy, :]: {s['bgbiomass'][iy, :]}")
        #print(f"the calculated bgbiomass is: {bgbiomass_calc}", flush = True)

        if twl > zToe:
            #print("Debug: TWL > zToe")
            x = s['x'][iy, :]
            zb = s['zb'][iy, :]
            eta = s['eta'][iy][0]
            R = s['R'][iy][0]
            sigma_s = s['sigma_s'][iy][0]

            # parameter set up
            dx = np.abs(x[1] - x[0])
            Bt = beta * 0.54  # dune toe slope trajectory
            Cs = p['Cs']
            dVResidual_prev = 0  # add this parameter in to be consistent with other codes  

            # find dune base location
            st0 = np.nanargmin(np.abs(zb - zToe))
            xToe = x[st0]
            #print(f"st0 = {st0}; xToe = {xToe}")

            # dune toe trajectory
            zbT = np.ones(len(x)) * np.nan
            zbT[st0:] = Bt * (x[st0:] - x[st0]) + zToe

            # correct toe trajectory that exceeds actual bed elevation
            ix = zbT > zb
            zbT[ix] = zb[ix]
        
            # calculate scaling coefficient for bg biomass
            #print("Entering the new bgCE calculation block")
            
            # Compute bgCE using min-max normalization formula
            bgCE = 0.8 - ((s['bgbiomass'][iy, :] - 42.380) / (2296 - 42.380)) * (0.8 - 0.2) # calculated with maximum-minimum normalization formula
            
            # Extract relevant bgCE values where zb < twl
            valid_bgCE = bgCE[st0:][zb[st0:] < twl]
            
            # Compute average belowground biomass coefficient, defaulting to 0.8 if no valid values exist
            if np.any(valid_bgCE):
                avg_bgCE = np.mean(valid_bgCE)
                #print(f"Debug: bgCE_relevant values where zb < twl: {valid_bgCE}")
                #print(f"Debug: Average bgCE where zb < twl: {avg_bgCE}")
            else:
                avg_bgCE = 0.8
                #print("bg biomass does NOT exist below TWLs")
            
            #Adjusted Cs based on the mean bgCE
            Cs_adjusted = Cs * avg_bgCE

            # initial volume calcs
            Vc = np.cumsum(dx * (zb[st0:] - zbT[st0:]))
            Vc = Vc - Vc[0]

            # collision calcs
            p_collision = 1 - norm.cdf(zToe, eta + wl, sigma_s)
            Nc = p_collision * dt / Tp

            # volume change calcs
            dV = 4 * Cs_adjusted * (np.max(twl - zToe, 0)) ** 2 * Nc ## Cs_adjusted has to be a scalar somehow!!
            dVT = dV - dVResidual_prev

            if dVT < 0:
                ds = 0
            else:
                ds = np.nanargmin(np.abs(Vc - dVT))

            st = st0 + ds
            # x_increment = x[st] - xToe
            dVResidual = Vc[ds] - dVT

            # lets add this residual back to the dune toe so have mass conservation
            dz = -dVResidual / dx
            numcells = np.size(np.arange(st0, st))

            # update morphology
            zb_new = zb.copy()
            zb_new[st0:st] = zbT[st0:st]

            #approach to redistribute residual sediment to the lower portion of the dune. needs to be tested
            #if numcells <= 1:
            #    zb_new[st] = zb_new[st] + dz
            #elif numcells < 3:
            #    zb_new[st - 1] = zb_new[st - 1] + dz
            #else:
            #    zb_new[(st0 + 1):st] = zb_new[(st0 + 1):st] + dz / [numcells - 1]
            
            # Compute erosion before updating the bed elevation
            dz_erosion = zb - zb_new

            
            # Update the spatial grid with the new bed levels
            s['zb'][iy,:] = zb_new
            s['dzbveg_waves'][iy, :] = dz_erosion
            #s['dVT'] = dVT
            
            # NEW erosion condition for vegetation removal
            # # --- Vegetation removal logic based on erosion ---
            # erosion_thresh = 0.3  # meters
            # ix_erosion = dz_erosion >= erosion_thresh  # boolean mask of erosion ≥ 0.3 m

            # if np.any(ix_erosion):  # if erosion exceeds threshold anywhere
            #     i_cutoff = np.where(ix_erosion)[0][-1]  # most landward erosion point (largest index)
    
            #     # Remove vegetation seaward (lower indices) of this point
            #     s['rhoveg'][iy, :i_cutoff+1] = 0.0
            #     s['hveg'][iy, :i_cutoff+1] = 0.0
            #     s['vegetated'][iy, :i_cutoff+1] = False
                
            #     print(f"row {iy}: Removed vegetation from x-indices 0 to {i_cutoff} due to erosion ≥ {erosion_thresh}m")
            
            
            # Check if the erosion exceeds -30 cm (-0.3 meters)
            ix_erosion = dz_erosion >= 0.3  # Erosion exceeds 30 cm
            
            if np.any(ix_erosion):  # Only act if erosion condition is met
                # Debug print before removal
                #print(f"[iy={iy}] rhoveg BEFORE erosion:")
                #print(s['rhoveg'][iy, :])
                
                # Remove vegetation from spatial grid if erosion condition is met
                i_start = np.where(ix_erosion)[0][0]    # Most landward eroded cell (highest index)
                s['rhoveg'][iy, :i_start+1] = 0.        # Remove vegetation seaward (lower indices)
                s['hveg'][iy, :i_start+1] = 0.
                s['vegetated'][iy, :i_start+1] = False
                s['lateral'][iy, :i_start+1] = False
                
                #print(f"EROSION CONDITION MET: Removing vegetation seaward of x={i_start} where erosion <= -0.3 m.")
                
                # Debug print after removal
                #print(f"[iy={iy}] rhoveg AFTER erosion:")
                #print(s['rhoveg'][iy, :])
                
                # s['rhoveg'][iy, ix_erosion] = 0  # Remove vegetation (set density to 0)
                # s['hveg'][iy, ix_erosion] = 0   # Remove vegetation height
                # s['vegetated'][iy, ix_erosion] = False 

    return s
