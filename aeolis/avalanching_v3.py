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
import numpy as np

# package modules
from aeolis.utils import *

# initialize logger
logger = logging.getLogger(__name__)

def angele_of_repose(s,p):
    '''Determine the dynamic and static angle of repose.
    
    Both the critical dynamic and static angle of repose are spatial varying
    and depend on surface moisture content and roots of present vegetation
    and .... 
        
    Parameters
    ----------
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
        
    Returns
    -------
    dict
        Spatial grids
        
    '''
        
    # comment Lisa: dependence on moisture content is not yet implemented 
    # Can we do something with theta dependent on vegetation cover (larger rhoveg = larger theta?)    
        
    theta_stat = p['theta_stat']
    theta_dyn  = p['theta_dyn']
    
    s['theta_stat'] = theta_stat
    s['theta_dyn'] = theta_dyn
        
    return s

def avalanche(s, p):
    '''Avalanching occurs if bed slopes exceed critical slopes.'''

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage import binary_dilation

    if not p['process_avalanche']:
        return s

    nx = p['nx'] + 1
    ny = p['ny'] + 1

    tan_stat = np.tan(np.deg2rad(s['theta_stat']))
    tan_dyn  = np.tan(np.deg2rad(s['theta_dyn']))
    E = 1.0  # ⬅️ Increased from 0.2 for stronger sediment redistribution
    max_iter_ava = p['max_iter_ava']

    max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s['zne'])
    s['gradh'] = grad_h.copy()

    unstable_ix = grad_h > tan_stat
    n_unstable = np.sum(unstable_ix)

    if n_unstable == 0:
        return s

    print(f"Avalanching triggered in {n_unstable} grid cell(s).")

    # === Plot setup ===
    DEBUG_PLOT = True
    if DEBUG_PLOT and not hasattr(avalanche, "fig"):
        plt.ion()
        avalanche.fig, avalanche.ax = plt.subplots(figsize=(6, 5))
        avalanche.cax = avalanche.ax.pcolormesh(s['zb'], shading='auto', cmap='terrain')
        avalanche.fig.colorbar(avalanche.cax, ax=avalanche.ax, label='Elevation (zb)')
        avalanche.ax.set_title('Avalanching — Elevation (zb)')
        avalanche.ax.set_xlabel('Cross-shore index')
        avalanche.ax.set_ylabel('Alongshore index')
        plt.tight_layout()

    for i in range(max_iter_ava):

        max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s['zne'])

        if max_grad_h < tan_dyn:
            print(f"  Iteration {i+1}: max slope below dynamic threshold — stopping.")
            break

        unstable_ix = grad_h > tan_dyn
        unstable_ix = binary_dilation(unstable_ix, iterations=1)  # Expand to neighbors

        slope_diff = np.zeros((ny, nx))
        slope_diff[unstable_ix] = np.tanh(grad_h[unstable_ix]) - np.tanh(0.9 * tan_dyn)

        flux_down = np.zeros((ny, nx, 4))
        flux_down[:, :, 0][unstable_ix] = slope_diff[unstable_ix] * grad_h_down[:, :, 0][unstable_ix] / grad_h[unstable_ix]
        flux_down[:, :, 2][unstable_ix] = slope_diff[unstable_ix] * grad_h_down[:, :, 2][unstable_ix] / grad_h[unstable_ix]

        q_out = 0.5 * (np.abs(flux_down[:, :, 0]) + np.abs(flux_down[:, :, 2]))
        q_in = np.zeros_like(q_out)

        # Cross-shore (x-direction) redistribution only
        q_in[:, 1:-1] = 0.5 * (
            np.maximum(flux_down[:, :-2, 0], 0.) - np.minimum(flux_down[:, 2:, 0], 0.) +
            np.maximum(flux_down[:, 2:, 2], 0.) - np.minimum(flux_down[:, :-2, 2], 0.)
        )

        delta_zb = E * (q_in - q_out)
        s['zb'] += delta_zb

        print(f"  Iteration {i+1}: Δzb max = {np.max(delta_zb):.6f}, adjusted cells = {np.sum(delta_zb != 0)}")

        # === Animate
        if DEBUG_PLOT:
            if np.isnan(s['zb']).all():
                print("Warning: s['zb'] contains only NaNs — skipping plot update.")
                continue

            for coll in avalanche.ax.collections:
                coll.remove()
    
            # Base plot: elevation
            zb_plot = np.ma.masked_invalid(s['zb'])
            avalanche.cax = avalanche.ax.pcolormesh(zb_plot, shading='auto', cmap='terrain')
        
            # Overlay unstable cells in red
            unstable_overlay = np.zeros_like(s['zb'], dtype=float)
            unstable_overlay[unstable_ix] = 1.0
            avalanche.ax.pcolormesh(
                np.ma.masked_where(unstable_overlay == 0, unstable_overlay),
                shading='auto', cmap='Reds', alpha=0.4
            )
        
            avalanche.ax.set_title(f'Avalanching — Elevation + Unstable Cells (iter {i+1})')
            avalanche.fig.canvas.draw()
            avalanche.fig.canvas.flush_events()
            plt.pause(0.1)

    # === Update water level
    s['zs'] = s['SWL'].copy()
    ix = s['zb'] > s['zs']
    s['zs'][ix] = s['zb'][ix]

    return s

# def avalanche(s, p):
#     '''Avalanching occurs if bed slopes exceed critical slopes.
    
#     Simulates the process of avalanching that is triggered by the exceedence
#     of a critical static slope ``theta_stat`` by the bed slope. The iteration
#     stops if the bed slope does not exceed the dynamic critical slope
#     ``theta_dyn``.
    
#     Parameters
#     ----------
#     s : dict
#         Spatial grids
#     p : dict
#         Model configuration parameters
        
#     Returns
#     -------
#     dict
#         Spatial grids

#     '''

#     if p['process_avalanche']:

#         nx = p['nx']+1
#         ny = p['ny']+1

#         # parameters

#         tan_stat = np.tan(np.deg2rad(s['theta_stat']))
#         tan_dyn = np.tan(np.deg2rad(s['theta_dyn']))
        

#         E = 0.2

#         grad_h_down = np.zeros((ny,nx,4))
#         flux_down = np.zeros((ny,nx,4))
#         slope_diff = np.zeros((ny,nx))
#         grad_h = np.zeros((ny,nx))

#         max_iter_ava = p['max_iter_ava']
        
#         max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s['zne'])
        
#         s['gradh'] = grad_h.copy()
        
        
#         import matplotlib.pyplot as plt


#         #initiate_avalanche = (max_grad_h > tan_stat) 
        
#         unstable_cells = grad_h > tan_stat
#         n_unstable = np.sum(unstable_cells)

#         if n_unstable > 0:
#             print(f"Avalanching triggered in {n_unstable} grid cells.")
            
#             for i in range(max_iter_ava):

#                 max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s['zne'])
            
#                 if max_grad_h < tan_dyn:
#                     print(f"  Avalanche iteration {i+1}: max slope below dynamic threshold → stopping.")
#                     break
            
#                 grad_h_nonerod = (s['zb'] - s['zne']) / s['ds']
#                 slope_diff = np.zeros((ny, nx))
            
#                 # Flag unstable cells only — no need to check grad_h_nonerod
#                 unstable_ix = grad_h > tan_dyn
            
#                 # Use full unstable_ix mask for all flux computations
#                 slope_diff[unstable_ix] = np.tanh(grad_h[unstable_ix]) - np.tanh(0.9 * tan_dyn)
            
#                 # Handle directional redistribution
#                 flux_down = np.zeros((ny, nx, 4))
            
#                 for d in range(4):
#                     flux_down[:, :, d][unstable_ix] = (
#                         slope_diff[unstable_ix] * grad_h_down[:, :, d][unstable_ix] / grad_h[unstable_ix]
#                     )
            
#                 # Compute q_out
#                 q_out = 0.5 * np.sum(np.abs(flux_down), axis=2)
            
#                 # Compute q_in from surrounding cells
#                 q_in = np.zeros((ny, nx))
            
#                 # Direction 0: +x (left to right)
#                 q_in[:, 1:-1] += 0.5 * (
#                     np.maximum(flux_down[:, :-2, 0], 0.) - np.minimum(flux_down[:, 2:, 0], 0.)
#                 )
#                 # Direction 1: +y (bottom to top)
#                 q_in[1:-1, :] += 0.5 * (
#                     np.maximum(flux_down[:-2, :, 1], 0.) - np.minimum(flux_down[2:, :, 1], 0.)
#                 )
#                 # Direction 2: -x (right to left)
#                 q_in[:, 1:-1] += 0.5 * (
#                     np.maximum(flux_down[:, 2:, 2], 0.) - np.minimum(flux_down[:, :-2, 2], 0.)
#                 )
#                 # Direction 3: -y (top to bottom)
#                 q_in[1:-1, :] += 0.5 * (
#                     np.maximum(flux_down[2:, :, 3], 0.) - np.minimum(flux_down[:-2, :, 3], 0.)
#                 )
            
#                 delta_zb = E * (q_in - q_out)
#                 s['zb'] += delta_zb
            
#                 print(f"  Avalanche iteration {i+1}: Δzb max = {np.max(delta_zb):.6f}, cells adjusted = {np.sum(delta_zb != 0)}")

#                 import matplotlib.pyplot as plt

#                 DEBUG_PLOT = True  # Toggle animation

#                 if DEBUG_PLOT:
#                     if not hasattr(avalanche, "fig"):
#                         plt.ion()
#                         avalanche.fig, avalanche.ax = plt.subplots(figsize=(6, 5))
#                         avalanche.cax = avalanche.ax.pcolormesh(s['zb'], shading='auto', cmap='terrain')
#                         avalanche.fig.colorbar(avalanche.cax, ax=avalanche.ax, label='Elevation (zb)')
#                         avalanche.ax.set_title('Avalanching — Elevation (zb)')
#                         avalanche.ax.set_xlabel('Cross-shore index')
#                         avalanche.ax.set_ylabel('Alongshore index')
#                         plt.tight_layout()
#                     else:
#                         for coll in avalanche.ax.collections:
#                             coll.remove()  # Clear previous mesh
#                         avalanche.cax = avalanche.ax.pcolormesh(s['zb'], shading='auto', cmap='terrain')
#                         avalanche.fig.canvas.draw()
#                         avalanche.fig.canvas.flush_events()
#                         plt.pause(0.1)
            
#         #if initiate_avalanche:
#             #for i in range(0,max_iter_ava):

#               #  grad_h_down *= 0.
#                # flux_down *= 0.
#                 #slope_diff *= 0.
#                 #grad_h *= 0.

#                 #max_grad_h, grad_h, grad_h_down = calc_gradients(s['zb'], nx, ny, s['ds'], s['dn'], s[ 'zne'])

#                 #if max_grad_h < tan_dyn:
#                  #   break

#                 ## Calculation of flux

#                 #grad_h_nonerod = (s['zb'] - s['zne']) / s['ds'] # HAS TO BE ADJUSTED!    
# 				
#                # ix = np.logical_and(grad_h > tan_dyn, grad_h_nonerod > 0)
#                 #slope_diff[ix] = np.tanh(grad_h[ix]) - np.tanh(0.9*tan_dyn)    
                
#                # ix = grad_h_nonerod < grad_h - tan_dyn 
#                 # slope_diff[ix] = np.tanh(grad_h_nonerod[ix])				                    

#                 # ix = grad_h != 0
                
                
#                 # if ny == 1:
#                 #     #1D interpretation
#                 #     flux_down[:,:,0][ix] = slope_diff[ix] * grad_h_down[:,:,0][ix] / grad_h[ix]
#                 #     flux_down[:,:,2][ix] = slope_diff[ix] * grad_h_down[:,:,2][ix] / grad_h[ix]
                    
#                 #     # Calculation of change in bed level
#                 #     q_in = np.zeros((ny,nx))
                    
#                 #     q_out = 0.5*np.abs(flux_down[:,:,0]) + 0.5*np.abs(flux_down[:,:,2])
                    
#                 #     q_in[0,1:-1] =   0.5*(np.maximum(flux_down[0,:-2,0],0.) \
#                 #                         - np.minimum(flux_down[0,2:,0],0.) \
#                 #                         + np.maximum(flux_down[0,2:,2],0.) \
#                 #                         - np.minimum(flux_down[0,:-2,2],0.))
#                 # else:
#                 #     # 2D interpretation
#                 #     flux_down[:,:,0][ix] = slope_diff[ix] * grad_h_down[:,:,0][ix] / grad_h[ix]
#                 #     flux_down[:,:,1][ix] = slope_diff[ix] * grad_h_down[:,:,1][ix] / grad_h[ix]
#                 #     flux_down[:,:,2][ix] = slope_diff[ix] * grad_h_down[:,:,2][ix] / grad_h[ix]
#                 #     flux_down[:,:,3][ix] = slope_diff[ix] * grad_h_down[:,:,3][ix] / grad_h[ix]

#                 #     # Calculation of change in bed level
#                 #     q_in = np.zeros((ny,nx))

#                 #     q_out = 0.5*np.abs(flux_down[:,:,0]) + 0.5* np.abs(flux_down[:,:,1]) + 0.5*np.abs(flux_down[:,:,2]) + 0.5* np.abs(flux_down[:,:,3])

#                 #     q_in[1:-1,1:-1] =   0.5*(np.maximum(flux_down[1:-1,:-2,0],0.) \
#                 #                         - np.minimum(flux_down[1:-1,2:,0],0.) \
#                 #                         + np.maximum(flux_down[:-2,1:-1,1],0.) \
#                 #                         - np.minimum(flux_down[2:,1:-1,1],0.) \

#                 #                         + np.maximum(flux_down[1:-1,2:,2],0.) \
#                 #                         - np.minimum(flux_down[1:-1,:-2,2],0.) \
#                 #                         + np.maximum(flux_down[2:,1:-1,3],0.) \
#                 #                         - np.minimum(flux_down[:-2,1:-1,3],0.))

#                 # s['zb'] += E * (q_in - q_out)
#                 # delta_zb = E * (q_in - q_out)
#                 # print(f"Max Δzb: {np.max(delta_zb):.6f}, Min Δzb: {np.min(delta_zb):.6f}")
                

#         # Ensure water level is up-to-date with bed level
#         s['zs'] = s['SWL'].copy()
#         ix = (s['zb'] > s['zs'])
#         s['zs'][ix] = s['zb'][ix]
        
#     return s	    


def calc_gradients(zb, nx, ny, ds, dn, zne):
    '''Calculates the downslope gradients in the bed that are needed for
    avalanching module
 
    Parameters
    ----------
        
        
    Returns
    -------
    np.ndarray
        Downslope gradients in 4 different directions (nx*ny, 4)
    '''
    
    grad_h_down = np.zeros((ny,nx,4))

    # Calculation of slope (positive x-direction)
    grad_h_down[:,1:-1,0] = zb[:,1:-1] - zb[:,2:] 
    ix = zb[:,2:] > zb[:,:-2]
    grad_h_down[:,1:-1,0][ix] = - (zb[:,1:-1][ix] - zb[:,:-2][ix])    
    ix = np.logical_and(zb[:,2:]>zb[:,1:-1], zb[:,:-2]>zb[:,1:-1])
    grad_h_down[:,1:-1,0][ix] = 0.

    # Calculation of slope (positive y-direction)
    grad_h_down[1:-1,:,1] = zb[1:-1,:] - zb[2:,:]    
    ix = zb[2:,:] > zb[:-2,:]
    grad_h_down[1:-1,:,1][ix] = - (zb[1:-1,:][ix] - zb[:-2,:][ix])    
    ix = np.logical_and(zb[2:,:]>zb[1:-1,:], zb[:-2,:]>zb[1:-1,:])
    grad_h_down[1:-1,:,1][ix] = 0.

    # Calculation of slope (negative x-direction)
    grad_h_down[:,1:-1,2] = zb[:,1:-1] - zb[:,:-2]    
    ix = zb[:,:-2] > zb[:,2:]
    grad_h_down[:,1:-1,2][ix] = - (zb[:,1:-1][ix] - zb[:,2:][ix])    
    ix = np.logical_and(zb[:,:-2]>zb[:,1:-1], zb[:,2:]>zb[:,1:-1])
    grad_h_down[:,1:-1,2][ix] = 0.

    # Calculation of slope (negative y-direction)
    grad_h_down[1:-1,:,3] = zb[1:-1,:] - zb[:-2,:]    
    ix = zb[:-2,:] > zb[2:,:]
    grad_h_down[1:-1,:,3][ix] = - (zb[1:-1,:][ix] - zb[2:,:][ix])    
    ix = np.logical_and(zb[:-2,:]>zb[1:-1,:], zb[2:,:]>zb[1:-1,:])
    grad_h_down[1:-1,:,3][ix] = 0.
  
    if ny == 1:
        #1D interpretation
        grad_h_down[:,0,:]  = 0
        grad_h_down[:,-1,:] = 0

    else:
        # 2D interpretation
        grad_h_down[:,0,:]  = 0
        grad_h_down[:,-1,:] = 0
        grad_h_down[0,:,:]  = 0
        grad_h_down[-1,:,:] = 0
        
    grad_h_down[:,:,0] /= ds
    grad_h_down[:,:,1] /= dn
    grad_h_down[:,:,2] /= ds
    grad_h_down[:,:,3] /= dn

    #grad_h2 = 0.5*grad_h_down[:,:,0]**2 + 0.5*grad_h_down[:,:,1]**2 + 0.5*grad_h_down[:,:,2]**2 + 0.5*grad_h_down[:,:,3]**2
    grad_h2 = np.max(np.abs(grad_h_down), axis=2)  # test replacement 
    
    if 0: #Sierd_com; to be changed in future release
        ix = zb < zne + 0.005
        grad_h2[ix] = 0.
    
    grad_h = grad_h2

    #grad_h = np.sqrt(grad_h2)
    grad_h = np.max(np.abs(grad_h_down), axis=2)  # test replacement 

    max_grad_h = np.max(grad_h)
    

    return max_grad_h, grad_h, grad_h_down
    

