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

import os
import re
import logging
from datetime import datetime
import numpy as np

# package modules
from aeolis.utils import *


# initialize logger
logger = logging.getLogger(__name__)


# check if netCDF4 is available
try:
    import netCDF4
    HAVE_NETCDF = True
except ImportError:
    HAVE_NETCDF = False
    logger.warning('No netCDF4 available, output is disabled')


def initialize(outputfile, outputvars, s, p, dimensions):
    '''Create empty CF-compatible netCDF4 output file

    Parameters
    ----------
    outputfile : str
        Name of netCDF4 output file
    outputvars : dictionary
        Spatial grids to be written to netCDF4 output file
    s : dict
        Spatial grids
    p : dict
        Model configuration parameters
    dimensions : dict
        Dictionary that specifies a tuple with the named dimensions
        for each spatial grid (e.g. ('ny', 'nx', 'nfractions'))

    Examples
    --------
    >>> netcdf.initialize('aeolis.nc',
    ...                   ['Ct', 'Cu', 'zb'],
    ...                   ['avg', 'max'],
    ...                   s, p, {'Ct':('ny','nx','nfractions'),
    ...                          'Cu':('ny','nx','nfractions'),
    ...                          'zb':('ny','nx')})

    '''

    # abort if netCDF4 is not available
    if not HAVE_NETCDF:
        return
        
    with netCDF4.Dataset(outputfile, 'w') as nc:

        # add dimensions
        nc.createDimension('s', p['nx']+1)  

        # this is where a quasi 2D model is stored in 1D    
        if p['ny']+1 == 3:
            nc.createDimension('n', 1)
        else:
            nc.createDimension('n', p['ny']+1)
        
        nc.createDimension('time', 0)
        nc.createDimension('nv', 2)
        nc.createDimension('nv2', 4)
        nc.createDimension('layers', p['nlayers'])
        nc.createDimension('fractions', p['nfractions'])
            
        # add global attributes
        # see http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/formats/DataDiscoveryAttConvention.html
        nc.Conventions = 'CF-1.6'
        nc.Metadata_Conventions = 'Unidata Dataset Discovery v1.0'
        #nc.featureType = 'grid'
        #nc.cdm_data_type = 'grid'
        nc.standard_name_vocabulary = 'CF-1.6'
        nc.title = ''
        nc.summary = ''
        nc.source = 'AeoLiS'
        nc.id = ''
        nc.naming_authority = ''
        nc.time_coverage_start = ''
        nc.time_coverage_end = ''
        nc.time_coverage_resolution = ''
        nc.geospatial_lat_min = 0
        nc.geospatial_lat_max = 0
        nc.geospatial_lat_units = 'degrees_north'
        nc.geospatial_lat_resolution = ''
        nc.geospatial_lon_min = 0
        nc.geospatial_lon_max = 0
        nc.geospatial_lon_units = 'degrees_east'
        nc.geospatial_lon_resolution = ''
        nc.geospatial_vertical_min = 0
        nc.geospatial_vertical_max = 0
        nc.geospatial_vertical_units = ''
        nc.geospatial_vertical_resolution = ''
        nc.geospatial_vertical_positive = ''
        nc.institution = ''
        nc.creator_name = ''
        nc.creator_url = ''
        nc.creator_email = ''
        nc.project = ''
        nc.processing_level = ''
        nc.references = ''
        nc.keywords_vocabulary = 'NASA/GCMD Earth Science Keywords. Version 6.0'
        nc.keywords = ''
        nc.acknowledgment = ''
        nc.comment = ''
        nc.contributor_name = ''
        nc.contributor_role = ''
        nc.date_created = datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%MZ')
        nc.date_modified = datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%MZ')
        nc.date_issued = datetime.strftime(datetime.utcnow(), '%Y-%m-%dT%H:%MZ')
        nc.publisher_name = ''
        nc.publisher_email = ''
        nc.publisher_url = ''
        nc.history = ''
        nc.license = ''
        nc.metadata_link = '0'
            
        # add variables
        nc.createVariable('s', 'float32', (u's'))
        nc.variables['s'].long_name = 's-coordinate'
        nc.variables['s'].units = '1'
        nc.variables['s'].valid_min = -np.inf
        nc.variables['s'].valid_max = np.inf
        
        nc.createVariable('n', 'float32', (u'n'))
        nc.variables['n'].long_name = 'n-coordinate'
        nc.variables['n'].units = '1'
        nc.variables['n'].valid_min = -np.inf
        nc.variables['n'].valid_max = np.inf
        
        nc.createVariable('x', 'float32', (u'n', u's'))
        nc.variables['x'].long_name = 'x-coordinate'
        nc.variables['x'].standard_name = 'projection_x_coordinate'
        nc.variables['x'].units = 'm'
        nc.variables['x'].axis = 'X'
        nc.variables['x'].valid_min = -np.inf
        nc.variables['x'].valid_max = np.inf
        nc.variables['x'].bounds = 'x_bounds'
        nc.variables['x'].grid_mapping = 'crs'
        
        nc.createVariable('y', 'float32', (u'n', u's'))
        nc.variables['y'].long_name = 'y-coordinate'
        nc.variables['y'].standard_name = 'projection_y_coordinate'
        nc.variables['y'].units = 'm'
        nc.variables['y'].axis = 'Y'
        nc.variables['y'].valid_min = -np.inf
        nc.variables['y'].valid_max = np.inf
        nc.variables['y'].bounds = 'y_bounds'
        nc.variables['y'].grid_mapping = 'crs'
        
        nc.createVariable('layers', 'float32', (u'layers',))
        nc.variables['layers'].long_name = 'bed layers'
        nc.variables['layers'].units = '1'
        nc.variables['layers'].valid_min = 0
        nc.variables['layers'].valid_max = np.inf
        
        nc.createVariable('fractions', 'float32', (u'fractions',))
        nc.variables['fractions'].long_name = 'sediment fractions'
        nc.variables['fractions'].units = 'm'
        nc.variables['fractions'].valid_min = 0
        nc.variables['fractions'].valid_max = np.inf
        
        nc.createVariable('lat', 'float32', (u'n', u's'))
        nc.variables['lat'].long_name = 'latitude'
        nc.variables['lat'].standard_name = 'latitude'
        nc.variables['lat'].units = 'degrees_north'
        nc.variables['lat'].valid_min = -np.inf
        nc.variables['lat'].valid_max = np.inf
        nc.variables['lat'].bounds = 'lat_bounds'
        nc.variables['lat'].ancillary_variables = ''
        
        nc.createVariable('lon', 'float32', (u'n', u's'))
        nc.variables['lon'].long_name = 'longitude'
        nc.variables['lon'].standard_name = 'longitude'
        nc.variables['lon'].units = 'degrees_east'
        nc.variables['lon'].valid_min = -np.inf
        nc.variables['lon'].valid_max = np.inf
        nc.variables['lon'].bounds = 'lon_bounds'
        nc.variables['lon'].ancillary_variables = ''
        
        nc.createVariable('time', 'float64', (u'time',))
        nc.variables['time'].long_name = 'time'
        nc.variables['time'].standard_name = 'time'
        nc.variables['time'].units = 'seconds since %s' % p['refdate']
        nc.variables['time'].calendar = 'julian'
        nc.variables['time'].axis = 'T'
        nc.variables['time'].bounds = 'time_bounds'
        #nc.variables['time'].ancillary_variables = ''
            
        nc.createVariable('x_bounds', 'float32', (u's', u'n', u'nv'))
        nc.variables['x_bounds'].units = 'm'
        nc.variables['x_bounds'].comment = 'x-coordinate values at the upper and lower bounds of each pixel.'
            
        nc.createVariable('y_bounds', 'float32', (u's', u'n', u'nv'))
        nc.variables['y_bounds'].units = 'm'
        nc.variables['y_bounds'].comment = 'y-coordinate values at the left and right bounds of each pixel.'
            
        nc.createVariable('lat_bounds', 'float32', (u's', u'n', u'nv2'))
        nc.variables['lat_bounds'].units = 'degrees_north'
        nc.variables['lat_bounds'].comment = 'latitude values at the north and south bounds of each pixel.'
            
        nc.createVariable('lon_bounds', 'float32', (u's', u'n', u'nv2'))
        nc.variables['lon_bounds'].units = 'degrees_east'
        nc.variables['lon_bounds'].comment = 'longitude values at the west and east bounds of each pixel.'
            
        nc.createVariable('time_bounds', 'float32', (u'time', u'nv'))
        nc.variables['time_bounds'].units = 'seconds since %s' % p['refdate']
        nc.variables['time_bounds'].comment = 'time bounds for each time value'

        meta = parse_metadata(outputvars)
        for var0, exts in outputvars.items():

            if var0 not in s:
                continue

            if var0 not in dimensions:
                continue

            dims = ['time'] + [d[1:] for d in dimensions[var0]]
            dims = ['s' if d == 'x' else d for d in dims]
            dims = ['n' if d == 'y' else d for d in dims]

            for ext in exts:
                if ext is None:
                    var = var0
                else:
                    var = '%s_%s' % (var0, ext)
            
                nc.createVariable(var, 'float32', dims)
                nc.variables[var].long_name = var
                #nc.variables[var].standard_name = 'sea_surface_height_above_mean_sea_level'
                nc.variables[var].scale_factor = 1.0
                nc.variables[var].add_offset = 0.0
                nc.variables[var].valid_min = -np.inf
                nc.variables[var].valid_max = np.inf
                nc.variables[var].coordinates = ' '.join(dims)
                nc.variables[var].grid_mapping = 'crs'
                nc.variables[var].source = ''
                nc.variables[var].references = ''
                #nc.variables[var].cell_methods = ''
                #nc.variables[var].ancillary_variables = ''
                #nc.variables[var].comment = ''

                if meta[var0]['units']:
                    nc.variables[var].units = meta[var0]['units']
            
        nc.createVariable('crs', 'int32', ())
        nc.variables['crs'].grid_mapping_name = 'stereographic'
        nc.variables['crs'].epsg_code = 'EPSG:28992'
        nc.variables['crs'].semi_major_axis = 6377397.155
        nc.variables['crs'].semi_minor_axis = 6356078.96282
        nc.variables['crs'].inverse_flattening = 299.1528128
        nc.variables['crs'].latitude_of_projection_origin = 52.0922178
        nc.variables['crs'].longitude_of_projection_origin = 5.23155
        nc.variables['crs'].scale_factor_at_projection_origin = 0.9999079
        nc.variables['crs'].false_easting = 155000.0
        nc.variables['crs'].false_northing = 463000.0
        nc.variables['crs'].proj4_params = '+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +towgs84=565.4174,50.3319,465.5542,-0.398957388243134,0.343987817378283,-1.87740163998045,4.0725 +no_defs'

        # store static data
        nc.variables['s'][:] = np.arange(p['nx']+1)

        # this is where a quasi 2D model is stored in 1D
        if p['ny']+1 == 3:
            nc.variables['n'][:] = np.arange(1)
            nc.variables['x'][:,:] = s['x'][0,:]
            nc.variables['y'][:,:] = s['y'][0,:]
        else:
            nc.variables['n'][:] = np.arange(p['ny']+1)
            nc.variables['x'][:,:] = s['x']
            nc.variables['y'][:,:] = s['y']

        nc.variables['layers'][:] = np.arange(p['nlayers'])
        nc.variables['fractions'][:] = p['grain_size']

        nc.variables['lat'][:,:] = 0.
        nc.variables['lon'][:,:] = 0.
        nc.variables['x_bounds'][:,:] = 0.
        nc.variables['y_bounds'][:,:] = 0.
        nc.variables['lat_bounds'][:,:] = 0.
        nc.variables['lon_bounds'][:,:] = 0.
        
        
        # store model settings writing arrays as attributes is not supported anymore
        if 0:
            for k, v in p.items():
                if k.startswith('_'):
                    continue
                k = 'par_%s' % k
                if v is None:
                    nc.setncattr(k, -1)
                elif isinstance(v, bool):
                    nc.setncattr(k, int(v))
                else:
                    nc.setncattr(k, np.real(v))
        
        
        # Liz Debugging
        #print("Variables passed to append:", nc.variables) # CONFIRMED
        if 'bgbiomass' in nc.variables:
            print("Type of bgbiomass in variables:", type(nc.variables['bgbiomass'])) # bg biomass is passed to append as a netcdf variable as bgbiomass(time, n, s)

def append(outputfile, variables):
    '''Append variables to existing netCDF4 output file

    Increments the time axis length with one and appends the provided
    spatial grids along the time axis. The ``variables`` dictionary
    should at least have the ``time`` field indicating the current
    simulation time. The CF time bounds are updated accordingly.

    Parameters
    ----------
    outputfile : str
        Name of netCDF4 output file
    variables : dict
        Dictionary with spatial grids and time

    Examples
    --------
    >>> netcdf.append('aeolis.nc', {'time', 3600.,
    ...                             'Ct', np.array([[0.,0., ... ,0.]]),
    ...                             'Cu', np.array([[1.,1., ... ,1.]]))

    See Also
    --------
    set_bounds

    '''

    # abort if netCDF4 is not available
    if not HAVE_NETCDF:
        return
    
    #with netCDF4.Dataset(outputfile, 'a') as nc:
       # Check and extract bgbiomass if it exists in the file
       # if isinstance(nc.variables['bgbiomass'][:], np.ma.MaskedArray):
       #     print(nc.variables['bgbiomass'])
       #     variables['bgbiomass'] = nc.variables['bgbiomass'][:].filled(0)  # Replace masked values with 0
       #     print("Before extraction:", type(nc.variables['bgbiomass']))
       #     print("After extraction:", type(variables['bgbiomass']))
       #     print("Shape of bgbiomass data:", variables['bgbiomass'].shape)
       #     #print(variables['bgbiomass'])
       # else:
       #     variables['bgbiomass'] = np.array(nc.variables['bgbiomass'][:])
       
       # Liz Debugging
       # if 'bgbiomass' in nc.variables:
       #     bgbiomass_var = nc.variables['bgbiomass']
       #     print("Type of bgbiomass:", type(bgbiomass_var))
       #     print("initial shape of bgbiomass:", bgbiomass_var.shape)
       #     # if bgbiomass_var.shape[0] == 0:  # If the first dimension is 0
       #     #     print("bgbiomass is empty. Initializing with zeros.")
       #     #     variables['bgbiomass'] = np.zeros((1, *bgbiomass_var.shape[1:]))
       #     # else:
       #     #     variables['bgbiomass'] = bgbiomass_var[:].filled(0)  # Replace masked values with zeros
       #     # print("bgbiomass shape after extraction:", variables['bgbiomass'].shape)
       # else:
       #     print("bgbiomass not found in the NetCDF file")
       
       
       # if 'bgbiomass' in nc.variables:
       #     print("Before extraction:", type(nc.variables['bgbiomass']))
       #     variables['bgbiomass'] = np.array(nc.variables['bgbiomass'][:])
       #     print("After extraction:", type(variables['bgbiomass']))
       #     print("Shape of bgbiomass data:", variables['bgbiomass'].shape)
       #     print(variables['bgbiomass'])
       # else:
       #     print("bgbiomass not found in the NetCDF file")
       
    
    #with netCDF4.Dataset(outputfile, 'a') as nc:
        # # Check if bgbiomass exists & handle correctly
        # if 'bgbiomass' in nc.variables: # should be 
        #     bgbiomass_var = nc.variables['bgbiomass']
        #     print("bgbiomass before extraction:", type(bgbiomass_var))
        #     print("Initial shape of bgbiomass:", bgbiomass_var.shape)
            
        #     # test to see if rhoveg is different... it is NOT!
        #     rhoveg_var = nc.variables['rhoveg']
        #     print("rhoveg before extraction:", type(rhoveg_var)) # comes out as netCDF4 variable
        #     print("Initial shape of rhoveg:", rhoveg_var.shape)  # comes out as shape of (0, 6, 228)
            
        #     # Print rhoveg and bgbiomass
        #     print("contents of bgbiomass_var:")
        #     print(bgbiomass_var[:])
            
        #     print("contents of rhoveg_var:")
        #     print(rhoveg_var[:])            
            
        #     # Assign bgbiomass values
        #     variables['bgbiomass'] = np.array(bgbiomass_var[:])
            
        #     # Debug output
        #     print("After extraction:", type(variables['bgbiomass']))
        #     print("Shape of bgbiomass data:", variables['bgbiomass'].shape)
        #     print(variables['bgbiomass'])
        # else:
        #     print("bgbiomass not found in the NetCDF")
        
    with netCDF4.Dataset(outputfile, 'a') as nc:
        # Get the current time index and add the time data
        i = nc.variables['time'].shape[0]
        print(f"Current time index: {i}")
        
        # Check if time needs to be appended as an array (if it's scalar, wrap it in an array)
        time_to_append = np.array([variables['time']]) if np.isscalar(variables['time']) else variables['time']
        
        #nc.variables['time'][i] = variables['time']
        nc.variables['time'][:] = np.append(nc.variables['time'][:], time_to_append)
        
        # # Check if bgbiomass exists & handle correctly
        # if 'bgbiomass' in nc.variables: # should be 
        #     bgbiomass_var = nc.variables['bgbiomass']
        #     print("1. bgbiomass before extraction:", type(bgbiomass_var)) # comes out as netCDF variable
        #     print("2. Initial shape of bgbiomass:", bgbiomass_var.shape) # comes out as shape of (1, 6, 228)
            
        #     # Assign bgbiomass values
        #     variables['bgbiomass'] = np.array(bgbiomass_var[:])
            
        #     # Debug output
        #     print("3. Data type after extraction:", type(variables['bgbiomass']))
        #     print("4. Shape of bgbiomass data:", variables['bgbiomass'].shape)
        #     #print(variables['bgbiomass'])
        # else:
        #     print("bgbiomass not found in the NetCDF")
        
        # Append variables, including bgbiomass
        for k, v in variables.items():
            print(f"Processing variable: {k}")
            print(f"Type of variable {k}: {type(v)}")
            
            if k == 'time':
                print("Time variable has already been updated")
                #nc.variables['time'][:] = np.append(nc.variables['time'][:], v)
                continue    
        
            if k == 'bgbiomass': # handle bgbiomass specifically
                if 'bgbiomass' in nc.variables:
                    bgbiomass_var = nc.variables['bgbiomass']
                    #rhoveg_var = nc.variables['rhoveg']
                    #print("1. bgbiomass before extraction:", type(bgbiomass_var)) # comes out as netCDF variable
                    #print("2. Initial shape of bgbiomass:", bgbiomass_var.shape) # comes out as shape of (1, 6, 228)
                    
                    # Debug the netcdf writing issues
                    print(f"1. Fill Value {getattr(bgbiomass_var, '_FillValue', 'None')}")
                    print(f" Shape: {bgbiomass_var.shape}")
                    
                    # Validate and debug bgbiomass data before writing
                    print("2. Data type of bgbiomass being written:", type(v))
                    print("3. Shape of bgbiomass being written:", v.shape)
                    print("4. Min value of bgbiomass:", np.min(v))
                    print("5. Max value of bgbiomass:", np.max(v))
                    print("6. Any NaN in bgbiomass:", np.isnan(v).any())
                    print("7. Any inf in bgbiomass:", np.isinf(v).any())
                    
                    # Write the bgbiomass data.... TEST
                    nc.variables['bgbiomass'][i, ...] = v
                    
                    # Extract & assign bgbiomass values
                    # v = np.array(bgbiomass_var[:])  # Extract bgbiomass data... THESE RUN, COMMENTING OUT TO SEE IF I CAN FIX A BUG... DO NOT DELETE!!
                    # variables['bgbiomass'] = v
                    #print("3. Data type after extraction:", type(v))
                    #print("4. Shape of bgbiomass data:", v.shape)
                else:
                    print("bgbiomass not found in the NetCDF")
                continue
            
            
            # this is where a quasi 2D model is stored in 1D
            if v.shape[0]==3:
                nc.variables[k][i,...] = v[1,:]
            # this is where a 2D model is stored in 2D
            else: 
                nc.variables[k][i,...] = v
            
            

    set_bounds(outputfile)


def set_bounds(outputfile):
    '''Sets CF time bounds

    Parameters
    ----------
    outputfile : str
        Name of netCDF4 output file

    '''
    
    # abort if netCDF4 is not available
    if not HAVE_NETCDF:
        return

    with netCDF4.Dataset(outputfile, 'a') as nc:
        i = nc.variables['time'].shape[0] - 1
        nc.variables['time_bounds'][i,0] = 0 if i == 0 else nc.variables['time'][i-1]
        nc.variables['time_bounds'][i,1] = nc.variables['time'][i]


def dump(outputfile, dumpfile, var='mass', ix=-1):
    '''Dumps time slice from netCDF4 output file to ASCII file

    This function can be used to use a specific time slice from a
    netCDF4 output file as input file for another AeoLiS model
    run. For example, the bed composition from a spinup run can be
    used as initial composition for other runs reducing the spinup
    time.

    Parameters
    ----------
    outputfile : str
        Name of netCDF4 output file
    dumpfile : str
        Name of ASCII dump file
    var : str, optional
        Name of spatial grid to be dumped (default: mass)
    ix : int
        Time slice index to be dumped (default: -1)

    Examples
    --------
    >>> # use bedcomp_file = bedcomp.txt in model configuration file
    ... netcdf.dump('aeolis.nc', 'bedcomp.txt', var='mass')

    '''

    # abort if netCDF4 is not available
    if not HAVE_NETCDF:
        return

    with netCDF4.Dataset(outputfile, 'r') as ds:
        m = ds.variables[var][ix,...]
        np.savetxt(dumpfile, m.reshape((-1, m.shape[1])))


def parse_metadata(outputvars):
    '''Parse metadata from constants.py

    Parses the Python comments in constants.py to extract meta data,
    like units, for the model state variables that can be used as
    netCDF4 meta data.

    Parameters
    ----------
    outputvars : dictionary
        Spatial grids to be written to netCDF4 output file

    Returns
    -------
    meta : dict
        Dictionary with meta data for the output variables

    '''

    pyfile = os.path.join(os.path.split(__file__)[0], 'constants.py')
    meta = {var:{'units':None} for var in outputvars.keys()}

    if os.path.exists(pyfile):
        with open(pyfile, 'r') as fp:
            for line in fp:
                m = re.match('^\s*\'(.*)\',\s*#\s*\[(.*)\]', line)
                if m:
                    var, units = m.groups()
                    if var in meta.keys():
                        if units == '-':
                            meta[var]['units'] = '1'
                        else:
                            meta[var]['units'] = units

    return meta
