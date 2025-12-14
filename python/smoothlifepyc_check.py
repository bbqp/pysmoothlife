import ctypes
import os

import numpy as np


class Tfunc_State(ctypes.Structure):
    _fields_ = [
        ('alpham', ctypes.c_float),
        ('alpham', ctypes.c_float),
        ('bint', ctypes.c_float * 2),
        ('dint', ctypes.c_float * 2),
    ]
    

class Grid(ctypes.Structure):
    _fields_ = [
        ('data', ctypes.POINTER(ctypes.c_float)),
        ('distsq_closest', ctypes.POINTER(ctypes.c_float)),
        ('distsq_furthest', ctypes.POINTER(ctypes.c_float)),
        ('dist_midpoint', ctypes.POINTER(ctypes.c_float)),
        ('weights', ctypes.POINTER(ctypes.c_float)),
        ('domain_indices', ctypes.POINTER(ctypes.c_int)),
        
        ('tfs', Tfunc_State),
        
        ('x0', ctypes.c_float),
        ('xn', ctypes.c_float),
        ('y0', ctypes.c_float),
        ('yn', ctypes.c_float),
        
        ('numx', ctypes.c_int),
        ('numy', ctypes.c_int),
        
        ('dx', ctypes.c_float),
        ('dy', ctypes.c_float),
        ('dA', ctypes.c_float),
        
        ('ri', ctypes.c_float),
        ('ri_ratio', ctypes.c_float),
        ('ro', ctypes.c_float),
        ('b', ctypes.c_float),
        
        ('cindexi_start', ctypes.c_int),
        ('cindexi_end', ctypes.c_int),
        ('cindexj_start', ctypes.c_int),
        ('cindexj_end', ctypes.c_int),
        ('cindexk_start', ctypes.c_int),
        ('cindexk_end', ctypes.c_int),
        
        ('aindexi_start', ctypes.c_int),
        ('aindexi_end', ctypes.c_int),
        ('aindexj_start', ctypes.c_int),
        ('aindexj_end', ctypes.c_int),
        ('aindexk_start', ctypes.c_int),
        ('aindexk_end', ctypes.c_int),
        
        ('cspani', ctypes.c_int),
        ('cspano', ctypes.c_int),
        ('cdiami', ctypes.c_int),
        ('cdiamo', ctypes.c_int),

        ('cstart', ctypes.c_int),
        ('cend', ctypes.c_int),
        ('astart', ctypes.c_int),
        ('aend', ctypes.c_int),
        ('fstart', ctypes.c_int),
        ('fend', ctypes.c_int),

        ('circle_area', ctypes.c_float),
        ('annulus_area', ctypes.c_float)
    ]
    

def print_grid(grid):
    if grid.data:
        print('grid.data:           ', hex(ctypes.addressof(grid.data.contents)))
    else:
        print('grid.data:            NULL')

    if grid.distsq_closest:
        print('grid.distsq_closest: ', hex(ctypes.addressof(grid.distsq_closest.contents)))
    else:
        print('grid.distsq_closest:  NULL')
        
    if grid.distsq_furthest:
        print('grid.distsq_furthest:', hex(ctypes.addressof(grid.distsq_furthest.contents)))
    else:
        print('grid.distsq_furthest: NULL')

    if grid.dist_midpoint:
        print('grid.dist_midpoint:  ', hex(ctypes.addressof(grid.dist_midpoint.contents)))
    else:
        print('grid.dist_midpoint:   NULL')

    if grid.weights:
        print('grid.weights:        ', hex(ctypes.addressof(grid.weights.contents)))
    else:
        print('grid.weights:         NULL')

    if grid.domain_indices:
        print('grid.domain_indices: ', hex(ctypes.addressof(grid.domain_indices.contents)))
    else:
        print('grid.domain_indices:  NULL')
    print('')

    print('grid.x0:  ', grid.x0)
    print('grid.xn:  ', grid.xn)
    print('grid.numx:', grid.numx)
    print('')

    print('grid.y0:  ', grid.y0)
    print('grid.y0:  ', grid.yn)
    print('grid.numy:', grid.numy)
    print('')

    print('grid.dx:', grid.dx)
    print('grid.dy:', grid.dy)
    print('grid.dA:', grid.dA)
    print('')
    
    print('grid.ri:      ', grid.ri)
    print('grid.ri_ratio:', grid.ri_ratio)
    print('grid.ro:      ', grid.ro)
    print('grid.b:      ', grid.b)
    print('')
        
    print('grid.cindexi_start:', grid.cindexi_start)
    print('grid.cindexi_end:  ', grid.cindexi_end)
    print('grid.cindexj_start:', grid.cindexj_start)
    print('grid.cindexj_end:  ', grid.cindexj_end)
    print('grid.cindexk_start:', grid.cindexk_start)
    print('grid.cindexk_end:  ', grid.cindexk_end)
    print('')
    
    print('grid.aindexi_start:', grid.aindexi_start)
    print('grid.aindexi_end:  ', grid.aindexi_end)
    print('grid.aindexj_start:', grid.aindexj_start)
    print('grid.aindexj_end:  ', grid.aindexj_end)
    print('grid.aindexk_start:', grid.aindexk_start)
    print('grid.aindexk_end:  ', grid.aindexk_end)
    print('')
    
    print('grid.cspani:', grid.cspani)
    print('grid.cspano:', grid.cspano)
    print('grid.cdiami:', grid.cdiami)
    print('grid.cdiamo:', grid.cdiamo)
    print('')

    print('grid.cstart:', grid.cstart)
    print('grid.cend:  ', grid.cend)
    print('grid.astart:', grid.astart)
    print('grid.aend:  ', grid.aend)
    print('grid.fstart:', grid.fstart)
    print('grid.fend:  ', grid.fend)
    print('')

    print('grid.circle_area: ', grid.circle_area)
    print('grid.annulus_area:', grid.annulus_area)
    print('')

def convert_pointer_to_array(ptr, dtype, size):
    return np.fromiter(iter=ptr, dtype=dtype, count=size)

if __name__ == '__main__':
    np.set_printoptions(linewidth=200)
    libpath = os.path.abspath('.') + os.sep + 'smoothlife.so'
    slib = ctypes.cdll.LoadLibrary(libpath)

    grid = Grid()
    
    x0 = ctypes.c_float(-1.0)
    xn = ctypes.c_float(1.0)
    numx = ctypes.c_int(129)
    
    y0 = ctypes.c_float(-1.0)
    yn = ctypes.c_float(1.0)
    numy = ctypes.c_int(129)
    
    ri_ratio = ctypes.c_float(0.025)
    ri = 2 * 0.025
    ro = 3 * ri
    b = ctypes.c_float(0.25 * (ro - ri))

    slib.grid_init(ctypes.byref(grid), x0, xn, numx, y0, yn, numy, ri_ratio, b)
    print_grid(grid)
    
    data = convert_pointer_to_array(grid.data, np.float32, grid.fend)
    weights = convert_pointer_to_array(grid.weights, np.float32, grid.aend)
    
    cweights = np.reshape(weights[:grid.cend], shape=(grid.cdiami, grid.cdiami), order='F')
    #aweights = np.reshape(weights[grid.astart:grid.astart], shape=(grid.cdiamo, grid.cdiamo), order='F')
    
    print(cweights)
    
    slib.grid_clear(ctypes.byref(grid))
    print_grid(grid)
    
    