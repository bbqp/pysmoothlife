import ctypes
import numpy as np
import os


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


class SmoothlifeState:
    def __init__(self, libpath=None, x0=-1, xn=1, numx=129, y0=-1, yn=1, numy=129, ri_ratio=0.05, b_ratio=0.25):
        if libpath is None:
            libpath = os.path.abspath('../../c/lib/smoothlife.so')
        
        self.libpath = libpath
        self.slib = ctypes.cdll.LoadLibrary(libpath)
    
        self.grid = Grid()
        
        self.x0 = ctypes.c_float(x0)
        self.xn = ctypes.c_float(xn)
        self.numx = ctypes.c_int(numx)
        
        self.y0 = ctypes.c_float(y0)
        self.yn = ctypes.c_float(yn)
        self.numy = ctypes.c_int(numy)
        
        self.ri_ratio = ctypes.c_float(ri_ratio)
        ri = 2 * ri_ratio
        ro = 3 * ri
        self.b = ctypes.c_float(b_ratio * (ro - ri))

        status = self.slib.grid_init(ctypes.byref(self.grid), self.x0, self.xn,
                                     self.numx, self.y0, self.yn, self.numy,
                                     self.ri_ratio, self.b)
        
        if status == 0:
            raise MemoryError('Unable to allocate memory for new grid object.')

    def resize(self, numx, numy):
        if (numx != self.numx.value) or (numy != self.numy.value):
            new_grid = Grid()
            status = slib.grid_init(ctypes.byref(new_grid), self.x0, self.xn,
                                    ctypes.c_float(numx), self.y0, self.yn,
                                    ctypes.c_float(numy), self.ri_ratio,
                                    self.b)
            
            if status != 0:
                self.clear()
                self.grid = new_grid
            else:
                raise MemoryError('Unable to allocate memory for new grid object.')

    def update(self):
        self.slib.grid_update_state_function(ctypes.byref(self.grid))
        
    def get_state_data_copy(self):
        ptr = self.grid.data
        dtype = np.float32
        count = self.grid.fend
        data = np.fromiter(iter=ptr, dtype=dtype, count=count)
        
        return np.reshape(data[self.grid.fstart:self.grid.fend],
                          shape=(self.grid.numx - 1, self.grid.numy - 1),
                          order='F')

    def set_value(self, val=0):
        self.slib.grid_set_value(ctypes.byref(self.grid), ctypes.c_float(val))
        
    def set_random(self):
        self.slib.grid_set_random(ctypes.byref(self.grid))
        
    def update(self):
        self.slib.grid_update_state_function(ctypes.byref(self.grid))
    
    def __del__(self):
        self.slib.grid_clear(ctypes.byref(self.grid))

    def print_grid(grid):
        if grid.data:
            print('grid.data:           ', hex(ctypes.addressof(self.grid.data.contents)))
        else:
            print('grid.data:            NULL')

        if grid.distsq_closest:
            print('grid.distsq_closest: ', hex(ctypes.addressof(self.grid.distsq_closest.contents)))
        else:
            print('grid.distsq_closest:  NULL')
            
        if grid.distsq_furthest:
            print('grid.distsq_furthest:', hex(ctypes.addressof(self.grid.distsq_furthest.contents)))
        else:
            print('grid.distsq_furthest: NULL')

        if grid.dist_midpoint:
            print('grid.dist_midpoint:  ', hex(ctypes.addressof(self.grid.dist_midpoint.contents)))
        else:
            print('grid.dist_midpoint:   NULL')

        if grid.weights:
            print('grid.weights:        ', hex(ctypes.addressof(self.grid.weights.contents)))
        else:
            print('grid.weights:         NULL')

        if grid.domain_indices:
            print('grid.domain_indices: ', hex(ctypes.addressof(self.grid.domain_indices.contents)))
        else:
            print('grid.domain_indices:  NULL')
        print('')

        print('grid.x0:  ', self.grid.x0)
        print('grid.xn:  ', self.grid.xn)
        print('grid.numx:', self.grid.numx)
        print('')

        print('grid.y0:  ', self.grid.y0)
        print('grid.y0:  ', self.grid.yn)
        print('grid.numy:', self.grid.numy)
        print('')

        print('grid.dx:', self.grid.dx)
        print('grid.dy:', self.grid.dy)
        print('grid.dA:', self.grid.dA)
        print('')
        
        print('grid.ri:      ', self.grid.ri)
        print('grid.ri_ratio:', self.grid.ri_ratio)
        print('grid.ro:      ', self.grid.ro)
        print('grid.b:      ', self.grid.b)
        print('')
            
        print('grid.cindexi_start:', self.grid.cindexi_start)
        print('grid.cindexi_end:  ', self.grid.cindexi_end)
        print('grid.cindexj_start:', self.grid.cindexj_start)
        print('grid.cindexj_end:  ', self.grid.cindexj_end)
        print('grid.cindexk_start:', self.grid.cindexk_start)
        print('grid.cindexk_end:  ', self.grid.cindexk_end)
        print('')
        
        print('grid.aindexi_start:', self.grid.aindexi_start)
        print('grid.aindexi_end:  ', self.grid.aindexi_end)
        print('grid.aindexj_start:', self.grid.aindexj_start)
        print('grid.aindexj_end:  ', self.grid.aindexj_end)
        print('grid.aindexk_start:', self.grid.aindexk_start)
        print('grid.aindexk_end:  ', self.grid.aindexk_end)
        print('')
        
        print('grid.cspani:', self.grid.cspani)
        print('grid.cspano:', self.grid.cspano)
        print('grid.cdiami:', self.grid.cdiami)
        print('grid.cdiamo:', self.grid.cdiamo)
        print('')

        print('grid.cstart:', self.grid.cstart)
        print('grid.cend:  ', self.grid.cend)
        print('grid.astart:', self.grid.astart)
        print('grid.aend:  ', self.grid.aend)
        print('grid.fstart:', self.grid.fstart)
        print('grid.fend:  ', self.grid.fend)
        print('')

        print('grid.circle_area: ', self.grid.circle_area)
        print('grid.annulus_area:', self.grid.annulus_area)
        print('')