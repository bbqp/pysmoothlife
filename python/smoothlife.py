# smoothlife.py
#
# An implementation of Conway's game of life using floating point
# numerics. Life and death is determined based on birth and death
# intervals.

class SmoothLife:

    def __init__(x0=-10, xn=10, numx=101, y0=-10, yn=10, numy=101,
                 b0=0.278, bn=0.365, d0=0.267, dn=0.445, r0=1, rn=3,
                 aa_width=1, sig1=None, sig2=None,
                 alpham=0.147, alphan=0.028):

        self.x0 = x0
        self.xn = xn
        self.numx = numx
        
        self.y0 = y0
        self.yn = yn
        self.numy = numy
        
        self.b0 = b0
        self.bn = bn
        
        self.d0 = d0
        self.dn = dn
        
        self.r0 = r0
        self.rn = rn
        
        self.inner_areas = np.zeros(shape=(numx, numy), dtype=np.float32)
        self.outer_areas = np.zeros(shape=(numx, numy), dtype=np.float32)
        self.annls_areas = np.zeros(shape=(numx, numy), dtype=np.float32)
        
        # The inner, outer, and annulus areas.
        self.inner_area = np.pi * self.r0**2
        self.outer_area = np.pi * self.rn**2
        self.annls_area = self.outer_area - self.inner_area
        
        # Grid points for the indices.
        self.J, self.I = np.meshgrid(np.arange(numy), np.arange(numx))
        
        # Set the anti-aliasing width.
        self.aa_width = aa_width
        
        # Set the step for the transition functions.
        self.alpham = alpham
        self.alphan = alphan
        
        # Set the transition functions.
        if sig1 is None:
            sig1 = lambda x, a, alpha: (1 + exp(-4 * (x - a) / alpha))**(-1)

        if sig2 is None:
            sig2 = lambda x, a, b, alpha: sig1(x, a, alpha) \
                 * (1 - sig1(x, b, alpha))
            
        self.sig1 = sig1
        self.sig2 = sig2
        self.sigm = lambda x, y, m, alpha: \
                    x * (1 - sig1(m, 0.5, alpha)) + y * sig1(m, 0.5, alpha)
        self.s = lambda n, m: \
                 sig1(n, sigm(self.b0, self.d0, m, alpham), alphan) \
               * (1 - sig1(n, sigm(self.b1, self.d1, m, alpham), alphan))

    def cell_area(self, i, j):
        # Scaled analogues of the inner and outer radii.
        grid_r0 = (self.r0 / max(self.xn, self.yn))
        grid_rn = (self.rn / max(self.xn, self.yn))
    
        # Get the bounding box of the inner area, starting from the lower
        # left corner and moving counterclockwise. Note that the cells we
        # grab necessarily underestimate the area. The rationale is that
        # we're going for speed for now, so any improvements will need to
        # keep this in mind.
        bboxi = (int(np.ceiling(i - self.r0)) % numx,
                 int(np.floor(i + self.r0)) % numx)
        bboxj = (int(np.ceiling(j - self.r0)) % numy,
                 int(np.floor(j + self.r0)) % numy)
        
        # The inner, outer, and annulus areas for this cell.
        ci_area = 0.0
        co_area = 0.0
        an_area = 0.0
        
        # The box of gridpoints.
        grid_ia = self.inner_area[bboxi[
        
        # Now find the indices will the cells inside of the inner circle.
        inside_I = np.where(

    def compute_cell_areas(self):
        inner = self.inner_areas.copy()
        outer = self.outer_areas.copy()
        annls = outer.copy()

        for i in range(self.numx):
            for j in range(self.numy):
                ia, oa = self.cell_area(i, j)
                inner[i, j] = ia
                outer[i, j] = oa
                annls[i, j] = oa - ia

        self.inner_areas = inner
        self.outer_areas = outer
        self.annls_areas = annls
    
    def update_cells(self):
        alive = np.where((self.b0 < self.annls_area) & (self.annls_areas < self.bn))
        dead = np.where((self.d0 < self.annls_area) & (self.annls_areas < self.dn))

        self.cells[dead]

if __name__ == '__main__':
    print('Under construction.')