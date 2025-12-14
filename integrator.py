"""
integrator.py

A class for numerically integrating a function f(x,y)
on a circular grid.
"""


import numpy as np


class Integrator:

    def __init__(self, X, Y):
        """
        Integrator.__init__()

        Creates an instance of an Integrator object.

        Inputs:

            r0 - The inner radius of the annulus.
            r1 - The outer radius of the annulus.
            b  - The antialiasing parameter.

        Outputs:

        An instance of an Integrator object.
        """

        self.X = X
        self.Y = Y
        self.Z = np.zeros(X.shape)
        self.L = np.zeros(X.shape)
        self.hx = X[0,1] - X[0,0]
        self.hy = Y[1,0] - Y[0,0]

    def get_inner_antialiasing_indices(self, c, r0, b):
        self.Z = (self.X - c[0])**2 + (self.Y - c[1])**2
        ind = np.where((self.Z > (r0 - b / 2)**2) &
                       (self.Z <= (r0 + b / 2)**2))
        inner_antialiasing_ind = set((i, j) for i, j in zip(*ind))

        return inner_antialiasing_ind

    def get_center_antialiasing_indices(self, c, r0, r1, b):
        self.Z = (self.X - c[0])**2 + (self.Y - c[1])**2
        ind = np.where((self.Z > (r0 + b / 2)**2) & \
                       (self.Z <= (r1 - b / 2)**2))
        center_antialiasing_ind = set((i, j) for i, j in zip(*ind))

        return center_antialiasing_ind

    def get_outer_antialiasing_indices(self, c, r1, b):
        self.Z = (self.X - c[0])**2 + (self.Y - c[1])**2
        ind = np.where((self.Z > (r1 - b / 2)**2) &
                       (self.Z <= (r1 + b / 2)**2))
        outer_antialiasing_ind = set((i, j) for i, j in zip(*ind))

        return outer_antialiasing_ind

    def integrate(self, c, r0, r1, b, F):
        # Compute a placeholder for all of the squared distances from the
        # center c.
        self.Z = (self.X - c[0])**2 + (self.Y - c[1])**2
        self.L = np.sqrt(self.Z)
        
        # Compute the normalizing factors for the areas of the small disk of
        # radius r0 and the large annulus with radii (r0 - b/2) and (r1 + b/2).
        N = np.pi * r0**2
        M = np.pi * (r1**2 - r0**2)

        # A placeholder for the area of each cell on the grid.
        cell_area = self.hx * self.hy

        # Get the inner, center, and outer indices.
        circ = np.where(self.Z < (r0 - b / 2)**2)
        inner = np.where((self.Z > (r0 - b / 2)**2) & (self.Z <= (r0 + b / 2)**2))
        center = np.where((self.Z > (r0 + b / 2)**2) & (self.Z <= (r1 - b / 2)**2))
        outer = np.where((self.Z > (r1 - b / 2)**2) & (self.Z <= (r1 + b / 2)**2))

        # Now compute the inner and outer weights. We omit the center weights,
        # which are all just 1.0.
        inner_weights = (r0 + b/2 - self.L[inner[0], inner[1]]) / b
        outer_weights = (r1 + b/2 - self.L[outer[0], outer[1]]) / b

        # Compute the approximate integral for the annulus.
        annulus_integral = np.dot(inner_weights, F[inner[0], inner[1]]) * cell_area
        annulus_integral += np.sum(F[center[0], center[1]]) * cell_area
        annulus_integral += np.dot(outer_weights, F[outer[0], outer[1]]) * cell_area
        annulus_integral /= M
        
        # Compute the approximate integral of the inner circle.
        inner_circle_integral = np.sum(F[circ[0], circ[1]]) * cell_area
        inner_circle_integral /= N

        return annulus_integral, inner_circle_integral


if __name__ == '__main__':
    # Bounds and point counts for each dimension.
    x0 = 0
    xn = 10
    y0 = 0
    yn = 10
    
    # Parameters for the annulus.
    c = x0 + 0.5 * (xn - x0), y0 + 0.5 * (yn - y0)
    r0 = 1.0
    r1 = 3.0
    b = 0.25

    # The dimenions, if each cell was a pixel.
    M = 40
    N = 40
    
    # The number of times we have to compute an integral (one for each cell).
    nruns = M * N
        
    # Create the grid.
    x = np.linspace(x0, xn, num=M-1, endpoint=False) + (xn - x0) / (2 * (M - 1))
    y = np.linspace(y0, yn, num=N-1, endpoint=False) + (yn - y0) / (2 * (N - 1))
    X, Y = np.meshgrid(x, y)

    # A function whose discrete values are given at the midpoint of each cell.
    F = np.ones(X.shape)
    
    # Create an integrator object.
    itgr = Integrator(X, Y)

    for run in range(nruns):
        print('Iteration {0:4d} of {1:4d}\r'.format(run + 1, nruns), end='')
        
        # Compute the approximate integral of F over the region.
        integral = itgr.integrate(c, r0, r1, b, F)
        
        # Display the integral's approximate value to the user.
        #print('Approximate integral over grid of size {0:d}: {1:f}'.format((N-1)**2, integral[0]))
        
        # For fun, pump out the exact answer.
        #print('Exact integral: {0:f}'.format(1))
        #print('')
    print('')