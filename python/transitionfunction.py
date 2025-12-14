import numpy as np
import matplotlib
import matplotlib.pyplot as plt


class TransitionFunction:

    def __init__(self, alpham=0.147, alphan=0.028, bint=[0.278, 0.365], dint=[0.267, 0.445], dtype=np.float64):
        self.alpham = dtype(alpham)
        self.alphan = dtype(alphan)
        self.bint = dtype(bint)
        self.dint = dtype(dint)

    def sig1(self, x, a, alpha):
        return 1 / (1 + np.exp(-4*(x - a) / alpha))

    def sig2(self, x, a, b, alpha):
        return self.sig1(x, a, alpha) * (1 - self.sig1(x, b, alpha))

    def sig3(self, x, y, m, alpha):
        sig1m = self.sig1(m, 0.5, alpha)

        return x * (1 - sig1m) + y * sig1m

    def __call__(self, n, m):
        b1, b2 = self.bint
        d1, d2 = self.dint

        a = self.sig3(b1, d1, m, self.alpham)
        b = self.sig3(b2, d2, m, self.alpham)

        return self.sig2(n, a, b, self.alphan)

    def plot(self, npts=100, show=True, save=False, savefile='transition_function.png'):
        p0, p1 = dtype([0, 1])
        pts = np.linspace(p0, p1, npts)
        XX, YY = np.meshgrid(pts, pts)
        ZZ = self.__call__(XX, YY)

        fig, ax = plt.subplots(1, subplot_kw={'projection': '3d'})
        surf = ax.plot_surface(XX, YY, ZZ, cmap=matplotlib.cm.coolwarm, linewidth=0.5, antialiased=False)

        # Customize the z axis.
        ax.set_xlabel('Annular Average')
        ax.set_ylabel('Circular Average')

        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)

        if show:
            plt.show()

        if save:
            fig.savefig(savefile)
