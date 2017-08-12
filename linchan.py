import numpy as np
import scipy.linalg

class LinearChannel(object):
    """Streamwise constant linearized channel flow
    """

    def __init__(self, ny, nz, re):
        self._ny = ny
        self._nz = nz
        self._re = re

    @property
    def ny(self):
        """Degree of Chebyshev polynomials used in y-direction"""
        return self._ny

    @property
    def nz(self):
        """Number of Fourier modes used in z-direction"""
        return self._nz

    @property
    def re(self):
        """Reynolds number"""
        return self._re


def cheb(n):
    """Compute differentiation matrix and grid for Chebyshev polynomials

    Collocation points consist of n + 1 points between -1 and 1

    Returns a tuple (D,x), where D is a differentiation matrix (a numpy array)
    and x is a numpy array containing the collocation points

    From Trefethen, Spectral Methods in Matlab, p. 54
    """
    if n == 0:
        D = 0
        x = 1
        return D, x

    x = np.cos(np.linspace(0, np.pi, n+1))
    # print(x)
    c = np.ones(n+1)
    c[0] = 2
    c[-1] = 2
    c[1::2] *= -1
    # print(c)
    X = np.outer(x, np.ones(n+1))
    dX = X - X.T
    D = np.outer(c, 1./c) / (dX + np.identity(n+1))
    D = D - np.diag(np.sum(D,1))
    # print(D)
    return D, x

def fourier(n):
    if n & 1 != 0:
        raise ValueError("n = %d is not even" % n)

    if n == 0:
        D = 0
        x = 1
        return D, x

    h = 2 * np.pi / n
    x = h * np.arange(n)
    column = np.zeros(n)
    column[1:] = 0.5 / np.tan(0.5 * x[1:])
    column[1::2] *= -1
    # print(column)
    D = scipy.linalg.toeplitz(column, -column)
    # print(D)
    return D, x


if __name__ == "__main__":
    ny = 16
    nz = 8
    reynolds = 100
    linchan = LinearChannel(ny, nz, reynolds)
    print("Linear channel: ny = %d, nz = %d, Re = %.0f" %
          (linchan.ny, linchan.nz, linchan.re))
    cheb(4)
    cheb(5)
    fourier(4)
    fourier(6)
