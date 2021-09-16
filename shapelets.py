import math
import numpy
import matplotlib.pyplot as plt
from numpy.polynomial.hermite import hermval
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import scipy.integrate as integrate

box = 20
xc = box/2.0
nmax = 15

beta = int(raw_input('beta='))
n = int(raw_input('n='))

def B1(x, n, beta):
      x = x - xc
      coef = numpy.zeros(n+1)
      coef[n] = 1
      fi_x = (2.0**n*math.sqrt(math.pi)*math.factorial(n))**(-0.5)*numpy.exp(-x**2.0/beta**2.0/2.0)*hermval(x/beta,coef)
      return beta**(-0.5)*fi_x

x=numpy.arange(0,box,0.5)
plt.plot(x,B1(x,n,beta))
plt.xlabel(r'$x$')
plt.ylabel(r'$fi_n(x)$')
plt.title(r'shapelets')
plt.show()

m = int(raw_input('m='))

from scipy import integrate
func = lambda x: B1(x,n,beta)*B1(x,m,beta)
print integrate.quad(func, -numpy.inf, numpy.inf)[0]
