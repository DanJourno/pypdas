'''
Description : A bound constrained QP class
Author      : Zheng Han
Date        : 04/07/2013
Usage       : to be filled
'''
import numpy as np
from scipy import sparse as sp
from scipy.io import loadmat
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import randint, randn

class BQP:
    # Class constructor
    def __init__(self,filename = None, hessian = [], linear = [], ub = [], A = None, x0 = None, z0 = None):
        # If no data file specified, read the data sperately
        if filename is None:
            self.H = hessian     # Hessian
            self.c = linear      # Linear coefficient
            self.u = ub          # Upper bounds
        else:
            # Load Matlab data
            mlb = loadmat(filename)
            self.H = mlb["H"]
            self.c = mlb["c"]
            self.u = mlb["u"]

        self.n = len(self.c)     # Problem size
        self.k = 0               # Iteration counter

        # Preprocess optional parameters
        if A is None:
            # Randomized initial A
            self.A = list(set(randint(self.n, size = self.n)))
        else:
            self.A = A           # Active-set estimate (optional)
        self. I = np.setdiff1d(range(self.n),self.A)

        if x0 is None:
            # Randomized initial x
            self.x = randn(self.n,1)    
        else:
            self.x = x0          # Initial point (optional)
        if z0 is None:
            # Randomized initial z
            self.z = randn(self.n,1)    
        else:
            self.z = z0          # Initial point (optional)

            
    # Subspace minimization (exact)
    def ssm(self, A = None):
        if A is None:
            A = self.A

        # Fix active primal variables (x[A])
        self.x[A] = self.u[A]
        I = np.setdiff1d(range(self.n),A)
        # Fix inactive dual variables (z[I])
        self.z[I] = 0

        # Solve inactive primal variables (x[I])
        if len(I) == 0:
            Hii = np.zeros((0,0))
            Hia = np.zeros((0,self.n))
        elif len(A) == 0:
            Hii = self.H
            Hia = np.zeros((self.n,0))
        else:
            Hia = self.H[[I],[A]]
            Hii = self.H[[I],[I]]

        if np.shape(Hii)[0] == 1:
            if len(A) == 0:
                rhs = self.c[I].ravel()
            else:
                rhs = (self.c[I]+Hia*self.u[A]).ravel()
            self.x[I] = -spsolve(Hii.data, rhs)[:,np.newaxis]
        elif np.shape(Hii)[0] != 0:
            # self.x[I] = -np.atleast_2d(spsolve(Hii,self.c[I]+Hia*self.u[A])).T
            if len(A) == 0:
                self.x = -spsolve(self.H,self.c)[:,np.newaxis]
            else:
                self.x[I] = -spsolve(Hii,self.c[I]+Hia*self.u[A])[:,np.newaxis]

        # Compute active dual variables (z[A])
        if len(A) != 0:
            self.z[A] = -self.H[A,:]*self.x - self.c[A]

    # Update partition
    def newp(self, x = None, z = None):
        if x is None:
            x = self.x
        if z is None:
            z = self.z

        # Violated x
        Vx = np.where(x > self.u)[0]
        # Violated z
        Vz = np.where(z < 0)[0]

        # Obtain new A and I
        self.A = np.union1d(Vx, np.setdiff1d(self.A, Vz) )
        self.I = np.union1d(Vz, np.setdiff1d(self.I, Vx) )

    # Print iteration
    def print_iter(self, k = None):
        if k is not None:
            count = k
        else:
            count = self.k

        print '{0:3d} {1:4d} {2:4d}   {3:+.3e}  {4:.3e}'.format(count,len(self.A),self.n-len(self.A),self.obj(),float(self.kkt_error()) )

    # Print title
    def print_title(self,rep = None):
        if rep is None:
            rep = 60
        print '='*rep
        print 'Iter |A|   |I|     obj        res'
        print '='*rep
    # Print frameline
    def print_line(self,rep = None):
        if rep is None:
            rep = 60
        print '='*rep

    # Compute the residual norm
    def kkt_error(self):
        dx = self.u - self.x
        norm_dx = norm(dx[np.where(dx < 0)])
        norm_dz = norm(self.z[np.where(self.z < 0)])
        return norm([norm_dx, norm_dz])

    # Compute the objective function value
    def obj(self):
        return np.dot(self.x.T, self.H*self.x/2 + self.c)[0,0]

# Call as function
if __name__ == "__main__":
    print "This is a BQP class"

