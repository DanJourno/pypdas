'''
Description : A class that can apply conjugate-gradient iterations on BQP class
Author      : Zheng Han
Usage       : To be filled
'''
import numpy as np
from numpy.linalg import solve, norm

class CG:
    # Class constructor
    def __init__(self, bqp = None):
        self.k = 0
        self.r = None
        self.p = None
        if bqp is None:
            pass
        else:
            self.bqp = bqp
            if len(bqp.I) == 0:
                self.A = np.zeros((0,0))
                self.b = np.zeros((0,0))
            elif len(bqp.A) == 0:
                self.A = bqp.H
                self.b = -bqp.c.ravel()
            else:
                self.A = bqp.H[[bqp.I],[bqp.I]]
#                self.b = -bqp.c[bqp.I] - bqp.H[ [bqp.I],[bqp.A]]*bqp.u[bqp.A]
                self.b = -bqp.c[bqp.I] - bqp.H[bqp.I,:][:,bqp.A]*bqp.u[bqp.A]

    # Apply CG for some number of iterations
    # clsbqp.fix() should be called before it
    def applycg(self, rep = 1, tol = 1e-16):
        # Initial point
        if self.r is None:
            self.r = self.A*self.bqp.x[self.bqp.I] - self.b
            self.p = -self.r.copy()
        if len(self.b) == 0:
            return

        # Apply CG iteration(s) on clsbqp
        for i in range(rep):
            Ap = self.A*self.p
            alfa = np.dot(self.r.T,self.r)[0,0]/np.dot(self.p.T, Ap)[0,0]
            self.bqp.x[self.bqp.I] += alfa*self.p
            r1 = self.r + alfa*Ap
            beta = np.dot(r1.T,r1)[0,0]/np.dot(self.r.T,self.r)[0,0]
            self.p = -r1 + beta*self.p
            self.r = r1
            self.k += 1
            if norm(self.r,np.inf) < tol:
                break

        # Compute active dual variables (z[A])
        if len(self.bqp.A) != 0:
            self.bqp.z[self.bqp.A] = -self.bqp.H[self.bqp.A,:]*self.bqp.x - self.bqp.c[self.bqp.A]

    # Reset CG to initial state
    def reset(self):
        self.__init__(self.bqp)

    # Print iteration
    def print_iter(self, k = None):
        if k is not None:
            count = k
        else:
            count = self.k

        print '{0:3d} {1:4d} {2:4d}   {3:+.3e}  {4:.3e}    {5:4d}     {6:.3e}'.format(self.bqp.k,len(self.bqp.A),self.bqp.n-len(self.bqp.A),self.bqp.obj(),float(self.bqp.kkt_error()),count,norm(self.r,np.inf) )

    # Print title
    def print_title(self,rep = None):
        if rep is None:
            rep = 60
        print '='*rep
        print 'Iter |A|   |I|     obj        res       CG-iter    CG-res'
        print '='*rep
