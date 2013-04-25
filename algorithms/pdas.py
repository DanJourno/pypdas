'''
Description : Primal-dual active-set methods
Author      : Zheng Han
Usage       : To be filled
'''
from pypdas.algorithms.clscg import CG
import time
import numpy as np
from numpy.linalg import norm
from scipy.sparse.linalg import cg,spsolve

def exactopt(bqp = None, limiter = 1000):
    # Initialize np.information collector
    collector = dict()

    if bqp is None:
        print "Empty BQP class."
        pass
    else:
        # Print title
        bqp.print_title()
        collector['Time'] = -time.time()
        collector['State'] = bqp.state = 'Suboptimal'
        # Algorithm loop
        for i in range(limiter):
            # Subspace minimization
            bqp.ssm()
            # Print iteration
            bqp.print_iter()
            # Increase counter
            bqp.k += 1
            # Check optimality
            if bqp.kkt_error() < 1e-12:
                bqp.state = 'Optimal'
                break
                return
            # New partition
            bqp.newp()

        collector['Time'] += time.time()
        collector['Iter'] = bqp.k
        collector['State'] = bqp.state
        # Print endline
        bqp.print_line()

        # Print collected np.info
        for key,val in collector.items():
            print key.ljust(10)+':', val
        return collector    

# Apply CG to solve the subproblem, should replicate exactopt
def cgopt(bqp = None,limiter = 1000):
    # Initialize np.information collector
    collector = dict()
    collector['Tcg'] = 0
    if bqp is None:
        print "Empty BQP class."
        pass
    else:
        collector['Time'] = -time.time()
        collector['State'] = bqp.state = 'Suboptimal'
        clscg = CG(bqp)
        # Print title
        clscg.print_title()
        # Algorithm loop
        for i in range(limiter):
            # Subspace minimization with CG iterations
            bqp.fix()
            clscg.applycg(rep = 1000)

            # Print iteration
            clscg.print_iter()
            # Increase counter
            bqp.k += 1
            collector['Tcg'] += clscg.k
            # Check optimality
            if bqp.kkt_error() < 1e-12:
                bqp.state = 'Optimal'
                break
                return
            # New partition
            bqp.newp()
            clscg.reset()

        collector['Time'] += time.time()
        collector['Iter'] = bqp.k

        collector['State'] = bqp.state
        # Print endline
        bqp.print_line()

        # Print collected np.info
        for key,val in collector.items():
            print key.ljust(10)+':', val
        return collector    

# Apply inexact method to obtain an exact update
def exupdate(bqp = None, freq = 1, limiter = 1000):
    # Initialize np.information collector
    collector = dict()
    collector['Tcg'] = 0
    # Initialize partition monitor
    pmonitor = dict()
    pmonitor['NumChange'] = 1
    if bqp is None:
        print "Empty BQP class."
        pass
    else:
        collector['Time'] = -time.time()
        collector['State'] = bqp.state = 'Suboptimal'
        clscg = CG(bqp)
        # Print title
        clscg.print_title()
        # Algorithm loop
        for i in range(limiter):
            # Fix known variables
            bqp.fix()


            # Check residual's infty norm |r_I|
#            nmHii = max(np.sum(abs(bqp.H[ [bqp.I],[bqp.I]]),axis=1))
#            nmHai = max(np.sum(abs(bqp.H[ [bqp.A],[bqp.I]]),axis=1))
            nmHii = norm(bqp.H[ [bqp.I],[bqp.I]].todense(),np.inf)
            nmHai = norm(bqp.H[ [bqp.A],[bqp.I]].todense(),np.inf)

            # Run CG until r is sufficiently small
            while True:
                if pmonitor['NumChange'] == 0:
                    freq = 1000
                    clscg.applycg(rep = freq)
                    bqp.k -= 1
                    break
                clscg.applycg(rep = freq)
                th1 = norm(bqp.u[bqp.I] - bqp.x[bqp.I],-np.inf)
                th2 = norm(bqp.z[bqp.A],-np.inf)
                th1 /= nmHii
                if nmHai != 0.0:
                    th2 /= nmHii*nmHai
                else:
                    th2 = np.inf
                the = min(th1,th2)
                if norm(clscg.r,np.inf) < the:
                    break

            # Print iteration
            clscg.print_iter()

                
            # Increase counter
            bqp.k += 1
            collector['Tcg'] += clscg.k
            # Check optimality
            if bqp.kkt_error() < 1e-12:
                bqp.state = 'Optimal'
                break
                return
            # New partition
            pmonitor['NumChange'] = bqp.newp()

            clscg.reset()

        collector['Time'] += time.time()
        collector['Iter'] = bqp.k
        collector['State'] = bqp.state
        # Print endline
        bqp.print_line()

        # Print collected np.info
        for key,val in collector.items():
            print key.ljust(10)+':', val
        return collector    

# Apply inexact method to obtain an inexact update
def inexupdate(bqp = None, freq = 1, limiter = 1000):
    # Initialize np.information collector
    collector = dict()
    collector['Tcg'] = 0
    if bqp is None:
        print "Empty BQP class."
        pass
    else:
        collector['Time'] = -time.time()
        collector['State'] = bqp.state = 'Suboptimal'
        # Initialize partition monitor
        pmonitor = dict()
        pmonitor['NumChange'] = 1

        clscg = CG(bqp)
        # Print title
        clscg.print_title()
        # Algorithm loop
        for i in range(limiter):
            # Fix known variables
            bqp.fix()


            # Run CG until a new partition is obtained or r is sufficiently small
            while True:
                clscg.applycg(rep = freq)
                if pmonitor['NumChange'] == 0:
                    freq = 1000
                    clscg.applycg(rep = freq)

                if len(clscg.r) > 1:
                    tmp1 = spsolve(clscg.A,clscg.r)
                elif len(clscg.r) == 1:
                    tmp1 = clscg.r/clscg.A.data
                else:
                    tmp = 0
                tmp2 = bqp.H[[bqp.A],[bqp.I]]*tmp1
                tmp = np.concatenate((tmp1,tmp2))[:,np.newaxis]
                # Violated x
                Vx = np.where((bqp.u - bqp.x < 0) & (bqp.u-bqp.x + tmp <0))[0]
                # Violated z
                Vz = np.where((bqp.z < 0) & (bqp.z +tmp<0))[0]
#                Vz = np.where(bqp.z < 0 & bqp.z +tmp<0)
                # Identified new A and I
                pmonitor['NumChange'] = len(Vx) + len(Vz)
                if pmonitor['NumChange'] > 0:
                    break

                if norm(clscg.r,np.inf) < 1.0e-16:
                    break

            # Print iteration
            clscg.print_iter()
            bqp.A = np.union1d(Vx, np.setdiff1d(bqp.A, Vz) )
            bqp.I = np.union1d(Vz, np.setdiff1d(bqp.I, Vx) )

                
            # Increase counter
            bqp.k += 1
            collector['Tcg'] += clscg.k
            # Check optimality
            if bqp.kkt_error() < 1e-12:
                bqp.state = 'Optimal'
                break
                return

            clscg.reset()

        collector['Time'] += time.time()
        collector['Iter'] = bqp.k
        collector['State'] = bqp.state
        # Print endline
        bqp.print_line()

        # Print collected np.info
        for key,val in collector.items():
            print key.ljust(10)+':', val
        return collector    
