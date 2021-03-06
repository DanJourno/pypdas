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
from scipy.sparse import identity
import copy

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
            NumChange = bqp.newp()
            if NumChange == 0:
                bqp.state = 'Optimal'
                break
                return

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
        collector['ResRatio'] = 0.0
        clscg = CG(bqp)
        # Print title
        clscg.print_title()
        # Algorithm loop
        for i in range(limiter):
            # Subspace minimization with CG iterations
            bqp.fix()
            # Record initial residual
            clscg.applycg(rep = 0)
            if len(clscg.r) > 0:
                r0 = norm(clscg.r.ravel(), np.inf)

            clscg.applycg(rep = 1000)
            if len(clscg.r) > 0 and r0 > 0:
                ratio = norm(clscg.r.ravel(), np.inf)/r0
                collector['ResRatio'] += ratio
            else:
                ratio = 0.0
            # Print iteration
            clscg.print_iter(rt=ratio)
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
        collector['ResRatio'] /= collector['Iter']
        # Print endline
        bqp.print_line(rep=75)

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
        collector['ResRatio'] = 0.0
        clscg = CG(bqp)
        # Print title
        clscg.print_title()
        # Algorithm loop
        for i in range(limiter):
            # Fix known variables
            bqp.fix()

            # Record initial residual
            clscg.applycg(rep = 0)
            if len(clscg.r) > 0:
                r0 = norm(clscg.r.ravel(), np.inf)
            # Check residual's infty norm |r_I|
#            nmHii = max(np.sum(abs(bqp.H[ [bqp.I],[bqp.I]]),axis=1))
#            nmHai = max(np.sum(abs(bqp.H[ [bqp.A],[bqp.I]]),axis=1))
            if len(bqp.I) != 0:
                nmHii = norm(np.linalg.inv(bqp.H[ [bqp.I],[bqp.I]].todense()),np.inf)
                nmHai = norm(bqp.H[ [bqp.A],[bqp.I]].todense(),np.inf)

    
            # Run CG until r is sufficiently small
            while True:
                if pmonitor['NumChange'] == 0:
                    clscg.applycg(rep = 1000)
                    bqp.k -= 1
                    ratio = norm(clscg.r.ravel(), np.inf)/r0
                    collector['ResRatio'] += ratio                
                    clscg.print_iter(rt=ratio)
                    break
                clscg.applycg(rep = freq)
                if len(bqp.I) == 0:
                    ratio = 0
                    break
                th1 = norm(bqp.u[bqp.I] - bqp.x[bqp.I],-np.inf)
                th2 = norm(bqp.z[bqp.A],-np.inf)
                th1 /= nmHii
                if nmHai != 0.0:
                    th2 /= nmHii*nmHai
                else:
                    th2 = np.inf
                the = min(th1,th2)
                if norm(clscg.r.ravel(),np.inf) < max(the, 1.0e-17):
                    if norm(clscg.r,np.inf) > 0:
                        ratio = norm(clscg.r.ravel(), np.inf)/r0
                    else:
                        ratio = 0.0
                    break

            # Try updating the partition to get NumChange
            tmpbqp = copy.copy(bqp)
            pmonitor['NumChange'] = tmpbqp.newp()

            # Print iteration if NumChange > 0 (suboptimal)
            if pmonitor['NumChange'] > 0:
                clscg.print_iter(rt=ratio)
                collector['ResRatio'] += ratio

            if bqp.kkt_error() < 1e-12 and clscg.k < 1:
                clscg.print_iter(rt=ratio)

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

            if pmonitor['NumChange'] > 0:
                clscg.reset()

        collector['Time'] += time.time()
        collector['Iter'] = bqp.k
        collector['State'] = bqp.state
        collector['ResRatio'] /= collector['Iter']
        # Print endline
        bqp.print_line(rep=75)

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
        collector['ResRatio'] = 0.0
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
            # Record initial residual
            clscg.applycg(rep = 0)
            if len(clscg.r) > 0:
                r0 = norm(clscg.r.ravel(), np.inf)
            else:
                r0 = 0.0
            # Run CG until a new partition is obtained or r is sufficiently small
            while True:
                clscg.applycg(rep = freq)

                if len(clscg.r) > 1:
                    tmp1 = spsolve(clscg.A,clscg.r)
                    if len(bqp.A) !=0:
                        tmp2 = bqp.H[[bqp.A],[bqp.I]]*tmp1
                        tmp = np.concatenate((tmp1,tmp2))[:,np.newaxis]
                    else:
                        tmp = np.concatenate((tmp1,[]))[:,np.newaxis]
                elif len(clscg.r) == 1:
                    tmp1 = clscg.r/clscg.A.data
                    tmp2 = bqp.H[[bqp.A],[bqp.I]]*tmp1
                    tmp = np.concatenate((tmp1,tmp2))[:,np.newaxis]
                else:
                    tmp = np.zeros((bqp.n,1))
                

                # Violated x
                Vx = np.where((bqp.u - bqp.x < 0) & (bqp.u-bqp.x + tmp <0))[0]
                # Violated z
                Vz = np.where((bqp.z < 0) & (bqp.z +tmp<0))[0]
#                Vz = np.where(bqp.z < 0 & bqp.z +tmp<0)
                # Identified new A and I
                pmonitor['NumChange'] = len(Vx) + len(Vz)
                if pmonitor['NumChange'] > 0:
                    if len(clscg.r) > 0 and r0 > 0:
                        ratio = norm(clscg.r.ravel(), np.inf)/r0
                    else:
                        ratio = 0.0
                    break
                elif pmonitor['NumChange'] == 0:
                    clscg.applycg(rep = 1000)
                    ratio = norm(clscg.r.ravel(), np.inf)/r0
                    break

                if norm(clscg.r.ravel(),np.inf) < 1.0e-16:
                    if len(clscg.r) > 0 and r0 > 0:
                        ratio = norm(clscg.r.ravel(), np.inf)/r0
                    else:
                        ratio = 0.0
                    break

            # Print iteration
            clscg.print_iter(rt=ratio)
            collector['ResRatio'] += ratio
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
        collector['ResRatio'] /= collector['Iter']
        # Print endline
        bqp.print_line(75)

        # Print collected np.info
        for key,val in collector.items():
            print key.ljust(10)+':', val
        return collector    


# Apply inexact method to obtain an inexact update by bounding on inv(Hii)
def inexupdate2(bqp = None, freq = 1, limiter = 1000, compB = True):
    # Initialize np.information collector
    collector = dict()
    collector['Tcg'] = 0
    if bqp is None:
        print "Empty BQP class."
        pass
    else:
        collector['Time'] = -time.time()
        collector['State'] = bqp.state = 'Suboptimal'
        collector['ResRatio'] = 0.0
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

            # Check if I is empty
            clscg.applycg(rep = 0)
            if len(clscg.r) == 0:
                ratio = 0.0
                clscg.applycg(rep = 0)
            else: # I is not empty
                # Initial residual norm
                r0 = norm(clscg.r.ravel(), np.inf)
                # Bound on the norm of inv(Hii)
                if compB:
                    B = norm(np.linalg.inv(clscg.A.todense()),1)
                else:
                    s = norm(bqp.H.todense(), np.inf)+1
                    B = np.sqrt(len(clscg.r))/( s - norm((s*identity(len(clscg.r)) - clscg.A).todense(), np.inf)/np.sqrt(len(clscg.r)) )
                    assert B > 0, "B should > 0"

                # Run CG until a new partition is obtained or r is sufficiently small
                while True:
                    # Run rep CG steps
                    clscg.applycg(rep = freq)
                    # Bound on the pertubation
                    pIL = B*np.min(clscg.r)
                    pIU = B*np.max(clscg.r)
                    assert pIL <= pIU, 'pIL should <= pIU'
                    # Obtain lower and upper bound on xI
                    xIL = copy.copy(bqp.x)
                    xIU = copy.copy(bqp.x)
                    xIL[bqp.I] -= pIU
                    xIU[bqp.I] -= pIL

                    # Obtain lower and upper bound on zA
                    if len(bqp.A) > 0:
                        zAU = copy.copy(bqp.z)
                        pAU = -pIL*np.dot((np.sign(-bqp.H[bqp.A,:][:,bqp.I].todense())),np.ones(len(bqp.I)))
                        zAU[bqp.A] += np.array(pAU).ravel()[:,np.newaxis]
                    else:
                        zAU = bqp.z
                
                    # Violated x
                    Vx = np.where(bqp.u - xIL <0)[0]
                    #assert set(Vx).issubset(set(np.where(bqp.u-bqp.x<0)[0])), 'subset Vx'

                    # Violated z
                    Vz = np.where(zAU < 0)[0]
                    #assert set(Vz).issubset(set(np.where(bqp.z<0)[0])),'subset Vz'

                    pmonitor['NumChange'] = len(Vx) + len(Vz)

                    if pmonitor['NumChange'] > 0:
                        if len(clscg.r) > 0 and r0 > 0:
                            ratio = norm(clscg.r.ravel(), np.inf)/r0
                        else:
                            ratio = 0.0
                        break
                        
                    # Check if optimality is reached
                    if norm(clscg.r.ravel(),np.inf) < 1.0e-16:
                        if len(clscg.r) > 0 and r0 > 0:
                            ratio = norm(clscg.r.ravel(), np.inf)/r0
                        else:
                            ratio = 0.0
                        break

                    
            # Print iteration
            clscg.print_iter(rt=ratio)
            if len(clscg.r) == 0:
                #import pdb;pdb.set_trace()
                bqp.newp()

            collector['ResRatio'] += ratio
            # Update partition
            try:
                bqp.A = np.union1d(Vx, np.setdiff1d(bqp.A, Vz) )
                bqp.I = np.union1d(Vz, np.setdiff1d(bqp.I, Vx) )
            except NameError:
                pass
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
        collector['ResRatio'] /= collector['Iter']
        # Print endline
        bqp.print_line(75)

        # Print collected np.info
        for key,val in collector.items():
            print key.ljust(10)+':', val
        return collector    
