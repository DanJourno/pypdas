'''
Description : Primal-dual active-set method
Author      : Zheng Han
Usage       : To be filled
'''
from pypdas.algorithms.clscg import CG
import time

def exactopt(bqp = None, limiter = 1000):
    # Initialize information collector
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

        # Print collected info
        for key,val in collector.items():
            print key.ljust(10)+':', val

# Apply CG to solve the subproblem, should replicate exactopt
def cgopt(bqp = None):
    # Initialize information collector
    collector = dict()

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
        while True:
            # Subspace minimization with CG iterations
            bqp.fix()
            clscg.applycg(rep = 1000)

            # Print iteration
            clscg.print_iter()
            # Increase counter
            bqp.k += 1
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

        # Print collected info
        for key,val in collector.items():
            print key.ljust(10)+':', val
