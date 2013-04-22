'''
Description : Primal-dual active-set method
Author      : Zheng Han
Usage       : To be filled
'''
from pypdas.algorithms.clscg import CG

def exactopt(bqp = None):
    if bqp is None:
        pass
    else:
        # Print title
        bqp.print_title()

        # Algorithm loop
        while True:
            # Subspace minimization
            bqp.ssm()
            # Print iteration
            bqp.print_iter()
            # Increase counter
            bqp.k += 1
            # Check optimality
            if bqp.kkt_error() < 1e-12:
                break
                return
            # New partition
            bqp.newp()

        # Print endline
        bqp.print_line()

# Apply CG to solve the subproblem, should replicate exactopt
def cgopt(bqp = None):
    if bqp is None:
        pass
    else:
        # Print title
        bqp.print_title()

        clscg = CG(bqp)
        # Algorithm loop
        while True:
            # Subspace minimization with CG iterations
            bqp.fix()
            clscg.applycg(rep = 1000)

            # Print iteration
            bqp.print_iter()
            # Increase counter
            bqp.k += 1
            # Check optimality
            if bqp.kkt_error() < 1e-12:
                break
                return
            # New partition
            bqp.newp()
            clscg.reset()
        # Print endline
        bqp.print_line()
