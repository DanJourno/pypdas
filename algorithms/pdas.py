'''
Description : Primal-dual active-set method
Author : Zheng Han
Usage : To be filled
'''

def optimize(bqp = None):
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
