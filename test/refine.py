#! /home/zhh210/package/sage-5.8/sage
from glob import glob
import numpy as np

cgop = []
exop = []
exup = []
ineup= []
ineup2= []
for size in ['b','q','w']:
    for cond in ['b','w','z']:
        cgop.append(size+'-cgop-'+cond)
        exop.append(size+'-exop-'+cond)
        exup.append(size+'-exup-'+cond)
        ineup.append(size+'-ineup-'+cond)
        ineup2.append(size+'-ineup2-'+cond)

print '|---cgop----'
print '| Call CG to solve subproblems exactly\n'
print '------------------------------------------------------------'
print ' size     cond      iter     cg-iter     time      res-ratio'
print '------------------------------------------------------------'

for i in cgop:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}  {5:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5])

print '|---exop----'
print '| Call default to solve subproblems exactly\n'
for i in exop:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3])

    

print '|---exup----'
print '| Call CG to solve subproblems inexactly but obtain exact update\n'
print '------------------------------------------------------------'
print ' size     cond      iter     cg-iter     time      res-ratio'
print '------------------------------------------------------------'
for i in exup:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}  {5:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5])


print '|---ineup----'
print '| Call CG to solve subproblems inexactly and obtain inexact update\n'
print '------------------------------------------------------------'
print ' size     cond      iter     cg-iter     time      res-ratio'
print '------------------------------------------------------------'
for i in ineup:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}  {5:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5])

print '|---ineup2----'
print '| Call CG to solve subproblems inexactly and obtain inexact update(version 2)\n'
print '------------------------------------------------------------'
print ' size     cond      iter     cg-iter     time      res-ratio'
print '------------------------------------------------------------'
for i in ineup2:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}  {5:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5])
