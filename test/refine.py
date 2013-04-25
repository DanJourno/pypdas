#! /home/zhh210/package/sage-5.8/sage
from glob import glob
import numpy as np

cgop = []
exop = []
exup = []
ineup= []
for size in ['b','q','w']:
    for cond in ['b','w','z']:
        cgop.append(size+'-cgop-'+cond)
        exop.append(size+'-exop-'+cond)
        exup.append(size+'-exup-'+cond)
        ineup.append(size+'-ineup-'+cond)

print '---cgop----\n'
for i in cgop:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4])

print '---exop----\n'
for i in exop:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3])

    
print '---exup----\n'
for i in exup:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4])

print '---ineup----\n'
for i in ineup:
    tmp = np.sum(np.genfromtxt(i)[:,1:-1],axis=0)/50
    print '{0:.1e}  {1:.1e}  {2:.2e}  {3:.2e}  {4:.4e}'.format(tmp[0],tmp[1],tmp[2],tmp[3],tmp[4])
