#! /home/zhh210/package/sage-5.8/sage
'''
A script to test the pdas algorithms on randexamples/.
'''
import sys
import copy
import io
sys.path.append("/home/zhh210/workspace")
from pypdas.algorithms.pdas import *
from pypdas.algorithms.clscg import *
from pypdas.bqp.clsbqp import *

dic = dict()
dic['b'] = 1.0e+2
dic['q'] = 1.0e+3
dic['w'] = 1.0e+4
dic['z'] = 1.0e+6

size = sys.argv[1][0]
cond = sys.argv[1][1]
rept = sys.argv[1][2:]

inputname = 'randexample/'+sys.argv[1]
ouputname = 'randresults/'+size+cond+rept+'.out'

a = BQP(inputname)
b = copy.copy(a)
c = copy.copy(a)
d = copy.copy(a)

size = dic[size]
if size == 1e+4:
    size /= 2
cond = dic[cond]
rept = int(rept)

print "---Test exactopt():---"
ca = exactopt(a)
fa = io.open(ouputname+'.exop','wb')
sa = '{0:3d}  {1:.1e}  {2:.1e}  {3:4d}   {4:.4e}  {5}'.format(rept,size,cond,int(ca['Iter']),int(ca['Time']),ca['State'].rjust(6))

fa.write(sa+'\n')
fa.close()

print "---Test cgopt():---"
cb = cgopt(b)
fb = io.open(ouputname+'.cgop','wb')
sb = '{0:3d}  {1:.1e}  {2:.1e}  {3:4d}  {4:5d}  {5:.4e}  {6}'.format(rept,size,cond,int(cb['Iter']),int(cb['Tcg']),float(cb['Time']),cb['State'].rjust(6))
fb.write(sb+'\n')
fb.close()

print "---Test exupdate():---"
cc = exupdate(c)
fc = io.open(ouputname+'.exup','wb')
sc = '{0:3d}  {1:.1e}  {2:.1e}  {3:4d}  {4:5d}  {5:.4e}  {6}'.format(rept,size,cond,int(cc['Iter']),int(cc['Tcg']),float(cc['Time']),cc['State'].rjust(6))
fc.write(sc+'\n')
fc.close()

print "---Test inexupdate():---"
cd = inexupdate(d)
fd = io.open(ouputname+'.ineup','wb')
sd = '{0:3d}  {1:.1e}  {2:.1e}  {3:4d}  {4:5d}  {5:.4e}  {6}'.format(rept,size,cond,int(cd['Iter']),int(cd['Tcg']),float(cd['Time']),cd['State'].rjust(6))
fd.write(sd+'\n')
fd.close()
