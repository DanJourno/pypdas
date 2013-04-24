#! /home/zhh210/package/sage-5.8/sage
'''
A script to test the algorithms
from pypdas.
'''
import sys
import copy
sys.path.append("/home/zhh210/workspace")
from pypdas.algorithms.pdas import *
from pypdas.algorithms.clscg import *
from pypdas.bqp.clsbqp import *


#a = BQP('test')
a = BQP('randexample/qz9')
print type(a.H),type(a.c),type(a.u),type(a.A),type(a.I)
b = copy.copy(a)
c = copy.copy(a)
d = copy.copy(a)
print "---Test exactopt():---"
exactopt(a)
print "---Test cgopt():---"
cgopt(b)
print "---Test exupdate():---"
exupdate(c)
print "---Test inexupdate():---"
inexupdate(d)
