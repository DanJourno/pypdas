#! /home/zheng/Program/sage-5.8/sage
'''
A script to test the algorithms
from pypdas.
'''
import sys
import copy
sys.path.append("/home/zheng/Documents/workspace/Python")
from pypdas.algorithms.pdas import *
from pypdas.algorithms.clscg import *
from pypdas.bqp.clsbqp import *


a = BQP('test')
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
