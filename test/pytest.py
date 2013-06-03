#! /usr/bin/python

# Load necessary modules
from pypdas.utility.random import sprandsym
from pypdas.bqp.clsbqp import BQP
from pypdas.algorithms import pdas
from copy import copy

# Generate a matlab data file containing a BQP with size=10, sparsity=0.5, cond=100 and with name 'test.mat'
sprandsym(100,0.5,0.01,1,'test')
a = BQP('test')
b = copy(a); c = copy(a); d = copy(a)

# Run different algorithms to solve the BQP
pdas.cgopt(a)
pdas.exupdate(b)
pdas.inexupdate(c)
pdas.inexupdate2(d)
