# cython: language_level=3
import numpy as np
import adios2 as ad2
from tqdm import tqdm

import os
import subprocess

from math import sqrt, floor, exp
import torch

import sys
sys.path.insert(0, '')

from XGC import hello

cdef public void call_hello(int x):
    hello(x)