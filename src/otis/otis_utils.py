#============================================================#
#         Welcome to the Orbiting TImescale for Sinking library
#============================================================
#         O  T  I  S
#         /\     /\     Created by Alice Damiano
#        /  \___/  \     contact: alice.damiano@inaf.it
#       (    o o    )
#       /     ^     \
#===============================================================
        
#===============================================================


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pandas as pd
from scipy.integrate import quad, odeint
from scipy.interpolate import RegularGridInterpolator, griddata
import math
import sympy as sp
from tqdm import tqdm
from scipy.integrate import solve_ivp


