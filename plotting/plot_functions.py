from scipy.special import erf
from scipy.special import dawsn
from scipy.stats import gamma
from scipy.stats import moment
from scipy.stats import binom
from scipy.stats import rv_continuous
from scipy.stats import hmean
from scipy.integrate import quad
from scipy.integrate import trapezoid as trapz
from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess
import os
import fnmatch
import pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.cm as cm
# %matplotlib inline
import seaborn as sns
import numpy as np
import math
import random

# AA simulation replicates
parallel_AA = 49

# define functions for analytical approximations
N = 2000
U = 0.025
def pi(a, x): # the fixation probability under stationary stabilizing selection
    return (erf(a / 2.0) - erf(a / 2.0 * (1 - 2 * x))) / (2 * erf(a / 2.0))

def h_plus(a, x):
    return erf(a / 2.0) + erf(a / 2.0 * (1 - 2 * x))
def h_minus(a, x):
    return erf(a / 2.0) - erf(a / 2.0 * (1 - 2 * x))
def h(a, x):
    return math.sqrt(math.pi) / (2.0 * a) * math.exp(a * a / 4.0 * (1 - 2 * x) * (1 - 2 * x)) / (erf(a / 2.0) * x * (1 - x)) * h_minus(a, x) * h_plus(a, x)

def tau(a, x, p): # the sojourn time under stationary stabilizing selection
    assert (x > 0.0 or x == 0.0) and (x < 1.0 or x == 1.0)
    if x < p:
        return 2 * N * h(a, x) * h_plus(a, p) / h_plus(a, x)
    else:
        return 2 * N * h(a, x) * h_minus(a, p) / h_minus(a, x)
    
def tau_M(a, x): # the folded sojourn time tau_M(a, x) := tau(a, x, 1 / (2 * N)) + tau(a, 1 - x, 1 / (2 * N))
    assert (x > 0.0 or x == 0.0) and (x < 0.5 or x == 0.5)
    if x < 1.0 / (2 * N):
        return 2 * N * x * 2 * math.exp(-a * a * x * (1.0 - x)) / (x * (1.0 - x))
    else:
        return 2 * math.exp(-a * a * x * (1.0 - x)) / (x * (1.0 - x))
    
def pi_(a): # the fixation probability of the stationary distribution under stationary stabilizing selection
    numerator = quad(lambda x: pi(a, x) * tau_M(a, x), 0.0, 0.5, points = [1/(2*N)])[0]
    denominator = quad(lambda x: tau_M(a, x), 0.0, 0.5, points = [1/(2*N)])[0]
    return np.divide(numerator, denominator)

def heterozygosity_(a): # the heterozygosity under stationary stabilizing selection
    return quad(lambda x: 2.0 * x * (1.0 - x) * tau(a, x, 1/(2*N)), 0.0, 1.0, points = [1/(2*N), 1.0 - 1/(2*N)])[0]

def heterozygosity_truncating_rare_alleles_(a, maf_threshold): # the heterozygosity under stationary stabilizing selection with alleles whose MAFs are below a threshold, namely, maf_threshold, truncated (for the error-corrected McDonald-Kreitman test)
    return quad(lambda x: 2.0 * x * (1.0 - x) * tau(a, x, 1/(2*N)), maf_threshold, 1.0 - maf_threshold, points=[1/(2*N)])[0]





def v_star(a, x): # the allelic contribution to the genetic variance
    return 2 * a * a * x * (1.0 - x)

def v(a, x): # the density of variance per unit mutational input under stationary stabilizing selection
    assert (x > 0.0 or x == 0.0) and (x < 0.5 or x == 0.5)
    return v_star(a, x) * tau_M(a, x)

def v_(a): # the marginal density of variance per unit mutational input under stationary stabilizing selection
    return quad(lambda x: v(a, x), 0.0, 0.5, points = [1/(2*N)])[0]

def V_A(E2Ns): # the genetic variance under stationary stabilizing selection
    V2Ns = E2Ns ** 2
    S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0., scale=float(V2Ns) / float(E2Ns))
    return float(2 * N * U) * quad(lambda ss: 4.0 * np.sqrt(np.abs(ss)) * dawsn(
        np.sqrt(np.abs(ss)) / 2.0) * S_dist.pdf(ss), 0.0, S_dist.ppf(0.9999999999))[0]









## recurrent shifts
# linear Lande
def E_Delta_x_stab_sel(a, x): # the expected change in allele frequency under stationary stabilizing selection, which is the same as under recurrent shifts using the linear Lande approximation for the change in allele frequency due to a shift
    return -a * a / (2 * N) * x * (1.0 - x) * (0.5 - x)

def jump(a, x, shift_s0, sigma_0_del): # the linear Lande approximation for the change in allele frequency due to a shift
    return a * x * (1.0 - x) * shift_s0 / sigma_0_del

def V_Delta_x_drift_and_shifts(a, x, shift_s0, sigma_0_del, p): # the variance of change in allele frequency under recurrent shifts, using the linear Lande approximation
    return x * (1.0 - x) / (2 * N) + p * jump(a, x, shift_s0, sigma_0_del) ** 2

def psi_shifts(a, y, shift_s0, sigma_0_del, p):
    return np.exp(-2 * quad(lambda z: E_Delta_x_stab_sel(a, z) / V_Delta_x_drift_and_shifts(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

def pi_recurrent_shifts(a, x, shift_s0, sigma_0_del, p): # the fixation probability under recurrent shifts, using the linear Lande approximation
    numerator = quad(lambda y: psi_shifts(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)

# def jump_stabilizing_selection(a, x, shift_s0, sigma_0_del):
#     t_1 = 2 * N / sigma_0_del ** 2 * math.log(abs(shift_s0) * sigma_0_del) # V_S = 2 * N
#     return a * x * (1.0 - x) * shift_s0 / sigma_0_del - a ** 2 / (2 * N) * x * (1.0 - x) * (0.5 - x) * (t_1 - (shift_s0 / sigma_0_del) ** 2 * 2 * N / t_1) # V_S = 2 * N

# def V_Delta_x_drift_and_shifts_stabilizing_selection(a, x, shift_s0, sigma_0_del, p):
#     return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_stabilizing_selection(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_stabilizing_selection(-a, x, shift_s0, sigma_0_del) ** 2

# def psi_shifts_stabilizing_selection(a, y, shift_s0, sigma_0_del, p):
#     return np.exp(-2 * quad(lambda z: E_Delta_x_stab_sel(a, z) / V_Delta_x_drift_and_shifts_stabilizing_selection(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

# def pi_recurrent_shifts_stabilizing_selection(a, x, shift_s0, sigma_0_del, p):
#     numerator = quad(lambda y: psi_shifts_stabilizing_selection(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
#     denominator = quad(lambda y: psi_shifts_stabilizing_selection(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
#     return np.divide(numerator, denominator)

# def jump_adjusted(a, x, shift_s0, sigma_0_del):
#     if x + jump_stabilizing_selection(a, x, shift_s0, sigma_0_del) < 0.5 / (2 * N):
#         return -x + jump(a, x, shift_s0, sigma_0_del) - jump_stabilizing_selection(a, x, shift_s0, sigma_0_del)
#     if x + jump_stabilizing_selection(a, x, shift_s0, sigma_0_del) > 1.0 - 0.5 / (2 * N):
#         return 1.0 - x + jump(a, x, shift_s0, sigma_0_del) - jump_stabilizing_selection(a, x, shift_s0, sigma_0_del)
#     return jump(a, x, shift_s0, sigma_0_del)

# def V_Delta_x_drift_and_shifts_adjusted(a, x, shift_s0, sigma_0_del, p):
#     return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_adjusted(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_adjusted(a, x, -shift_s0, sigma_0_del) ** 2 - p ** 2 * (0.5 * (jump_adjusted(a, x, shift_s0, sigma_0_del) + jump_adjusted(a, x, -shift_s0, sigma_0_del))) ** 2

# def psi_shifts_adjusted(a, y, shift_s0, sigma_0_del, p):
#     return math.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_adjusted(a, z, shift_s0, sigma_0_del) + 0.5 * p * jump_adjusted(a, z, -shift_s0, sigma_0_del)) / V_Delta_x_drift_and_shifts_adjusted(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

# def pi_recurrent_shifts_adjusted(a, x, shift_s0, sigma_0_del, p):
#     numerator = quad(lambda y: psi_shifts_adjusted(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
#     denominator = quad(lambda y: psi_shifts_adjusted(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
#     return np.divide(numerator, denominator)

# def V_Delta_x_drift_and_shifts_t_1(a, x, shift_s0, sigma_0_del, p):
#     t_1 = 2 * N / sigma_0_del ** 2 * np.log(abs(shift_s0) * sigma_0_del)
#     return t_1 * x * (1.0 - x) / (2 * N) + p * jump(a, x, shift_s0, sigma_0_del) ** 2

# def V_Delta_x_drift_and_shifts_t_1_stabilizing_selection(a, x, shift_s0, sigma_0_del, p):
#     t_1 = 2 * N / sigma_0_del ** 2 * np.log(abs(shift_s0) * sigma_0_del)
#     return t_1 * x * (1.0 - x) / (2 * N) + 0.5 * p * jump_stabilizing_selection(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_stabilizing_selection(a, x, -shift_s0, sigma_0_del) ** 2 - ((1.0 - p) * t_1 * E_Delta_x_stab_sel(a, x) + 0.5 * p * jump_stabilizing_selection(a, x, shift_s0, sigma_0_del) + 0.5 * p * jump_stabilizing_selection(a, x, shift_s0, sigma_0_del)) ** 2




# nonlinear Lande
def F_d(a, x, shift_s0, sigma_0_del):
    return 1 / a * (np.exp(a * shift_s0 / sigma_0_del) - 1) / (1 + x * (np.exp(a * shift_s0 / sigma_0_del) - 1))

def jump_nonlinear(a, x, shift_s0, sigma_0_del): # the nonlinear Lande approximation for the change in allele frequency due to a shift
    return a * x * (1.0 - x) * F_d(a, x, shift_s0, sigma_0_del)

def V_Delta_x_drift_and_shifts_nonlinear(a, x, shift_s0, sigma_0_del, p): # the variance of change in allele frequency under recurrent shifts, using the nonlinear Lande approximation
    return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_nonlinear(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_nonlinear(a, x, -shift_s0, sigma_0_del) ** 2 - p ** 2 * (0.5 * (jump_nonlinear(a, x, shift_s0, sigma_0_del) + jump_nonlinear(a, x, -shift_s0, sigma_0_del))) ** 2

def psi_shifts_nonlinear(a, y, shift_s0, sigma_0_del, p):
    return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_nonlinear(a, z, shift_s0, sigma_0_del) + 0.5 * p * jump_nonlinear(a, z, -shift_s0, sigma_0_del)) / V_Delta_x_drift_and_shifts_nonlinear(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

def pi_recurrent_shifts_nonlinear(a, x, shift_s0, sigma_0_del, p): # the fixation probability under recurrent shifts, using the nonlinear Lande approximation
    numerator = quad(lambda y: psi_shifts_nonlinear(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonlinear(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)

# def barD_L(shift_s0, sigma_0_del, t):
#     return shift_s0 / sigma_0_del * 2 * N / t

# def bars_D(a, shift_s0, sigma_0_del, t):
#     return a * barD_L(shift_s0, sigma_0_del, t) / (2 * N)

# def barD_L_squared(shift_s0, sigma_0_del, t):
#     return 1 / 2 * shift_s0 ** 2 * 2 * N / t

# def I_D(shift_s0, sigma_0_del, t):
#     return shift_s0 / sigma_0_del * (1 - barD_L(shift_s0, sigma_0_del, t) / (shift_s0 * sigma_0_del) * (1 + 1 / 6 * (shift_s0 * sigma_0_del) ** 2 / (2 * N)))

# def I_S(shift_s0, sigma_0_del, t):
#     return 1 / 2 * t / (2 * 2 * N) * (1 - barD_L_squared(shift_s0, sigma_0_del, t) / (2 * N) * (2 - barD_L_squared(shift_s0, sigma_0_del, t) / (2 * N)))

# def bars_S_nonlinear(a, x, shift_s0, sigma_0_del, t, sign_I_D):
#     return -a ** 2 / (2 * N) * ((1 / 2 - x) * (1 - barD_L_squared(shift_s0, sigma_0_del, t) / (2 * N)) - a * x * (1 - x) * (sign_I_D * I_D(shift_s0, sigma_0_del, t) - a * (1 / 2 - x) * I_S(shift_s0, sigma_0_del, t)))

# def F_t(a, x, shift_s0, sigma_0_del, t, sign_I_D): # F_d + stabilizing selection
#     return 1 / a * (math.exp((bars_D(a, shift_s0, sigma_0_del, t) + bars_S_nonlinear(a, x, shift_s0, sigma_0_del, t, sign_I_D)) * t) - 1) / (1 + x * (math.exp((bars_D(a, shift_s0, sigma_0_del, t) + bars_S_nonlinear(a, x, shift_s0, sigma_0_del, t, sign_I_D)) * t) - 1))

# def jump_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del):
#     sign_I_D = -shift_s0 / abs(shift_s0)
#     t_1 = 2 * N / sigma_0_del ** 2 * math.log(abs(shift_s0) * sigma_0_del)
#     return a * x * (1 - x) * F_t(a, x, shift_s0, sigma_0_del, t_1, sign_I_D)

# def V_Delta_x_drift_and_shifts_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del, p):
#     return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_stabilizing_selection_nonlinear(a, x, -shift_s0, sigma_0_del) ** 2 - p ** 2 * (0.5 * (jump_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del) + jump_stabilizing_selection_nonlinear(a, x, -shift_s0, sigma_0_del))) ** 2

# def psi_shifts_stabilizing_selection_nonlinear(a, y, shift_s0, sigma_0_del, p):
#     t_1 = 2 * N / sigma_0_del ** 2 * math.log(abs(shift_s0) * sigma_0_del)
#     return np.exp(-2 * quad(lambda z: ((1.0 - p * round(t_1)) * E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_stabilizing_selection_nonlinear(a, z, shift_s0, sigma_0_del) + 0.5 * p * jump_stabilizing_selection_nonlinear(a, z, -shift_s0, sigma_0_del)) / V_Delta_x_drift_and_shifts_stabilizing_selection_nonlinear(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

# def pi_recurrent_shifts_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del, p):
#     numerator = quad(lambda y: psi_shifts_stabilizing_selection_nonlinear(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
#     denominator = quad(lambda y: psi_shifts_stabilizing_selection_nonlinear(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
#     return np.divide(numerator, denominator)

# def V_Delta_x_drift_and_shifts_t_1_nonlinear(a, x, shift_s0, sigma_0_del, p):
#     t_1 = 2 * N / sigma_0_del ** 2 * np.log(abs(shift_s0) * sigma_0_del)
#     return t_1 * x * (1.0 - x) / (2 * N) + 0.5 * p * jump_nonlinear(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_nonlinear(a, x, -shift_s0, sigma_0_del) ** 2 - (0.5 * p * jump_nonlinear(a, x, shift_s0, sigma_0_del) + 0.5 * p * jump_nonlinear(a, x, -shift_s0, sigma_0_del)) ** 2

# def psi_shifts_t_1_nonlinear(a, y, shift_s0, sigma_0_del, p):
#     return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_nonlinear(a, z, shift_s0, sigma_0_del) + 0.5 * p * jump_nonlinear(a, z, -shift_s0, sigma_0_del)) / V_Delta_x_drift_and_shifts_t_1_nonlinear(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

# def pi_recurrent_shifts_t_1_nonlinear(a, x, shift_s0, sigma_0_del, p):
#     numerator = quad(lambda y: psi_shifts_t_1_nonlinear(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
#     denominator = quad(lambda y: psi_shifts_t_1_nonlinear(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
#     return np.divide(numerator, denominator)

# def V_Delta_x_drift_and_shifts_t_1_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del, p):
#     t_1 = 2 * N / sigma_0_del ** 2 * np.log(abs(shift_s0) * sigma_0_del)
#     return t_1 * x * (1.0 - x) / (2 * N) + 0.5 * p * jump_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_stabilizing_selection_nonlinear(a, x, -shift_s0, sigma_0_del) ** 2 - ((1.0 - p) * t_1 * E_Delta_x_stab_sel(a, x) + 0.5 * p * jump_stabilizing_selection_nonlinear(a, x, shift_s0, sigma_0_del) + 0.5 * p * jump_stabilizing_selection_nonlinear(a, x, -shift_s0, sigma_0_del)) ** 2



# linear non-Lande
def f(a):
    if a == 0:
        return 0.0
    return 2 * a * a * a * math.exp(-a * a / 4) / (math.sqrt(math.pi) * erf(a / 2.0))

def g(a, E2Ns): # the distribution of the effect size squared
    return 1.0 / E2Ns * math.exp(-1.0 / E2Ns * a * a) * abs(a)

def C(E2Ns): # the amplification factor
    return quad(lambda a: v_(a) * g(a, E2Ns), 0, np.inf)[0] / quad(lambda a: f(a) * g(a, E2Ns), 0, np.inf)[0] - 1

def A(population_directory, E2Ns):
    eta = read_d2ax_over_shift(population_directory)
    return eta * (1 + C(E2Ns)) - 1

def jump_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value): # the linear non-Lande approximation for the change in allele frequency due to a shift
    if A_value is not None: # then A_value overwrites A(population_directory, E2Ns)
        return a * x * (1.0 - x) * (1 + A_value) * shift_s0 / sigma_0_del
    return a * x * (1.0 - x) * (1 + A(population_directory, E2Ns)) * shift_s0 / sigma_0_del

def V_Delta_x_drift_and_shifts_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value): # the variance of change in allele frequency under recurrent shifts, using the linear non-Lande approximation
    return x * (1.0 - x) / (2 * N) + p * jump_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) ** 2

def psi_shifts_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value):
    return np.exp(-2 * quad(lambda z: E_Delta_x_stab_sel(a, z) / V_Delta_x_drift_and_shifts_nonLande(a, z, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.1, y)[0])

def pi_recurrent_shifts_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value): # the fixation probability under recurrent shifts, using the linear non-Lande approximation
    numerator = quad(lambda y: psi_shifts_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)

# def jump_stabilizing_selection_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value):
#     t_1 = 2 * N / sigma_0_del ** 2 * math.log(shift_s0 * sigma_0_del) # V_S = 2 * N
#     if A_value is not None: # then A_value overwrites population_directory and E2Ns when calculating the value of A
#         return a * x * (1.0 - x) * (1 + A_value) * shift_s0 / sigma_0_del - a ** 2 / (2 * N) * x * (1.0 - x) * (0.5 - x) * (t_1 - ((1 + A_value) * shift_s0 / sigma_0_del) ** 2 * 2 * N / t_1) # V_S = 2 * N
#     return a * x * (1.0 - x) * (1 + A(population_directory, E2Ns)) * shift_s0 / sigma_0_del - a ** 2 / (2 * N) * x * (1.0 - x) * (0.5 - x) * (t_1 - ((1 + A(population_directory, E2Ns)) * shift_s0 / sigma_0_del) ** 2 * 2 * N / t_1) # V_S = 2 * N

# def V_Delta_x_drift_and_shifts_stabilizing_selection_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value):
#     return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_stabilizing_selection_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) ** 2 + 0.5 * p * jump_stabilizing_selection_nonLande(-a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) ** 2

# def psi_shifts_stabilizing_selection_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value):
#     return np.exp(-2 * quad(lambda z: E_Delta_x_stab_sel(a, z) / V_Delta_x_drift_and_shifts_stabilizing_selection_nonLande(a, z, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.1, y)[0])

# def pi_recurrent_shifts_stabilizing_selection_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value):
#     numerator = quad(lambda y: psi_shifts_stabilizing_selection_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, x)[0]
#     denominator = quad(lambda y: psi_shifts_stabilizing_selection_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, 1.0)[0]
#     return np.divide(numerator, denominator)



# nonlinear non-Lande
def F_d_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value):
    if A_value is not None: # then A_value overwrites A(population_directory, E2Ns)
        return 1 / a * (np.exp((1 + A_value) * a * shift_s0 / sigma_0_del) - 1) / (1 + x * (np.exp((1 + A_value) * a * shift_s0 / sigma_0_del) - 1))
    return 1 / a * (np.exp((1 + A(population_directory, E2Ns)) * a * shift_s0 / sigma_0_del) - 1) / (1 + x * (np.exp((1 + A(population_directory, E2Ns)) * a * shift_s0 / sigma_0_del) - 1))

def jump_nonlinear_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value): # the nonlinear non-Lande approximation for the change in allele frequency due to a shift
    return a * x * (1.0 - x) * F_d_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value)

def V_Delta_x_drift_and_shifts_nonlinear_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value): # the variance of change in allele frequency under recurrent shifts, using the nonlinear non-Lande approximation
    return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_nonlinear_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) ** 2 + 0.5 * p * jump_nonlinear_nonLande(a, x, -shift_s0, sigma_0_del, population_directory, E2Ns, A_value) ** 2 - p ** 2 * (0.5 * (jump_nonlinear_nonLande(a, x, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) + jump_nonlinear_nonLande(a, x, -shift_s0, sigma_0_del, population_directory, E2Ns, A_value))) ** 2

def psi_shifts_nonlinear_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value):
    return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_nonlinear_nonLande(a, z, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) + 0.5 * p * jump_nonlinear_nonLande(a, z, -shift_s0, sigma_0_del, population_directory, E2Ns, A_value)) / V_Delta_x_drift_and_shifts_nonlinear_nonLande(a, z, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.1, y)[0])

def pi_recurrent_shifts_nonlinear_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value): # the fixation probability under recurrent shifts, using the nonlinear non-Lande approximation
    numerator = quad(lambda y: psi_shifts_nonlinear_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonlinear_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)



# dissect the role of each part
def psi_shifts_nonlinear_drift(a, y, shift_s0, sigma_0_del, p):
    return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_nonlinear(a, z, shift_s0, sigma_0_del) + 0.5 * p * jump_nonlinear(a, z, -shift_s0, sigma_0_del)) / (z * (1.0 - z) / (2 * N)), 0.1, y)[0])

def pi_recurrent_shifts_nonlinear_drift(a, x, shift_s0, sigma_0_del, p): # the fixation probability under recurrent shifts, with the expected change in allele frequency using the nonlinear Lande approximation and the variance of change in allele frequency using the drift term
    numerator = quad(lambda y: psi_shifts_nonlinear_drift(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonlinear_drift(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)

def psi_shifts_nonlinear_linear_V_Delta_x(a, y, shift_s0, sigma_0_del, p):
    return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_nonlinear(a, z, shift_s0, sigma_0_del) + 0.5 * p * jump_nonlinear(a, z, -shift_s0, sigma_0_del)) / V_Delta_x_drift_and_shifts(a, z, shift_s0, sigma_0_del, p), 0.1, y)[0])

def pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a, x, shift_s0, sigma_0_del, p): # the fixation probability under recurrent shifts, with the expected change in allele frequency using the nonlinear Lande approximation and with the variance of change in allele frequency using the linear Lande approximation
    numerator = quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x(a, y, shift_s0, sigma_0_del, p), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)

def psi_shifts_nonlinear_linear_V_Delta_x_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value):
    return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * jump_nonlinear_nonLande(a, z, shift_s0, sigma_0_del, population_directory, E2Ns, A_value) + 0.5 * p * jump_nonlinear_nonLande(a, z, -shift_s0, sigma_0_del, population_directory, E2Ns, A_value)) / V_Delta_x_drift_and_shifts_nonLande(a, z, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.1, y)[0])

def pi_recurrent_shifts_nonlinear_linear_V_Delta_x_nonLande(a, x, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value): # the fixation probability under recurrent shifts, with the expected change in allele frequency using the nonlinear non-Lande approximation and with the variance of change in allele frequency using the linear non-Lande approximation
    numerator = quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x_nonLande(a, y, shift_s0, sigma_0_del, p, population_directory, E2Ns, A_value), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)

def V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a, x, shift_s0, sigma_0_del, p): # the nonlinear approximation for the variance of change in allele frequency, without the asymmetric term
    return x * (1.0 - x) / (2 * N) + 0.5 * p * jump_nonlinear(a, x, shift_s0, sigma_0_del) ** 2 + 0.5 * p * jump_nonlinear(a, x, -shift_s0, sigma_0_del) ** 2



def tau_shifts(a, x, p_initial, shift_s0, sigma_0_del, p): # the sojourn time under recurrent shifts
    assert (x > 0.0 or x == 0.0) and (x < 1.0 or x == 1.0)
    if x < p_initial:
        return 2 * (1 - pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a, p_initial, shift_s0, sigma_0_del, p)) / (V_Delta_x_drift_and_shifts(a, x, shift_s0, sigma_0_del, p) * psi_shifts_nonlinear_linear_V_Delta_x(a, x, shift_s0, sigma_0_del, p)) * quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x(a, y, shift_s0, sigma_0_del, p), 0.0, x)[0]
    else:
        return 2 * pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a, p_initial, shift_s0, sigma_0_del, p) / (V_Delta_x_drift_and_shifts(a, x, shift_s0, sigma_0_del, p) * psi_shifts_nonlinear_linear_V_Delta_x(a, x, shift_s0, sigma_0_del, p)) * quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x(a, y, shift_s0, sigma_0_del, p), x, 1.0)[0]

def tau_shifts_M(a, x, shift_s0, sigma_0_del, p): # the folded sojourn time under recurrent shifts
    assert (x > 0.0 or x == 0.0) and (x < 0.5 or x == 0.5)
    return (tau_shifts(a, x, 1.0 / (2 * N), shift_s0, sigma_0_del, p) + tau_shifts(a, 1.0 - x, 1.0 / (2 * N), shift_s0, sigma_0_del, p))

def tau_shifts_M_conditional_on_fixation(a, x, shift_s0, sigma_0_del, p): # the folded sojourn time conditional on fixation under recurrent shifts
    return tau_shifts_M(a, x, shift_s0, sigma_0_del, p) * pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a, x, shift_s0, sigma_0_del, p) / pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a, 1.0 / (2 * N), shift_s0, sigma_0_del, p)

def fixation_time_shifts(a, shift_s0, sigma_0_del, p): # the fixation time under recurrent shifts
    return quad(lambda x: tau_shifts_M_conditional_on_fixation(a, x, shift_s0, sigma_0_del, p), 0.0, 0.5, points=[1/(2*N)])[0]

def heterozygosity_shifts_(a, shift_s0, sigma_0_del, p): # the heterozygosity under recurrent shifts
    return quad(lambda x: 2.0 * x * (1.0 - x) * tau_shifts(a, x, p_initial=1/(2*N), shift_s0=shift_s0, sigma_0_del=sigma_0_del, p=p), 0.0, 1.0, points = [1/(2*N)])[0]

def heterozygosity_shifts_truncating_rare_alleles_(a, shift_s0, sigma_0_del, p, maf_threshold): # the heterozygosity under recurrent shifts with alleles whose MAFs are below a threshold, namely, maf_threshold, truncated (for the error-corrected McDonald-Kreitman test)
    return quad(lambda x: 2.0 * x * (1.0 - x) * tau_shifts(a, x, p_initial=1/(2*N), shift_s0=shift_s0, sigma_0_del=sigma_0_del, p=p), maf_threshold, 1.0 - maf_threshold, points=[1/(2*N)])[0]



# a distribution of the shift sizes
def psi_shifts_nonlinear_linear_V_Delta_x_gamma_distributed_shift_sizes(a, y, sigma_0_del, p):
    return np.exp(-2 * quad(lambda z: (E_Delta_x_stab_sel(a, z) + 0.5 * p * quad(lambda shift_s0: np.exp(-shift_s0) * jump_nonlinear(a, z, shift_s0, sigma_0_del), 0.0, np.inf)[0] + 0.5 * p * quad(lambda shift_s0: np.exp(-shift_s0) * jump_nonlinear(a, z, -shift_s0, sigma_0_del), 0.0, np.inf)[0]) / V_Delta_x_drift_and_shifts(a, z, np.sqrt(2), sigma_0_del, p), 0.1, y)[0])

def pi_recurrent_shifts_nonlinear_linear_V_Delta_x_gamma_distributed_shift_sizes(a, x, sigma_0_del, p): # the fixation probability under recurrent shifts, where the shift size follows a gamma distribution, for now with shape=1 and scale=1
    numerator = quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x_gamma_distributed_shift_sizes(a, y, sigma_0_del, p), 0.0, x)[0]
    denominator = quad(lambda y: psi_shifts_nonlinear_linear_V_Delta_x_gamma_distributed_shift_sizes(a, y, sigma_0_del, p), 0.0, 1.0)[0]
    return np.divide(numerator, denominator)





def v_shifts(a, x, shift_s0, sigma_0_del, p): # the density of variance per unit mutational input under recurrent shifts
    assert (x > 0.0 or x == 0.0) and (x < 0.5 or x == 0.5)
    return v_star(a, x) * tau_shifts_M(a, x, shift_s0, sigma_0_del, p)

def v_shifts_(a, shift_s0, sigma_0_del, p): # the marginal density of variance per unit mutational input under recurrent shifts
    return quad(lambda x: v_shifts(a, x, shift_s0, sigma_0_del, p), 0.0, 0.5, points=[1/(2*N)])[0]

def V_A_shifts(E2Ns, shift_s0, sigma_0_del, p): # the genetic variance under recurrent shifts
    V2Ns = E2Ns ** 2
    S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0., scale=float(V2Ns) / float(E2Ns))
    return 2 * N * U * quad(lambda ss: v_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=sigma_0_del, p=p) * S_dist.pdf(ss), 0.0, S_dist.ppf(0.9999999999))[0]

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ...
allows you to iterarte over all possible pairs in a list or set
EG. for v, w in pairwise([5,4,6])
    print v, w"""
    import itertools
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

# obtain the fraction of contribution to change in mean phenotype by standing variation, i.e., \eta, to calculate A
def read_d2ax_over_shift(population_directory):
    import os
    result_directory = os.path.join(population_directory, 'U0_frozen_d2ax_over_shift_efs_bins_H_sD')
    with open(os.path.join(result_directory, os.listdir(result_directory)[0]), 'rb') as f:
        import pickle
        x = pickle.load(f)
    return sum(x[1][list(x[1].keys())[-1]])

# obtain the number of fixed mutations in my AA simulation results, to calculate the fixation probability
def my_AA_results_num_fixed_nm(my_AA_directory, particular_2Ns, is_particular_2Ns):
    import os
    import pickle
    with open(os.path.join(my_AA_directory, 'fixation_probability'), 'rb') as f:
        num_fixed_nm = pickle.load(f)[0]['num_fixed_nm']
        
    if particular_2Ns is None:
        return [sum(x) for x in zip(num_fixed_nm, num_fixed_nm[::-1])][int(0.5 * len(num_fixed_nm) + 1):-1]
    if not is_particular_2Ns:
        return [sum(x) for x in zip(num_fixed_nm[:-2 * len(particular_2Ns)], num_fixed_nm[-1 - 2 * len(particular_2Ns)::-1])][int(0.5 * (len(num_fixed_nm) - 2 * len(particular_2Ns)) + 1):-1]
    return [sum(x) for x in zip(num_fixed_nm[-2 * len(particular_2Ns):], num_fixed_nm[:-1 - 2 * len(particular_2Ns):-1])][len(particular_2Ns):]

# obtain the number of new mutations in my AA simulation results, to calculate the fixation probability
def my_AA_results_num_new_mutations(my_AA_directory, particular_2Ns, is_particular_2Ns):
    import os
    import pickle
    with open(os.path.join(my_AA_directory, 'fixation_probability'), 'rb') as f:
        num_new_mutations = pickle.load(f)[1]
        
    if particular_2Ns is None:
        return [sum(x) for x in zip(num_new_mutations, num_new_mutations[::-1])][int(0.5 * len(num_new_mutations) + 1):-1]
    if not is_particular_2Ns:
        return [sum(x) for x in zip(num_new_mutations[:-2 * len(particular_2Ns)], num_new_mutations[-1 - 2 * len(particular_2Ns)::-1])][int(0.5 * (len(num_new_mutations) - 2 * len(particular_2Ns)) + 1):-1]
    return [sum(x) for x in zip(num_new_mutations[-2 * len(particular_2Ns):], num_new_mutations[:-1 - 2 * len(particular_2Ns):-1])][len(particular_2Ns):]

# obtain the integral heterozygosity in my AA simulation results
def my_AA_results_integral_heterozygosity(my_AA_directory):
    with open(os.path.join(my_AA_directory, 'fixation_probability'), 'rb') as f:
        integral_heterozygosity = pickle.load(f)[0]['integral_heterozygosity']
    
    return integral_heterozygosity[int(0.5 * len(integral_heterozygosity) + 1):-1]

# obtain the relative contribution to short-term phenotypic adaptation in my AA simulation results
def my_AA_results_d2ax_t_1(my_AA_directory, trait_increasing):
    """
    if trait_increasing, then convert to trait increasing
    """
    with open(os.path.join(my_AA_directory, 'fixation_probability'), 'rb') as f:
        d2ax_t_1 = pickle.load(f)[0]['d2ax_t_1']
        
    if not trait_increasing:
        return [sum(x) for x in zip(d2ax_t_1, d2ax_t_1[::-1])][int(0.5 * len(d2ax_t_1) + 1):-1]
    else:
        if sum(d2ax_t_1) >= 0:
            return [sum(x) for x in zip(d2ax_t_1, d2ax_t_1[::-1])][int(0.5 * len(d2ax_t_1) + 1):-1]
        else:
            return [-sum(x) for x in zip(d2ax_t_1, d2ax_t_1[::-1])][int(0.5 * len(d2ax_t_1) + 1):-1]

### calculate the quantities of interest, as a function of the effect size, the shift size, or the rate of shifts
def calculate_fixation_probability_mean_and_se_across_effect_sizes_from_num_fixed_and_num_new_mutations(my_AA_directories, indices_ss, num_new_mutations, nE):
    num_fixed = np.sum([read_pickle(os.path.join(my_AA_directory, 'num_fixed_nE_%d' % nE)) for my_AA_directory in my_AA_directories], axis=0)
    fixation_probability_mean = [num_fixed[index_ss + 1] / (num_new_mutations / (nE + 1) * len(my_AA_directories)) for index_ss in indices_ss]
    fixation_probability_se = np.sqrt(np.divide(np.multiply(fixation_probability_mean, np.subtract(1, fixation_probability_mean)), num_new_mutations / (nE + 1) * len(my_AA_directories)))
    return fixation_probability_mean, fixation_probability_se

def calculate_fixation_probability_mean_and_se_across_shifts_from_num_fixed_and_num_new_mutations(my_AA_directories_across_shifts, index_ss, num_new_mutations, nE):
    assert len(my_AA_directories_across_shifts[0]) == parallel_AA
    num_fixed_wrt_shift_s0 = [np.sum([read_pickle(os.path.join(my_AA_directory, 'num_fixed_nE_%d' % nE))[index_ss] for my_AA_directory in my_AA_directories_across_shifts[index_shift_s0]]) for index_shift_s0 in range(len(my_AA_directories_across_shifts))]
    fixation_probability_mean_wrt_shift_s0 = np.divide(num_fixed_wrt_shift_s0, num_new_mutations / (nE + 1) * len(my_AA_directories_across_shifts[0]))
    fixation_probability_se_wrt_shift_s0 = np.sqrt(np.divide(np.multiply(fixation_probability_mean_wrt_shift_s0, np.subtract(1, fixation_probability_mean_wrt_shift_s0)), num_new_mutations / (nE + 1) * len(my_AA_directories_across_shifts[0])))
    return fixation_probability_mean_wrt_shift_s0, fixation_probability_se_wrt_shift_s0

def calculate_heterozygosity_mean_and_se_across_effect_sizes(my_AA_directories, indices_ss, nE):
    heterozygosity_each_generation = [[h for my_AA_directory in my_AA_directories for h in read_pickle(os.path.join(my_AA_directory, 'heterozygosity_nE_%d' % nE))[index_ss + 1]] for index_ss in indices_ss]
    heterozygosity_mean = np.mean(heterozygosity_each_generation, axis=1)
    heterozygosity_se = np.divide(np.std(heterozygosity_each_generation, axis=1), np.sqrt(np.size(heterozygosity_each_generation, axis=1)))
    return heterozygosity_mean, heterozygosity_se

def calculate_heterozygosity_with_frequency_cutoff_mean_and_se_across_effect_sizes(my_AA_directories, indices_ss, nE, maf_threshold):
    heterozygosity_each_generation = [[h for my_AA_directory in my_AA_directories for h in read_pickle(os.path.join(my_AA_directory, 'heterozygosity_nE_%d_maf_threshold_' % nE + str(maf_threshold).replace('.', '_')))[index_ss + 1]] for index_ss in indices_ss]
    heterozygosity_mean = np.mean(heterozygosity_each_generation, axis=1)
    heterozygosity_se = np.divide(np.std(heterozygosity_each_generation, axis=1), np.sqrt(np.size(heterozygosity_each_generation, axis=1)))
    return heterozygosity_mean, heterozygosity_se

def calculate_heterozygosity_mean_and_se_across_shifts(my_AA_directories_across_shifts, index_ss, nE):
    heterozygosity_each_generation_wrt_shift_s0 = [[h for my_AA_directory in my_AA_directories_across_shifts[index_shift_s0] for h in read_pickle(os.path.join(my_AA_directory, 'heterozygosity_nE_%d' % nE))[index_ss + 1]] for index_shift_s0 in range(len(my_AA_directories_across_shifts))]
    heterozygosity_mean_wrt_shift_s0 = [np.mean(heterozygosity_each_generation) for heterozygosity_each_generation in heterozygosity_each_generation_wrt_shift_s0]
    heterozygosity_se_wrt_shift_s0 = [np.divide(np.std(heterozygosity_each_generation), np.sqrt(np.size(heterozygosity_each_generation))) for heterozygosity_each_generation in heterozygosity_each_generation_wrt_shift_s0]
    return heterozygosity_mean_wrt_shift_s0, heterozygosity_se_wrt_shift_s0

def calculate_heterozygosity_with_frequency_cutoff_mean_and_se_across_shifts(my_AA_directories_across_shifts, index_ss, nE, maf_threshold):
    heterozygosity_each_generation_wrt_shift_s0 = [[h for my_AA_directory in my_AA_directories_across_shifts[index_shift_s0] for h in read_pickle(os.path.join(my_AA_directory, 'heterozygosity_nE_%d_maf_threshold_' % nE + str(maf_threshold).replace('.', '_')))[index_ss + 1]] for index_shift_s0 in range(len(my_AA_directories_across_shifts))]
    heterozygosity_mean_wrt_shift_s0 = [np.mean(heterozygosity_each_generation) for heterozygosity_each_generation in heterozygosity_each_generation_wrt_shift_s0]
    heterozygosity_se_wrt_shift_s0 = [np.divide(np.std(heterozygosity_each_generation), np.sqrt(np.size(heterozygosity_each_generation))) for heterozygosity_each_generation in heterozygosity_each_generation_wrt_shift_s0]
    return heterozygosity_mean_wrt_shift_s0, heterozygosity_se_wrt_shift_s0

def calculate_contribution_to_change_mean_and_se_across_effect_sizes(my_AA_directories, indices_ss, E2Ns, nE):
    V2Ns = E2Ns ** 2
    S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
    a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]
    contribution_to_change_each_shift = dict()
    for index_ss in indices_ss:
        contribution_to_change_each_shift[index_ss] = list()
    for my_AA_directory in my_AA_directories:
        with open(os.path.join(my_AA_directory, 'fixations_and_extinctions'), 'rb') as f:
            d2ax_t_1_scaled = pickle.load(f)[-1]
        for d2ax_t_1_scaled_each_shift in d2ax_t_1_scaled:
            for index_ss in indices_ss:
                contribution_to_change_each_shift[index_ss].append(0.0)
            for mut in [mut for mut in d2ax_t_1_scaled_each_shift if np.searchsorted(a_list_pos, abs(mut[0])) - 1 in indices_ss]:
                contribution_to_change_each_shift[np.searchsorted(a_list_pos, abs(mut[0])) - 1][-1] += mut[1]
    contribution_to_change_mean = [np.mean(contribution_to_change_each_shift[index_ss]) for index_ss in indices_ss]
    contribution_to_change_se = [np.divide(np.std(contribution_to_change_each_shift[index_ss]), np.sqrt(np.size(contribution_to_change_each_shift[index_ss]))) for index_ss in indices_ss]
    return contribution_to_change_mean, contribution_to_change_se

def calculate_contribution_to_change_mean_and_se_across_shifts(my_AA_directories_across_shifts, index_ss, E2Ns, nE):
    V2Ns = E2Ns ** 2
    S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
    a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]
    contribution_to_change_mean_wrt_shift_s0 = list()
    contribution_to_change_se_wrt_shift_s0 = list()
    for index_shift_s0 in range(len(my_AA_directories_across_shifts)):
        contribution_to_change_each_shift = list()
        for my_AA_directory in my_AA_directories_across_shifts[index_shift_s0]:
            with open(os.path.join(my_AA_directory, 'fixations_and_extinctions'), 'rb') as f:
                d2ax_t_1_scaled = pickle.load(f)[-1]
            for d2ax_t_1_scaled_each_shift in d2ax_t_1_scaled:
                contribution_to_change_each_shift.append(0.0)
                for mut in [mut for mut in d2ax_t_1_scaled_each_shift if np.searchsorted(a_list_pos, abs(mut[0])) - 1 == index_ss]:
                    contribution_to_change_each_shift[-1] += mut[1]
        contribution_to_change_mean = np.mean(contribution_to_change_each_shift)
        contribution_to_change_se = np.divide(np.std(contribution_to_change_each_shift), np.sqrt(np.size(contribution_to_change_each_shift)))
        contribution_to_change_mean_wrt_shift_s0 += [contribution_to_change_mean]
        contribution_to_change_se_wrt_shift_s0 += [contribution_to_change_se]
    return contribution_to_change_mean_wrt_shift_s0, contribution_to_change_se_wrt_shift_s0

def read_pickle(file):
    with open(file, 'rb') as f:
        return pickle.load(f)
