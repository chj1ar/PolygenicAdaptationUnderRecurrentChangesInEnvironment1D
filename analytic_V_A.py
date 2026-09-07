import argparse
import os
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain analytic V_A by iteration')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the converged analytic V_A')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation. U=Lu')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-D_s0_min', '--shift_s0_min', type=float, help='The min shift size in units of phenotypic standard deviation')
parser.add_argument('-D_s0_max', '--shift_s0_max', type=float, help='The max shift size in units of phenotypic standard deviation')
parser.add_argument('--number_shift_s0', type=int, help='The number of shift sizes')
parser.add_argument('--index_shift_s0', type=int, help='The index of this shift size')
parser.add_argument('--rate_of_shifts_min', type=float, help='The min rate of shifts')
parser.add_argument('--rate_of_shifts_max', type=float, help='The max rate of shifts')
parser.add_argument('--number_rates_of_shifts', type=int, help='The number of rates of shifts')
parser.add_argument('--index_rate_of_shift', type=int, help='The index of this rate of shifts')

# save the arguments passed into from the terminal
args = parser.parse_args()

N = args.N
U = args.U
shift_s0 = np.geomspace(args.shift_s0_min, args.shift_s0_max, args.number_shift_s0, dtype=np.float32)[args.index_shift_s0]
rate_of_shift = np.geomspace(args.rate_of_shifts_min, args.rate_of_shifts_max, args.number_rates_of_shifts, dtype=np.float32)[args.index_rate_of_shift]
V_A_RHS = 2 * V_A(E2Ns=args.E2Ns)
V_A_LHS = V_A(E2Ns=args.E2Ns)
while abs(V_A_LHS / V_A_RHS - 1) > 0.01:
    V_A_RHS = V_A_LHS
    V_A_LHS = V_A_shifts(E2Ns=args.E2Ns, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A_RHS), p=rate_of_shift)
assert abs(V_A_LHS / V_A_RHS - 1) <= 0.01
with open(os.path.join(args.save_directory, 'analytic_V_A_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (args.index_shift_s0, args.index_rate_of_shift, round(args.E2Ns))), 'wb') as f:
    pickle.dump(V_A_LHS, f, protocol=pickle.HIGHEST_PROTOCOL)