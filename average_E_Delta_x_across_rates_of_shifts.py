import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain the average expected change in allele frequency across rates of shifts')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the average expected change in allele frequency across rates of shifts')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-D_s0', '--shift_s0', type=float, help='The shift size in units of phenotypic standard deviation')
parser.add_argument('--index_shift_s0', type=int, help='The index of this shift size')
parser.add_argument('--rate_of_shifts_min', type=float, help='The min rate of shifts')
parser.add_argument('--rate_of_shifts_max', type=float, help='The max rate of shifts')
parser.add_argument('--number_rates_of_shifts', type=int, help='The number of rates of shifts')
parser.add_argument('-ss', type=float, help='The effect size squared')

# save the arguments passed into from the terminal
args = parser.parse_args()

if args.ss.is_integer():
    args.ss = int(args.ss)

with open(os.path.join(args.save_directory, 'average_E_Delta_x_across_rates_of_shifts_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_E2Ns_%d_ss_' % (args.index_shift_s0, round(args.E2Ns)) + str(args.ss).replace('.', '_')), 'wb') as f:
    pickle.dump(calculate_average_E_Delta_x_across_rates_of_shifts(N=args.N, U=args.U, E2Ns=args.E2Ns, index_shift_s0=args.index_shift_s0, shift_s0=args.shift_s0, rate_of_shifts_min=args.rate_of_shifts_min, rate_of_shifts_max=args.rate_of_shifts_max, rate_of_shifts_partitioning=args.number_rates_of_shifts, ss=args.ss), f, protocol=pickle.HIGHEST_PROTOCOL)
