import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain the average expected change in allele frequency across shift sizes')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the average expected change in allele frequency across shift sizes')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-D_s0_min', '--shift_s0_min', type=float, help='The min shift size in units of phenotypic standard deviation')
parser.add_argument('-D_s0_max', '--shift_s0_max', type=float, help='The max shift size in units of phenotypic standard deviation')
parser.add_argument('--number_shift_s0', type=int, help='The number of shift sizes')
parser.add_argument('--rate_of_shift', type=float, help='The rate of shifts')
parser.add_argument('--index_rate_of_shift', type=int, help='The index of this rate of shifts')
parser.add_argument('-ss', type=float, help='The effect size squared')

# save the arguments passed into from the terminal
args = parser.parse_args()

if args.ss.is_integer():
    args.ss = int(args.ss)

with open(os.path.join(args.save_directory, 'average_E_Delta_x_across_shift_sizes_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_rate_of_shift_%d_E2Ns_%d_ss_' % (args.index_rate_of_shift, round(args.E2Ns)) + str(args.ss).replace('.', '_')), 'wb') as f:
    pickle.dump(calculate_average_E_Delta_x_across_shift_sizes(N=args.N, U=args.U, E2Ns=args.E2Ns, shift_s0_min=args.shift_s0_min, shift_s0_max=args.shift_s0_max, shift_s0_partitioning=args.number_shift_s0, index_rate_of_shift=args.index_rate_of_shift, rate_of_shift=args.rate_of_shift, ss=args.ss), f, protocol=pickle.HIGHEST_PROTOCOL)
