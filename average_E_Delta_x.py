import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain the average expected change in allele frequency')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the average expected change in allele frequency')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation. U=Lu')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-D_s0', '--shift_s0', type=float, help='The shift size in units of phenotypic standard deviation')
parser.add_argument('--index_shift_s0', type=int, help='The index of this shift size')
parser.add_argument('--rate_of_shift', type=float, help='The rate of shifts')
parser.add_argument('--index_rate_of_shift', type=int, help='The index of this rate of shifts')
parser.add_argument('--ss_min', type=float, help='The min effect size squared')
parser.add_argument('--ss_max', type=float, help='The max effect size squared')

# save the arguments passed into from the terminal
args = parser.parse_args()

with open(os.path.join(args.save_directory, 'average_E_Delta_x_across_effect_sizes_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (args.index_shift_s0, args.index_rate_of_shift, round(args.E2Ns))), 'wb') as f:
    pickle.dump(calculate_average_E_Delta_x_across_effect_sizes(N=args.N, U=args.U, E2Ns=args.E2Ns, index_shift_s0=args.index_shift_s0, shift_s0=args.shift_s0, index_rate_of_shift=args.index_rate_of_shift, rate_of_shift=args.rate_of_shift, ss_min=args.ss_min, ss_max=args.ss_max), f, protocol=pickle.HIGHEST_PROTOCOL)
