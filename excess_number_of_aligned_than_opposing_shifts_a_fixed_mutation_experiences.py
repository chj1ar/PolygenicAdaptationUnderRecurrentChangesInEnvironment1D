import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain the excess number of aligned than opposing shifts a fixed mutation experiences')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the excess number of aligned than opposing shifts a fixed mutation experiences')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-D_s0', '--shift_s0', type=float, help='The shift size in units of phenotypic standard deviation')
parser.add_argument('--index_shift_s0', type=int, help='The index of this shift size')
parser.add_argument('--rate_of_shift', type=float, help='The rate of shifts')
parser.add_argument('--index_rate_of_shift', type=int, help='The index of this rate of shifts')
parser.add_argument('-ss', type=float, help='The effect size squared')
parser.add_argument('--index_ss', type=int, help='The index of this effect size squared')
parser.add_argument('-nE', type=int, help='The number of effect size bins')

# save the arguments passed into from the terminal
args = parser.parse_args()

with open(os.path.join(args.save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (args.index_shift_s0, args.index_rate_of_shift, round(args.E2Ns), args.index_ss, args.nE)), 'wb') as f:
    pickle.dump(calculate_excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences(N=args.N, U=args.U, E2Ns=args.E2Ns, index_shift_s0=args.index_shift_s0, shift_s0=args.shift_s0, index_rate_of_shift=args.index_rate_of_shift, rate_of_shift=args.rate_of_shift, index_ss=args.index_ss, ss=args.ss, nE=args.nE), f, protocol=pickle.HIGHEST_PROTOCOL)
