import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='record contribution to change for representative parameter combinations')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the contribution to change')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation. U=Lu')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-nE', type=int, help='The number of effect size bins')
parser.add_argument('-D_s0', '--shift_s0', type=float, help='The shift size in units of phenotypic standard deviation')
parser.add_argument('--rate_of_shift', type=float, help='The rate of shifts')
parser.add_argument('--parallel_AA', type=int, help='The number of AA simulation results in parallel')
parser.add_argument('--indices_ss', nargs='*', type=int, help='The indices in the effect size squared bin')

# save the arguments passed into from the terminal
args = parser.parse_args()

contribution_to_change_mean, contribution_to_change_se = calculate_contribution_to_change_mean_and_se_across_effect_sizes(my_AA_directories=[os.path.join(args.save_directory, 'my_AA_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_shift_s0_' + str(args.shift_s0).replace('.', '_') + '_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(args.E2Ns), i)) for i in range(1, 1 + args.parallel_AA)], indices_ss=args.indices_ss, E2Ns=args.E2Ns, nE=args.nE)

with open(os.path.join(args.save_directory, 'contribution_to_change_mean_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_shift_s0_' + str(args.shift_s0).replace('.', '_') + '_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d' % round(args.E2Ns)), 'wb') as f:
    pickle.dump(contribution_to_change_mean, f, protocol=pickle.HIGHEST_PROTOCOL)

with open(os.path.join(args.save_directory, 'contribution_to_change_se_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_shift_s0_' + str(args.shift_s0).replace('.', '_') + '_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d' % round(args.E2Ns)), 'wb') as f:
    pickle.dump(contribution_to_change_se, f, protocol=pickle.HIGHEST_PROTOCOL)
