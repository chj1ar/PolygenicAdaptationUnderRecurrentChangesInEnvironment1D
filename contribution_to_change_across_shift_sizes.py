import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain contribution to change across shift sizes')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the contribution to change across shift sizes')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation. U=Lu')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('-D_s0_list', '--shift_s0_list', nargs='*', type=float, help='The list of shift sizes')
parser.add_argument('--rate_of_shift', type=float, help='The rate of shifts')
parser.add_argument('-nE', type=int, help='The number of effect size bins')
parser.add_argument('--index_ss', type=int, help='The index of the effect size squared')
parser.add_argument('--parallel_AA', type=int, help='The number of AA simulation replicates')

# save the arguments passed into from the terminal
args = parser.parse_args()

args.shift_s0_list = [round(shift_s0) if shift_s0 == round(shift_s0) else shift_s0 for shift_s0 in args.shift_s0_list]

contribution_to_change_mean_wrt_shift_s0, contribution_to_change_se_wrt_shift_s0 = calculate_contribution_to_change_mean_and_se_across_shifts(my_AA_directories_across_shifts=[[os.path.join(args.save_directory, 'my_AA_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(args.E2Ns), j)) for j in range(1, 1 + args.parallel_AA)] for shift_s0 in args.shift_s0_list], index_ss=args.index_ss, E2Ns=args.E2Ns, nE=args.nE)

with open(os.path.join(args.save_directory, 'contribution_to_change_mean_wrt_shift_s0_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d_nE_%d_index_ss_%d' % (round(args.E2Ns), args.nE, args.index_ss)), 'wb') as f:
    pickle.dump(contribution_to_change_mean_wrt_shift_s0, f, protocol=pickle.HIGHEST_PROTOCOL)

with open(os.path.join(args.save_directory, 'contribution_to_change_se_wrt_shift_s0_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d_nE_%d_index_ss_%d' % (round(args.E2Ns), args.nE, args.index_ss)), 'wb') as f:
    pickle.dump(contribution_to_change_se_wrt_shift_s0, f, protocol=pickle.HIGHEST_PROTOCOL)
