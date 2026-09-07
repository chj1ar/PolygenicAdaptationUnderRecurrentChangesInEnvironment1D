import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain analytic fixation probability where the size of each shift independently follows a (Gamma) distribution')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the analytic fixation probability')
parser.add_argument('-N', type=int, help='The population size')
parser.add_argument('-U', type=float, help='The mutation rate per gamete per generation. U=Lu')
parser.add_argument('-E2Ns', type=float, help='The expectation of the effect size distribution')
parser.add_argument('--index_shift_s0', type=int, help='The index of the shift size for the analytic V_A')
parser.add_argument('--index_rate_of_shift', type=int, help='The index of the rate of shifts')
parser.add_argument('--rate_of_shift', type=float, help='The rate of shifts')
parser.add_argument('--ss_min', type=float, help='The minimal effect size squared')
parser.add_argument('--ss_max', type=float, help='The maximal effect size squared')

# save the arguments passed into from the terminal
args = parser.parse_args()

with open(os.path.join(args.save_directory, 'analytic_V_A_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (args.index_shift_s0, args.index_rate_of_shift, round(args.E2Ns))), 'rb') as f:
    analytic_V_A = pickle.load(f)

with open(os.path.join(args.save_directory, 'analytic_fixation_probability_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_shape_1_scale_1_rate_of_shift_' + str(args.rate_of_shift).replace('.', '_') + '_E2Ns_%d' % round(args.E2Ns)), 'wb') as f:
    pickle.dump([pi_recurrent_shifts_nonlinear_linear_V_Delta_x_gamma_distributed_shift_sizes(a=np.sqrt(ss), x=1.0 / (2 * args.N), sigma_0_del=np.sqrt(analytic_V_A), p=args.rate_of_shift) * 2.0 * args.N for ss in np.geomspace(args.ss_min, args.ss_max)], f, protocol=pickle.HIGHEST_PROTOCOL)
