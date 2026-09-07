import argparse
from plot_functions import *

parser = argparse.ArgumentParser(description='obtain analytic heterozygosity with sufficiently rare alleles truncated')
parser.add_argument('-sD', '--save_directory', type=str, help='The directory of the file to save the analytic heterozygosity')
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
parser.add_argument('--maf_threshold', type=float, help='The threshold of MAF below which alleles are truncated')

# save the arguments passed into from the terminal
args = parser.parse_args()

shift_s0 = np.geomspace(args.shift_s0_min, args.shift_s0_max, args.number_shift_s0, dtype=np.float32)[args.index_shift_s0]
rate_of_shift = np.geomspace(args.rate_of_shifts_min, args.rate_of_shifts_max, args.number_rates_of_shifts, dtype=np.float32)[args.index_rate_of_shift]

with open(os.path.join(args.save_directory, 'analytic_V_A_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (args.index_shift_s0, args.index_rate_of_shift, round(args.E2Ns))), 'rb') as f:
    analytic_V_A_index_shift_s0_index_rate_of_shift = pickle.load(f)

with open(os.path.join(args.save_directory, 'analytic_heterozygosity_N_%d_U_' % args.N + str(args.U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_maf_threshold_' % (args.index_shift_s0, args.index_rate_of_shift, round(args.E2Ns)) + str(args.maf_threshold).replace('.', '_')), 'wb') as f:
    pickle.dump(heterozygosity_shifts_truncating_rare_alleles_(a=np.sqrt(args.E2Ns), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A_index_shift_s0_index_rate_of_shift), p=rate_of_shift, maf_threshold=args.maf_threshold), f, protocol=pickle.HIGHEST_PROTOCOL)
