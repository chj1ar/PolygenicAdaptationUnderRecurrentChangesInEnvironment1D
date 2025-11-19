"""
Simulate the AA simulation.
"""
from populations_functions_record_and_resume import SimulatePopulations
import argparse

# default parameters
default_params = dict()
default_params['Vs'] = -1
default_params['V2Ns'] = -1
default_params['nE'] = 6
default_params['nF'] = 11
default_params['bTN'] = 10
default_params['rate_of_shift'] = None
default_params['directories_for_empirical_V_A'] = None
default_params['sigma_0_del'] = None
default_params['particular_2Ns'] = []

parser = argparse.ArgumentParser(description="Simulate a population that experiences a single shift in environment.")
parser.add_argument('-sD', '--save_directory', type=str, help="The directory in which the simulation results are saved")
parser.add_argument('-N', type=int, help="The population size")
parser.add_argument('-U', type=float, help="The mutation rate per gamete per generation. U = Lu")
parser.add_argument('-Vs', type=float, default=default_params['Vs'], help="The fitness parameter")
parser.add_argument('-D_s0', '--shift_s0', type=float, help="\Lambda in units of sqrt(V_A(0))")
parser.add_argument('-E2Ns', type=float, help="E[a^2]. a^2 is Gamma distributed.")
parser.add_argument('-V2Ns', type=float, default=default_params['V2Ns'], help="V(a^2)")
parser.add_argument('-nE', type=int, default=default_params['nE'], help="The number of effect size bins. effect_size_bins_pos = np.geomspace(0.3, 7.0, nE)")
parser.add_argument('-nF', type=int, default=default_params['nF'], help="The number of frequency bins. freq_bins_minor_pos = np.linspace(0.0, 0.5, nF)[1:-1]")
parser.add_argument('-bTN', '--burn_time_N', type=int, default=default_params['bTN'], help="Burning time in units of N generations")
parser.add_argument('--rate_of_shift', type=float, default=default_params['rate_of_shift'], help="Rate of recurrent shifts. If None, then stabilizing selection.")
parser.add_argument('-nnM', '--number_new_mutations', type=int, help="Number of new mutations to empirically calculate the fixation probability")
parser.add_argument('--complete-jumps', dest='complete_jumps', action='store_true', default=False, help="No shifts during jumps if true")
parser.add_argument('--directories_for_empirical_V_A', type=str, default=default_params['directories_for_empirical_V_A'], help="The simulation result directories where the empirical average V_A is used in this simulation. We input only the directory with the last array index of the directories, but will use all the directories")
parser.add_argument('-sigma_0_del', type=float, default=default_params['sigma_0_del'], help="Square root of the phenotypic variance, should be the empirical average V_A from some previous iteration of this simulation")
parser.add_argument('--particular_2Ns', nargs='*', type=float, default=default_params['particular_2Ns'], help="Particular values of a^2 in addition to the effect size distribution")
args = parser.parse_args()

popSimulator = SimulatePopulations(N=args.N, U=args.U, Vs=args.Vs, shift_s0=args.shift_s0, E2Ns=args.E2Ns, V2Ns=args.V2Ns, nE=args.nE, nF=args.nF, burn_time_N=args.burn_time_N, rate_of_shift=args.rate_of_shift, directories_for_empirical_V_A=args.directories_for_empirical_V_A, sigma_0_del=args.sigma_0_del, nnm=args.number_new_mutations, particular_2Ns=args.particular_2Ns, save_directory=args.save_directory)

popSimulator.initiate_population()

if not args.complete_jumps:
    popSimulator.run_population_under_recurrent_shifts()
else:
    popSimulator.run_population_under_recurrent_shifts_complete_jumps()

popSimulator.save_stats()
