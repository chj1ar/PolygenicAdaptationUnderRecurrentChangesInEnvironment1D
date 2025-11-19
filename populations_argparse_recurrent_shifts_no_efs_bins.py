"""
Simulate the AA simulations without effect size bins
"""
from populations_functions_record_and_resume import SimulatePopulations
import argparse

def none_or_str(argument):
    if argument == 'None':
        return None
    return argument

# default parameters
default_parameters = dict()
default_parameters['Vs'] = -1
default_parameters['V2Ns'] = -1
default_parameters['shape'] = None
default_parameters['scale'] = None
default_parameters['shift_s0'] = None
default_parameters['sigma_0_del'] = None

parser = argparse.ArgumentParser(description="Simulate a population without effect size bins")
parser.add_argument('-sD', '--save_directory', type=str, help="The directory in which the simulation results are saved")
parser.add_argument('-N', type=int, help="The population size")
parser.add_argument('-U', type=float, help="The mutation rate per gamete per generation. U = Lu")
parser.add_argument('-Vs', type=float, default=default_parameters['Vs'], help="The fitness parameter")
parser.add_argument('-E2Ns', type=float, help="E[a^2]")
parser.add_argument('-V2Ns', type=float, default=default_parameters['V2Ns'], help="V(a^2)")
parser.add_argument('-bTN', '--burn_time_N', type=int, help="Burn time in units of N generations")
parser.add_argument('-shape', type=float, default=default_parameters['shape'], help="The shape parameter of the gamma distributed shift sizes in units of the phenotypic standard deviation")
parser.add_argument('-scale', type=float, default=default_parameters['scale'], help="The scale parameter of the gamma distributed shift sizes in units of the phenotypic standard deviation")
parser.add_argument('-D_s0', '--shift_s0', type=float, default=default_parameters['shift_s0'], help="The shift size in units of the phenotypic standard deviation")
parser.add_argument('--rate_of_shift', type=float, help="The rate of shifts")
parser.add_argument('-nnM', '--number_new_mutations', type=int, help="Number of new mutations to measure the fixation probability")
parser.add_argument('--directories_for_empirical_V_A', type=none_or_str, help="The simulation result directories where the empirical V_A is used in this simulation. We input only the directory with the last array index of the directories, but will use all the directories")
parser.add_argument('-sigma_0_del', type=float, default=default_parameters['sigma_0_del'], help="The phenotypic standard deviation")
args = parser.parse_args()

popSimulator = SimulatePopulations(N=args.N, U=args.U, Vs=args.Vs, E2Ns=args.E2Ns, V2Ns=args.V2Ns, nE=None, nF=None, burn_time_N=args.burn_time_N, shape=args.shape, scale=args.scale, shift_s0=args.shift_s0, rate_of_shift=args.rate_of_shift, directories_for_empirical_V_A=args.directories_for_empirical_V_A, sigma_0_del=args.sigma_0_del, nnm=args.number_new_mutations, particular_2Ns=list(), save_directory=args.save_directory)

popSimulator.initiate_population(no_efs_bins=True)

popSimulator.run_population_under_recurrent_shifts_no_efs_bins()

popSimulator.save_stats_no_efs_bins()
