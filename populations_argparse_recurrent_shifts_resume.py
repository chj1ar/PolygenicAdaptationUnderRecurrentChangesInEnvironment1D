"""
Resume an AA simulation from the intermediate state recorded
"""
from populations_functions_record_and_resume import SimulatePopulations
import argparse

parser = argparse.ArgumentParser(description="Resume an AA simulation from the intermediate folder.")
parser.add_argument('-iF', '--intermediate_Folder', type=str, help="The folder in which the intermediate state is recorded")
args = parser.parse_args()

popSimulator = SimulatePopulations(intermediate_Folder=args.intermediate_Folder)

popSimulator.initiate_population()

popSimulator.run_population_under_recurrent_shifts_resume()

popSimulator.save_stats()
