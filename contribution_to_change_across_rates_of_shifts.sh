#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=contribution_to_change_across_rates_of_shifts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=16
shift_s0=1.5
nE=11
index_ss=$1
parallel_AA=49

python3 ~/my_code/contribution_to_change_across_rates_of_shifts.py -N $N -U $U -E2Ns $E2Ns -D_s0 $shift_s0 --rates_of_shifts_list 0.0001875 0.000375 0.00075 0.0015 0.003 0.006 0.012 0.024 0.048 0.096 0.192 0.384 0.768 -nE $nE --index_ss $index_ss --parallel_AA $parallel_AA -sD /insomnia001/depts/pas_lab/users/jc5473
