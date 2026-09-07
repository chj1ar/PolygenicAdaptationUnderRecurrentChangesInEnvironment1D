#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=contribution_to_change_across_shift_sizes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=16
rate_of_shift=0.048
nE=11
index_ss=$1
parallel_AA=49

python3 ~/my_code/contribution_to_change_across_shift_sizes.py -N $N -U $U -E2Ns $E2Ns -D_s0_list 0.5 0.75 1 1.5 2 3 --rate_of_shift $rate_of_shift -nE $nE --index_ss $index_ss --parallel_AA $parallel_AA -sD /insomnia001/depts/pas_lab/users/jc5473
