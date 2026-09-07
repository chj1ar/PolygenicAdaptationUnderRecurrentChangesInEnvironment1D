#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=contribution_to_change
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=12:00:00

parallel_AA=49
N=2000
U=0.025
E2Ns=16
nE=11
shift_s0=$1
rate_of_shift=$2

python3 ~/my_code/contribution_to_change.py -N $N -U $U -E2Ns $E2Ns -nE $nE -D_s0 $shift_s0 --rate_of_shift $rate_of_shift --parallel_AA $parallel_AA --indices_ss 0 1 3 6 9 -sD /insomnia001/depts/pas_lab/users/jc5473
