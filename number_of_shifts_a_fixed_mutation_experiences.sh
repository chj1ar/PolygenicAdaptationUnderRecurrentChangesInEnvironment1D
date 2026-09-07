#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=number_of_shifts_a_fixed_mutation_experiences
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=1
index_shift_s0=$1
shift_s0=$2
index_rate_of_shift=$3
rate_of_shift=$4
index_ss=$5
ss=$6
nE=100

python3 ~/my_code/number_of_shifts_a_fixed_mutation_experiences.py -N $N -U $U -E2Ns $E2Ns -D_s0 $shift_s0 --index_shift_s0 $index_shift_s0 --rate_of_shift $rate_of_shift --index_rate_of_shift $index_rate_of_shift -ss $ss --index_ss $index_ss -nE $nE -sD /insomnia001/depts/pas_lab/users/jc5473
