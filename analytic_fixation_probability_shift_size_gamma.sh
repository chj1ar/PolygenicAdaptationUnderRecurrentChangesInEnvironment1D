#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=analytic_fixation_probability_shift_size_gamma
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512mb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=1
ss_min=0.1
ss_max=10
index_shift_s0=$1
index_rate_of_shift=$2
rate_of_shift=$3

python3 ~/my_code/analytic_fixation_probability_shift_size_gamma.py -N $N -U $U -E2Ns $E2Ns --index_shift_s0 $index_shift_s0 --index_rate_of_shift $index_rate_of_shift --rate_of_shift $rate_of_shift --ss_min $ss_min --ss_max $ss_max -sD /insomnia001/depts/pas_lab/users/jc5473
