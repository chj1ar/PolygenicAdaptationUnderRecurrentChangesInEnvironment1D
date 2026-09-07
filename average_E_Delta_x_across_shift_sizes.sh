#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=average_E_Delta_x_across_shift_sizes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512mb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=1
shift_s0_min=0.5
shift_s0_max=5
number_shift_s0=32
index_rate_of_shift=$1
rate_of_shift=$2
ss=$3

python3 ~/my_code/average_E_Delta_x_across_shift_sizes.py -N $N -U $U -E2Ns $E2Ns -D_s0_min $shift_s0_min -D_s0_max $shift_s0_max --number_shift_s0 $number_shift_s0 --rate_of_shift $rate_of_shift --index_rate_of_shift $index_rate_of_shift -ss $ss -sD /insomnia001/depts/pas_lab/users/jc5473
