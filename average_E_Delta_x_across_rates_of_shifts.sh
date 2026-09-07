#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=average_E_Delta_x_across_rates_of_shifts
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512mb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=1
index_shift_s0=$1
shift_s0=$2
rate_of_shifts_min=0.0001875
rate_of_shifts_max=0.768
number_rates_of_shifts=57
ss=$3

python3 ~/my_code/average_E_Delta_x_across_rates_of_shifts.py -N $N -U $U -E2Ns $E2Ns -D_s0 $shift_s0 --index_shift_s0 $index_shift_s0 --rate_of_shifts_min $rate_of_shifts_min --rate_of_shifts_max $rate_of_shifts_max --number_rates_of_shifts $number_rates_of_shifts -ss $ss -sD /insomnia001/depts/pas_lab/users/jc5473
