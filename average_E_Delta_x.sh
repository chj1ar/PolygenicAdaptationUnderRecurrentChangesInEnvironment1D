#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --job-name=average_E_Delta_x
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
index_rate_of_shift=$3
rate_of_shift=$4
ss_min=0.1
ss_max=100

python3 ~/my_code/average_E_Delta_x.py -N $N -U $U -E2Ns $E2Ns -D_s0 $shift_s0 --index_shift_s0 $index_shift_s0 --rate_of_shift $rate_of_shift --index_rate_of_shift $index_rate_of_shift --ss_min $ss_min --ss_max $ss_max -sD /insomnia001/depts/pas_lab/users/jc5473
