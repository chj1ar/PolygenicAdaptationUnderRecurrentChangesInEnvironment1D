#!/usr/bin/sh

#SBATCH --account=pas_lab
#SBATCH --partition=short
#SBATCH --job-name=analytic_heterozygosity_truncating_rare_alleles
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512mb
#SBATCH --time=12:00:00

N=2000
U=0.025
E2Ns=1
index_shift_s0=$1
index_rate_of_shift=$2

python3 ~/my_code/analytic_heterozygosity_truncating_rare_alleles.py -N $N -U $U -E2Ns $E2Ns -D_s0_min 0.5 -D_s0_max 5 --number_shift_s0 32 --index_shift_s0 $index_shift_s0 --rate_of_shifts_min 0.0001875 --rate_of_shifts_max 0.768 --number_rates_of_shifts 57 --index_rate_of_shift $index_rate_of_shift --maf_threshold 0.05 -sD /insomnia001/depts/pas_lab/users/jc5473
