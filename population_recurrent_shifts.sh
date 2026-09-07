#!/usr/bin/sh
#SBATCH --job-name=populations
#SBATCH --account=pas_lab
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --time=12:00:00
#SBATCH --array=1-49

N=2000
U=0.025
E2Ns=16
shape=1
scale=1
rate_of_shift=0.768
python3 ~/my_code/populations_argparse_recurrent_shifts_no_efs_bins.py -N $N -U $U -shape $shape -scale $scale --rate_of_shift $rate_of_shift -E2Ns $E2Ns -bTN 10 -nnM 10000000 -sD /insomnia001/depts/pas_lab/users/jc5473/my_AA_N_${N}_U_${U/./_}_shape_${shape/./_}_scale_${scale/./_}_rate_of_shift_${rate_of_shift/./_}_E2Ns_${E2Ns}_no_efs_bins_d2ax_t_1_scaled_analytic_V_A_$SLURM_ARRAY_TASK_ID
