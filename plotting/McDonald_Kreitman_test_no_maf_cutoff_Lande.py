"""
Panels A, B, and C plot the simulation results and the analytic approximations of the two qualitatively different behaviors of the proportion of fixations due to adaptation, where shifts have some effect, alongside the corresponding fixation probability and the corresponding heterozygosity. Panels D and E plot the simulation results and the analytic approximations of the proportion of fixations due to adaptation as a function of the shift size or the rate of shifts, where the other is fixed.
The simulation results are saved with the prefix 'my_AA_'. The analytic approximations are the diffusion approximations detailed in the manuscript, and are calculated, by definition, from the fixation probability and the heterozygosity.
"""

from plot_functions import *

fontsize = 8

num_new_mutations = 10000000

E2Ns = 1
nE = 100
indices_ss = [21, 38, 62, 86, 98]
shift_s0_partitioning = 32
rate_of_shifts_partitioning = 57

analytic_V_A = np.zeros((shift_s0_partitioning, rate_of_shifts_partitioning))
for index_shift_s0 in range(shift_s0_partitioning):
    for index_rate_of_shift in range(rate_of_shifts_partitioning):
        # each individual analytic V_A is obtained by iteratively applying V_A_shifts() onto the value from the previous iteration until convergence; the initial value is the analytic V_A under no shifts
        with open(os.path.join(save_directory, 'analytic_V_A_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
            analytic_V_A[index_shift_s0][index_rate_of_shift] = pickle.load(f)

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]
a_list_pos_between = [math.sqrt(a * a_next) for a, a_next in pairwise(a_list_pos)]

fig = plt.figure(figsize=(10, 6.5), dpi=300)
gs = fig.add_gridspec(2, 5, width_ratios=[25, 3, 25, 3, 25], height_ratios=[1, 1])
# just for spacing
plt.subplot(gs[0, 1]).axis('off')
plt.subplot(gs[1, 1]).axis('off')
plt.subplot(gs[0, 3]).axis('off')
plt.subplot(gs[1, 3]).axis('off')
# not used
plt.subplot(gs[1, 4]).axis('off')

shift_s0_list = [0.5, 0.75, 1, 1.5, 2, 3, 5]
rate_of_shifts_list = [0.0001875, 0.000375, 0.00075, 0.0015, 0.003, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384, 0.768]
shift_s0_partitioning = 32
rate_of_shifts_partitioning = 57

# make the other panels
# the fixation probability
ax = plt.subplot(gs[0, 0])
ax.text(0, 1.01, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold')

neutral, = ax.plot(np.geomspace(0.1, 10), [1 for _ in np.geomspace(0.1, 10)], 'k:')

no_shifts, = ax.plot(np.geomspace(0.1, 10), [pi(a=np.sqrt(ss), x=1.0 / (2 * N)) * 2.0 * N for ss in np.geomspace(0.1, 10)], 'k--')

ax.text(x=10, y=1.05, s='Neutral', horizontalalignment='right', color='k', fontsize=fontsize)

simulations = ax.errorbar([], [], yerr=[], fmt='o', color='k')
analytics, = ax.plot([], [], 'k-')

ax.set_xscale('log')
ax.set_xlabel('Effect size squared', fontsize=fontsize)
ax.set_ylabel('Fixation probability (relative to neutral)', fontsize=fontsize)
ax.legend(handles=[no_shifts, analytics, simulations], labels=['No shifts', 'Analytics', 'All-allele simulations'], bbox_to_anchor=(0., 1.), loc='upper left', fontsize=fontsize)

# the heterozygosity
ax = plt.subplot(gs[0, 2])
ax.text(0, 1.01, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold')

neutral, = ax.plot(np.geomspace(0.1, 10), [1 for _ in np.geomspace(0.1, 10)], 'k:')

no_shifts, = ax.plot(np.geomspace(0.1, 10), [heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) for ss in np.geomspace(0.1, 10)], 'k--')

ax.text(x=10, y=1.005, s='Neutral', horizontalalignment='right', color='k', fontsize=fontsize)

simulations = ax.errorbar([], [], yerr=[], fmt='o', color='k')
analytics, = ax.plot([], [], 'k-')

ax.set_xscale('log')
ax.set_xlabel('Effect size squared', fontsize=fontsize)
ax.set_ylabel('Heterozygosity (relative to neutral)', fontsize=fontsize)
ax.legend(handles=[no_shifts, analytics, simulations], labels=['No shifts', 'Analytics', 'All-allele simulations'], bbox_to_anchor=(0., 0.), loc='lower left', fontsize=fontsize)

# the McDonald-Kreitman test
ax = plt.subplot(gs[0, 4])
ax.text(0, 1.01, 'C', transform=ax.transAxes, fontsize=12, fontweight='bold')

neutral, = ax.plot(np.geomspace(0.1, 10), [0 for _ in np.geomspace(0.1, 10)], 'k:')

no_shifts, = ax.plot(np.geomspace(0.1, 10), [1 - heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) / pi(a=np.sqrt(ss), x=1 / (2 * N)) / (2 * N) for ss in np.geomspace(0.1, 10)], 'k--')

simulations = ax.errorbar([], [], yerr=[], fmt='o', color='k')
analytics, = ax.plot([], [], 'k-')

ax.set_xscale('log')
ax.set_xlabel('Effect size squared', fontsize=fontsize)
ax.set_ylabel(r'$\alpha := 1 - \frac{P_n}{P_s} / \frac{D_n}{D_s}$', fontsize=fontsize)
ax.legend(handles=[no_shifts, analytics, simulations], labels=['No shifts', 'Analytics', 'All-allele simulations'], bbox_to_anchor=(0., 0.), loc='lower left', fontsize=fontsize)

index_shift_s0, shift_s0 = 15, 1.5
for index_rate_of_shift, rate_of_shift, color in ((19, 0.003, 'g'), (37, 0.048, 'lightgreen')):
    # simulations
    ## obtain the fixation probability of alleles with each a^2
    fixation_probability_mean, fixation_probability_se = calculate_fixation_probability_mean_and_se_across_effect_sizes_from_num_fixed_and_num_new_mutations(my_AA_directories=[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), i)) for i in range(1, 1 + parallel_AA)], indices_ss=indices_ss, num_new_mutations=num_new_mutations, nE=nE)
    ax = plt.subplot(gs[0, 0])
    ax.errorbar(np.square(a_list_pos_between)[indices_ss], np.multiply(fixation_probability_mean, 2 * N), yerr=np.multiply(fixation_probability_se, 1.96 * 2 * N), fmt='o', color=color)
    ## obtain the heterozygosity of alleles with each a^2 in each generation, which are separated by 10N generations
    heterozygosity_mean, heterozygosity_se = calculate_heterozygosity_mean_and_se_across_effect_sizes(my_AA_directories=[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), i)) for i in range(1, 1 + parallel_AA)], indices_ss=indices_ss, nE=nE)
    ax = plt.subplot(gs[0, 2])
    ax.errorbar(np.square(a_list_pos_between)[indices_ss], np.divide(heterozygosity_mean, (2 * N * U / (nE + 1)) * heterozygosity_(a=0.1)), yerr=1.96 * np.divide(heterozygosity_se, (2 * N * U / (nE + 1)) * heterozygosity_(a=0.1)), fmt='o', color=color)
    ## plot the McDonald-Kreitman test
    ax = plt.subplot(gs[0, 4])
    ax.errorbar(np.square(a_list_pos_between)[indices_ss], np.subtract(1, np.divide(np.divide(heterozygosity_mean, (2 * N * U / (nE + 1)) * heterozygosity_(a=0.1) * 2 * N), fixation_probability_mean)), yerr=1.96 * np.divide(np.divide(heterozygosity_se, (2 * N * U / (nE + 1)) * heterozygosity_(a=0.1) * 2 * N), fixation_probability_mean), fmt='o', color=color)
    # analytics
    fixation_probability_relative_to_neutral = [pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a=np.sqrt(ss), x=1.0 / (2 * N), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) * 2.0 * N for ss in np.geomspace(0.1, 10)]
    heterozygosity_relative_to_neutral = [heterozygosity_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / heterozygosity_(a=0.1) for ss in np.geomspace(0.1, 10)]
    ## the fixation probability
    ax = plt.subplot(gs[0, 0])
    ax.plot(np.geomspace(0.1, 10), fixation_probability_relative_to_neutral, color=color)
    ## the heterozygosity
    ax = plt.subplot(gs[0, 2])
    ax.plot(np.geomspace(0.1, 10), heterozygosity_relative_to_neutral, color=color)
    ## the McDonald-Kreitman test
    ax = plt.subplot(gs[0, 4])
    ax.plot(np.geomspace(0.1, 10), np.subtract(1, np.divide(heterozygosity_relative_to_neutral, fixation_probability_relative_to_neutral)), color=color)

# Graph II for the McDonald-Kreitman test (w.r.t. the shift size or the rate of shifts; different effect sizes)

ax = plt.subplot(gs[1, 0])
ax.text(0, 1.01, 'D', transform=ax.transAxes, fontsize=12, fontweight='bold')

index_rate_of_shift, rate_of_shift = 28, 0.012

neutral, = ax.plot([0] + shift_s0_list, [0 for _ in [0] + shift_s0_list], 'k:')

simulations = [None, None, None]
no_shifts = [None, None, None]
analytics = [None, None, None]

for i, (index_ss, ss, color) in enumerate([[indices_ss[0], 0.25, 'b'], [indices_ss[2], 1, 'purple'], [indices_ss[4], 4.25, 'r']]):
    fixation_probability_mean_wrt_shift_s0, fixation_probability_se_wrt_shift_s0 = calculate_fixation_probability_mean_and_se_across_shifts_from_num_fixed_and_num_new_mutations(my_AA_directories_across_shifts=[[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), i)) for i in range(1, 1 + parallel_AA)] for shift_s0 in shift_s0_list], index_ss=index_ss, num_new_mutations=num_new_mutations, nE=nE)
    heterozygosity_mean_wrt_shift_s0, heterozygosity_se_wrt_shift_s0 = calculate_heterozygosity_mean_and_se_across_shifts(my_AA_directories_across_shifts=[[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), j)) for j in range(1, 1 + parallel_AA)] for shift_s0 in shift_s0_list], index_ss=index_ss, nE=nE)
    simulations[i] = ax.errorbar(shift_s0_list, np.subtract(1, np.divide(np.divide(heterozygosity_mean_wrt_shift_s0, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1) * 2 * N), fixation_probability_mean_wrt_shift_s0)), yerr=1.96 * np.divide(np.divide(heterozygosity_se_wrt_shift_s0, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1) * 2 * N), fixation_probability_mean_wrt_shift_s0), fmt='o', color=color)
    no_shifts[i], = ax.plot([0] + shift_s0_list, [1 - heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) / pi(a=np.sqrt(ss), x=1 / (2 * N)) / (2 * N) for _ in [0] + shift_s0_list], color=color, linestyle='--')
    analytics[i], = ax.plot([0] + list(np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning)), [1 - heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) / pi(a=np.sqrt(ss), x=1 / (2 * N)) / (2 * N)] + [1 - heterozygosity_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / heterozygosity_(a=0.1) / pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a=np.sqrt(ss), x=1 / (2 * N), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / (2 * N) for index_shift_s0, shift_s0 in enumerate(np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning))], color=color)

# ax.text(x=shift_s0_list[0], y=1.005, s='Neutral', color='k', fontsize=fontsize)
ax.text(x=3.3, y=0.135, s='$a^2=%g$' % 0.25, color='b', rotation=30, rotation_mode='anchor', fontsize=fontsize)
ax.text(x=3.3, y=0.37, s='$a^2=%g$' % 1, color='purple', rotation=45, rotation_mode='anchor', fontsize=fontsize)
ax.text(x=3.3, y=0.67, s='$a^2=%g$' % 4, color='r', rotation=45, rotation_mode='anchor', fontsize=fontsize)

simulation_mean = ax.errorbar([], [], yerr=[], fmt='o', color='k', label='All-allele simulations')
no_shifts_mean, = ax.plot([], [], 'k--', label='No shifts')
analytic_mean, = ax.plot([], [], 'k-', label='Analytics')

ax.set_xscale('symlog', linthresh=shift_s0_list[0], linscale=0.15)
ax.set_xlabel(r'Shift size ($\Lambda / \sqrt{V_A}$)', fontsize=fontsize)
ax.set_ylabel(r'$\alpha := 1 - \frac{P_n}{P_s} / \frac{D_n}{D_s}$', fontsize=fontsize)
ax.set_xticks(ticks=[0] + shift_s0_list, labels=[0] + [shift_s0 if index_shift_s0 % 2 == 0 else None for index_shift_s0, shift_s0 in enumerate(shift_s0_list)], fontsize=fontsize)
# ax.set_title('$E[a^2]=%d, 2Np=%g$' % (round(E2Ns), 2 * N * rate_of_shift), fontsize=fontsize)
ax.legend(bbox_to_anchor=(0., 1.), loc='upper left', fontsize=fontsize)
ax.minorticks_off()

ax = plt.subplot(gs[1, 2])
ax.text(0, 1.01, 'E', transform=ax.transAxes, fontsize=12, fontweight='bold')

index_shift_s0, shift_s0 = 15, 1.5

neutral, = ax.plot([0] + list(2 * N * np.array(rate_of_shifts_list)), [0 for _ in [0] + rate_of_shifts_list], 'k:')

simulations = [None, None, None]
no_shifts = [None, None, None]
analytics = [None, None, None]

for i, (index_ss, ss, color) in enumerate([[indices_ss[0], 0.25, 'b'], [indices_ss[2], 1, 'purple'], [indices_ss[4], 4.25, 'r']]):
    fixation_probability_mean_wrt_rate_of_shifts, fixation_probability_se_wrt_rate_of_shifts = calculate_fixation_probability_mean_and_se_across_shifts_from_num_fixed_and_num_new_mutations(my_AA_directories_across_shifts=[[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), i)) for i in range(1, 1 + parallel_AA)] for rate_of_shift in rate_of_shifts_list], index_ss=index_ss, num_new_mutations=num_new_mutations, nE=nE)
    heterozygosity_mean_wrt_rate_of_shifts, heterozygosity_se_wrt_rate_of_shifts = calculate_heterozygosity_mean_and_se_across_shifts(my_AA_directories_across_shifts=[[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), j)) for j in range(1, 1 + parallel_AA)] for rate_of_shift in rate_of_shifts_list], index_ss=index_ss, nE=nE)
    simulations[i] = ax.errorbar(2 * N * np.array(rate_of_shifts_list), np.subtract(1, np.divide(np.divide(heterozygosity_mean_wrt_rate_of_shifts, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1) * 2 * N), fixation_probability_mean_wrt_rate_of_shifts)), yerr=1.96 * np.divide(np.divide(heterozygosity_se_wrt_rate_of_shifts, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1) * 2 * N), fixation_probability_mean_wrt_rate_of_shifts), fmt='o', color=color)
    no_shifts[i], = ax.plot([0] + list(2 * N * np.array(rate_of_shifts_list)), [1 - heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) / pi(a=np.sqrt(ss), x=1 / (2 * N)) / (2 * N) for _ in [0] + rate_of_shifts_list], color=color, linestyle='--')
    analytics[i], = ax.plot([0] + list(2 * N * np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning)), [1 - heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) / pi(a=np.sqrt(ss), x=1 / (2 * N)) / (2 * N)] + [1 - heterozygosity_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / heterozygosity_(a=0.1) / pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a=np.sqrt(ss), x=1 / (2 * N), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / (2 * N) for index_rate_of_shift, rate_of_shift in enumerate(np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning))], color=color)

ax.set_xscale('symlog', linthresh=2 * N * rate_of_shifts_list[0], linscale=0.6)
ax.set_xlabel(r'Rate of shifts ($2Np$)', fontsize=fontsize)
ax.set_xticks(ticks=[0] + list(2 * N * np.array(rate_of_shifts_list)), labels=[0] + [2 * N * rate_of_shift if index_rate_of_shift % 2 == 0 and 2 * N * rate_of_shift != int(2 * N * rate_of_shift) else int(2 * N * rate_of_shift) if index_rate_of_shift % 2 == 0 else None for index_rate_of_shift, rate_of_shift in enumerate(rate_of_shifts_list)], fontsize=fontsize)
# ax.set_title('$E[a^2]=%d, \Lambda=%g\sqrt{V_A}$' % (round(E2Ns), shift_s0), fontsize=fontsize)
ax.minorticks_off()

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'McDonald_Kreitman_test_no_maf_cutoff_Lande.pdf'), bbox_inches='tight')
plt.show()
