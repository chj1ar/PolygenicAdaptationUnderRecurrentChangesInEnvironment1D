"""
The fixation probability where the size of each shift is drawn independently from a Gamma distribution, in the non-Lande case.
"""

from plot_functions import *

fontsize = 8

num_new_mutations = 10000000

E2Ns = 16
nE = 11
indices_ss = [0, 1, 3, 6, 9]
shift_s0_partitioning_nonLande = 29
rate_of_shifts_partitioning_nonLande = 66

analytic_V_A_nonLande = np.zeros((shift_s0_partitioning_nonLande, rate_of_shifts_partitioning_nonLande))
for index_shift_s0 in range(shift_s0_partitioning_nonLande):
    for index_rate_of_shift in range(rate_of_shifts_partitioning_nonLande):
        with open(os.path.join(save_directory, 'analytic_V_A_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
            analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift] = pickle.load(f)

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]
a_list_pos_between = [math.sqrt(a * a_next) for a, a_next in pairwise(a_list_pos)]

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

shift_s0_list = [0.5, 0.75, 1, 1.5, 2, 3]
rate_of_shifts_list = [0.0001875, 0.000375, 0.00075, 0.0015, 0.003, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384, 0.768]
shift_s0_partitioning = 29
rate_of_shifts_partitioning = 66

neutral, = ax.plot(np.geomspace(1, 100), [1 for _ in np.geomspace(1, 100)], 'k:')

no_shifts, = ax.plot(np.geomspace(1, 100), [pi(a=np.sqrt(ss), x=1.0 / (2 * N)) * 2.0 * N for ss in np.geomspace(1, 100)], 'k--')

# index_shift_s0, shift_s0 = 17, 1.5
for index_rate_of_shift, rate_of_shift, color in ((22, 0.003, 'r'), (44, 0.048, 'mistyrose'), (54, 0.192, 'lightgray')):
    fixation_probability_mean, fixation_probability_se = calculate_fixation_probability_mean_and_se_across_effect_sizes_from_num_fixed_and_num_new_mutations(my_AA_directories=[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shape_1_scale_1_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_analytic_V_A_%d' % (round(E2Ns), i)) for i in range(1, 1 + parallel_AA)], indices_ss=indices_ss, num_new_mutations=num_new_mutations, nE=nE)
    ax.errorbar(np.square(a_list_pos_between)[indices_ss], np.multiply(fixation_probability_mean, 2 * N), yerr=np.multiply(fixation_probability_se, 1.96 * 2 * N), fmt='o', color=color)
    # ax.plot(np.geomspace(1, 100), [pi_recurrent_shifts_nonlinear_linear_V_Delta_x_gamma_distributed_shift_sizes(a=np.sqrt(ss), x=1.0 / (2 * N), sigma_0_del=np.sqrt(analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) * 2.0 * N for ss in np.geomspace(1, 100)], color=color)
    with open(os.path.join(save_directory, 'analytic_fixation_probability_N_%d_U_' % N + str(U).replace('.', '_') + '_shape_1_scale_1_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d' % round(E2Ns)), 'rb') as f:
        ax.plot(np.geomspace(1, 100), pickle.load(f), color=color)

ax.text(x=100, y=1.05, s='Neutral', horizontalalignment='right', color='k', fontsize=fontsize)

simulations = ax.errorbar([], [], yerr=[], fmt='o', color='k')
analytics, = ax.plot([], [], 'k-')

ax.set_xscale('log')
ax.set_xlabel('Effect size squared', fontsize=fontsize)
ax.set_ylabel('Fixation probability (relative to neutral)', fontsize=fontsize)
ax.legend(handles=[no_shifts, analytics, simulations], labels=['No shifts', 'Analytics', 'All-allele simulations'], bbox_to_anchor=(0., 1.), loc='upper left', fontsize=fontsize)

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'fixation_probability_shift_size_gamma_nonLande.pdf'), bbox_inches='tight')
plt.show()
