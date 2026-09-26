"""
The heterozygosity in the non-Lande case.
"""

from plot_functions import *

fontsize = 8

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

fig = plt.figure(figsize=(8, 6.5), dpi=300)
gs = fig.add_gridspec(2, 5, width_ratios=[25, 3, 25, 1, 3], height_ratios=[1, 1])
# just for spacing
plt.subplot(gs[0, 1]).axis('off')
plt.subplot(gs[1, 1]).axis('off')
plt.subplot(gs[0, 3]).axis('off')
plt.subplot(gs[1, 3]).axis('off')
# not used
plt.subplot(gs[1, 4]).axis('off')
# make gs[0, 4] a colorbar
shift_s0_list = [0.5, 0.75, 1, 1.5, 2, 3]
rate_of_shifts_list = [0.0001875, 0.000375, 0.00075, 0.0015, 0.003, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384, 0.768]
shift_s0_partitioning = 29
rate_of_shifts_partitioning = 66

x = 2 * N * np.geomspace(rate_of_shifts_list[0] / (rate_of_shifts_list[-1] / rate_of_shifts_list[0]) ** (0.5 / (rate_of_shifts_partitioning - 1)), rate_of_shifts_list[-1] * (rate_of_shifts_list[-1] / rate_of_shifts_list[0]) ** (0.5 / (rate_of_shifts_partitioning - 1)), rate_of_shifts_partitioning + 1)
y = np.geomspace(shift_s0_list[0] / (shift_s0_list[-1] / shift_s0_list[0]) ** (0.5 / (shift_s0_partitioning - 1)), shift_s0_list[-1] * (shift_s0_list[-1] / shift_s0_list[0]) ** (0.5 / (shift_s0_partitioning - 1)), shift_s0_partitioning + 1)
X, Y = np.meshgrid(x, y)

Z = np.empty((shift_s0_partitioning, rate_of_shifts_partitioning))
for index_shift_s0 in range(shift_s0_partitioning):
    for index_rate_of_shift in range(rate_of_shifts_partitioning):
        with open(os.path.join(save_directory, 'analytic_heterozygosity_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
            Z[index_shift_s0, index_rate_of_shift] = pickle.load(f)
Z = np.divide(Z, heterozygosity_(a=np.sqrt(E2Ns)))
z_min, z_max = Z.min(), Z.max()

norm = colors.Normalize(vmin=z_min, vmax=z_max)
sm = plt.cm.ScalarMappable(cmap='Greens_r', norm=norm)
cbar = plt.colorbar(sm, cax=plt.subplot(gs[0, 4]), orientation='vertical')
cbar.set_label('Fold increase in heterozygosity', fontsize=fontsize)
cbar.ax.set_yticks([1.01, z_max])
cbar.ax.set_yticklabels([1.01, '%.2f' % z_max])
cbar.ax.axhline(y=1.01, color='lightgray', linestyle=':')
cbar.ax.axhline(y=1.1, color='lightgray', linestyle='-')
cbar.ax.minorticks_off()
# make the other panels
ax = plt.subplot(gs[0, 0])
ax.text(0, 1.01, 'A', transform=ax.transAxes, fontsize=12, fontweight='bold')

neutral, = ax.plot(np.geomspace(1, 100), [1 for _ in np.geomspace(1, 100)], 'k:')

no_shifts, = ax.plot(np.geomspace(1, 100), [heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) for ss in np.geomspace(1, 100)], 'k--')

index_shift_s0, shift_s0 = 17, 1.5
for index_rate_of_shift, rate_of_shift, color in ((22, 0.003, 'g'), (44, 0.048, 'limegreen')):
    heterozygosity_mean, heterozygosity_se = calculate_heterozygosity_mean_and_se_across_effect_sizes(my_AA_directories=[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), i)) for i in range(1, 1 + parallel_AA)], indices_ss=indices_ss, nE=nE)
    ax.errorbar(np.square(a_list_pos_between)[indices_ss], np.divide(heterozygosity_mean, (2 * N * U / (nE + 1)) * heterozygosity_(a=0.1)), yerr=1.96 * np.divide(heterozygosity_se, (2 * N * U / (nE + 1)) * heterozygosity_(a=0.1)), fmt='o', color=color)
    ax.plot(np.geomspace(1, 100), [heterozygosity_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / heterozygosity_(a=0.1) for ss in np.geomspace(1, 100)], color=color)

ax.text(x=100, y=1.005, s='Neutral', horizontalalignment='right', color='k', fontsize=fontsize)

simulations = ax.errorbar([], [], yerr=[], fmt='o', color='k')
analytics, = ax.plot([], [], 'k-')

ax.set_xscale('log')
ax.set_xlabel('Effect size squared', fontsize=fontsize)
ax.set_ylabel('Heterozygosity (relative to neutral)', fontsize=fontsize)
ax.legend(handles=[no_shifts, analytics, simulations], labels=['No shifts', 'Analytics', 'All-allele simulations'], bbox_to_anchor=(1., 0.95), loc='upper right', fontsize=fontsize)

ax = plt.subplot(gs[0, 2])
ax.text(0, 1.01, 'B', transform=ax.transAxes, fontsize=12, fontweight='bold')

norm = colors.Normalize(vmin=z_min, vmax=z_max)
ax.pcolormesh(X, Y, Z, cmap='Greens_r', norm=norm)
ax.plot([2 * N * np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning)[np.searchsorted(Z[index_shift_s0, :], 1.1)] if Z[index_shift_s0, -1] >= 1.1 else None for index_shift_s0, shift_s0 in enumerate(np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning))], np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning), color='lightgray', linestyle='-')
ax.plot([2 * N * np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning)[np.searchsorted(Z[index_shift_s0, :], 1.01)] if Z[index_shift_s0, -1] >= 1.01 else None for index_shift_s0, shift_s0 in enumerate(np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning))], np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning), color='lightgray', linestyle=':')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(2 * N * rate_of_shifts_list[0], 2 * N * rate_of_shifts_list[-1])
ax.set_ylim(shift_s0_list[0], shift_s0_list[-1])
ax.set_xlabel(r'Rate of shifts ($2Np$)', fontsize=fontsize)
ax.set_ylabel(r'Shift size ($\Lambda / \sqrt{V_A}$)', fontsize=fontsize)
assert len(rate_of_shifts_list) % 2 == 1
ax.set_xticks(ticks=list(2 * N * np.array(rate_of_shifts_list)), labels=[2 * N * rate_of_shift if index_rate_of_shift % 2 == 0 and 2 * N * rate_of_shift != int(2 * N * rate_of_shift) else int(2 * N * rate_of_shift) if index_rate_of_shift % 2 == 0 else None for index_rate_of_shift, rate_of_shift in enumerate(rate_of_shifts_list)], fontsize=fontsize)
ax.set_yticks(ticks=shift_s0_list, labels=shift_s0_list, fontsize=fontsize)
ax.minorticks_off()
ax.text(x=4, y=1, s='1%', color='k', rotation=-63, rotation_mode='anchor')
ax.text(x=48, y=1, s='10%', color='k', rotation=-63, rotation_mode='anchor')
ax.scatter(2 * N * np.array([0.003, 0.048]), [1.5, 1.5], color=['lightgray', 'k'])

# panels C and D: the heterozygosity as a function of the shift size (or the rate of shifts) for a given rate of shifts (or shift size), for a few effect sizes

ax = plt.subplot(gs[1, 0])
ax.text(0, 1.01, 'C', transform=ax.transAxes, fontsize=12, fontweight='bold')

index_rate_of_shift, rate_of_shift = 44, 0.048

neutral, = ax.plot([0] + shift_s0_list, [1 for _ in [0] + shift_s0_list], 'k:')

simulations = [None, None, None]
no_shifts = [None, None, None]
analytics = [None, None, None]

for i, (index_ss, ss, color) in enumerate([[0, 2.15, 'b'], [3, 7.55, 'purple'], [9, 34.1, 'r']]):
    heterozygosity_mean_wrt_shift_s0, heterozygosity_se_wrt_shift_s0 = calculate_heterozygosity_mean_and_se_across_shifts(my_AA_directories_across_shifts=[[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), j)) for j in range(1, 1 + parallel_AA)] for shift_s0 in shift_s0_list], index_ss=index_ss, nE=nE)
    simulations[i] = ax.errorbar(shift_s0_list, np.divide(heterozygosity_mean_wrt_shift_s0, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1)), yerr=1.96 * np.divide(heterozygosity_se_wrt_shift_s0, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1)), fmt='o', color=color)
    no_shifts[i], = ax.plot([0] + shift_s0_list, [heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) for _ in [0] + shift_s0_list], color=color, linestyle='--')
    analytics[i], = ax.plot([0] + list(np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning)), [heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1)] + [heterozygosity_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / heterozygosity_(a=0.1) for index_shift_s0, shift_s0 in enumerate(np.geomspace(shift_s0_list[0], shift_s0_list[-1], shift_s0_partitioning))], color=color)

ax.text(x=0, y=1.005, s='Neutral', color='k', fontsize=fontsize)
ax.text(x=shift_s0_list[0], y=0.75, s=r'$a^2=%g$' % 2, color='b', fontsize=fontsize)
ax.text(x=shift_s0_list[0], y=0.275, s=r'$a^2=%g$' % 8, color='purple', fontsize=fontsize)
ax.text(x=shift_s0_list[0], y=0.125, s=r'$a^2=%g$' % 32, color='r', fontsize=fontsize)

simulation_mean = ax.errorbar([], [], yerr=[], fmt='o', color='k', label='All-allele simulations')
no_shifts_mean, = ax.plot([], [], 'k--', label='No shifts')
analytic_mean, = ax.plot([], [], 'k-', label='Analytics')

ax.set_xscale('symlog', linthresh=shift_s0_list[0], linscale=0.15)
ax.set_xlabel(r'Shift size ($\Lambda / \sqrt{V_A}$)', fontsize=fontsize)
ax.set_ylabel('Heterozygosity (relative to neutral)', fontsize=fontsize)
ax.set_xticks(ticks=[0] + shift_s0_list, labels=[0] + shift_s0_list, fontsize=fontsize)
# ax.set_title('$E[a^2]=%d, 2Np=%g$' % (round(E2Ns), 2 * N * rate_of_shift), fontsize=fontsize)
ax.legend(bbox_to_anchor=(0., 0.67), loc='upper left', fontsize=fontsize)
ax.minorticks_off()

ax = plt.subplot(gs[1, 2])
ax.text(0, 1.01, 'D', transform=ax.transAxes, fontsize=12, fontweight='bold')

index_shift_s0, shift_s0 = 17, 1.5

neutral, = ax.plot([0] + list(2 * N * np.array(rate_of_shifts_list)), [1 for _ in [0] + rate_of_shifts_list], 'k:')

simulations = [None, None, None]
no_shifts = [None, None, None]
analytics = [None, None, None]

for i, (index_ss, ss, color) in enumerate([[0, 2.15, 'b'], [3, 7.55, 'purple'], [9, 34.1, 'r']]):
    heterozygosity_mean_wrt_rate_of_shifts, heterozygosity_se_wrt_rate_of_shifts = calculate_heterozygosity_mean_and_se_across_shifts(my_AA_directories_across_shifts=[[os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_d2ax_t_1_scaled_every_2N_generations_analytic_V_A_%d' % (round(E2Ns), j)) for j in range(1, 1 + parallel_AA)] for rate_of_shift in rate_of_shifts_list], index_ss=index_ss, nE=nE)
    simulations[i] = ax.errorbar(2 * N * np.array(rate_of_shifts_list), np.divide(heterozygosity_mean_wrt_rate_of_shifts, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1)), yerr=1.96 * np.divide(heterozygosity_se_wrt_rate_of_shifts, 2 * N * U / (nE + 1) * heterozygosity_(a=0.1)), fmt='o', color=color)
    no_shifts[i], = ax.plot([0] + list(2 * N * np.array(rate_of_shifts_list)), [heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1) for _ in [0] + rate_of_shifts_list], color=color, linestyle='--')
    analytics[i], = ax.plot([0] + list(2 * N * np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning)), [heterozygosity_(a=np.sqrt(ss)) / heterozygosity_(a=0.1)] + [heterozygosity_shifts_(a=np.sqrt(ss), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) / heterozygosity_(a=0.1) for index_rate_of_shift, rate_of_shift in enumerate(np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning))], color=color)

ax.set_xscale('symlog', linthresh=2 * N * rate_of_shifts_list[0], linscale=0.6)
ax.set_xlabel(r'Rate of shifts ($2Np$)', fontsize=fontsize)
ax.set_xticks(ticks=[0] + list(2 * N * np.array(rate_of_shifts_list)), labels=[0] + [2 * N * rate_of_shift if index_rate_of_shift % 2 == 0 and 2 * N * rate_of_shift != int(2 * N * rate_of_shift) else int(2 * N * rate_of_shift) if index_rate_of_shift % 2 == 0 else None for index_rate_of_shift, rate_of_shift in enumerate(rate_of_shifts_list)], fontsize=fontsize)
# ax.set_title('$E[a^2]=%d, \Lambda=%g\sqrt{V_A}$' % (round(E2Ns), shift_s0), fontsize=fontsize)
ax.minorticks_off()

plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'heterozygosity_nonLande.pdf'), bbox_inches='tight')
plt.show()
