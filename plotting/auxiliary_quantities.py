#!/usr/bin/env python
# coding: utf-8

# # Auxiliary quantities
# 
# The average expected change in allele frequency
# 
# The variance of change in allele frequency
# 
# The total number of shifts that a fixed mutation experiences
# 
# The excess number of aligned shifts than opposing shifts that a fixed mutation experiences

# In[1]:


from plot_functions import *


# In[2]:


# common parameter values
fontsize = 8
N = 2000
U = 0.025
parallel_AA = 49
ss_min = 0.1
ss_max = 100
save_directory = '/insomnia001/depts/pas_lab/users/jc5473/'


# In[3]:


with open(os.path.join(save_directory, 'empirical_average_V_A_Lande'), 'rb') as f:
    empirical_average_V_A = pickle.load(f)
    
with open(os.path.join(save_directory, 'empirical_average_V_A_nonLande'), 'rb') as f:
    empirical_average_V_A_nonLande = pickle.load(f)

E2Ns = 1
shift_s0_partitioning = 32
rate_of_shifts_partitioning = 57

analytic_V_A = np.zeros((shift_s0_partitioning, rate_of_shifts_partitioning))
for index_shift_s0 in range(shift_s0_partitioning):
    for index_rate_of_shift in range(rate_of_shifts_partitioning):
        with open(os.path.join(save_directory, 'analytic_V_A_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
            analytic_V_A[index_shift_s0][index_rate_of_shift] = pickle.load(f)

E2Ns = 16
shift_s0_partitioning_nonLande = 29
rate_of_shifts_partitioning_nonLande = 66

analytic_V_A_nonLande = np.zeros((shift_s0_partitioning_nonLande, rate_of_shifts_partitioning_nonLande))
for index_shift_s0 in range(shift_s0_partitioning_nonLande):
    for index_rate_of_shift in range(rate_of_shifts_partitioning_nonLande):
        with open(os.path.join(save_directory, 'analytic_V_A_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
            analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift] = pickle.load(f)


# ## The average expected change in allele frequency

# In[4]:


# the average expected change in allele frequency, in the Lande case
E2Ns = 1

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

ax.plot(np.geomspace(ss_min, ss_max), [quad(lambda x_0: E_Delta_x_stab_sel(a=np.sqrt(ss), x=x_0) * tau(a=np.sqrt(ss), x=x_0, p=1 / (2 * N)), 0.0, 1.0, points=[1/(2*N)])[0] / quad(lambda x_0: tau(a=np.sqrt(ss), x=x_0, p=1 / (2 * N)), 0.0, 1.0, points=[1/(2*N)])[0] for ss in np.geomspace(ss_min, ss_max)], color='lightgray', label='no shifts')

index_shift_s0, shift_s0 = 15, 1.5
index_rate_of_shift, rate_of_shift = 37, 0.048
# ax.plot(np.geomspace(ss_min, ss_max), calculate_average_E_Delta_x_across_effect_sizes(N=N, U=U, E2Ns=E2Ns, index_shift_s0=index_shift_s0, shift_s0=shift_s0, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss_min=ss_min, ss_max=ss_max), color='k', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))
with open(os.path.join(save_directory, 'average_E_Delta_x_across_effect_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
    ax.plot(np.geomspace(ss_min, ss_max), pickle.load(f), color='k', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))

shift_s0_list = [0.5, 0.75, 1, 1.5, 2, 3, 5]
rate_of_shifts_list = [0.0001875, 0.000375, 0.00075, 0.0015, 0.003, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384, 0.768]
shift_s0_partitioning = 32
rate_of_shifts_partitioning = 57
index_shift_s0, shift_s0 = 15, 1.5
index_rate_of_shift = np.searchsorted([pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a=1, x=1.0 / (2 * N), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) for index_rate_of_shift, rate_of_shift in enumerate(np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning))], 1 / (2 * N)) - 1
rate_of_shift = np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning)[index_rate_of_shift]
# ax.plot(np.geomspace(ss_min, ss_max), calculate_average_E_Delta_x_across_effect_sizes(N=N, U=U, E2Ns=E2Ns, index_shift_s0=index_shift_s0, shift_s0=shift_s0, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss_min=ss_min, ss_max=ss_max), color='k', linestyle=':', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))
with open(os.path.join(save_directory, 'average_E_Delta_x_across_effect_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
    ax.plot(np.geomspace(ss_min, ss_max), pickle.load(f), color='k', linestyle=':', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))

ax.set_xscale('log')
ax.set_ylim(-3e-5, 2e-4)
plt.xlabel('Effect size squared', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[5]:


E2Ns = 1

index_shift_s0, shift_s0 = 15, 1.5
rate_of_shifts_partitioning = 57
rate_of_shifts_min=0.0001875
rate_of_shifts_max=0.768

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
    # ax.plot(2 * N * np.geomspace(rate_of_shifts_min, rate_of_shifts_max, rate_of_shifts_partitioning), calculate_average_E_Delta_x_across_rates_of_shifts(N=N, U=U, E2Ns=E2Ns, index_shift_s0=index_shift_s0, shift_s0=shift_s0, rate_of_shifts_min=rate_of_shifts_min, rate_of_shifts_max=rate_of_shifts_max, rate_of_shifts_partitioning=rate_of_shifts_partitioning, ss=ss), color=color, label=r'$a^2=%g$' % ss)
    with open(os.path.join(save_directory, 'average_E_Delta_x_across_rates_of_shifts_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_E2Ns_%d_ss_' % (index_shift_s0, round(E2Ns)) + str(ss).replace('.', '_')), 'rb') as f:
        ax.plot(2 * N * np.geomspace(rate_of_shifts_min, rate_of_shifts_max, rate_of_shifts_partitioning), pickle.load(f), color=color, label=r'$a^2=%g$' % ss)

ax.set_xscale('log')
plt.xlabel(r'$2Np$', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}$' % shift_s0, fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_across_rates_of_shifts_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[6]:


E2Ns = 1

index_rate_of_shift, rate_of_shift = 28, 0.012
shift_s0_partitioning = 32
shift_s0_min = 0.5
shift_s0_max = 5

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
    # ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), calculate_average_E_Delta_x_across_shift_sizes(N=N, U=U, E2Ns=E2Ns, shift_s0_min=shift_s0_min, shift_s0_max=shift_s0_max, shift_s0_partitioning=shift_s0_partitioning, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss=ss), color=color, label=r'$a^2=%g$' % ss)
    with open(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_rate_of_shift_%d_E2Ns_%d_ss_' % (index_rate_of_shift, round(E2Ns)) + str(ss).replace('.', '_')), 'rb') as f:
        ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), pickle.load(f), color=color, label=r'$a^2=%g$' % ss)

ax.set_xscale('log')
plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(ticks=[0.5, 1, 2, 5], labels=[0.5, 1, 2, 5], fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.title(r'$2Np=%g$' % (2 * N * rate_of_shift), fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[7]:


E2Ns = 1

index_rate_of_shift, rate_of_shift = 37, 0.048
shift_s0_partitioning = 32
shift_s0_min = 0.5
shift_s0_max = 5

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
    # ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), calculate_average_E_Delta_x_across_shift_sizes(N=N, U=U, E2Ns=E2Ns, shift_s0_min=shift_s0_min, shift_s0_max=shift_s0_max, shift_s0_partitioning=shift_s0_partitioning, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss=ss), color=color, label=r'$a^2=%g$' % ss)
    with open(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_rate_of_shift_%d_E2Ns_%d_ss_' % (index_rate_of_shift, round(E2Ns)) + str(ss).replace('.', '_')), 'rb') as f:
        ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), pickle.load(f), color=color, label=r'$a^2=%g$' % ss)

ax.set_xscale('log')
plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(ticks=[0.5, 1, 2, 5], labels=[0.5, 1, 2, 5], fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.title(r'$2Np=%g$' % (2 * N * rate_of_shift), fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_frequent_shifts_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[39]:


# the average expected change in allele frequency, in the non-Lande case
E2Ns = 16

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

ax.plot(np.geomspace(ss_min, ss_max), [quad(lambda x_0: E_Delta_x_stab_sel(a=np.sqrt(ss), x=x_0) * tau(a=np.sqrt(ss), x=x_0, p=1 / (2 * N)), 0.0, 1.0, points=[1/(2*N)])[0] / quad(lambda x_0: tau(a=np.sqrt(ss), x=x_0, p=1 / (2 * N)), 0.0, 1.0, points=[1/(2*N)])[0] for ss in np.geomspace(ss_min, ss_max)], color='lightgray', label='no shifts')

index_shift_s0, shift_s0 = 17, 1.5
index_rate_of_shift, rate_of_shift = 54, 0.192
# ax.plot(np.geomspace(ss_min, ss_max), calculate_average_E_Delta_x_across_effect_sizes(N=N, U=U, E2Ns=E2Ns, index_shift_s0=index_shift_s0, shift_s0=shift_s0, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss_min=ss_min, ss_max=ss_max), color='k', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))
with open(os.path.join(save_directory, 'average_E_Delta_x_across_effect_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
    ax.plot(np.geomspace(ss_min, ss_max), pickle.load(f), color='k', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))

shift_s0_list = [0.5, 0.75, 1, 1.5, 2, 3]
rate_of_shifts_list = [0.0001875, 0.000375, 0.00075, 0.0015, 0.003, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.384, 0.768]
shift_s0_partitioning = 29
rate_of_shifts_partitioning = 66
index_shift_s0, shift_s0 = 17, 1.5
index_rate_of_shift = np.searchsorted([pi_recurrent_shifts_nonlinear_linear_V_Delta_x(a=1, x=1.0 / (2 * N), shift_s0=shift_s0, sigma_0_del=np.sqrt(analytic_V_A_nonLande[index_shift_s0][index_rate_of_shift]), p=rate_of_shift) for index_rate_of_shift, rate_of_shift in enumerate(np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning))], 1 / (2 * N))
rate_of_shift = np.geomspace(rate_of_shifts_list[0], rate_of_shifts_list[-1], rate_of_shifts_partitioning)[index_rate_of_shift]
# ax.plot(np.geomspace(ss_min, ss_max), calculate_average_E_Delta_x_across_effect_sizes(N=N, U=U, E2Ns=E2Ns, index_shift_s0=index_shift_s0, shift_s0=shift_s0, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss_min=ss_min, ss_max=ss_max), color='k', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))
with open(os.path.join(save_directory, 'average_E_Delta_x_across_effect_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns))), 'rb') as f:
    ax.plot(np.geomspace(ss_min, ss_max), pickle.load(f), color='k', linestyle=':', label=r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift))

ax.set_xscale('log')
ax.set_ylim(-3e-5, 6e-5)
plt.xlabel('Effect size squared', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[40]:


E2Ns = 16

index_shift_s0, shift_s0 = 17, 1.5
rate_of_shifts_partitioning = 66
rate_of_shifts_min = 0.0001875
rate_of_shifts_max = 0.768

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
    # ax.plot(2 * N * np.geomspace(rate_of_shifts_min, rate_of_shifts_max, rate_of_shifts_partitioning), calculate_average_E_Delta_x_across_rates_of_shifts(N=N, U=U, E2Ns=E2Ns, index_shift_s0=index_shift_s0, shift_s0=shift_s0, rate_of_shifts_min=rate_of_shifts_min, rate_of_shifts_max=rate_of_shifts_max, rate_of_shifts_partitioning=rate_of_shifts_partitioning, ss=ss), color=color, label=r'$a^2=%g$' % ss)
    with open(os.path.join(save_directory, 'average_E_Delta_x_across_rates_of_shifts_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_E2Ns_%d_ss_' % (index_shift_s0, round(E2Ns)) + str(ss).replace('.', '_')), 'rb') as f:
        ax.plot(2 * N * np.geomspace(rate_of_shifts_min, rate_of_shifts_max, rate_of_shifts_partitioning), pickle.load(f), color=color, label=r'$a^2=%g$' % ss)

ax.set_xscale('log')
plt.xlabel(r'$2Np$', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}$' % shift_s0, fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_across_rates_of_shifts_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[41]:


E2Ns = 16

index_rate_of_shift, rate_of_shift = 44, 0.048
shift_s0_partitioning = 29
shift_s0_min = 0.5
shift_s0_max = 3

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
    # ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), calculate_average_E_Delta_x_across_shift_sizes(N=N, U=U, E2Ns=E2Ns, shift_s0_min=shift_s0_min, shift_s0_max=shift_s0_max, shift_s0_partitioning=shift_s0_partitioning, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss=ss), color=color, label=r'$a^2=%g$' % ss)
    with open(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_rate_of_shift_%d_E2Ns_%d_ss_' % (index_rate_of_shift, round(E2Ns)) + str(ss).replace('.', '_')), 'rb') as f:
        ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), pickle.load(f), color=color, label=r'$a^2=%g$' % ss)

ax.set_xscale('log')
plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(ticks=[1, 3], labels=[1, 3], fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.title(r'$2Np=%g$' % (2 * N * rate_of_shift), fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[42]:


E2Ns = 16

index_rate_of_shift, rate_of_shift = 54, 0.192
shift_s0_partitioning = 29
shift_s0_min = 0.5
shift_s0_max = 3

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
    # ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), calculate_average_E_Delta_x_across_shift_sizes(N=N, U=U, E2Ns=E2Ns, shift_s0_min=shift_s0_min, shift_s0_max=shift_s0_max, shift_s0_partitioning=shift_s0_partitioning, index_rate_of_shift=index_rate_of_shift, rate_of_shift=rate_of_shift, ss=ss), color=color, label=r'$a^2=%g$' % ss)
    with open(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_N_%d_U_' % N + str(U).replace('.', '_') + '_index_rate_of_shift_%d_E2Ns_%d_ss_' % (index_rate_of_shift, round(E2Ns)) + str(ss).replace('.', '_')), 'rb') as f:
        ax.plot(np.geomspace(shift_s0_min, shift_s0_max, shift_s0_partitioning), pickle.load(f), color=color, label=r'$a^2=%g$' % ss)

ax.set_xscale('log')
plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
plt.ylabel(r'$\overline{E[\Delta x]}$', fontsize=fontsize)
plt.xticks(ticks=[1, 3], labels=[1, 3], fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.legend(fontsize=fontsize)
plt.title(r'$2Np=%g$' % (2 * N * rate_of_shift), fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'average_E_Delta_x_across_shift_sizes_frequent_shifts_nonLande.pdf'), bbox_inches='tight')
plt.show()


# ## The variance of change in allele frequency

# In[5]:


# the variance of change in allele frequency, in the Lande case
fontsize = 8

N = 2000
U = 0.025
E2Ns = 1

shift_s0 = 1.5

for x in (1 / (2 * N), 0.1, 0.5):
    for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = fig.add_subplot(1, 1, 1)
    
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [x * (1 - x) / (2 * N) for _ in np.geomspace(0.0001875, 0.768)], color=color, alpha=0.15, label=r'$a^2=%g$, drift' % ss)
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [V_Delta_x_drift_and_shifts(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for rate_of_shift in np.geomspace(0.0001875, 0.768)], color=color, linestyle='--', label=r'$a^2=%g$, linear' % ss)
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for rate_of_shift in np.geomspace(0.0001875, 0.768)], color=color, linestyle=':', label=r'$a^2=%g$, nonlinear w/o asymmetry' % ss)
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [V_Delta_x_drift_and_shifts_nonlinear(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for rate_of_shift in np.geomspace(0.0001875, 0.768)], color=color, label=r'$a^2=%g$, nonlinear' % ss)
    
        ax.set_xscale('log')
        plt.xlabel(r'$2Np$', fontsize=fontsize)
        plt.ylabel(r'$V(\Delta x)$', fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title(r'$\Lambda=%g\sqrt{V_A}, x_0=%g$' % (shift_s0, x), fontsize=fontsize)
        plt.savefig(os.path.join(save_directory, 'V_Delta_x_across_rates_of_shifts_x_' + str(x).replace('.', '_') + '_ss_' + str(ss).replace('.', '_') + '_E2Ns_%d.pdf' % round(E2Ns)), bbox_inches='tight')
        plt.show()


# In[6]:


fontsize = 8

N = 2000
U = 0.025
E2Ns = 1

rate_of_shift = 0.012

for x in (1 / (2 * N), 0.1, 0.5):
    for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = fig.add_subplot(1, 1, 1)
    
        ax.plot(np.geomspace(0.5, 5), [x * (1 - x) / (2 * N) for _ in np.geomspace(0.5, 5)], color=color, alpha=0.15, label=r'$a^2=%g$, drift' % ss)
        ax.plot(np.geomspace(0.5, 5), [V_Delta_x_drift_and_shifts(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 5)], color=color, linestyle='--', label=r'$a^2=%g$, linear' % ss)
        ax.plot(np.geomspace(0.5, 5), [V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 5)], color=color, linestyle=':', label=r'$a^2=%g$, nonlinear w/o asymmetry' % ss)
        ax.plot(np.geomspace(0.5, 5), [V_Delta_x_drift_and_shifts_nonlinear(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 5)], color=color, label=r'$a^2=%g$, nonlinear' % ss)
    
        ax.set_xscale('log')
        plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
        plt.ylabel(r'$V(\Delta x)$', fontsize=fontsize)
        plt.xticks(ticks=[0.5, 1, 2, 5], labels=[0.5, 1, 2, 5], fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title(r'$2Np=%g, x_0=%g$' % (2 * N * rate_of_shift, x), fontsize=fontsize)
        plt.savefig(os.path.join(save_directory, 'V_Delta_x_across_shift_sizes_x_' + str(x).replace('.', '_') + '_ss_' + str(ss).replace('.', '_') + '_E2Ns_%d.pdf' % round(E2Ns)), bbox_inches='tight')
        plt.show()


# In[7]:


fontsize = 8

N = 2000
U = 0.025
E2Ns = 1

rate_of_shift = 0.048

for x in (1 / (2 * N), 0.1, 0.5):
    for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = fig.add_subplot(1, 1, 1)
    
        ax.plot(np.geomspace(0.5, 5), [x * (1 - x) / (2 * N) for _ in np.geomspace(0.5, 5)], color=color, alpha=0.15, label=r'$a^2=%g$, drift' % ss)
        ax.plot(np.geomspace(0.5, 5), [V_Delta_x_drift_and_shifts(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 5)], color=color, linestyle='--', label=r'$a^2=%g$, linear' % ss)
        ax.plot(np.geomspace(0.5, 5), [V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 5)], color=color, linestyle=':', label=r'$a^2=%g$, nonlinear w/o asymmetry' % ss)
        ax.plot(np.geomspace(0.5, 5), [V_Delta_x_drift_and_shifts_nonlinear(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 5)], color=color, label=r'$a^2=%g$, nonlinear' % ss)
    
        ax.set_xscale('log')
        plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
        plt.ylabel(r'$V(\Delta x)$', fontsize=fontsize)
        plt.xticks(ticks=[0.5, 1, 2, 5], labels=[0.5, 1, 2, 5], fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title(r'$2Np=%g, x_0=%g$' % (2 * N * rate_of_shift, x), fontsize=fontsize)
        plt.savefig(os.path.join(save_directory, 'V_Delta_x_across_shift_sizes_frequent_shifts_x_' + str(x).replace('.', '_') + '_ss_' + str(ss).replace('.', '_') + '_E2Ns_%d.pdf' % round(E2Ns)), bbox_inches='tight')
        plt.show()


# In[8]:


# the variance of change in allele frequency, in the non-Lande case
fontsize = 8

N = 2000
U = 0.025
E2Ns = 16

shift_s0 = 1.5

for x in (1 / (2 * N), 0.1, 0.5):
    for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = fig.add_subplot(1, 1, 1)
    
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [x * (1 - x) / (2 * N) for _ in np.geomspace(0.0001875, 0.768)], color=color, alpha=0.15, label=r'$a^2=%g$, drift' % ss)
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [V_Delta_x_drift_and_shifts(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for rate_of_shift in np.geomspace(0.0001875, 0.768)], color=color, linestyle='--', label=r'$a^2=%g$, linear' % ss)
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for rate_of_shift in np.geomspace(0.0001875, 0.768)], color=color, linestyle=':', label=r'$a^2=%g$, nonlinear w/o asymmetry' % ss)
        ax.plot(2 * N * np.geomspace(0.0001875, 0.768), [V_Delta_x_drift_and_shifts_nonlinear(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for rate_of_shift in np.geomspace(0.0001875, 0.768)], color=color, label=r'$a^2=%g$, nonlinear' % ss)
    
        ax.set_xscale('log')
        plt.xlabel(r'$2Np$', fontsize=fontsize)
        plt.ylabel(r'$V(\Delta x)$', fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title(r'$\Lambda=%g\sqrt{V_A}, x_0=%g$' % (shift_s0, x), fontsize=fontsize)
        plt.savefig(os.path.join(save_directory, 'V_Delta_x_across_rates_of_shifts_x_' + str(x).replace('.', '_') + '_ss_' + str(ss).replace('.', '_') + '_E2Ns_%d.pdf' % round(E2Ns)), bbox_inches='tight')
        plt.show()


# In[9]:


fontsize = 8

N = 2000
U = 0.025
E2Ns = 16

rate_of_shift = 0.048

for x in (1 / (2 * N), 0.1, 0.5):
    for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = fig.add_subplot(1, 1, 1)
    
        ax.plot(np.geomspace(0.5, 3), [x * (1 - x) / (2 * N) for _ in np.geomspace(0.5, 3)], color=color, alpha=0.15, label=r'$a^2=%g$, drift' % ss)
        ax.plot(np.geomspace(0.5, 3), [V_Delta_x_drift_and_shifts(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 3)], color=color, linestyle='--', label=r'$a^2=%g$, linear' % ss)
        ax.plot(np.geomspace(0.5, 3), [V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 3)], color=color, linestyle=':', label=r'$a^2=%g$, nonlinear w/o asymmetry' % ss)
        ax.plot(np.geomspace(0.5, 3), [V_Delta_x_drift_and_shifts_nonlinear(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 3)], color=color, label=r'$a^2=%g$, nonlinear' % ss)
    
        ax.set_xscale('log')
        plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
        plt.ylabel(r'$V(\Delta x)$', fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title(r'$2Np=%g, x_0=%g$' % (2 * N * rate_of_shift, x), fontsize=fontsize)
        plt.savefig(os.path.join(save_directory, 'V_Delta_x_across_shift_sizes_x_' + str(x).replace('.', '_') + '_ss_' + str(ss).replace('.', '_') + '_E2Ns_%d.pdf' % round(E2Ns)), bbox_inches='tight')
        plt.show()


# In[10]:


fontsize = 8

N = 2000
U = 0.025
E2Ns = 16

rate_of_shift = 0.192

for x in (1 / (2 * N), 0.1, 0.5):
    for ss, color in ((0.5, 'r'), (5, 'g'), (35, 'b')):
        fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
        ax = fig.add_subplot(1, 1, 1)
    
        ax.plot(np.geomspace(0.5, 3), [x * (1 - x) / (2 * N) for _ in np.geomspace(0.5, 3)], color=color, alpha=0.15, label=r'$a^2=%g$, drift' % ss)
        ax.plot(np.geomspace(0.5, 3), [V_Delta_x_drift_and_shifts(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 3)], color=color, linestyle='--', label=r'$a^2=%g$, linear' % ss)
        ax.plot(np.geomspace(0.5, 3), [V_Delta_x_drift_and_shifts_nonlinear_without_asymmetry(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 3)], color=color, linestyle=':', label=r'$a^2=%g$, nonlinear w/o asymmetry' % ss)
        ax.plot(np.geomspace(0.5, 3), [V_Delta_x_drift_and_shifts_nonlinear(a=np.sqrt(ss), x=x, shift_s0=shift_s0, sigma_0_del=np.sqrt(V_A(E2Ns=E2Ns)), p=rate_of_shift) for shift_s0 in np.geomspace(0.5, 3)], color=color, label=r'$a^2=%g$, nonlinear' % ss)
    
        ax.set_xscale('log')
        plt.xlabel(r'$\Lambda / \sqrt{V_A}$', fontsize=fontsize)
        plt.ylabel(r'$V(\Delta x)$', fontsize=fontsize)
        plt.xticks(fontsize=fontsize)
        plt.yticks(fontsize=fontsize)
        plt.legend(fontsize=fontsize)
        plt.title(r'$2Np=%g, x_0=%g$' % (2 * N * rate_of_shift, x), fontsize=fontsize)
        plt.savefig(os.path.join(save_directory, 'V_Delta_x_across_shift_sizes_frequent_shifts_x_' + str(x).replace('.', '_') + '_ss_' + str(ss).replace('.', '_') + '_E2Ns_%d.pdf' % round(E2Ns)), bbox_inches='tight')
        plt.show()


# ## The total number of shifts that a fixed mutation experiences

# In[14]:


# the histogram of the total number of shifts that a fixed mutation experiences, in the Lande case
# linear scale
E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)
    
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_frequent_shifts_linear_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[15]:


# log scale
E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)
    
# ax.axvline(x=rate_of_shift * np.mean([fixation_timing[1] - fixation_timing[0] for fixation_timing in fixation_timings]), color='k', label='Naive')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_frequent_shifts_log_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[16]:


E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 6, 0.012

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)
    
# ax.axvline(x=rate_of_shift * np.mean([fixation_timing[1] - fixation_timing[0] for fixation_timing in fixation_timings]), color='k', label='Naive')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_linear_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[17]:


E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 6, 0.012

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)
    
# ax.axvline(x=rate_of_shift * np.mean([fixation_timing[1] - fixation_timing[0] for fixation_timing in fixation_timings]), color='k', label='Naive')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_log_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[18]:


# the histogram of the total number of shifts that a fixed mutation experiences, in the non-Lande case
# linear scale
E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 10, 0.192

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)

ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_frequent_shifts_linear_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[19]:


# log scale
E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 10, 0.192

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)

ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_frequent_shifts_log_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[20]:


E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)

ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_linear_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[21]:


E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_shifts_experienced = list()
    # naive_samples = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     shift_times = [tM for tM in range(len(dist_ess_over_time)) if abs(dist_ess_over_time[tM] - dist_ess_over_time[tM - 1]) >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             numbers_of_shifts_experienced.append(np.searchsorted(shift_times, fixation_timing[1]) - np.searchsorted(shift_times, fixation_timing[0]))
    #             naive_samples.append(rate_of_shift * (fixation_timing[1] - fixation_timing[0]))
    with open(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_shifts_experienced, naive_samples = pickle.load(f)
    
    ax.axvline(x=np.mean(naive_samples), color=color, linestyle='--')
    ax.axvspan(xmin=np.mean(numbers_of_shifts_experienced) - 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), xmax=np.mean(numbers_of_shifts_experienced) + 1.96 * np.divide(np.std(numbers_of_shifts_experienced), np.sqrt(np.size(numbers_of_shifts_experienced))), color=color, alpha=0.15)
    ax.axvline(x=np.mean(numbers_of_shifts_experienced), color=color)
    
    ax.hist(numbers_of_shifts_experienced, density=True, color=color, label=r'$a^2=%g$' % ss)

ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.plot([], [], color='k', linestyle='--', label='Naive')
ax.set_xscale('log')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('Number of shifts a fixed mutation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'number_of_shifts_a_fixed_mutation_experiences_log_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# ## The excess number of aligned shifts than opposing shifts that a fixed mutation experiences

# In[30]:


# the histogram of the excess number of aligned shifts than opposing shifts that a fixed mutation experiences, in the Lande case
# log scale
E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.set_xscale('symlog', linthresh=5, linscale=1)
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_frequent_shifts_symlog_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[31]:


# linear scale
E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_frequent_shifts_linear_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[32]:


E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 6, 0.012

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.set_xscale('symlog', linthresh=5, linscale=1)
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_symlog_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[33]:


E2Ns = 1
nE = 100

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 6, 0.012

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((38, 0.5, 'r'), (98, 4, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_linear_scale_Lande.pdf'), bbox_inches='tight')
plt.show()


# In[34]:


# the histogram of the excess number of aligned shifts than opposing shifts that a fixed mutation experiences, in the non-Lande case
# linear scale
E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 10, 0.192

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_frequent_shifts_linear_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[35]:


# log scale
E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 10, 0.192

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.set_xscale('symlog', linthresh=5, linscale=1)
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_frequent_shifts_symlog_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[36]:


E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_linear_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[37]:


E2Ns = 16
nE = 11

V2Ns = E2Ns ** 2
S_dist = gamma(float(E2Ns) ** 2 / float(V2Ns), loc=0.,
                                               scale=float(V2Ns) / float(E2Ns))
a_list_pos = [math.sqrt(S_dist.ppf((i + 1) / (nE + 1))) for i in range(nE)]

index_shift_s0, shift_s0 = 3, 1.5
index_rate_of_shift, rate_of_shift = 8, 0.048

fig = plt.figure(figsize=(3.5, 3.5), dpi=300)
ax = fig.add_subplot(1, 1, 1)

for index_ss, ss, color in ((2, 5, 'r'), (9, 35, 'g')):
    # numbers_of_aligned_shifts_experienced = list()
    # numbers_of_opposing_shifts_experienced = list()
    # parallel_AA = 49
    # for i in range(1, 1 + parallel_AA):
    #     with open(os.path.join(save_directory, 'my_AA_N_%d_U_' % N + str(U).replace('.', '_') + '_shift_s0_' + str(shift_s0).replace('.', '_') + '_rate_of_shift_' + str(rate_of_shift).replace('.', '_') + '_E2Ns_%d_no_efs_bins_time_arising_fixing_empirical_average_V_A_%d/fixations_and_extinctions' % (round(E2Ns), i)), 'rb') as f:
    #         fixations_and_extinctions, _, dist_ess_over_time, _, _, _, fixation_timings = pickle.load(f)
    #     trait_increasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] >= 0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     trait_decreasing_shift_times = [tM for tM in range(len(dist_ess_over_time)) if dist_ess_over_time[tM] - dist_ess_over_time[tM - 1] <= -0.5 * shift_s0 * np.sqrt(empirical_average_V_A_nonLande[index_shift_s0][index_rate_of_shift])]
    #     for fixation_timing, a in zip(fixation_timings, [mut[0] for mut in fixations_and_extinctions if mut[1]]):
    #         if np.searchsorted(a_list_pos, np.abs(a)) == index_ss + 1:
    #             number_of_trait_increasing_shifts_experienced = np.searchsorted(trait_increasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_increasing_shift_times, fixation_timing[0])
    #             number_of_trait_decreasing_shifts_experienced = np.searchsorted(trait_decreasing_shift_times, fixation_timing[1]) - np.searchsorted(trait_decreasing_shift_times, fixation_timing[0])
    #             if a > 0:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #             else:
    #                 numbers_of_aligned_shifts_experienced.append(number_of_trait_decreasing_shifts_experienced)
    #                 numbers_of_opposing_shifts_experienced.append(number_of_trait_increasing_shifts_experienced)
    with open(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_N_%d_U_' % N + str(U).replace('.', '_') + '_index_shift_s0_%d_index_rate_of_shift_%d_E2Ns_%d_index_ss_%d_nE_%d' % (index_shift_s0, index_rate_of_shift, round(E2Ns), index_ss, nE)), 'rb') as f:
        numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced = pickle.load(f)
    ax.hist(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced), density=True, color=color, label=r'$a^2=%g$' % ss)
    ax.axvline(x=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), color=color)
    ax.axvspan(xmin=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) - 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), xmax=np.mean(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)) + 1.96 * np.divide(np.std(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)), np.sqrt(np.size(np.subtract(numbers_of_aligned_shifts_experienced, numbers_of_opposing_shifts_experienced)))), color=color, alpha=0.15)

ax.axvline(x=0, color='k', linestyle='--', label='#aligned = #opposing')
ax.plot([], [], color='k', label=r'$mean\pm 1.96SEs$')
ax.set_xscale('symlog', linthresh=5, linscale=1)
ax.legend(bbox_to_anchor=(0., -.15), loc='upper left', fontsize=fontsize)
plt.xlabel('#aligned - #opposing shifts a fixation experiences', fontsize=fontsize)
plt.ylabel('Density', fontsize=fontsize)
plt.title(r'$\Lambda=%g\sqrt{V_A}, 2Np=%g$' % (shift_s0, 2 * N * rate_of_shift), fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.savefig(os.path.join(save_directory, 'excess_number_of_aligned_than_opposing_shifts_a_fixed_mutation_experiences_symlog_scale_nonLande.pdf'), bbox_inches='tight')
plt.show()


# In[ ]:




