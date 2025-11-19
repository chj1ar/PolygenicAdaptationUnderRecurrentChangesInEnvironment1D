import pickle
from mutation_functions import Mutation
import numpy as np
import math
import os
from scipy.special import dawsn
from scipy.integrate import quad
from scipy.stats import gamma



class SimulatePopulations(object):
    """
    Class to store one set of simulation parameters, evolve with this set of simulation parameters, and save the resulting stats under this set of simulation parameters.
    If intermediate_Folder is not None, then use the parameters there instead, and save the simulation results into intermediate_Folder.
    Attributes:
        N: (int) population size
        Vs: (float) fitness parameter
        U: (float) mutation rate (per gamete per generation)
        shape: (float) the shape parameter of gamma distributed shift sizes in units of the phenotypic standard deviation
        scale: (float) the scale parameter of gamma distributed shift sizes in units of the phenotypic standard deviation
        shift_s0: (float) shift size in units of the phenotypic standard deviation
        E2Ns: (float) the expectation of population-scaled selection coefficients
        V2Ns: (float) the variance of population-scaled selection coefficients
        nE: (int) number of effect size bins (currently for the fixation probability)
        nF: (int) number of frequency bins (currently unused)
        burn_time_N: (int) burn time in units of population size
        rate_of_shift: (float) rate of shifts
        directories_for_empirical_V_A: (str) the directories where the empirical average V_A is used here
        sigma_0_del: (float) specify this if you want to run a simulation using its square as the empirical V_A (usually from some previous iteration of this simulation)
        nnm: (int) number of new mutations to calculate the fixation probability
        particular_2Ns: (List[float]) particular values of population-scaled selection coefficients
        save_directory: (str) the directory in which you want your simulation results to be saved
    """
    def __init__(self, intermediate_Folder=None, N=5000, U=0.01, Vs=-1, shape=None, scale=None, shift_s0=None, E2Ns=2.0, V2Ns=-1, nE=6, nF=11, burn_time_N=10, rate_of_shift=None, directories_for_empirical_V_A=None, sigma_0_del=None, nnm=None, particular_2Ns=[], save_directory=None):
        if intermediate_Folder is not None:
            save_directory = intermediate_Folder
            with open(os.path.join(intermediate_Folder, 'identifiers.txt'), 'r') as f:
                N = int(f.readline().split('N: ')[1])
                Vs = float(f.readline().split('Vs: ')[1])
                U = float(f.readline().split('U: ')[1])
                shift_s0 = float(f.readline().split('shift_s0: ')[1])
                _ = float(f.readline().split('var_0: ')[1])
                sigma_0_del = float(f.readline().split('sigma_0_del: ')[1])
                _ = float(f.readline().split('shift: ')[1])
                E2Ns = float(f.readline().split('E2Ns: ')[1])
                V2Ns = float(f.readline().split('V2Ns: ')[1])
                nE = int(f.readline().split('nE: ')[1])
                nF = int(f.readline().split('nF: ')[1])
                _ = int(f.readline().split('burn_time: ')[1])
                burn_time_N = int(f.readline().split('burn_time_N: ')[1])
                rate_of_shift = float(f.readline().split('rate_of_shift: ')[1])
                nnm = int(f.readline().split('nnm: ')[1])

        self.N = N
        # if V_S < 0, set it to be 2N
        if Vs < 0:
            self.Vs = float(2 * self.N)
        else:
            self.Vs = Vs
        self.U = U
        self.E2Ns = E2Ns
        # if V2Ns < 0, set it to be E2Ns^2, and hence the gamma distribution becomes exponential distribution with mean E2Ns
        if V2Ns < 0:
            self.V2Ns = float(self.E2Ns ** 2)
        else:
            self.V2Ns = V2Ns
        self.shape = shape
        self.scale = scale
        self.shift_s0 = shift_s0
        self._VAR_0 = self._var_0()
        self.nE = nE
        self.nF = nF
        self.burn_time_N = burn_time_N
        self._BURN_TIME = int(self.burn_time_N * self.N)
        self.rate_of_shift = rate_of_shift
        self.directories_for_empirical_V_A = directories_for_empirical_V_A
        if self.directories_for_empirical_V_A is None: # if directories_for_empirical_V_A is not specified, then we seek sigma_0_del
            if sigma_0_del is None: # if sigma_0_del is also not specified, then we use the analytical V_A
                self.sigma_0_del = self._VAR_0 ** 0.5
            else:
                self.sigma_0_del = sigma_0_del
        else:
            mean_of_var_ess = list()
            parallel = int(self.directories_for_empirical_V_A.split('_')[-1])
            for i in range(1, 1 + parallel):
                if not os.path.exists(os.path.join('_'.join(self.directories_for_empirical_V_A.split('_')[:-1]) + '_' + str(i), 'fixation_probability')):
                    continue
                with open(os.path.join('_'.join(self.directories_for_empirical_V_A.split('_')[:-1]) + '_' + str(i), 'fixation_probability'), 'rb') as f:
                    var_ess_over_time = pickle.load(f)[3]
                mean_of_var_ess.append(np.mean(var_ess_over_time[self._BURN_TIME:]))
            self.sigma_0_del = np.sqrt(np.mean(mean_of_var_ess))
        self.nnm = nnm
        self.particular_2Ns = particular_2Ns
        self.save_directory = save_directory

    def initiate_population(self, no_efs_bins):
        self._population_class = _Population(N=self.N, U=self.U, Vs=self.Vs, E2Ns=self.E2Ns, V2Ns=self.V2Ns, nE=self.nE, nF=self.nF, particular_2Ns=self.particular_2Ns, no_efs_bins=no_efs_bins)

    def run_population_under_recurrent_shifts(self):
        """
        The population undergoes the stationary process of recurrent shifts throughout, so I combine burn-in and running
        Returns:
            /
        """
        self._BURN_TIME = int(self.burn_time_N * self.N)
        # Whether there are still segregating mutations that are the first `self.nnm` new mutations
        self._THERE_ARE_MUTANTS = True
        # Count the number of new mutations
        nnm = 0
        # Count the number of generations
        tM = 0
        # Whether we have frozen and recorded contribution to phenotypic change
        self._have_we_frozen_and_recorded = False

        if self.rate_of_shift is not None:
            # self._t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(self._shift)
            # time_to_record_jump_sizes_after_jump = -1
            frozen_time = -1 - self.Vs / self.sigma_0_del ** 2 * math.log(2 * math.sqrt(self.Vs))
            waiting_time = np.random.geometric(p=self.rate_of_shift)

        while self._THERE_ARE_MUTANTS:
            # recurrent shifts
            if self.rate_of_shift is not None:
                if tM == waiting_time:
                    direction_of_shift = np.random.binomial(n=1, p=0.5) * 2 - 1
                    if self.shape is None and self.scale is None:
                        assert self.shift_s0 is not None
                        self._population_class.shift_optimum(self.shift_s0 * self.sigma_0_del * direction_of_shift)
                        t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(self.shift_s0 * self.sigma_0_del)
                    else:
                        assert self.shape is not None and self.scale is not None
                        shift_s0 = np.random.gamma(shape=self.shape, scale=self.scale)
                        self._population_class.shift_optimum(shift_s0 * self.sigma_0_del * direction_of_shift)
                        t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(shift_s0 * self.sigma_0_del)
                    # freeze at the first complete jump after burning
                    if tM >= self._BURN_TIME and not self._have_we_frozen_and_recorded:
                        self._population_class.freeze_recurrent_shifts()
                        frozen_time = waiting_time

                    # # record jump sizes during main run
                    # if tM >= self._BURN_TIME:
                    #     self._population_class.record_jump_sizes_before_jump()
                    #     time_to_record_jump_sizes_after_jump = waiting_time + round(self._t_1)
                    waiting_time += np.random.geometric(p=self.rate_of_shift)

            # burn-in
            if tM < self._BURN_TIME:
                self._population_class.next_gen(is_new_mutation=False, is_later_new_mutation=False)
            # main run
            else:
                # if self.rate_of_shift is not None:
                #     if tM == time_to_record_jump_sizes_after_jump:
                #         self._population_class.record_jump_sizes_after_jump(self._shift * direction_of_shift)
                # record contribution to phenotypic change
                if self.rate_of_shift is not None:
                    if tM == frozen_time + round(t_1):
                        self._population_class.record_contribution_to_phenotypic_change_recurrent_shifts()
                        self._have_we_frozen_and_recorded = True
                # record integral heterozygosity in this generation
                self._population_class.record_integral_heterozygosity_current_generation()
                # record the state
                if tM % self.N == 0:
                    self._save_identifiers()
                    self._population_class.record_mutations_and_phenotypic_distribution(nnm=nnm, directory=self.save_directory)
                if nnm < self.nnm:
                    nnm += self._population_class.next_gen(is_new_mutation=True, is_later_new_mutation=False)[1]
                else:
                    self._THERE_ARE_MUTANTS = self._population_class.next_gen(is_new_mutation=True, is_later_new_mutation=True)

            tM += 1

        # Record fixation probability
        self._population_class.record_fixation_probability(nm=True)
        # Record integral heterozygosity
        self._population_class.record_integral_heterozygosity()

    def run_population_under_recurrent_shifts_resume(self):
        """
        Resume the simulation from the most recently recorded state
        Returns:
            /
        """
        self._THERE_ARE_MUTANTS = True
        nnm = self._population_class.load_muts(directory=self.save_directory)
        tM = 0

        if self.rate_of_shift is not None:
            self._t_1 = self.Vs / self._VAR_0 * math.log(self._shift)
            waiting_time = np.random.geometric(p=self.rate_of_shift)

        while self._THERE_ARE_MUTANTS:
            # recurrent shifts
            if self.rate_of_shift is not None:
                if tM == waiting_time:
                    direction_of_shift = np.random.binomial(n=1, p=0.5) * 2 - 1
                    self._population_class.shift_optimum(self._shift * direction_of_shift)
                    waiting_time += np.random.geometric(p=self.rate_of_shift)

            if nnm < self.nnm:
                nnm += self._population_class.next_gen(is_new_mutation=True, is_later_new_mutation=False)[1]
            else:
                self._THERE_ARE_MUTANTS = self._population_class.next_gen(is_new_mutation=True, is_later_new_mutation=True)

            tM += 1

        # Record fixation probability
        self._population_class.record_fixation_probability(nm=True)

    def run_population_under_recurrent_shifts_no_efs_bins(self):
        """
        No effect size bins in order to obtain more informative simulation results
        Returns:
            /
        """
        self._BURN_TIME = int(self.burn_time_N * self.N)
        # Whether there are still segregating mutations that are the first `self.nnm` new mutations
        self._THERE_ARE_MUTANTS = True
        # Count the number of new mutations
        nnm = 0
        # Count the number of generations
        tM = 0
        # # Whether we have frozen and recorded contribution to phenotypic change
        # self._have_we_frozen_and_recorded = False
        # make sure that the recorded d2ax's are independent
        time_to_record_d2ax = 0
        # # record jump sizes after burn-in
        # indices = list()

        if self.rate_of_shift is not None:
            # self._t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(self._shift)
            # time_to_record_jump_sizes_t_1_after_jump = -1
            frozen_time = -1 - self.Vs / self.sigma_0_del ** 2 * math.log(2 * math.sqrt(self.Vs))
            waiting_time = np.random.geometric(p=self.rate_of_shift)

        while self._THERE_ARE_MUTANTS:
            # recurrent shifts
            if self.rate_of_shift is not None:
                if tM == waiting_time:
                    direction_of_shift = np.random.binomial(n=1, p=0.5) * 2 - 1
                    if self.shape is None and self.scale is None:
                        assert self.shift_s0 is not None
                        self._population_class.shift_optimum(self.shift_s0 * self.sigma_0_del * direction_of_shift)
                        t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(self.shift_s0 * self.sigma_0_del)
                    else:
                        assert self.shape is not None and self.scale is not None
                        shift_s0 = np.random.gamma(shape=self.shape, scale=self.scale)
                        self._population_class.shift_optimum(shift_s0 * self.sigma_0_del * direction_of_shift)
                        t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(shift_s0 * self.sigma_0_del)
                    # freeze at complete jumps after burn-in
                    if tM >= self._BURN_TIME:
                        self._population_class.freeze_recurrent_shifts()
                        frozen_time = waiting_time

                    # # record jump sizes until the next shift after burn-in
                    # if tM >= self._BURN_TIME:
                    #     if len(indices) > 0:
                    #         self._population_class.record_jump_sizes_until_next_shift_after_jump_no_efs_bins(indices, waiting_time_increment)
                    #     self._population_class.remove_extinct_fixed_from_seg_no_efs_bins()
                    #     indices = self._population_class.record_jump_sizes_before_jump_no_efs_bins()
                    #     time_to_record_jump_sizes_t_1_after_jump = waiting_time + round(self._t_1)
                    # else:
                    #     self._population_class.remove_extinct_fixed_from_seg_no_efs_bins()

                    waiting_time_increment = np.random.geometric(p=self.rate_of_shift)
                    waiting_time += waiting_time_increment

            # burn-in
            if tM < self._BURN_TIME:
                self._population_class.next_gen_no_efs_bins(is_new_mutation=False, is_later_new_mutation=False, t=tM)
            # main run
            else:
                # # record jump sizes during the rapid phase
                # if tM == time_to_record_jump_sizes_t_1_after_jump:
                #     self._population_class.record_jump_sizes_t_1_after_jump_no_efs_bins(indices)
                # record relative allelic contribution to short-term phenotypic adaptation
                if self.rate_of_shift is not None:
                    if tM == frozen_time + round(t_1):
                        self._population_class.record_contribution_to_phenotypic_change_recurrent_shifts_no_efs_bins(t=tM)
                        time_to_record_d2ax = tM
                # record heterozygosity in this generation
                if tM % self._BURN_TIME == 0:
                    self._population_class.record_heterozygosity_current_generation_no_efs_bins()
                # # record the state
                # if tM % self.N == 0:
                #     self._save_identifiers()
                #     self._population_class.record_mutations_and_phenotypic_distribution(nnm=nnm, directory=self.save_directory)
                if nnm < self.nnm:
                    nnm += self._population_class.next_gen_no_efs_bins(is_new_mutation=True, is_later_new_mutation=False, t=tM)
                else:
                    self._THERE_ARE_MUTANTS = self._population_class.next_gen_no_efs_bins(is_new_mutation=True, is_later_new_mutation=True, t=tM)

            tM += 1

        # # Record fixation probability
        # self._population_class.record_fixation_probability(nm=True)
        # # Record integral heterozygosity
        # self._population_class.record_integral_heterozygosity()

    def run_population_under_recurrent_shifts_complete_jumps(self):
        """
        To test whether the discrepancies between all-allele simulation results and diffusion approximation are because
        of incomplete jumps, we force jumps to be complete
        Returns:
            /
        """
        self._BURN_TIME = int(self.burn_time_N * self.N)
        # Whether there are still segregating mutations that are the first `self.nnm` new mutations
        self._THERE_ARE_MUTANTS = True
        # Count the number of new mutations
        nnm = 0
        # Count the number of generations
        tM = 0

        if self.rate_of_shift is not None:
            if self.sigma_0_del is None: # if sigma_0_del is not specified, then we use the analytical V_A
                self._t_1 = self.Vs / self._VAR_0 * math.log(self._shift)
            else:
                self._t_1 = self.Vs / self.sigma_0_del ** 2 * math.log(self._shift)
            self._rate_of_shift_single_allele = 1.0 / (1.0 / self.rate_of_shift - self._t_1) # 1/p_s = 1/p - t_1
            # time_to_record_jump_sizes_after_jump = -1
            waiting_time = np.random.geometric(p=self.rate_of_shift)

        while self._THERE_ARE_MUTANTS:
            # recurrent shifts
            if self.rate_of_shift is not None:
                if tM == waiting_time:
                    # burn-in
                    if tM < self._BURN_TIME:
                        direction_of_shift = np.random.binomial(n=1, p=0.5) * 2 - 1
                        self._population_class.shift_optimum(self._shift * direction_of_shift)
                        waiting_time += np.random.geometric(p=self.rate_of_shift)
                    # main run
                    else:
                        direction_of_shift = np.random.binomial(n=1, p=0.5) * 2 - 1
                        self._population_class.shift_optimum(self._shift * direction_of_shift)
                        # # record jump sizes
                        # self._population_class.record_jump_sizes_before_jump()
                        # time_to_record_jump_sizes_after_jump = waiting_time + round(self._t_1)

                        waiting_time += round(self._t_1) + np.random.geometric(p=self._rate_of_shift_single_allele)


            # burn-in
            if tM < self._BURN_TIME:
                self._population_class.next_gen(is_new_mutation=False, is_later_new_mutation=False)
            # main run
            else:
                # # record jump sizes
                # if tM == time_to_record_jump_sizes_after_jump:
                #     self._population_class.record_jump_sizes_after_jump(self._shift * direction_of_shift)

                if nnm < self.nnm:
                    nnm += self._population_class.next_gen(is_new_mutation=True, is_later_new_mutation=False)[1]
                else:
                    self._THERE_ARE_MUTANTS = self._population_class.next_gen(is_new_mutation=True, is_later_new_mutation=True)

            tM += 1

        # Record fixation probability
        self._population_class.record_fixation_probability(nm=True)



    def save_stats(self):
        self._save_identifiers()
        with open(os.path.join(self.save_directory, 'fixation_probability'), 'wb') as f:
            pickle.dump(self._population_class.stats(), f, protocol=pickle.HIGHEST_PROTOCOL)

    def save_stats_no_efs_bins(self):
        self._save_identifiers()
        with open(os.path.join(self.save_directory, 'fixations_and_extinctions'), 'wb') as f:
            pickle.dump(self._population_class.stats_no_efs_bins(), f, protocol=pickle.HIGHEST_PROTOCOL)

    def _save_identifiers(self):
        if not os.path.exists(self.save_directory):
            os.makedirs(self.save_directory)
        with open(os.path.join(self.save_directory, 'identifiers.txt'), 'w') as f:
            f.write('N: ' + str(self.N) + '\n')
            f.write('Vs: ' + str(self.Vs) + '\n')
            f.write('U: ' + str(self.U) + '\n')
            f.write('shape: ' + str(self.shape) + '\n')
            f.write('scale: ' + str(self.scale) + '\n')
            f.write('shift_s0: ' + str(self.shift_s0) + '\n')
            f.write('var_0: ' + str(self._VAR_0) + '\n')
            f.write('sigma_0_del: ' + str(self.sigma_0_del) + '\n')
            f.write('E2Ns: ' + str(self.E2Ns) + '\n')
            f.write('V2Ns: ' + str(self.V2Ns) + '\n')
            f.write('nE: ' + str(self.nE) + '\n')
            f.write('nF: ' + str(self.nF) + '\n')
            f.write('burn_time: ' + str(self._BURN_TIME) + '\n')
            f.write('burn_time_N: ' + str(self.burn_time_N) + '\n')
            f.write('rate_of_shift: ' + str(self.rate_of_shift) + '\n')
            f.write('nnm: ' + str(self.nnm) + '\n')

    def _var_0(self):
        """
        V_A(0) = \int_0^\infty 2NU v(a)g(a)da, where v(a) = 4a D_+(a/2)
        Returns:
            V_A(0)
        """
        S_dist = gamma(float(self.E2Ns) ** 2 / float(self.V2Ns), loc=0., scale=float(self.V2Ns) / float(self.E2Ns))
        return float(2 * self.N * self.U) * quad(lambda ss: 4.0 * np.sqrt(np.abs(ss)) * dawsn(
            np.sqrt(np.abs(ss)) / 2.0) * S_dist.pdf(ss), 0.0, S_dist.ppf(0.9999999999))[0]



class _Population(object):
    """
    The container class for the information about the population at every generation

    Attributes:
        N: (int) population size
        Vs: (float) fitness parameter
        _U: (float) mutation rate (per gamete per generation)
        _segregating: (set) set of mutations currently segregating in the population
        _frozen: (set) set of mutations segregating at the frozen time
        _FITNESS_OPTIMUM: (float) the phenotypic value that is optimal in fitness
        _effect_size_squared_dist: (distribution) the distribution of a^2, which is a Gamma distribution
        _dist_ess: (float) _FITNESS_OPTIMUM - mean phenotype. _ess means essential to record
        _var_ess: (float) phenotypic variance
        _mu3: (float) third central moment
        _mean_fixed_ess: (float) the contribution of fixed alleles to the mean phenotype
        _dist_ess_over_time: (List[float]) _dist_ess at every generation
        _var_ess_over_time: (List[float]) _var_ess at every generation
        _mu3_over_time: (List[float]) _mu3 at every generation
        _mean_fixed_ess_over_time: (List[float]) _mean_fixed_ess at every generation
        _index: (int) (currently unused)
    """
    def __init__(self, N, U, Vs, E2Ns, V2Ns, nE, nF, particular_2Ns, no_efs_bins):
        self.N = N
        # if V_S < 0, set it to be 2N
        if Vs < 0:
            self.Vs = float(2 * self.N)
        else:
            self.Vs = Vs
        self._U = U
        self.E2Ns = E2Ns
        # if V2Ns < 0, set it to be E2Ns^2, and hence the gamma distribution becomes exponential distribution with mean E2Ns
        if V2Ns < 0:
            self.V2Ns = float(self.E2Ns ** 2)
        else:
            self.V2Ns = V2Ns
        self._effect_size_squared_dist = gamma(float(self.E2Ns) ** 2 / float(self.V2Ns), loc=0.,
                                               scale=float(self.V2Ns) / float(self.E2Ns))
        self.nE = nE
        self.nF = nF
        self.particular_2Ns = particular_2Ns

        self._segregating = set()

        self._frozen = set()
        self._frozen_dist_ess = 0.0
        self._frozen_var_ess = 0.0

        self._FITNESS_OPTIMUM = 0.0

        self._dist_ess = 0.0
        self._var_ess = 0.0
        self._mu3 = 0.0
        self._mean_fixed_ess = 0.0
        self._dist_ess_over_time = list()
        self._var_ess_over_time = list()
        self._mu3_over_time = list()
        self._mean_fixed_ess_over_time = list()

        if not no_efs_bins:
            self._set_up_histograms()
        else:
            self._set_up_no_histograms()
            self._index = 0

    def _set_up_histograms(self):
        """
        Set up bins for histograms. Called on initialization.
        Returns:
            /
        """
        self._set_up_histograms_effects()
        # record jump sizes
        self._set_up_histograms_frequencies()

    # record jump sizes
    def _set_up_histograms_frequencies(self):
        """
        Set up frequency bins for histograms. Called on initialization. (currently unused)
        Returns:
            /
        """
        self.fbins = [i / self.nF for i in range(1, self.nF)]
        self._hist_fbins = dict()
        self._hist_fbins['jump_sizes'] = [[] for _ in range(self.nF)]

    def _set_up_histograms_effects(self):
        """
        Set up effect size bins for histograms. Called on initialization.
        Returns:
            /
        """
        self.effect_size_bins_pos = [math.sqrt(self._effect_size_squared_dist.ppf((i + 1) / (self.nE + 1))) for i in range(self.nE)]
        self.effect_size_bins = sorted([-a for a in self.effect_size_bins_pos] + [0.0] + self.effect_size_bins_pos)
        self.num_efs_bins = len(self.effect_size_bins) + 1
        self._hist_efs_bins = dict()
        self._hist_efs_bins['num_fixed'] = [0 for _ in range(self.num_efs_bins)]
        self._hist_efs_bins['frac_fixed'] = [0.0 for _ in range(self.num_efs_bins)]
        self._num_standing_variants_efs_bins = [0 for _ in range(self.num_efs_bins)]
        self._hist_efs_bins['num_fixed_nm'] = [0 for _ in range(self.num_efs_bins + 2 * len(self.particular_2Ns))]
        self._hist_efs_bins['frac_fixed_nm'] = [0.0 for _ in range(self.num_efs_bins + 2 * len(self.particular_2Ns))]
        self._num_new_mutations_efs_bins = [0 for _ in range(self.num_efs_bins + 2 * len(self.particular_2Ns))]
        self._hist_efs_bins['integral_heterozygosity'] = [[] for _ in range(self.num_efs_bins)] # in each generation
        self._hist_efs_bins['d2ax_t_1'] = [0.0 for _ in range(self.num_efs_bins)]

        # record jump sizes
        self._hist_efs_bins['jump_sizes'] = [[] for _ in range(self.num_efs_bins)]

        # to test the fixations of ups_and_downs mutations, all we need is mut.update_ups_and_downs, mut.is_ups_and_downs and mut.which_ups_and_downs
        # self._UPS_AND_DOWNS = ['111', '110', '101', '100', '011', '010', '001', '000']
        # self._hist_efs_bins['num_fixed_nm_ups_and_downs'] = dict() # for these 2, the keys are the 8 possibilities
        # self._num_ups_and_downs_new_mutations_efs_bins = dict()
        # for ups_and_downs in self._UPS_AND_DOWNS:
        #     self._hist_efs_bins['num_fixed_nm_ups_and_downs'][ups_and_downs] = [0 for _ in range(self.num_efs_bins)]
        #     self._num_ups_and_downs_new_mutations_efs_bins[ups_and_downs] = [0 for _ in range(self.num_efs_bins)]
        # self._NUPS = [5, 7]  # just for shift_s0=1.0 and 0.75
        # self._hist_efs_bins['num_fixed_nm_more_ups'] = dict() # for these 2, the keys are the numbers of ups that are of interest
        # self._num_more_ups_new_mutations_efs_bins = dict()
        # for nups in self._NUPS:
        #     self._hist_efs_bins['num_fixed_nm_more_ups'][nups] = [0 for _ in range(self.num_efs_bins)]
        #     self._num_more_ups_new_mutations_efs_bins[nups] = [0 for _ in range(self.num_efs_bins)]

        # self._hist_efs_bins['d2ax'] = [0.0 for _ in range(self.num_efs_bins)]

        # for '*_per_unit_mut_input' stats
        self.one_over_mutational_input_efs_bins_pos = [
            1.0 / float(2 * self.N * self._U * quad(self._effect_size_squared_dist.pdf, efsbinni, efsbinni_next)[0]) for
            efsbinni, efsbinni_next
            in
            zip([0.0] + self.effect_size_bins_pos,
                self.effect_size_bins_pos + [math.sqrt(self._effect_size_squared_dist.ppf(0.99999999999))])]
        self.one_over_mutational_input_efs_bins = [one_over for one_over in self.one_over_mutational_input_efs_bins_pos[
            ::-1]] + self.one_over_mutational_input_efs_bins_pos

    def _set_up_no_histograms(self):
        """
        Set up non-histograms. Called on initialization.
        _fixations_and_extinctions: fixed and extinct mutations
        _heterozygosity: the heterozygosity at every heterozygosity-recording generation
        _fixation_timings: the generations at which fixed mutations arise and fix
        _d2ax_t_1_scaled: the relative allelic contribution to short-term phenotypic change after a shift divided by (shift size / phenotypic variance)
        Returns:
            /
        """
        self._fixations_and_extinctions = list()
        self._heterozygosity = list()
        self._d2ax_t_1_scaled = list()
        self._fixation_timings = list()
        # self._jump_sizes_until_next_shift = list()
        # self._jump_sizes_t_1 = list()

    def _var_0(self):
        """
        V_A(0) = \int_0^\infty 2NU v(a)g(a)da, where v(a) = 4a D_+(a/2)
        Returns:
            V_A(0)
        """
        return float(2 * self.N * self._U) * quad(lambda ss: 4.0 * np.sqrt(np.abs(ss)) * dawsn(
            np.sqrt(np.abs(ss)) / 2.0) * self._effect_size_squared_dist.pdf(ss), 0.0, self._effect_size_squared_dist.ppf(0.9999999999999))[0]

    def var(self):
        return self._var_ess

    def freeze_standing_variants(self):
        """
        Record the number of standing variants in each bin of effect size and the frequency of each mutation at the time of the shift
        Returns:
            /
        """
        for mut in self._segregating:
            mut.record_frozen_freq()
            efsbinni = np.searchsorted(self.effect_size_bins, mut.a(minor=True, frozen=True))
            self._num_standing_variants_efs_bins[efsbinni] += 1

    def freeze_recurrent_shifts(self):
        """
        Record the frequency of each mutation segregating at the frozen time
        Returns:
            /
        """
        self._frozen_dist_ess = self._dist_ess
        self._frozen_var_ess = self._var_ess
        self._frozen = set()
        for mut in self._segregating:
            mut.record_frozen_freq_recurrent_shifts()
            self._frozen.add(mut)


    def record_mutations_and_phenotypic_distribution(self, nnm, directory):
        """
        Record the state of the simulation from time to time, in case the program crashes
        Args:
            directory: the directory to record
        Returns:
            /
        """
        if not os.path.exists(directory):
            os.makedirs(directory)
        with open(os.path.join(directory, 'record_mutations_and_phenotypic_distribution'), 'wb') as f:
            pickle.dump([(nnm, self._FITNESS_OPTIMUM, self._mean_fixed_ess), (self._hist_efs_bins['num_fixed_nm'],
            self._num_new_mutations_efs_bins)] + list(self._segregating), f, protocol=pickle.HIGHEST_PROTOCOL)

    def load_muts(self, directory):
        """
        Load the recorded state of the simulation
        Args:
            directory: the directory in which the state is recorded

        Returns:
            nnm: the nnm in run_population_under_recurrent_shifts
        """
        with open(os.path.join(directory, 'record_mutations_and_phenotypic_distribution'), 'rb') as f:
            muts = pickle.load(f)

        nnm, self._FITNESS_OPTIMUM, self._mean_fixed_ess = muts[0]
        self._hist_efs_bins['num_fixed_nm'], self._num_new_mutations_efs_bins = muts[1]
        for mut in muts[2:]:
            self._segregating.add(mut)
        self._update_essential_moments()
        return nnm


    def record_fixation_probability(self, nm):
        """
        Record the fixation probability of standing variants if nm is False, and of the first `self.nnm` new mutations if nm is True.
        Returns:
            /
        """
        if not nm:
            for efsbinni in range(self.num_efs_bins):
                num = float(self._num_standing_variants_efs_bins[efsbinni])
                if num > 0:
                    self._hist_efs_bins['frac_fixed'][efsbinni] = float(self._hist_efs_bins['num_fixed'][efsbinni]) / num
                else:
                    self._hist_efs_bins['frac_fixed'][efsbinni] = 0
        else:
            for efsbinni in range(self.num_efs_bins + 2 * len(self.particular_2Ns)):
                num = float(self._num_new_mutations_efs_bins[efsbinni])
                if num > 0:
                    self._hist_efs_bins['frac_fixed_nm'][efsbinni] = float(self._hist_efs_bins['num_fixed_nm'][efsbinni]) / num
                else:
                    self._hist_efs_bins['frac_fixed_nm'][efsbinni] = 0

    def record_contribution_to_phenotypic_change(self):
        """
        Only if a single shift, record the contribution to phenotypic change of standing variants.
        Returns:
            /
        """
        self._hist_efs_bins['d2ax_per_unit_mut_input'] = [
            self._hist_efs_bins['d2ax'][efsbinni] * self.one_over_mutational_input_efs_bins[efsbinni] for efsbinni in
            range(self.num_efs_bins)]
        self._hist_efs_bins['d2ax_scaled_per_unit_mut_input'] = [
            self._hist_efs_bins['d2ax_per_unit_mut_input'][efsbinni] * self._SPECIAL_SCALING_FACTOR for efsbinni in
            range(self.num_efs_bins)]

    def record_contribution_to_phenotypic_change_recurrent_shifts(self):
        """
        record the contribution to phenotypic change of mutations segregating at the frozen time
        Returns:
            /
        """
        for mut in self._frozen:
            efsbinni = np.searchsorted(self.effect_size_bins, mut.a())
            d2ax = mut.contribution_to_phenotypic_change_recurrent_shifts(frozen=True)
            self._hist_efs_bins['d2ax_t_1'][efsbinni] += d2ax

    def record_contribution_to_phenotypic_change_recurrent_shifts_no_efs_bins(self, t):
        """
        No effect size bins
        Args:
            t: (int) the current generation, in order to distinguish frozen times

        Returns:
            /
        """
        self._d2ax_t_1_scaled.append([])
        for mut in self._frozen:
            d2ax = mut.contribution_to_phenotypic_change_recurrent_shifts(frozen=True)
            self._d2ax_t_1_scaled[-1].append((mut.a(), d2ax / (self._frozen_dist_ess / self._frozen_var_ess), t))

    def record_integral_heterozygosity_current_generation(self):
        """
        record the integral heterozygosity in the current generation
        Returns:
            /
        """
        for efsbinni in range(self.num_efs_bins):
            assert len(self._hist_efs_bins['integral_heterozygosity'][efsbinni]) == len(
                self._hist_efs_bins['integral_heterozygosity'][0]) # which equals the current generation
        for efsbinni in range(self.num_efs_bins):
            self._hist_efs_bins['integral_heterozygosity'][efsbinni].append(0.0)
        for mut in self._segregating:
            efsbinni = np.searchsorted(self.effect_size_bins, mut.a())
            self._hist_efs_bins['integral_heterozygosity'][efsbinni][-1] += 2.0 * mut.x() * (1.0 - mut.x())

    def record_integral_heterozygosity(self):
        """
        calculate the mean and standard error of integral heterozygosity over generations after burning
        Returns:
            /
        """
        for efsbinni in range(self.num_efs_bins):
            assert len(self._hist_efs_bins['integral_heterozygosity'][efsbinni]) == len(
                self._hist_efs_bins['integral_heterozygosity'][0]) # which equals the total number of generations after burning
        for efsbinni in range(self.num_efs_bins):
            integral_heterozygosity_mean = np.mean(self._hist_efs_bins['integral_heterozygosity'][efsbinni])
            integral_heterozygosity_se = np.std(self._hist_efs_bins['integral_heterozygosity'][efsbinni]) / np.sqrt(
                np.size(self._hist_efs_bins['integral_heterozygosity'][efsbinni]))
            self._hist_efs_bins['integral_heterozygosity'][efsbinni] = [integral_heterozygosity_mean, integral_heterozygosity_se]

    def record_heterozygosity_current_generation_no_efs_bins(self):
        """
        No effect size bins
        Returns:
            /
        """
        self._heterozygosity.append(tuple([(mut.a(), mut.x()) for mut in self._segregating]))

    def record_jump_sizes_before_jump(self):
        for mut in self._segregating:
            mut.x_for_jump_sizes.append([mut.a(), mut.x(), self._dist_ess])

    def record_jump_sizes_after_jump(self, shift):
        for mut in self._segregating:
            if len(mut.x_for_jump_sizes) > 0:
                assert len(mut.x_for_jump_sizes[-1]) == 3
                mut.x_for_jump_sizes[-1].append(mut.x())
                # D(t_0) == \Lambda only
                if abs(mut.x_for_jump_sizes[-1][2] - shift) < 1.0: # \delta
                    efsbinni = np.searchsorted(self.effect_size_bins, mut.a())
                    self._hist_efs_bins['jump_sizes'][efsbinni].append(mut.x_for_jump_sizes[-1])
                    fbinni = np.searchsorted(self.fbins, mut.x_for_jump_sizes[-1][1])
                    self._hist_fbins['jump_sizes'][fbinni].append(mut.x_for_jump_sizes[-1])

    def record_jump_sizes_before_jump_no_efs_bins(self):
        for mut in self._segregating:
            mut.jump_sizes_until_next_shift.append([mut.a(), mut.x(), self._dist_ess])
            mut.jump_sizes_t_1.append([mut.a(), mut.x(), self._dist_ess])
        return [mut.index() for mut in self._segregating]

    def record_jump_sizes_until_next_shift_after_jump_no_efs_bins(self, indices, duration):
        for mut in self._segregating:
            if mut.index() in indices:
                if len(mut.jump_sizes_until_next_shift) > 0:
                    assert len(mut.jump_sizes_until_next_shift[-1]) == 3
                    mut.jump_sizes_until_next_shift[-1].extend([mut.x(), duration])
                    self._jump_sizes_until_next_shift.append(mut.jump_sizes_until_next_shift[-1])

    def record_jump_sizes_t_1_after_jump_no_efs_bins(self, indices):
        for mut in self._segregating:
            if mut.index() in indices:
                if len(mut.jump_sizes_t_1) > 0:
                    assert len(mut.jump_sizes_t_1[-1]) == 3
                    mut.jump_sizes_t_1[-1].append(mut.x())
                    self._jump_sizes_t_1.append(mut.jump_sizes_t_1[-1])

    def stats(self):
        """
        The stats to record after the simulation finishes running
        Returns:
            _hist_efs_bins: fixation probability of new mutations and number of fixed new mutations in each effect size bin
            _num_new_mutations_efs_bins: number of new mutations in each effect size bin
            _dist_ess_over_time: _dist_ess at every generation
            _var_ess_over_time: _var_ess at every generation
            _mu3_over_time: _mu3 at every generation
            _hist_fbins: (currently unused)
        """
        return self._hist_efs_bins, self._num_new_mutations_efs_bins, self._dist_ess_over_time, self._var_ess_over_time, self._mu3_over_time, self._hist_fbins

    def stats_no_efs_bins(self):
        """
        No effect size bins
        Returns:
            _fixations_and_extinctions: fixed and extinct mutations
            _heterozygosity: the heterozygosity at every heterozygosity-recording generation
            _dist_ess_over_time: _dist_ess at every generation
            _var_ess_over_time: _var_ess at every generation
            _mu3_over_time: _mu3 at every generation
            _mean_fixed_ess_over_time: _mean_fixed_ess at every generation
            _fixation_timings: the generations at which fixed mutations arise and fix
            _d2ax_t_1_scaled: the relative allelic contribution to short-term phenotypic change after a shift divided by (shift size / phenotypic variance)
        """
        return self._fixations_and_extinctions, self._heterozygosity, self._dist_ess_over_time, self._var_ess_over_time, self._mu3_over_time, self._mean_fixed_ess_over_time, self._fixation_timings, self._d2ax_t_1_scaled

    def shift_optimum(self, shift):
        """
        Shift the fitness optimum by shift
        Args:
            shift: shift size
        Returns:
            /
        """
        self._FITNESS_OPTIMUM += shift
        self._SPECIAL_SCALING_FACTOR = self._var_0() / float(shift)
        self._update_essential_moments()
        self._update_aligned_mutations(direction_of_shift=np.sign(shift))

    def _update_aligned_mutations(self, direction_of_shift):
        """
        Update the 'up' status of each mutation when it experiences its first 3 shifts.
        Args:
            direction_of_shift: (int, 1 or -1) the direction of this shift, i.e., the first shift

        Returns:
            /
        """
        for mut in self._segregating:
            if mut.is_new_mutation() and not mut.is_later_new_mutation():
                mut.update_ups_and_downs(update_ups_and_downs=(mut.a() * direction_of_shift > 0))



    def next_gen(self, is_new_mutation, is_later_new_mutation):
        """
        Advance a generation
        Args:
            is_new_mutation (boolean): whether the mutations arise after burn-in, i.e., new mutations.
            is_later_new_mutation (boolean): whether there are already no less than `self.nnm` new mutations.
        Returns (possibly):
            (boolean) whether there are still segregating mutations that are standing variants
            (int) the number of new mutations that arise in this generation
            (boolean) whether there are still segregating mutations that are the first `self.nnm` new mutations
        """
        self._wright_fisher()
        # count the number nnm_this_gen of mutations that arise in this generation, which is used if is_new_mutation is True and is_later_new_mutation is False.
        nmut = len(self._segregating)
        self._de_novo_mutations(is_new_mutation=is_new_mutation, is_later_new_mutation=is_later_new_mutation)
        nnm_this_gen = len(self._segregating) - nmut
        self._remove_extinct_fixed_from_seg(is_new_mutation=is_new_mutation)
        self._update_essential_moments()

        # returns
        # during the burning period, there is no return.
        if not is_new_mutation:
            return
        # then is_new_mutation is True.
        # if not is_later_new_mutation means either we are interested in standing variants or there are less than `self.nnm` new mutations, and so we return 2 variables: whether there are still segregating mutations that are standing variants, and the number of new mutations that arise in this generation.
        if not is_later_new_mutation:
            for mut in self._segregating:
                if not mut.is_new_mutation():
                    return True, nnm_this_gen
            return False, nnm_this_gen
        # then is_new_mutation is True and is_later_new_mutation is True. Return whether there are segregating mutations that are the first `self.nnm` new mutations.
        for mut in self._segregating:
            if mut.is_new_mutation() and not mut.is_later_new_mutation():
                return True
        return False

    def _de_novo_mutations(self, is_new_mutation, is_later_new_mutation):
        """
        Add de novo mutations into self._segregating. #DNM ~ Poisson(2NU)
        Returns:
            /
        """
        for _ in range(np.random.poisson(2 * self.N * self._U)):
            mut = self._get_mutation(is_new_mutation=is_new_mutation, is_later_new_mutation=is_later_new_mutation)
            self._segregating.add(mut)

    def _get_mutation(self, is_new_mutation, is_later_new_mutation):
        """
        Return a de novo mutation with effect size sampled from a Gamma distribution and simple random effect direction.
        Also, with a certain probability to get particular effect sizes.
        Returns:
            mut: an instance of class Mutation
        """
        p = np.random.random()
        if int(p // 0.1) < len(self.particular_2Ns):
            abs_a = np.sqrt(self.particular_2Ns[int(p // 0.1)])
        else:
            abs_a = np.sqrt(np.random.gamma(shape=float(self.E2Ns ** 2) / float(self.V2Ns), scale=float(self.V2Ns) / float(self.E2Ns)))
        derived_sign = np.random.binomial(n=1, p=0.5) * 2 - 1
        a = abs_a * derived_sign
        mut = Mutation(a=a, x=1.0 / float(2 * self.N), is_new_mutation=is_new_mutation, is_later_new_mutation=is_later_new_mutation, N=self.N)
        if is_new_mutation and not is_later_new_mutation:
            if int(p // 0.1) < len(self.particular_2Ns):
                if derived_sign > 0:
                    efsbinni = self.num_efs_bins + len(self.particular_2Ns) + int(p // 0.1)
                else:
                    efsbinni = self.num_efs_bins + len(self.particular_2Ns) - 1 - int(p // 0.1)
            else:
                efsbinni = np.searchsorted(self.effect_size_bins, mut.a())
            self._num_new_mutations_efs_bins[efsbinni] += 1
        return mut

    def next_gen_no_efs_bins(self, is_new_mutation, is_later_new_mutation, t):
        """
        No effect size bins
        Args:
            t: (int) the current generation, in order to record self._fixation_timings
        """
        self._wright_fisher()
        # count the number nnm_this_gen of mutations that arise in this generation. Used if is_new_mutation and not is_later_new_mutation
        nmut = len(self._segregating)
        self._de_novo_mutations_no_efs_bins(is_new_mutation=is_new_mutation, is_later_new_mutation=is_later_new_mutation, time_arising=t)
        nnm_this_gen = len(self._segregating) - nmut
        self._remove_extinct_fixed_from_seg_no_efs_bins(time_fixing=t)
        self._update_essential_moments()

        # returns
        if not is_new_mutation:
            return
        if not is_later_new_mutation:
            return nnm_this_gen
        for mut in self._segregating:
            if mut.is_new_mutation() and not mut.is_later_new_mutation():
                return True
        return False

    def _de_novo_mutations_no_efs_bins(self, is_new_mutation, is_later_new_mutation, time_arising):
        """
        No effect size bins
        Args:
            time_arising: (int) the generation at which these mutations arise, in order to record self._fixation_timings
        """
        for _ in range(np.random.poisson(2 * self.N * self._U)):
            mut = self._get_mutation_no_efs_bins(is_new_mutation=is_new_mutation, is_later_new_mutation=is_later_new_mutation, time_arising=time_arising)
            self._segregating.add(mut)

    def _get_mutation_no_efs_bins(self, is_new_mutation, is_later_new_mutation, time_arising):
        """
        No effect size bins
        Args:
            time_arising: (int) the generation at which this mutation arises, in order to record self._fixation_timings
        """
        abs_a = np.sqrt(np.random.gamma(shape=float(self.E2Ns ** 2) / float(self.V2Ns), scale=float(self.V2Ns) / float(self.E2Ns)))
        derived_sign = np.random.binomial(n=1, p=0.5) * 2 - 1
        a = abs_a * derived_sign
        mut = Mutation(a=a, x=1.0 / float(2 * self.N), is_new_mutation=is_new_mutation, is_later_new_mutation=is_later_new_mutation, N=self.N, index=self._index, time_arising=time_arising)
        self._index += 1
        return mut

    def _remove_extinct_fixed_from_seg(self, is_new_mutation):
        """
        Remove fixed and extinct mutations from self._segregating.
        Returns:
            /
        """
        fixed_and_extinct_mutants = set()
        MINOR = True
        FROZEN = True

        for mut in self._segregating:
            if mut.is_fixed():
                self._mean_fixed_ess += 2 * mut.a()
            # if is_new_mutation and not mut.is_new_mutation() and mut.is_fixed(minor=MINOR, frozen=FROZEN):
            #     efsbinni = np.searchsorted(self.effect_size_bins, mut.a(minor=MINOR, frozen=FROZEN))
            #     self._hist_efs_bins['num_fixed'][efsbinni] += 1
            if mut.is_fixed() or mut.is_extinct():
                fixed_and_extinct_mutants.add(mut)
                # if is_new_mutation and not mut.is_new_mutation():
                #     efsbinni = np.searchsorted(self.effect_size_bins, mut.a(minor=MINOR, frozen=FROZEN))
                #     d2ax = mut.contribution_to_phenotypic_change(frozen=FROZEN)
                #     self._hist_efs_bins['d2ax'][efsbinni] += d2ax
            if mut.is_new_mutation() and not mut.is_later_new_mutation() and mut.is_fixed():
                if abs(mut.a()) in [np.sqrt(s) for s in self.particular_2Ns]:
                    if mut.a() > 0:
                        efsbinni = self.num_efs_bins + len(self.particular_2Ns) + [np.sqrt(s) for s in self.particular_2Ns].index(abs(mut.a()))
                    else:
                        efsbinni = self.num_efs_bins + len(self.particular_2Ns) - 1 - [np.sqrt(s) for s in self.particular_2Ns].index(abs(mut.a()))
                else:
                    efsbinni = np.searchsorted(self.effect_size_bins, mut.a())
                self._hist_efs_bins['num_fixed_nm'][efsbinni] += 1
                # if mut.is_ups_and_downs():
                #     self._hist_efs_bins['num_fixed_nm_ups_and_downs'][mut.which_ups_and_downs()][efsbinni] += 1
                # for nups in self._NUPS:
                #     if mut.is_all_ups(nups):
                #         self._hist_efs_bins['num_fixed_nm_more_ups'][nups][efsbinni] += 1

        for mut in fixed_and_extinct_mutants:
            # ups and downs
            # if mut.is_new_mutation() and not mut.is_later_new_mutation():
            #     if mut.is_ups_and_downs():
            #         efsbinni = np.searchsorted(self.effect_size_bins, mut.a())
            #         self._num_ups_and_downs_new_mutations_efs_bins[mut.which_ups_and_downs()][efsbinni] += 1
            #     for nups in self._NUPS:
            #         if mut.is_all_ups(nups):
            #             self._num_more_ups_new_mutations_efs_bins[nups][efsbinni] += 1

            self._segregating.remove(mut)

    def _remove_extinct_fixed_from_seg_no_efs_bins(self, time_fixing):
        """
        No effect size bins
        Args:
            time_fixing: (int) the generation at which fixed mutations fix, in order to record self._fixation_timings
        """
        fixed_and_extinct_mutants = set()
        for mut in self._segregating:
            if mut.is_fixed():
                self._mean_fixed_ess += 2 * mut.a()
                # record the time of arising and the time of fixing
                if mut.is_new_mutation() and not mut.is_later_new_mutation():
                    self._fixation_timings.append((mut.time_arising(), time_fixing))
            if mut.is_fixed() or mut.is_extinct():
                fixed_and_extinct_mutants.add(mut)
                # record the fixation status of every mutation
                if mut.is_new_mutation() and not mut.is_later_new_mutation():
                    self._fixations_and_extinctions.append((mut.a(), mut.is_fixed()))

        for mut in fixed_and_extinct_mutants:
            self._segregating.remove(mut)

    def _update_essential_moments(self):
        """
        Update self._dist_ess, self._var_ess, and self._mu3
        mean phenotype = \Sigma 2ax = \Sigma_segregating 2ax + \Sigma_fixed 2a
        phenotypic variance = \Sigma 2a^2 x(1 - x) = \Sigma_segregating 2a^2 x(1 - x)
        Returns:
            /
        """
        mean_seg = 0.0
        var = 0.0
        mu3 = 0.0
        for mut in self._segregating:
            mean_seg += 2 * mut.a() * mut.x()
            var += 2 * mut.a() ** 2 * mut.x() * (1.0 - mut.x())
            mu3 += 4 * mut.a() ** 3 * mut.x() * (1.0 - mut.x()) * (0.5 - mut.x())
        mean = mean_seg + self._mean_fixed_ess
        self._dist_ess = self._FITNESS_OPTIMUM - mean
        self._var_ess = var
        self._mu3 = mu3
        self._dist_ess_over_time.append(self._dist_ess)
        self._var_ess_over_time.append(self._var_ess)
        self._mu3_over_time.append(self._mu3)
        self._mean_fixed_ess_over_time.append(self._mean_fixed_ess)

    def _wright_fisher(self):
        """
        Update frequencies of currently segregating mutations according to Wright-Fisher
        Returns:
            /
        """
        for mut in self._segregating:
            # the average fitness of individuals with i (= 0, 1 or 2) copies of this mutation
            w_0 = math.exp(-(-self._dist_ess - 2.0 * mut.a() * mut.x()) ** 2 / (2.0 * self.Vs))
            w_1 = math.exp(-(-self._dist_ess - 2.0 * mut.a() * mut.x() + mut.a()) ** 2 / (2.0 * self.Vs))
            w_2 = math.exp(-(-self._dist_ess - 2.0 * mut.a() * mut.x() + 2.0 * mut.a()) ** 2 / (2.0 * self.Vs))
            # the overall average fitness
            w = mut.x() ** 2 * w_2 + 2.0 * mut.x() * (1.0 - mut.x()) * w_1 + (1.0 - mut.x()) ** 2 * w_0

            # update the frequency of this mutation
            expected_new_freq = (mut.x() ** 2 * w_2 + mut.x() * (1.0 - mut.x()) * w_1) / w
            mut.update_freq(np.random.binomial(n=2 * self.N, p=expected_new_freq) / float(2 * self.N))
