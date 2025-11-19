"""
Class for mutations
"""

class Mutation(object):
    """
    Attributes:
        _x: (float) the frequency
        _a: (float) the effect size
        _is_new_mutation: (boolean) whether the mutation is a new mutation, i.e., arises after burn-in
        _is_later_new_mutation: (boolean) whether there are already sufficient number of new mutations
        _ups_and_downs: (a list of 3 boolean variables) only for _is_new_mutation and not _is_later_new_mutation. each element is True if the 1st/2nd/3rd shift that the mutation experiences is aligned with its effect direction; False if opposing or no shifts (currently unused)
        _ups_and_downs_updated: (a list of 3 boolean variables) whether each element of _3ups has been updated
        _N: (int) the population size
        _frozen_freq: (float) the frequency at the time of the shift
        _frozen_freq_recurrent_shifts: (float) the frequency at the frozen time
        lifetime: (int) the number of generations since the mutation arose
    """
    def __init__(self, x, a, is_new_mutation, is_later_new_mutation, N, index, time_arising):
        self._x = x
        self._a = a
        self._is_new_mutation = is_new_mutation
        self._is_later_new_mutation = is_later_new_mutation
        self._ups_and_downs = [False, False, False]  # each element turns into True when the corresponding condition is satisfied
        self._ups_and_downs_updated = [False, False, False]
        self._more_ups = []
        self._N = N
        self.jump_sizes_until_next_shift = []
        self.jump_sizes_t_1 = []

        self.lifetime = 0

        self._index = index

        self._time_arising = time_arising

    def x(self, minor=False, frozen=False):
        """
        If minor, then return the current frequency of the minor allele at freeze time if frozen is True or at present if frozen is False; if not minor, then return the frequency of the allele itself.
        Args:
            minor and frozen: (boolean) whether we are interested in the minor allele at freeze time

        Returns:
            (float) the current frequency of either the minor allele at freeze time or the allele itself
        """
        if minor:
            if frozen:
                if self._frozen_freq > 0.5:
                    return 1.0 - self._x
            else:
                if self._x > 0.5:
                    return 1.0 - self._x
        return self._x

    def a(self, minor=False, frozen=False):
        """
        if minor and frozen, then return the effect size of the minor allele at freeze time; if not minor or frozen, then return the effect size of the allele itself.
        Args:
            minor: (boolean) whether we are interested in the minor allele
            frozen: (boolean) whether we are interested in the frozen allele

        Returns:
            (float) the effect size of either the minor allele at freeze time or the allele itself
        """
        if minor and frozen:
            if self._frozen_freq > 0.5:
                return -self._a
        return self._a

    def record_frozen_freq(self):
        self._frozen_freq = self._x

    def record_frozen_freq_recurrent_shifts(self):
        self._frozen_freq_recurrent_shifts = self._x

    def update_freq(self, x):
        self._x = x
        self.lifetime += 1

    def update_ups_and_downs(self, update_ups_and_downs):
        if not self._ups_and_downs_updated[0]:
            self._ups_and_downs_updated[0] = True
            if update_ups_and_downs:
                self._ups_and_downs[0] = True
        elif not self._ups_and_downs_updated[1]:
            self._ups_and_downs_updated[1] = True
            if update_ups_and_downs:
                self._ups_and_downs[1] = True
        elif not self._ups_and_downs_updated[2]:
            self._ups_and_downs_updated[2] = True
            if update_ups_and_downs:
                self._ups_and_downs[2] = True

        self._more_ups.append(update_ups_and_downs)

    def delta_x(self, minor=False, frozen=False):
        """
        Return the frequency change of the previous or frozen minor allele if minor is True or of the allele itself if minor is False; since the frozen time if frozen is True or since the previous generation if frozen is False.
        Args:
            minor: (boolean) whether we are interested in the minor allele
            frozen: (boolean) whether we are interested in the frozen allele

        Returns:
            the frequency change
        """
        if frozen:
            if not minor:
                return self._x - self._frozen_freq

    def delta_x_recurrent_shifts(self, minor, frozen):
        """
        Return the frequency change of the (minor) allele since the frozen time if frozen is True
        Args:
            minor:
            frozen:

        Returns:
            the frequency change
        """
        if frozen:
            if not minor:
                return self._x - self._frozen_freq_recurrent_shifts

    def is_new_mutation(self):
        return self._is_new_mutation

    def is_later_new_mutation(self):
        return self._is_later_new_mutation

    def index(self):
        return self._index

    def time_arising(self):
        return self._time_arising

    def is_ups_and_downs(self):
        return self._ups_and_downs_updated == [True, True, True]

    def is_all_ups(self, nups):
        if len(self._more_ups) < nups:
            return False
        for up in self._more_ups[:nups]:
            if not up:
                return False
        return True

    def which_ups_and_downs(self):
        assert self.is_ups_and_downs()
        if self._ups_and_downs == [True, True, True]:
            return '111'
        elif self._ups_and_downs == [True, True, False]:
            return '110'
        elif self._ups_and_downs == [True, False, True]:
            return '101'
        elif self._ups_and_downs == [True, False, False]:
            return '100'
        elif self._ups_and_downs == [False, True, True]:
            return '011'
        elif self._ups_and_downs == [False, True, False]:
            return '010'
        elif self._ups_and_downs == [False, False, True]:
            return '001'
        elif self._ups_and_downs == [False, False, False]:
            return '000'

    def is_fixed(self, minor=False, frozen=False):
        """
        Whether this mutation fixes
        Args:
            minor: (boolean) whether we are interested in the minor allele
            frozen: (boolean) whether we are interested in the frozen allele

        Returns:
            (boolean) if minor and frozen, then True if the minor allele at freeze time fixes; if not minor or frozen, then True if the allele itself fixes
        """
        if minor and frozen:
            if self._frozen_freq > 0.5:
                if self._x == 0.0:
                    return True
                else:
                    return False
        if self._x == 1.0:
            return True
        else:
            return False

    def is_extinct(self, minor=False, frozen=False):
        """
        Whether this mutation goes extinct
        Args:
            minor: (boolean) whether we are interested in the minor allele
            frozen: (boolean) whether we are interested in the frozen allele

        Returns:
            (boolean) if minor and frozen, then True if the minor allele at freeze time goes extinct; if not minor or frozen, then True if the allele itself goes extinct
        """
        if minor and frozen:
            if self._frozen_freq > 0.5:
                if self._x == 1.0:
                    return True
                else:
                    return False
        if self._x == 0.0:
            return True
        else:
            return False

    def contribution_to_phenotypic_change(self, frozen=False):
        """
        Return the contribution to phenotypic change since the previous generation if frozen is False, or the frozen time if frozen is True.
        Args:
            frozen: (boolean) if frozen, since the frozen time; if not frozen, since the previous generation

        Returns:
            contribution to phenotypic change
        """
        return self.delta_x(minor=False, frozen=frozen) * 2.0 * self._a

    def contribution_to_phenotypic_change_recurrent_shifts(self, frozen):
        """
        Return the contribution to phenotypic change since the frozen time if frozen is True.
        Args:
            frozen:

        Returns:
            contribution to phenotypic change
        """
        return self.delta_x_recurrent_shifts(minor=False, frozen=frozen) * 2.0 * self._a
