import numpy as np
from scipy import special

from rna_voila.constants import MINVAL
from rna_voila.vlsv import get_expected_dpsi, matrix_area, get_expected_psi
from rna_voila.voila_log import voila_log



class EventDescription:

    @classmethod
    def a5ss(cls, evt_desc):
        """
        Using lsv type, does this lsv have 5 prime splice sites.
        :return: boolean
        """
        return 'A5SS' in [cls.reference_exon_ss(evt_desc), cls.other_exons_ss(evt_desc)]

    @classmethod
    def a3ss(cls, evt_desc):
        """
        Using lsv type, does this lsv have 3 prime splice sites.
        :return: boolean
        """
        return 'A3SS' in [cls.reference_exon_ss(evt_desc), cls.other_exons_ss(evt_desc)]

    @classmethod
    def reference_exon_ss(cls, evt_desc):
        """
        Check for 3 prime or 5 prime splice sites in reference exon.
        :return: list of strings
        """
        try:

            ss = filter(lambda x: x != 'i', evt_desc.split('|')[1:])
            ss = map(lambda x: x.split('.')[0].split('e')[0], ss)

            if len(set(ss)) > 1:
                if evt_desc[0] == 's':
                    return 'A5SS'
                else:
                    return 'A3SS'

        except IndexError:
            return 'na'

    @classmethod
    def other_exons_ss(cls, evt_desc):
        """
        Find 3 prime or 5 prime splice sites in exons that aren't the reference exon.
        :return: List of strings
        """
        try:

            ss = filter(lambda lt: lt != 'i', evt_desc.split('|')[1:])
            exons = {}
            for x in ss:
                exon = x.split('.')[0].split('e')[1]
                ss = x.split('.')[1]
                try:
                    exons[exon].add(ss)
                except KeyError:
                    exons[exon] = {ss}

            if any(len(values) > 1 for values in exons.values()):
                if evt_desc[0] == 's':
                    return 'A3SS'
                else:
                    return 'A5SS'

        except IndexError:
            return 'na'

    @classmethod
    def exon_skipping(cls, evt_desc):
        """
        Using lsv type, does this lsv have exon skipping.
        :return: boolean
        """
        try:
            return cls.exon_count(evt_desc) > 2
        except TypeError:
            return False

    @classmethod
    def exon_count(cls, evt_desc):
        """
        Using lsv type, how many exons are in this lsv.
        :return: integer
        """
        try:
            exons = filter(lambda x: x != 'i', evt_desc.split('|')[1:])
            exons = map(lambda x: x.split('.')[0].split('e')[1], exons)
            return len(set(exons)) + 1
        except IndexError:
            return 'na'


def unpack_means(value):
    """
    In the case where lsv is binary, we need to generate the other junction data.
    :param value: means
    :return: list
    """
    if np.size(value, 0) == 1:
        value = np.append(value, np.array(1 - value[0]))
    return np.nan_to_num(value).tolist()


def unpack_bins(value):
    """
    In the case where lsv is binary, we need to generate the other junction data.
    :param value: bins
    :return: list
    """
    if np.size(value, 0) == 1:
        value = np.append(value, [np.flip(value[-1], 0)], axis=0)
    return np.nan_to_num(value).tolist()


def generate_excl_incl(means):
    """
    Create exclusion and inclusion values for plots.
    :param means: lsv means data
    :return: list
    """
    l = []
    for mean in means:
        mean = np.nan_to_num(mean)
        if mean < 0:
            l.append((-mean, 0))
        else:
            l.append((0, mean))
    return l


def generate_means(bins):
    """
    Create means where not available in voila file.
    :param bins: bins data
    :return: list
    """
    m = []
    for b in bins:
        m.append(np.nan_to_num(get_expected_dpsi(b)))
    return m


def generate_prior_removed_expected_dpsi(ir, prior, bins=None, pr_removed_bins=None):
    """
    Similar input / output use case to the function below, but gets expected dpsi instead of matrix area (prob)
    """
    x = []

    for pr_removed_bin in generate_bins_prior_removed(ir, prior, bins) if bins is not None else pr_removed_bins:
        x.append(abs(get_expected_dpsi(pr_removed_bin)))

    return x

def generate_high_probability_non_changing(ir, prior, non_changing_threshold, bins=None, pr_removed_bins=None):
    """
    Calculate the probability of non changing lsv junctions.
    :param ir: Does this lsv have intron retention.
    :param prior: prior matrix from rna_voila file.
    :param non_changing_threshold: non-changing threshold set by user.
    :param bins: bins data from rna_voila file.
    :return: list
    """
    x = []

    for pr_removed_bin in generate_bins_prior_removed(ir, prior, bins) if bins is not None else pr_removed_bins:
        x.append(matrix_area(pr_removed_bin, non_changing_threshold, non_changing=True))

    return x

def generate_bins_prior_removed(ir, prior, bins):
    new_bins = []
    prior = prior[1 if ir else 0]

    for bin in bins:
        bin = np.array(bin)
        bin += MINVAL
        bin /= bin.sum()
        A = np.log(bin) - prior
        R = np.exp(A - special.logsumexp(A))
        new_bins.append(R)

    return new_bins


def generate_standard_deviations(bins):
    """
    Calulate standard deviations for lsv junctions.
    :param bins: bins data from rna_voila file.
    :return: list
    """
    v = np.sqrt(generate_variances(bins))
    return v

def generate_variances(bins):
    """
    Calulate variances for lsv junctions.
    :param bins: bins data from rna_voila file.
    :return: list
    """
    v = []
    for b in bins:
        epsi = get_expected_psi(b)
        # b used to be a numpy array and now it's a list...
        step_bins = 1.0 / len(b)
        projection_prod = b * np.arange(step_bins / 2, 1, step_bins) ** 2
        var = np.sum(projection_prod) - (epsi ** 2)
        if var < 0:
            v.append(0.0)
            if abs(var) < step_bins:
                # if var is < 0 but also of less magnitude then the bin size, we assume a rounding error and
                # pass along a zero variance
                pass
            else:
                # otherwise, this should not happen, so we send a warning
                log = voila_log()
                log.warning("Negative variance found when calculating over %s" % str(b))
        else:
            v.append(var)


    return v
