import operator
from abc import ABC
from functools import reduce
from itertools import chain

import numpy as np
import scipy.stats

import rna_voila.config
#from rna_voila.api.matrix_hdf5 import DeltaPsi, Psi, Heterogen, MatrixType
from rna_voila.api.matrix_utils import unpack_means, unpack_bins, generate_excl_incl, generate_means, \
    generate_high_probability_non_changing, generate_variances, generate_standard_deviations

from rna_voila.exceptions import LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile, LsvIdNotFoundInAnyVoilaFile
from rna_voila.vlsv import is_lsv_changing, matrix_area, get_expected_psi
from rna_voila.api.matrix_utils import EventDescription
from rna_voila import constants
from multiprocessing import Pool
from itertools import combinations
import rna_majiq as nm
from collections import namedtuple
import numpy as np

def preconfig_group_names_psi(cov_file):
    return nm.PsiCoverage.from_zarr(cov_file).prefixes

def preconfig_group_names_dpsi(cov_file):
    return nm.DeltaPsiDataset.from_zarr(cov_file).prefixes

def preconfig_group_names_het(cov_file):
    return nm.HeterogenDataset.from_zarr(cov_file).prefixes

def get_lsvid2lsvidx(sg_zarr, cov_zarr, append_to):
    lsvid2lsvidx = append_to
    if type(cov_zarr) is dict:
        for k, v in cov_zarr.items():
            lsvid2lsvidx[k] = {}
            events = v.get_events(sg_zarr.introns, sg_zarr.junctions)
            lsv_ids = sg_zarr.exon_connections.event_id(events.ref_exon_idx, events.event_type)
            for lsv_idx, lsv_id in enumerate(lsv_ids):
                lsvid2lsvidx[k][lsv_id] = lsv_idx
    else:
        events = cov_zarr.get_events(sg_zarr.introns, sg_zarr.junctions)
        lsv_ids = sg_zarr.exon_connections.event_id(events.ref_exon_idx, events.event_type)
        for lsv_idx, lsv_id in enumerate(lsv_ids):
            lsvid2lsvidx[lsv_id] = lsv_idx

    return lsvid2lsvidx

def get_lsvidx2lsvid(sg_zarr, cov_zarr, append_to):
    # TODO this does not append at the moment ; check whether there will ever be any discrepancy in
    # different cov files in practice
    events = cov_zarr.get_events(sg_zarr.introns, sg_zarr.junctions)
    lsv_ids = sg_zarr.exon_connections.event_id(events.ref_exon_idx, events.event_type)
    return lsv_ids


LSVTypeData = namedtuple('LSVTypeData', 'target binary a5ss a3ss exon_skipping')
def get_lsvtype_cache(sg_zarr, cov_zarr):
    res = []
    if type(cov_zarr) is dict:
        cov_zarr = next(iter(cov_zarr.values()))
    events = cov_zarr.get_events(sg_zarr.introns, sg_zarr.junctions)
    event_descriptions = sg_zarr.exon_connections.event_description(events.ref_exon_idx, events.event_type)
    for e_idx, event_description in enumerate(event_descriptions):
        td = LSVTypeData(
            events.event_type[e_idx] == b"t",
            events.event_size[e_idx] == 2,
            EventDescription.a5ss(event_description),
            EventDescription.a3ss(event_description),
            EventDescription.exon_skipping(event_description)
        )
        res.append(td)

    return res


def get_matrix_format_str():
    config = rna_voila.config.ViewConfig()

    if config.voila_file is not None:
        return 'h'
    elif config.majiq_file is not None:
        return 'z'
    raise NotImplementedError("Invalid Matrix File Format")


"""
View matrix using zarr also implicitly requires splicegraph connection!
"""

class ViewMatrix(ABC):




    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def lsv(self, lsv_id):
        raise NotImplementedError()

    def lsvs(self, gene_id=None):
        """
        Get a generator for all the lsvs.  If gene id is specified, then only return lsv for that gene.
        :param gene_id: gene id
        :return: generator of lsv objects
        """

        if gene_id:
            for lsv_id in self.lsv_ids(gene_ids=[gene_id]):

                yield self.lsv(lsv_id)
        else:
            for lsv_id in self.lsv_ids():
                yield self.lsv(lsv_id)

classNameTypeMap = {
    'PsiCoverage': constants.ANALYSIS_PSI,
    'DeltaPsiDataset': constants.ANALYSIS_DELTAPSI,
    'HeterogenDataset': constants.ANALYSIS_HETEROGEN
}

class ViewMatrixType(ViewMatrix):
    def __init__(self, cov_object):
        if type(cov_object) is str:
            cov_object = rna_voila.config.ViewConfig().cov_zarr[cov_object]
        #super().__init__(matrix_hdf5, lsv_id, fields)
        self.q = cov_object
        self.sg = rna_voila.config.ViewConfig().sg_zarr
        self._lsvs = self.sg.exon_connections.lsvs()

    @property
    def analysis_type(self):
        return classNameTypeMap[type(self.q).__name__]

    def lsv_ids(self, gene_ids=None):
        """
        Get a set of lsv ids from all voila files for specified gene ids. If gene ids is None, then get all lsv ids.
        ENSMUSG00000032735:t:61814198-61814332
        :param gene_ids: list of gene ids
        :return:
        """

        events = self.q.get_events(self.sg.introns, self.sg.junctions)

        if not gene_ids:
            yield from events.ec_idx
        else:
            for gene_id in gene_ids:
                yield from events.ec_idx[events.slice_for_gene(self.sg.genes[gene_id])]

    @property
    def means(self):
        raise NotImplementedError()

    @property
    def bins(self):
        """
        Get bins data from rna_voila file.
        :return: list
        """
        return unpack_bins(self.get('bins'))

    @property
    def junction_count(self):
        """
        Get count of all junctions in this lsv.
        :return: integer
        """
        return len(tuple(self.means))



    @property
    def splice_graph_experiment_names(self):
        """
        experiment names for the splice graph drop down.
        :return: List
        """

        #TODO need to get groups, not just experiments
        config = rna_voila.config.ViewConfig()
        groups = []
        for group in config.sgc_zarr.prefixes:
            groups.append([group])
        return groups

    @property
    def experiment_names(self):
        """
        Group names for this set of het voila files.
        :return: list
        """
        #TODO need to get groups, not just experiments

        # config = rna_voila.config.ViewConfig()
        # group_names = []
        # for group in config.sgc_zarr.prefixes:
        #     group_names.append(group)
        # return group_names
        if hasattr(self.q, 'prefixes'):
            return self.q.prefixes
        group_names = []
        for comparison in self.q.comparisons:
            for group_name in comparison:
                group_names.append(str(group_name))
        return group_names

    @property
    def group_names(self):
        """
        Group names for this set of het voila files.
        :return: list
        """
        #TODO need to get groups, not just experiments
        #return ['combined']

        # config = rna_voila.config.ViewConfig()
        # group_names = []
        # for group in config.sgc_zarr.prefixes:
        #     group_names.append(group)
        # return group_names
        if hasattr(self.q, 'prefixes'):
            return self.q.prefixes
        group_names = []
        for comparison in self.q.comparisons:
            for group_name in comparison:
                group_names.append(str(group_name))
        return group_names

    @property
    def gene_ids(self):
        gene_idx_with_lsv = np.unique(self._lsvs.connection_gene_idx())
        return self.sg.genes.df.gene_id[gene_idx_with_lsv].values.tolist()

    def psi(self, lsv_id):
        obj = ViewPsi(cov_object=self.q)
        return obj.lsv(lsv_id)


class LSV_common:

    @property
    def lsv_type(self):
        ref_exon_idx = self._lsvs.ref_exon_idx[self.lsv_id]
        event_type = self._lsvs.event_type[self.lsv_id]
        return self.sg.exon_connections.event_description((ref_exon_idx,), (event_type,))[0]

    @property
    def reference_exon(self):
        ref_exon_idx = self._lsvs.ref_exon_idx[self.lsv_id]
        return (self.sg.exons.start[ref_exon_idx], self.sg.exons.end[ref_exon_idx])

    @property
    def gene_id(self):
        gene_idx = self._lsvs.connection_gene_idx(self.ec_idx_s.start)
        return self.sg.genes.gene_id[gene_idx]

    @property
    def _lsv_id(self):
        return rna_voila.config.ViewConfig().lsvidx2lsvid[self.lsv_id]

    @property
    def junctions(self):
        """
        Including intron retention, a 2d list of junction coordinates. Intron retention is the last set of coordinates.
        :return: numpy matrix
        """

        #ec_idx = np.arange(self._lsvs.ec_idx_start[self.e_idx_s], self._lsvs.ec_idx_end[self.e_idx_s])
        ec_idx = self._lsvs.select_eidx_to_select_ecidx(self.lsv_id)
        junctions = list(zip(self._lsvs.connection_start(ec_idx).tolist(), self._lsvs.connection_end(ec_idx).tolist()))
        # if lsvs.has_intron(e_idx):
        #     introns = junctions[-1:]
        #     junctions = junctions[:-1]
        return junctions

    @property
    def target(self):
        ref_exon_idx = self._lsvs.ref_exon_idx[self.lsv_id]
        event_type = self._lsvs.event_type[self.lsv_id]
        return self.sg.exon_connections.is_target_LSV(ref_exon_idx, event_type)

    @property
    def source(self):
        ref_exon_idx = self._lsvs.ref_exon_idx[self.lsv_id]
        event_type = self._lsvs.event_type[self.lsv_id]
        return self.sg.exon_connections.is_source_LSV(ref_exon_idx, event_type)

    @property
    def binary(self):
        return self.get_attr('binary')

    @property
    def complex(self):
        return self.get_attr('complex')



    @property
    def a5ss(self):
        return self.get_attr('a5ss')

    @property
    def a3ss(self):
        return self.get_attr('a3ss')

    @property
    def exon_skipping(self):
        return self.get_attr('exon_skipping')

    @property
    def exon_count(self):
        return self._lsvs.event_size[self.lsv_id]


    @property
    def intron_retention(self):
        ref_exon_idx = self._lsvs.ref_exon_idx[self.lsv_id]
        event_type = self._lsvs.event_type[self.lsv_id]
        return self.sg.exon_connections.has_intron(ref_exon_idx, event_type)


class ViewPsis(ViewMatrixType):

    def __init__(self, cov_files=None, cov_object=None, group_order_override=None):
        if group_order_override:
            group_order_override = group_order_override.copy()
        self.group_order_override = group_order_override
        if cov_files is None and cov_object is None:
            cached = rna_voila.config.ViewConfig().cov_zarr
            if cached:
                cov_object = cached
        if cov_files is None:
            cov_files = rna_voila.config.ViewConfig().cov_files
        self.cov_files = cov_files
        if cov_object is None:
            cov_object = nm.PsiCoverage.from_zarr(cov_files)
        self.cov_object = cov_object
        super().__init__(cov_object)

    @property
    def group_names(self):

        if self.group_order_override:
            return self.group_order_override
        elif type(rna_voila.config.ViewConfig().cov_zarr) is dict and rna_voila.config.ViewConfig().psicov_grouping_file:
            return list(rna_voila.config.ViewConfig().cov_zarr.keys())

        return self.q.prefixes

    @property
    def experiment_names(self):
        if type(rna_voila.config.ViewConfig().cov_zarr) is dict:
            return rna_voila.config.ViewConfig().cov_zarr_combined.prefixes
            # exp_names = []
            # for cov_obj in rna_voila.config.ViewConfig().cov_zarr.values():
            #     exp_names += cov_obj.prefixes
            # return exp_names
        return self.q.prefixes

    def lsv_ids(self, gene_ids=None):
        events = rna_voila.config.ViewConfig().cov_zarr_combined.get_events(self.sg.introns, self.sg.junctions)

        if not gene_ids:
            yield from events.ec_idx
        else:
            for gene_id in gene_ids:
                yield from events.ec_idx[events.slice_for_gene(self.sg.genes[gene_id])]

    class PsiLSV(ViewMatrixType, LSV_common):

        def __init__(self, cov_files, lsv_id, group_names, experiment_names, cov_object=None,):
            if cov_object is None:
                cov_object = nm.PsiCoverage.from_zarr(cov_files)
            super().__init__(cov_object)
            try:
                self.lsv_id = int(lsv_id)
            except ValueError:
                lsvidx = lsv_id
                if type(lsvidx) is bytes:
                    lsvidx = lsvidx.decode()
                self.lsv_id = rna_voila.config.ViewConfig().lsvid2lsvidx[lsvidx]


            self.ec_idx_s = self._lsvs.connections_slice_for_event(self.lsv_id)
            self._group_names = group_names
            self._experiment_names = experiment_names



        @property
        def means(self):
            """
            Get means data from rna_voila file.
            :return: list
            """
            return self.means_packed

        @property
        def means_packed(self):
            """
            Get means data from rna_voila file.
            :return: list
            """
            return self.q.bootstrap_psi_mean[self.ec_idx_s].to_numpy().T[0].tolist()

        @property
        def bins(self):
            bins = self.q.approximate_discretized_pmf(ec_idx=self.ec_idx_s, nbins=40, midpoint_approximation=True).to_numpy()
            if bins.shape[1] > 1:
                print("trying to get bins and more than one group returned, this should not happen")
                assert False
            return bins[:, 0, :]

        @property
        def group_bins(self):
            """
            Get bins in a dictionary where the key in the name of the group it belongs to.
            :return: generator of key, value
            """

            out = {}
            if rna_voila.config.ViewConfig().psicov_grouping_file:
                for group, cov in self.q.items():
                    bins = cov.approximate_discretized_pmf(ec_idx=self.ec_idx_s, nbins=40, midpoint_approximation=True).mean("prefix").to_numpy()
                    bins = np.nan_to_num(bins).tolist()
                    out[group] = bins
            else:
                bins = self.q.approximate_discretized_pmf(ec_idx=self.ec_idx_s, nbins=40, midpoint_approximation=True).to_numpy()
                for i, group in enumerate(self.experiment_names):
                    _bins = np.nan_to_num(bins[:, i]).tolist()
                    out[group] = _bins

            return out


        @property
        def all_group_means(self):
            return self.group_means

        @property
        def group_means(self):
            """
            Get means data from rna_voila file.
            :return: generator
            """
            out = {}
            if rna_voila.config.ViewConfig().psicov_grouping_file:
                for group, cov in self.q.items():
                    means = cov.bootstrap_psi_mean[self.ec_idx_s].mean("prefix").to_numpy()
                    means = np.nan_to_num(means).tolist()
                    out[group] = means
            else:
                means = self.q.bootstrap_psi_mean[self.ec_idx_s].to_numpy()
                for i, group in enumerate(self.experiment_names):
                    _means = np.nan_to_num(means[:, i]).tolist()
                    out[group] = _means

            return out


        @property
        def variances(self):
            """
            Create variance data of bins data.
            :return: list
            """
            return generate_variances(self.bins)

        @property
        def standard_deviations(self):
            """
            Create variance data of bins data.
            :return: list
            """
            return generate_standard_deviations(self.bins)



    def lsv(self, lsv_id):
        """
        Get lsv object by lsv id.
        :param lsv_id: lsv id
        :return: lsv object
        """
        return self.PsiLSV(self.cov_files, lsv_id, self.group_names, self.experiment_names, cov_object=self.cov_object)

class ViewPsi(ViewPsis):
    def __init__(self, cov_file=None, cov_object=None):
        if not cov_file:
            cov_file = rna_voila.config.ViewConfig().cov_file
        self.cov_file = cov_file
        super().__init__(cov_files=[self.cov_file], cov_object=cov_object)

    def psi(self, lsv_id):
        return self.lsv(lsv_id)



class ViewDeltaPsi(ViewMatrixType):
    def __init__(self, cov_file=None):
        """
        View for delta psi matrix.  This is used in creation of tsv and html files.
        """
        self.config = rna_voila.config.ViewConfig()


        if cov_file is None:
            self.cov_object = self.config.cov_zarr
        else:
            self.cov_object = nm.DeltaPsiDataset.from_zarr(cov_file)
        super().__init__(self.cov_object)

    @property
    def group_names(self):
        """
        Group names for this set of het voila files.
        :return: list
        """
        #TODO verify ??
        return self.q.comparisons[0]

    class DeltaPsiLSV(ViewMatrixType, LSV_common):

        def __init__(self, cov_object, lsv_id):
            super().__init__(cov_object)
            if type(lsv_id) is str:
                self.lsv_id = rna_voila.config.ViewConfig().lsvid2lsvidx[lsv_id]
            else:
                self.lsv_id = lsv_id

            self.ec_idx_s = self._lsvs.connections_slice_for_event(self.lsv_id)

        @property
        def means(self):
            """
            Create mean data from bins data.
            :return: list
            """

            return self.q.bootstrap_posterior.mean[0, self.ec_idx_s]

        @property
        def bins(self):
            """
            Create mean data from bins data.
            :return: list
            """

            return self.q.bootstrap_posterior.p[0, self.ec_idx_s]
            #return

        @property
        def group_bins(self):
            """
            Get dictionary of bins by group name.
            :return: generator of key, value
            """
            group_names = self.q.comparisons[0]
            bins = self.q.groups.approximate_discretized_pmf(nbins=80, ec_idx=self.ec_idx_s).to_numpy()

            return {g: p for g, p in zip(group_names, unpack_bins(bins))}




        @property
        def group_means(self):
            """
            Get dictionary of mean by group name.
            :return: generator of key, value
            """
            group_names = self.q.comparisons[0]
            means = self.q.groups.bootstrap_psi_mean[:, self.ec_idx_s].to_numpy().T
            means = np.nan_to_num(means)

            #print({g: p for g, p in zip(group_names, means)})
            return {g: p for g, p in zip(group_names, means)}


        @property
        def excl_incl(self):
            """
            Using means data, create exclude/include list.
            :return: list
            """
            return generate_excl_incl(self.means.values.tolist())

        def is_lsv_changing(self, threshold):
            """
            Is lsv changing based on threshold supplied by user.
            :param threshold: threshold value, usually 0.2
            :return: boolean
            """
            return is_lsv_changing(self.means, threshold)

        def probability_threshold(self):
            """
            Get probability that the junction in an LSV are changing.

            :return: list of p values
            """
            args = self.matrix_hdf5.args
            probability_threshold = args.probability_threshold
            threshold = args.probability_threshold
            return any(matrix_area(b, threshold=threshold) >= probability_threshold for b in self.bins)

        def changing(self, changing_threshold, probability_changing_threshold, junc_i=None):
            """
            Yet another changing function wrapper, because I can't figure out how the legacy ones were even being used
            This is used for creating the changing column in voila classifier. It functions similar in API to the
            non changing function below

            """
            for _junc_i in range(len(self.bins)) if junc_i is None else [junc_i]:
                junc_bins = self.bins[_junc_i]
                if matrix_area(junc_bins, threshold=changing_threshold) >= probability_changing_threshold:
                    return True

            return False

        def high_probability_non_changing(self, non_changing_threshold=None, junc_i=None):
            """
            Get probability that junctions in an lsv aren't changing.
            :return: list
            """
            if non_changing_threshold is None:
                non_changing_threshold = self.config.non_changing_between_group_dpsi

            bins = self.bins if junc_i is None else [self.bins[junc_i]]

            result = generate_high_probability_non_changing(self.intron_retention, self.matrix_hdf5.prior,
                                                          non_changing_threshold, bins)

            return result if junc_i is None else result[0]

    def lsv(self, lsv_id):
        """
        Get delta psi object by lsv id.
        :param lsv_id: lsv id
        :return: delta psi object
        """
        return self.DeltaPsiLSV(self.cov_object, lsv_id)


class ViewHeterogens(ViewMatrixType):

    def __init__(self, cov_file=None, group_order_override=None):
        """
        View for delta psi matrix.  This is used in creation of tsv and html files.
        """
        self.group_order_override = group_order_override
        self.config = rna_voila.config.ViewConfig()


        if cov_file is None:
            self.cov_object = nm.HeterogenDataset.from_zarr(self.config.cov_file)
        else:
            self.cov_object = nm.HeterogenDataset.from_zarr(cov_file)
        super().__init__(self.cov_object)

    @property
    def group_names(self):
        if self.group_order_override:
            return self.group_order_override
        return list(self.q.groups.keys())

    @property
    def experiment_names(self):
        return [self.q.groups[grp].prefixes for grp in self.q.groups]

    class HeterogenLSV(ViewMatrixType, LSV_common):

        def __init__(self, cov_object, lsv_id, group_names):
            super().__init__(cov_object)
            self._group_names = group_names
            try:
                self.lsv_id = int(lsv_id)
            except:
                self.lsv_id = rna_voila.config.ViewConfig().lsvid2lsvidx[lsv_id]


            self.ec_idx_s = self._lsvs.connections_slice_for_event(self.lsv_id)


        @property
        def groups_quantified(self):

            # mu_psi has shape (njunc, ngroups, max(nexperiments))
            mu_psi = np.array(self.mu_psi)[0]
            num_quant = (mu_psi >= 0).sum(axis=-1)
            return num_quant

        def get_attr(self, attr):
            """
            For attributes that exist is each het file, this will get all values and confirm they're all equal.
            :param attr: attribute found in het voila file.
            :return: attribute value
            """
            config = rna_voila.config.ViewConfig()
            if config.strict_indexing:
                s = set()
                for f in config.voila_files:
                    with ViewHeterogen(f) as m:
                        try:
                            het = m.lsv(self.lsv_id)
                            s.add(getattr(het, attr))
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                            pass
                assert len(s) == 1, s
                return s.pop()
            else:
                for f in config.voila_files:
                    with ViewHeterogen(f) as m:
                        try:
                            het = m.lsv(self.lsv_id)
                            return getattr(het, attr)
                        except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                            pass

        @property
        def group_bins(self):
            """
            Associate means values with it's experiment group name.
            :return: Generator key/value
            """

            mean_psi = self.mean_psi
            mean_psi = np.array(mean_psi)
            mean_psi = mean_psi.transpose((1, 0, 2))
            for group_name, mean in zip(self._group_names, mean_psi):
                yield group_name, mean.tolist()

        def junction_heat_map(self, stat_name, junc_idx):
            hets_grps = list(self.q.groups.keys())
            hets_grps_len = len(hets_grps)
            s = np.empty((hets_grps_len, hets_grps_len))

            s.fill(-1)


            stat_idx = list(self.q.stats.to_series()).index(stat_name)
            group_names = self._group_names

            for comp_idx, (group_name1, group_name2) in enumerate(self.q.comparisons):

                grp1median = self.q.groups[group_name1].raw_psi_mean_population_median[self.ec_idx_s]
                grp2median = self.q.groups[group_name2].raw_psi_mean_population_median[self.ec_idx_s]


                stat_value = self.q.approximate_pvalue[comp_idx, int(self.ec_idx_s.start)+int(junc_idx), stat_idx]
                dpsi_value = (grp1median - grp2median)[junc_idx]

                x = group_names.index(group_name1)
                y = group_names.index(group_name2)

                if x - y > 0:
                    s[x][y] = stat_value
                    s[y][x] = dpsi_value
                else:
                    s[x][y] = dpsi_value
                    s[y][x] = stat_value



            return np.nan_to_num(s).tolist()

        @property
        def median_psi(self):
            return self._median_psi()

        def quantile_psi(self, quantile):
            return self._median_psi(quantile=quantile)

        def _median_psi(self, quantile=None):
            """ Find median_psi per group in all het voila files

            Returns
            -------
            np.array(shape=(num_groups, num_junctions))
                The median PSI per junction for each group. -1 if not quantified
            """
            group_names = list(self.q.groups.keys())
            juncs_len = self.ec_idx_s.stop - self.ec_idx_s.start
            grps_len = len(group_names)
            # fill result
            median_psi = np.empty((grps_len, juncs_len))
            median_psi.fill(-1)

            for grp_idx, group_name in enumerate(self.q.groups):
                if quantile is None:
                    medians = self.q.groups[group_name].raw_psi_mean_population_median[self.ec_idx_s]
                else:
                    medians = self.q.groups[group_name].raw_psi_mean_population_quantile([quantile], ec_idx=self.ec_idx_s)[:,0]
                #mus = np.nanmedian(self.q.psi_pmf(group_name, self.ec_idx_s), axis=(1,2))
                median_psi[grp_idx, :] = medians
            #
            #
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         try:
            #             het = m.lsv(self.lsv_id)
            #
            #             if quantile is None:
            #                 medians = het.median_psi()
            #             else:
            #                 medians = het.quantile_psi(quantile)
            #             # get index into medians (second axis) per group
            #             for ndx, grp in enumerate(m.group_names):
            #                 idx = group_names.index(grp)  # index to median_psi
            #                 median_psi[idx, :] = medians[:, ndx]  # fill median_psi
            #
            #         except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
            #             pass
            return median_psi

        @property
        def mu_psi(self):
            """
            Find mu_psi in all het voila files and create a matrix that matches the new unified set of group/experiment
            names.  In the case where the experiments are all the same number, this will fill the empty space in the
            matrix with -1.
            :return: List
            """

            group_names = list(self.q.groups.keys())
            experiments = list(self.q.groups.values())
            exps_len = max(e.df.sizes['prefix'] for e in experiments)
            juncs_len = self.ec_idx_s.stop - self.ec_idx_s.start
            grps_len = len(group_names)
            mu_psi = np.empty((grps_len, juncs_len, exps_len))
            mu_psi.fill(-1)

            for grp_idx, group_name in enumerate(self.q.groups):
                mus = self.q.groups[group_name].raw_psi_mean[self.ec_idx_s].to_numpy()
                mu_psi[grp_idx][0:mus.shape[0], 0:mus.shape[1]] = mus

            mu_psi = mu_psi.transpose((1, 0, 2))
            #
            # voila_files = rna_voila.config.ViewConfig().voila_files
            # group_names = self.matrix_hdf5.group_names
            # experiment_names = self.matrix_hdf5.experiment_names
            # exps_len = max(len(e) for e in experiment_names)
            # juncs_len = len(self.junctions)
            # grps_len = len(group_names)
            # mu_psi = np.empty((grps_len, juncs_len, exps_len))
            # mu_psi.fill(-1)
            #
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         try:
            #
            #             het = m.lsv(self.lsv_id)
            #
            #             mus = het.mu_psi
            #             mus = mus.transpose((1, 0, 2))
            #
            #             for mu, grp in zip(mus, m.group_names):
            #                 idx = group_names.index(grp)
            #                 mu_shp = mu.shape
            #                 mu_psi[idx][0:mu_shp[0], 0:mu_shp[1]] = mu
            #
            #         except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
            #             pass
            #
            # mu_psi = mu_psi.transpose((1, 0, 2))
            # ret shape is num_groups, num juncs, num experiments
            return mu_psi.tolist()

        @property
        def mean_psi(self):
            """
             Find mean_psi in all het voila files and create a matrix that matches the new unified set of group/experiment
            names.  In the case where the experiments are all the same number, this will fill the empty space in the
            matrix with -1.
            :return: List
            """
            group_names = list(self.q.groups.keys())
            juncs_len = self.ec_idx_s.stop - self.ec_idx_s.start
            grps_len = len(group_names)
            mean_psi = np.empty((grps_len, juncs_len, 40))
            mean_psi.fill(-1)

            for grp_idx, group_name in enumerate(self.q.groups):
                mus = np.mean(self.q.psi_pmf(group_name, self.ec_idx_s), axis=1).to_numpy()
                mean_psi[grp_idx][0:mus.shape[0], 0:mus.shape[1]] = mus

            mean_psi = mean_psi.transpose((1, 0, 2))
            #
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         try:
            #             het = m.lsv(self.lsv_id)
            #
            #             means = het.mean_psi
            #             means = means.transpose((1, 0, 2))
            #
            #             for mn, grp in zip(means, m.group_names):
            #                 idx = group_names.index(grp)
            #                 mn_shp = mn.shape
            #                 mean_psi[idx][0:mn_shp[0], 0:mn_shp[1]] = mn
            #         except (LsvIdNotFoundInVoilaFile, GeneIdNotFoundInVoilaFile):
            #             pass
            #
            # mean_psi = mean_psi.transpose((1, 0, 2))
            return mean_psi.tolist()

        @property
        def junction_psisamples_stats(self):
            """
            Get stats from psisamples quantiles with values
            :return: generator key/value
            """

            for stat_name in list(self.q.stats.to_series()):
                stat_idx = list(self.q.stats.to_series()).index(stat_name)
                stat_value = self.q.approximate_pvalue_quantiles[0, int(self.ec_idx_s.start):int(self.ec_idx_s.stop),
                             stat_idx].to_numpy()[:,0]

                yield stat_name, stat_value

            # config = rna_voila.config.ViewConfig()
            # voila_files = config.voila_files
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         het = m.lsv(self.lsv_id)
            #         groups = '-'.join(m.group_names)
            #         stat_names = m.stat_names
            #         try:
            #             for name, stat in zip(stat_names, het.junction_psisamples_stats.T):
            #                 if len(voila_files) == 1:
            #                     yield f"{name}_quantile", stat
            #                 else:
            #                     yield f"{groups} {name}_quantile", stat
            #         except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
            #             pass

        @property
        def junction_stats(self):
            """
            This gets associates stat test names with their values.
            :return: generator key/value
            """
            # hets_grps = list(self.q.groups.keys())
            # hets_grps_len = len(hets_grps)
            # s = np.empty((hets_grps_len, hets_grps_len))
            #
            # s.fill(-1)
            #
            # stat_idx = list(self.q.stats.to_series()).index(stat_name)
            # group_names = self._group_names
            #
            # for comp_idx, (group_name1, group_name2) in enumerate(self.q.comparisons):
            #     grp1median = self.q.groups[group_name1].raw_psi_mean_population_median[self.ec_idx_s]
            #     grp2median = self.q.groups[group_name2].raw_psi_mean_population_median[self.ec_idx_s]
            #
            #     #stat_value = self.q.approximate_pvalue[comp_idx, int(self.ec_idx_s.start) + int(junc_idx), stat_idx]
            #     stat_value = self.q.approximate_pvalue[comp_idx, int(self.ec_idx_s.start):int(self.ec_idx_s.end), stat_idx]
            #

            for stat_name in list(self.q.stats.to_series()):
                stat_idx = list(self.q.stats.to_series()).index(stat_name)
                stat_value = self.q.approximate_pvalue[0, int(self.ec_idx_s.start):int(self.ec_idx_s.stop),
                             stat_idx].to_numpy()
                yield stat_name, stat_value
            #
            # config = rna_voila.config.ViewConfig()
            # voila_files = config.voila_files
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         het = m.lsv(self.lsv_id)
            #         groups = '-'.join(m.group_names)
            #         stat_names = m.stat_names
            #         try:
            #             for name, stat in zip(stat_names, het.junction_stats.T):
            #                 if len(voila_files) == 1:
            #                     yield name, stat
            #                 else:
            #                     yield groups + ' ' + name, stat
            #         except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
            #             pass

        def changing(
            self,
            pvalue_threshold: float = 0.05,
            between_group_dpsi: float = 0.05
        ):
            """ Boolean of heuristic for changing heterogeneous events

            Parameters
            ----------
            pvalue_threshold: float
                Minimum p-value for which an LSV/junction can return true. Uses
                minimum p-value from all tests provided
            between_group_dpsi: float
                Maximum absolute difference in median values of PSI for which
                an LSV/junction can return true

            Returns
            -------
            Generator yielding group names and boolean array per junction
            """


            yield "changing", 'N/A'
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         het = m.lsv(self.lsv_id)
            #         groups = '-'.join(m.group_names)
            #         try:
            #             changing = het.changing(
            #                 pvalue_threshold=pvalue_threshold,
            #                 between_group_dpsi=between_group_dpsi
            #             )
            #             if len(voila_files) == 1:
            #                 yield "changing", changing
            #             else:
            #                 yield f"{groups} changing", changing
            #         except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
            #             pass

        def nonchanging(
            self,
            pvalue_threshold: float = 0.05,
            within_group_iqr: float = 0.10,
            between_group_dpsi: float = 0.05
        ):
            """ Boolean of heuristic for nonchanging heterogeneous events

            We define highly-confident non-changing events from MAJIQ Heterogen
            as being (1) above a nominal p-value threshold, (2) within-group
            variance is sufficiently low as measured by IQR, (3) between-group
            dPSI is sufficiently low as measured by difference in medians. We
            accept that between-group dPSI threshold may be redundant in
            combination with the other two thresholds.

            Parameters
            ----------
            pvalue_threshold: float
                Minimum p-value for which an LSV/junction can return true. Uses
                minimum p-value from all tests provided
            within_group_iqr: float
                Maximum IQR within a group for which an LSV/junction can return
                true
            between_group_dpsi: float
                Maximum absolute difference in median values of PSI for which
                an LSV/junction can return true

            Returns
            -------
            Generator yielding group names and boolean array per junction
            """

            yield "nonchanging", 'N/A'
            #
            # config = rna_voila.config.ViewConfig()
            # voila_files = config.voila_files
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         het = m.lsv(self.lsv_id)
            #         groups = '-'.join(m.group_names)
            #         try:
            #             nonchanging = het.nonchanging(
            #                 pvalue_threshold=pvalue_threshold,
            #                 within_group_iqr=within_group_iqr,
            #                 between_group_dpsi=between_group_dpsi
            #             )
            #             if len(voila_files) == 1:
            #                 yield "nonchanging", nonchanging
            #             else:
            #                 yield f"{groups} nonchanging", nonchanging
            #         except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
            #             pass

        @property
        def junction_scores(self):
            """
            This gets associates stat test score names with their values. (if any exist)
            :return: generator key/value
            """
            list(self.q.stats.to_series())
            score_names = ()
            for score_name in score_names:
                yield score_name, 0

            #
            # config = rna_voila.config.ViewConfig()
            # voila_files = config.voila_files
            # for f in voila_files:
            #     with ViewHeterogen(f) as m:
            #         groups = '-'.join(m.group_names)
            #         score_names = ('tnom_score',) if 'TNOM' in m.stat_names else ()
            #         try:
            #             try:
            #                 score_vals = m.get(self.lsv_id, 'tnom_score').T
            #             except KeyError:
            #                 score_names = ()
            #
            #             for score_name in score_names:
            #                 if len(voila_files) == 1:
            #                     yield score_name, score_vals
            #                 else:
            #                     yield groups + ' ' + score_name, score_vals
            #
            #         except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
            #             pass



    @property
    def stat_names(self):
        """
        Gets list of stat test names.
        :return: List
        """
        return list(self.cov_object.stats.to_series())

    @property
    def psi_samples_summary(self):
        """
        Summarize number of psi_samples used for posterior samples statistics

        :return: str indicating:
            + unknown: when all samples are missing this value in their metadata
            + inconsistent: different values found (or some missing/not)
            + value shared/consistent for all comparisons
        """
        return "no metadata"

        # values = set()  # type: Set[str]
        # voila_files = rna_voila.config.ViewConfig().voila_files
        # for f in voila_files:
        #     with ViewHeterogen(f) as m:
        #         try:
        #             values.add(str(m.psi_samples))
        #         except KeyError:
        #             values.add(
        #                 "missing metadata (older version of MAJIQ)"
        #             )
        # if len(values) > 1:
        #     return f"inconsistent ({sorted(values)})"
        # elif values:  # only 1
        #     return values.pop()
        # else:
        #     return f"no metadata (input files: {voila_files})"

    @property
    def test_percentile_summary(self):
        """
        Summarize percentile taken over psi-sample statistics

        :return: str indicating:
            + unknown: when all samples are missing this value in their metadata
            + inconsistent: different values found (or some missing/not)
            + value shared/consistent for all comparisons
        """
        return "no metadata"

        # values = set()  # type: Set[str]
        # voila_files = rna_voila.config.ViewConfig().voila_files
        # for f in voila_files:
        #     with ViewHeterogen(f) as m:
        #         try:
        #             values.add(str(m.test_percentile))
        #         except KeyError:
        #             values.add(
        #                 "missing metadata (older version of MAJIQ)"
        #             )
        # if len(values) > 1:
        #     return f"inconsistent ({sorted(values)})"
        # elif values:  # only 1
        #     return values.pop()
        # else:
        #     return f"no metadata (input files: {voila_files})"

    @property
    def junction_psisamples_stats_column_names(self):
        return (f"{x}_quantile" for x in self.junction_stats_column_names)

    @property
    def junction_stats_column_names(self):
        """
        Stat column names for tsv output.
        :return: generator
        """

        for name in self.stat_names:
            yield name

        # voila_files = rna_voila.config.ViewConfig().voila_files
        #
        # for f in voila_files:
        #     with ViewHeterogen(f) as m:
        #         groups = '-'.join(m.group_names)
        #         for name in m.stat_names:
        #             if len(voila_files) == 1:
        #                 yield name
        #             else:
        #                 yield groups + ' ' + name

    def _per_group_column_names(self, prefix):
        yield prefix

        # voila_files = rna_voila.config.ViewConfig().voila_files
        #
        # if len(voila_files) == 1:
        #     yield prefix
        # else:
        #     for f in voila_files:
        #         with ViewHeterogen(f) as m:
        #             groups = '-'.join(m.group_names)
        #             yield f"{groups} {prefix}"

    @property
    def changing_column_names(self):
        """ Column names associated with changing values
        """
        return self._per_group_column_names('changing')

    @property
    def nonchanging_column_names(self):
        """ Column names associated with nonchanging values
        """
        return self._per_group_column_names('nonchanging')

    @property
    def junction_scores_column_names(self):
        """
        Stat column names for tsv output.
        :return: generator
        """
        score_names = ('tnom_score',) if 'TNOM' in self.stat_names else ()
        for name in score_names:
            yield name
        #
        # voila_files = rna_voila.config.ViewConfig().voila_files
        #
        # for f in voila_files:
        #     with ViewHeterogen(f) as m:
        #         groups = '-'.join(m.group_names)
        #         score_names = ('tnom_score',) if 'TNOM' in m.stat_names else ()
        #         for name in score_names:
        #             if len(voila_files) == 1:
        #                 yield name
        #             else:
        #                 yield groups + ' ' + name

    # @property
    # def experiment_names(self):
    #     """
    #     Experiment names for this set of het voila files.
    #     :return: List
    #     """
    #     config = rna_voila.config.ViewConfig()
    #     exp_names = {}
    #     for f in config.voila_files:
    #         with ViewHeterogen(f) as m:
    #             for exp, grp in zip(m.experiment_names, m.group_names):
    #                 exp_names[grp] = exp
    #
    #     return [exp_names[grp] for grp in self.group_names]

    # @property
    # def group_names(self):
    #     """
    #     Group names for this set of het voila files.
    #     :return: list
    #     """
    #     if self.group_order_override:
    #         return self.group_order_override
    #     config = rna_voila.config.ViewConfig()
    #     grp_names = []
    #     for f in config.voila_files:
    #         with ViewHeterogen(f) as m:
    #             for grp in m.group_names:
    #                 if not grp in grp_names:
    #                     grp_names.append(grp)
    #
    #     return grp_names
    #
    # @property
    # def splice_graph_experiment_names(self):
    #     """
    #     experiment names for the splice graph drop down.
    #     :return: List
    #     """
    #
    #     config = rna_voila.config.ViewConfig()
    #     exp_names = {}
    #     for f in config.voila_files:
    #         with ViewHeterogen(f) as m:
    #             for exp, grp in zip(m.splice_graph_experiment_names, m.group_names):
    #                 exp_names[grp] = exp
    #
    #     return [exp_names[grp] for grp in self.group_names]
    #
    # @property
    # def gene_ids(self):
    #     """
    #     Get a set of gene ids from all het voila files.
    #     :return: generator
    #     """
    #
    #     voila_files = rna_voila.config.ViewConfig().voila_files
    #     vhs = [ViewHeterogen(f) for f in voila_files]
    #     yield from set(chain(*(v.gene_ids for v in vhs)))
    #     for v in vhs:
    #         v.close()

    @staticmethod
    def pair_merge(pairs):
        if len(pairs) == 1:
            if type(pairs[0]) is set:
                return pairs[0]
            else:
                with ViewHeterogen(pairs[0]) as v:
                    return set(v.lsv_ids())
        else:
            if type(pairs[0]) is set:
                s1 = pairs[0]
            else:
                with ViewHeterogen(pairs[0]) as v:
                    s1 = set(v.lsv_ids())
            if type(pairs[1]) is set:
                s2 = pairs[1]
            else:
                with ViewHeterogen(pairs[1]) as v:
                    s2 = set(v.lsv_ids())
            return s1 | s2


    # def lsv_ids(self, gene_ids=None):
    #     """
    #     Get a set of lsv ids from all voila files for specified gene ids. If gene ids is None, then get all lsv ids.
    #     :param gene_ids: list of gene ids
    #     :return:
    #     """
    #     assert NotImplementedError()
    #     if gene_ids or len(rna_voila.config.ViewConfig().voila_files) == 1:
    #         voila_files = rna_voila.config.ViewConfig().voila_files
    #         vhs = [ViewHeterogen(f) for f in voila_files]
    #         yield from set(chain(*(v.lsv_ids(gene_ids) for v in vhs)))
    #         for v in vhs:
    #             v.close()
    #     else:
    #         vhs = rna_voila.config.ViewConfig().voila_files
    #         p = Pool(rna_voila.config.ViewConfig().nproc)
    #         while len(vhs) > 1:
    #             vhs = [vhs[i:i + 2] for i in range(0, len(vhs), 2)]
    #             vhs = p.map(self.pair_merge, vhs)
    #
    #         yield from vhs[0]

    def lsv(self, lsv_id):
        """
        Get delta psi object by lsv id.
        :param lsv_id: lsv id
        :return: delta psi object
        """
        return self.HeterogenLSV(self.cov_object, lsv_id, self.group_names)

class ViewHeterogen(ViewMatrix):
    def __init__(self, voila_file):
        """
        This represents a single het voila file.  ViewHeterogens uses this class to retrieve data from the individual
        files.
        :param voila_file: voila file name
        """
        super().__init__(voila_file)

    class _ViewHeterogen(ViewMatrixType):
        def __init__(self, matrix_hdf5, lsv_id):
            super().__init__(matrix_hdf5, lsv_id)

        @property
        def dpsi(self):
            """
            Calculated the absolute difference in psi for heat map.
            :return: list
            """
            return [abs(reduce(operator.__sub__, (get_expected_psi(b) for b in bs))) for bs in self.mean_psi]

        @property
        def dpsi_signed(self):
            """
            Calculated the difference in psi for heat map. (with negative values possible)
            :return: list
            """
            return [reduce(operator.__sub__, (get_expected_psi(b) for b in bs)) for bs in self.mean_psi]

        @property
        def dpsi_median_signed(self):
            """
            Calculated the difference in median psi for heat map. (with negative values possible)
            :return: list
            """
            return [reduce(operator.__sub__, bs) for bs in self.median_psi()]

        @property
        def mean_psi(self):
            return self.get('mean_psi')

        @property
        def mu_psi(self):
            return self.get('mu_psi')

        @property
        def mu_psi_nanmasked(self):
            """ Mask missing/unquantified values in mu_psi with nan
            """
            mu_psi = self.mu_psi  # load mu psi ndarray into memory
            return np.where(mu_psi >= 0, mu_psi, np.nan)  # mask unquantified

        @property
        def junction_stats(self):
            """ Junction statistics computed on posterior means
            """
            return self.get('junction_stats')

        @property
        def junction_psisamples_stats(self):
            """ Quantile from junction statistics on posterior samples
            """
            return self.get('junction_psisamples_stats')

        def median_psi(self, mu_psi=None):
            """ Get group medians of psi_mean (per junction if present)

            Parameters
            ----------
            mu_psi: Optional[np.array(shape=(..., 2, num_experiments))]
                Values to compute median per group, with unquantified values
                masked as nan. If not specified, use self.mu_psi_nanmasked.
                Preceding axes typically correspond to junctions.

            Returns
            -------
            np.array(shape=(..., 2))
                Median value of quantified values of PSI per group/junction
            """
            if mu_psi is None:
                mu_psi = self.mu_psi_nanmasked
            return np.nanmedian(mu_psi, axis=-1)

        def quantile_psi(self, quantile):
            return np.nanquantile(self.mu_psi_nanmasked, quantile, axis=-1)

        def iqr_psi(self, mu_psi=None):
            """ Get group IQRs of psi_mean (per junction if present)

            Parameters
            ----------
            mu_psi: Optional[np.array(shape=(..., 2, num_experiments))]
                Values to compute median per group, with unquantified values
                masked as nan. If not specified, use self.mu_psi_nanmasked.
                Preceding axes typically correspond to junctions.

            Returns
            -------
            np.array(shape=(..., 2))
                IQR of PSI per group/junction
            """
            if mu_psi is None:
                mu_psi = self.mu_psi_nanmasked
            return scipy.stats.iqr(mu_psi, axis=-1, nan_policy="omit")


        def dpsi_median(self,
                         junc_i: int = None):

            if junc_i is None:
                junc_i = slice(None)  # vectorize over all junctions

            # get masked mean values of psi per group/experiment
            mu_psi = self.mu_psi_nanmasked[junc_i]
            # use to compute median psi per group
            median_psi = self.median_psi(mu_psi)
            # did difference in medians pass threshold?
            dpsi_median = np.abs(median_psi[..., 1] - median_psi[..., 0])
            return dpsi_median

        def changing(
            self,
            pvalue_threshold: float = 0.05,
            between_group_dpsi: float = 0.2,
            junc_i: int = None
        ):


            if junc_i is None:
                junc_i = slice(None)  # vectorize over all junctions

            # all statistics must be less than p-value threshold
            pvalue_passed = np.nanmax(self.junction_stats[junc_i], axis=-1) <= pvalue_threshold

            dpsi_passed = self.dpsi_median(junc_i) >= between_group_dpsi

            # pvalue and dpsi thresholds must all pass

            return pvalue_passed & dpsi_passed

        def nonchanging(
            self,
            pvalue_threshold: float = 0.05,
            within_group_iqr: float = 0.10,
            between_group_dpsi: float = 0.05,
            junc_i: int = None
        ):
            """ Boolean of heuristic for nonchanging heterogeneous events

            We define highly-confident non-changing events from MAJIQ Heterogen
            as being (1) above a nominal p-value threshold, (2) within-group
            variance is sufficiently low as measured by IQR, (3) between-group
            dPSI is sufficiently low as measured by difference in medians. We
            accept that between-group dPSI threshold may be redundant in
            combination with the other two thresholds.

            Parameters
            ----------
            pvalue_threshold: float
                Minimum p-value for which an LSV/junction can return true. Uses
                minimum p-value from all tests provided
            within_group_iqr: float
                Maximum IQR within a group for which an LSV/junction can return
                true
            between_group_dpsi: float
                Maximum absolute difference in median values of PSI for which
                an LSV/junction can return true
            junc_i: int
                Only calculate result for one junction from the LSV at the specified index
                None (default) returns a list of
            """
            if junc_i is None:
                junc_i = slice(None)  # vectorize over all junctions

            # all statistics must be greater than p-value threshold
            pvalue_passed = np.nanmin(self.junction_stats[junc_i], axis=-1) >= pvalue_threshold
            # mean values of PSI per group/experiment for requested junctions
            # missing quantifications are masked as nan
            mu_psi = self.mu_psi_nanmasked[junc_i]
            # use to compute IQR and median
            iqr_psi = self.iqr_psi(mu_psi)
            median_psi = self.median_psi(mu_psi)
            # did IQR pass threshold?
            iqr_passed = (iqr_psi <= within_group_iqr).all(axis=-1)  # both groups
            # did difference in medians pass threshold?
            dpsi_passed = np.abs(median_psi[..., 1] - median_psi[..., 0]) <= between_group_dpsi

            # pvalue, iqr, and dpsi thresholds must all pass
            return pvalue_passed & iqr_passed & dpsi_passed

    def lsv(self, lsv_id):
        return self._ViewHeterogen(self, lsv_id)



class ViewMulti:
    """
    View for set of Voila  files.  This is used in creation of tsv and html files.
    Base class
    """
    def __init__(self, view_class):
        """
        :param view_class: class for view of single item (ex, ViewPsi, ViewHeterogen)
        """
        self._group_names = None
        self.view_class = view_class
        config = rna_voila.config.ViewConfig()
        if self.view_class is ViewPsi:
            self.qmulti = nm.PsiCoverage.from_zarr(config.cov_files)
        elif self.view_class is ViewDeltaPsi:
            self.qmulti = nm.DeltaPsiDataset.from_zarr(config.cov_files)
        else:
            self.qmulti = nm.HeterogenDataset.from_zarr(config.cov_files)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    @property
    def experiment_names(self):
        """
        Experiment names for this set of het voila files.
        :return: List
        """
        config = rna_voila.config.ViewConfig()
        exp_names = {}
        for f in config.voila_files:
            with self.view_class(f) as m:
                for exp, grp in zip(m.experiment_names, m.group_names):
                    exp_names[grp] = exp

        return [exp_names[grp] for grp in self.group_names]

    @property
    def group_names(self):
        """
        Group names for this set of het voila files.
        :return: list
        """
        config = rna_voila.config.ViewConfig()
        grp_names = []
        for f in config.voila_files:
            with self.view_class(f) as m:
                for grp in m.group_names:
                    if not grp in grp_names:
                        grp_names.append(grp)

        return grp_names

    @property
    def splice_graph_experiment_names(self):
        """
        experiment names for the splice graph drop down.
        :return: List
        """

        config = rna_voila.config.ViewConfig()
        exp_names = {}
        for f in config.voila_files:
            with self.view_class(f) as m:
                for exp, grp in zip(m.splice_graph_experiment_names, m.group_names):
                    exp_names[grp] = exp

        return [exp_names[grp] for grp in self.group_names]

    @property
    def gene_ids(self):
        """
        Get a set of gene ids from all het voila files.
        :return: generator
        """

        voila_files = rna_voila.config.ViewConfig().voila_files
        vhs = [self.view_class(f) for f in voila_files]
        yield from set(chain(*(v.gene_ids for v in vhs)))
        for v in vhs:
            v.close()

    def lsv_ids_(self, gene_ids=None):
        """
        Get a set of lsv ids from all voila files for specified gene ids. If gene ids is None, then get all lsv ids.
        ENSMUSG00000032735:t:61814198-61814332
        :param gene_ids: list of gene ids
        :return:
        """

        events = self.qmulti.get_events(self.sg.introns, self.sg.junctions)


        if not gene_ids:
            ref_exons = events.ref_exon_idx
            event_types = events.event_type
        elif len(gene_ids) == 1:
            ref_exons = events.ref_exon_idx[events.slice_for_gene(self.sg.genes.gene_id.index(gene_ids[0]))]
            event_types = events.event_type[events.slice_for_gene(self.sg.genes.gene_id.index(gene_ids[0]))]
        else:
            ref_exons = np.concatenate(
                [
                    events.ref_exon_idx[events.slice_for_gene(self.sg.genes[gene_id])] for gene_id in gene_ids
                ])
            event_types = np.concatenate(
                [
                    events.event_type[events.slice_for_gene(self.sg.genes[gene_id])] for gene_id in gene_ids
                ])
            #events_slice = Ki[events.slice_for_gene(self.sg.genes.gene_id.index(gene_id)) for gene_id in gene_ids]]

        yield from self.sg.exon_connections.event_id(ref_exons, event_types)



    def lsv(self, lsv_id):
        raise NotImplementedError()

    def lsvs(self, gene_id=None):
        """
        Get all lsvs for set of het voila files.
        :param gene_id: gene id
        :return:
        """

        if gene_id:
            gene_ids = [gene_id]
        else:
            gene_ids = None

        for lsv_id in self.lsv_ids(gene_ids):
            yield self.lsv(lsv_id)


class _ViewMulti:
    """
    Base class for .lsv() call return of other multi objects
    """
    def __init__(self, matrix_hdf5, lsv_id, view_class):
        """
        :param view_class: class for view of single item (ex, ViewPsi, ViewHeterogen)
        """
        self.matrix_hdf5 = matrix_hdf5
        self.lsv_id = lsv_id
        self.view_class = view_class

    def _get_prop(self, prop, cast=None):
        """
        Look for the first input voila file with the property which exists. This does NOT validate the property
        is the same across all files where lsv id exists.

        cast should be specified when trying to retrieve generator style properties, because the generator
        will be processed inside this function, and attempting to iterate the generator outside of
        will cause an unrelated error because you have exited the "with" block.
        :return: property value
        """
        config = rna_voila.config.ViewConfig()
        propval = None
        for f in config.cov_files:
            with ViewPsi(f) as m:
                if propval is None:
                    try:

                        propval = getattr(m.lsv(self.lsv_id), prop)
                        if cast:
                            propval = cast(propval)

                    except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                        pass
                if propval is not None:
                    return propval

    def _get_prop_multi(self, prop):
        """
        Look for some attribute in all input files, and return a group:property dict
        it is assummed that the property itself will be read as a group:property dict
        from each individual file.
        :return: property value
        """
        config = rna_voila.config.ViewConfig()
        propval = None
        groups_to_props = {}
        for f in config.voila_files:
            with ViewPsi(f) as m:
                try:
                    propval = dict(getattr(m.lsv(self.lsv_id), prop))

                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass
                if propval is not None:
                    for key in propval:
                        groups_to_props[key] = propval[key]

        return groups_to_props


    def get_attr(self, attr):
        """
        For attributes that exist is each het file, this will get all values and confirm they're all equal.
        :param attr: attribute found in het voila file.
        :return: attribute value
        """

        voila_files = rna_voila.config.ViewConfig().voila_files
        s = set()
        for f in voila_files:
            with self.view_class(f) as m:
                try:
                    inner = m.lsv(self.lsv_id)
                    s.add(getattr(inner, attr))
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    pass
        if len(s) == 0:
            raise LsvIdNotFoundInAnyVoilaFile
        assert len(s) == 1, s
        return s.pop()

    # @property
    # def reference_exon(self):
    #     self.sg
    #     return self.get_attr('reference_exon')

    @property
    def target(self):
        return self.get_attr('target')

    @property
    def source(self):
        return self.get_attr('source')

    @property
    def binary(self):
        return self.get_attr('binary')

    @property
    def complex(self):
        return self.get_attr('complex')

    @property
    def gene_id(self):
        return self.get_attr('gene_id')

    @property
    def a5ss(self):
        return self.get_attr('a5ss')

    @property
    def a3ss(self):
        return self.get_attr('a3ss')

    @property
    def exon_skipping(self):
        return self.get_attr('exon_skipping')

    @property
    def exon_count(self):
        return self.get_attr('exon_count')

    @property
    def lsv_type(self):
        return self.get_attr('lsv_type')

    @property
    def intron_retention(self):
        return self.get_attr('intron_retention')


    @property
    def junctions(self):
        """
        Finds first file with junctions list for this specific lsv id. This does NOT validate the junctions list
        is the same across all files where lsv id exists.
        :return: numpy array
        """
        return self._get_prop('junctions')
