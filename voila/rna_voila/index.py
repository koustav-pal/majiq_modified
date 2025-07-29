import json
from itertools import chain

import h5py
import numpy as np

from rna_voila import constants
from rna_voila.api import ViewHeterogens, ViewDeltaPsi, ViewPsi, ViewPsis, ViewHeterogen
from rna_voila.api.matrix_hdf5 import MatrixHdf5
from rna_voila.api.view_splice_graph import ViewSpliceGraph
from rna_voila.config import ViewConfig
from rna_voila.exceptions import UnknownAnalysisType, IndexNotFound, UnknownIndexFieldType, GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile
from rna_voila.api.matrix_utils import unpack_means, unpack_bins, generate_excl_incl, generate_means, \
    generate_high_probability_non_changing, generate_variances, generate_standard_deviations
from rna_voila.vlsv import matrix_area
from rna_voila.voila_log import voila_log
import multiprocessing

import time
import hashlib
import os
from tqdm import tqdm
import pickle

lsv_filters = ['a5ss', 'a3ss', 'exon_skipping', 'target', 'source', 'binary', 'complex', 'intron_retention']
psi_keys = ['lsv_id', 'gene_id', 'gene_name'] + lsv_filters
dpsi_keys = ['lsv_id', 'gene_id', 'gene_name', 'excl_incl', 'dpsi_threshold', 'confidence_threshold'] + lsv_filters
het_keys = ['lsv_id', 'gene_id', 'gene_name', 'dpsi_threshold', 'stat_threshold'] + lsv_filters

skip_strict_indexing = False

class Index:
    def __init__(self, force_create=False, voila_files=None):
        """
        Factory class to generate the index for the supplied analysis type.
        """

        self.config = ViewConfig()
        analysis_type = self.config.analysis_type
        self.force_index = force_create or self.config.force_index
        self.index_voila_files = voila_files

        global skip_strict_indexing
        skip_strict_indexing = not self.config.enable_type_indexing

        # The case where there's no analysis type, we're assuming this is splice graph only.
        if analysis_type:

            if analysis_type == constants.ANALYSIS_PSI:
                index = self._psi
            elif analysis_type == constants.ANALYSIS_DELTAPSI:
                index = self._deltapsi
            elif analysis_type == constants.ANALYSIS_HETEROGEN:
                index = self._heterogen
            else:
                raise UnknownAnalysisType(analysis_type)

            index()

    @staticmethod
    def _index_in_voila(voila_file, remove_index=False):
        """
        Check if index has already been created in voila file. If remove_index has been set, then attempt to use the
        h5py api to remove the index dataset from the file.

        If the case of multiple voila files, we hash all of them and check that the hash of that group of hashes
        matches what is stored in the voila file. Here we are assuming that the 'order' of the files is not important
        (for example, if we check for a match in the first file but the previous index was stored in a different
        file, the index will be rebuilt again for the same data)
        :param voila_file:
        :param remove_index:
        :return:
        """

        if os.path.exists(voila_file):
            with MatrixHdf5(voila_file, 'r', pre_config=True) as m:
                return m.has_index()
        return False


    @staticmethod
    def _get_files_hash(voila_files):
        if ViewConfig().index_file:
            # if we use a separate file for indexing, we can verify the hash of all inputs in the verification
            h = hashlib.sha1()
            for filename in sorted(voila_files):
                with open(filename, 'rb') as f:
                    h.update(f.read())
            return h.hexdigest()
        else:
            # for now we get a hash based on the combined name of all files
            # (this is because we can not easily do it based on the content because we change the content for this process)
            return hashlib.sha1(''.join(voila_files).encode('utf-8')).hexdigest()

    @staticmethod
    def _write_index(voila_file, voila_index, dtype):
        """
        Helper method to write voila index to voila file using specific numpy data type.

        :param voila_file: location and name of voila file.
        :param voila_index: Array of index data.
        :param dtype: numpy data type string.
        :return: None
        """
        voila_index = np.array(voila_index, dtype=np.dtype(dtype))
        voila_files = ViewConfig().voila_files

        with MatrixHdf5(voila_file, 'a', pre_config=True) as m:
            hashval = Index._get_files_hash(voila_files)
            m.write_index(voila_index, hashval)

    @staticmethod
    def _get_voila_index_file():
        c = ViewConfig()
        if c.index_file:
            return c.index_file
        else:
            return c.voila_file

    @staticmethod
    def _create_dtype(voila_index):
        first_row = voila_index[0]
        dtype = []
        for field in first_row:
            if isinstance(field, str):
                dtype.append(len(field))
            elif isinstance(field, bool):
                dtype.append('?')
            elif isinstance(field, np.float64):
                dtype.append('f4')
            else:
                raise UnknownIndexFieldType(field)

        for row in voila_index:
            for idx, field in enumerate(row):
                if dtype[idx] not in ('?', 'f4'):
                    dtype[idx] = max(len(field), dtype[idx])

        dtype = [d if d in ('?', 'f4') else 'S' + str(d) for d in dtype]
        dtype = ','.join(dtype)

        return dtype

    @staticmethod
    def _heterogen_pool_add_index(args):
        """
        Multithread inner function for each iteration of _heterogen loop below
        """
        lsv_id, g, vfs, q = args
        config = ViewConfig()

        g_dpsi_thresh = np.array([-1])

        with ViewHeterogens() as m:
            het = m.lsv(lsv_id)
            gene_id = het.gene_id
            gene_name = g[gene_id]
            g_stats_thresh = np.array([2 for _ in m.stat_names])

            if skip_strict_indexing:
                lsv_f = [True for _ in lsv_filters]
            else:
                lsv_f = [bool(getattr(het, f)) for f in lsv_filters]

        for vf in config.voila_files:
            if vfs and vf not in vfs:
                continue
            with ViewHeterogen(vf) as m:
                try:
                    het = m.lsv(lsv_id)

                    dpsi_thresh = het.dpsi_median()
                    g_dpsi_thresh = np.maximum(np.max(np.abs(dpsi_thresh)), g_dpsi_thresh)

                    # get min val for each stat
                    stats_thresh = het.junction_stats
                    g_stats_thresh = np.minimum(np.min(stats_thresh, axis=0), g_stats_thresh)
                except (GeneIdNotFoundInVoilaFile, LsvIdNotFoundInVoilaFile):
                    continue

        g_dpsi_thresh = g_dpsi_thresh.tolist()[0]
        g_dpsi_thresh = json.dumps(g_dpsi_thresh)

        g_stats_thresh = g_stats_thresh.tolist()
        g_stats_thresh = json.dumps(g_stats_thresh)

        row = (lsv_id, gene_id, gene_name, g_dpsi_thresh, g_stats_thresh)

        # For some reason, numpy needs these in tuples.
        row = tuple(chain(row, lsv_f))

        if q:
            q.put(row)
        return row

    def _heterogen(self):
        """
        Create index for heterogen analysis type.  This is an odd case as there can be more then one voila file.  We
        create the index in the first voila file.  First is determined by ordering the file name and location.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = self.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            with ViewSpliceGraph() as sg:
                g = sg.gene_ids2gene_names

            if config.nproc > 1:
                manager = multiprocessing.Manager()
                q = manager.Queue()

                with ViewHeterogens() as m:
                    lsv_ids = [(x, g, self.index_voila_files, q) for x in m.lsv_ids()]


                p = multiprocessing.Pool(config.nproc)
                work_size = len(lsv_ids)

                # voila_index = p.map(self._heterogen_pool_add_index, zip(lsv_ids, range(work_size), repeat(work_size)))
                voila_index = p.map_async(self._heterogen_pool_add_index, lsv_ids, chunksize=config.parallel_chunksize if config.parallel_chunksize > 0 else None)

                # monitor loop
                while True:
                    if voila_index.ready():
                        break
                    else:
                        size = q.qsize()
                        print("Indexing LSV IDs: %d / %d" % (size, work_size))
                        time.sleep(2)


                voila_index = voila_index.get()

            else:
                voila_index = []
                with ViewHeterogens() as m:
                    for x in m.lsv_ids():
                        args = (x, g, self.index_voila_files, None)
                        voila_index.append(self._heterogen_pool_add_index(args))

            dtype = self._create_dtype(voila_index)
            log.info('Writing index: ' + voila_file)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _deltapsi_pool_add_index(args):
        """
        Multithread inner function for each iteration of _deltapsi loop below
        """
        lsv_id, g, q = args

        with ViewDeltaPsi() as m:
            dpsi = m.lsv(lsv_id)
            excl_incl = max(abs(a - b) for a, b in dpsi.excl_incl)
            gene_id = dpsi.gene_id

            dpsi_thresh = dpsi.means
            dpsi_thresh = np.abs(dpsi_thresh)
            dpsi_thresh = dpsi_thresh.tolist()
            dpsi_thresh = json.dumps(dpsi_thresh)

            # Calculate a list of confidence for 10 (0 thru 1) values of threshold. This is done because we
            # don't want to store the bins data in the index.  Although, this might be a good idea once the
            # index gets large enough.
            bins = dpsi.bins
            confidence_thresh = list(max(matrix_area(b, x) for b in bins) for x in np.linspace(0, 1, 10))
            confidence_thresh = json.dumps(confidence_thresh)

            gene_name = g[gene_id]

            row = (lsv_id, gene_id, gene_name, excl_incl, dpsi_thresh, confidence_thresh)

            if skip_strict_indexing:
                lsv_f = [True for _ in lsv_filters]
            else:
                lsv_f = [bool(getattr(dpsi, f)) for f in lsv_filters]


        # For some reason, numpy needs these in tuples.
        row = tuple(chain(row, lsv_f))
        if q:
            q.put(row)
        return row

    def _deltapsi(self):
        """
        Generates index for delta psi analysis type.

        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = self.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            with ViewSpliceGraph() as sg:
                g = sg.gene_ids2gene_names

            if config.nproc > 1:

                manager = multiprocessing.Manager()
                q = manager.Queue()

                with ViewDeltaPsi() as m:
                    lsv_ids = [(x, g, q) for x in m.lsv_ids()]
                p = multiprocessing.Pool(config.nproc)
                work_size = len(lsv_ids)

                voila_index = p.map_async(self._deltapsi_pool_add_index, lsv_ids, chunksize=config.parallel_chunksize if config.parallel_chunksize > 0 else None)

                # monitor loop
                while True:
                    if voila_index.ready():
                        break
                    else:
                        size = q.qsize()
                        print("Indexing LSV IDs: %d / %d" % (size, work_size))
                        time.sleep(2)


                voila_index = voila_index.get()

            else:
                voila_index = []
                with ViewDeltaPsi() as m:
                    for x in m.lsv_ids():
                        args = (x, g, None)
                        voila_index.append(self._deltapsi_pool_add_index(args))

            dtype = self._create_dtype(voila_index)
            log.info('Writing index: ' + voila_file)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _psi_pool_add_index(args):
        """
        Multithread inner function for each iteration of _psi loop below
        """
        lsv_id, g, m, q = args
        lsv = m.lsv(lsv_id)
        gene_id = lsv.gene_id

        gene_name = g[gene_id]

        row = (lsv_id, gene_id, gene_name)

        if skip_strict_indexing:
            lsv_f = [True for _ in lsv_filters]
        else:
            lsv_f = [bool(getattr(lsv, f)) for f in lsv_filters]

        # For some reason, numpy needs these in tuples.
        row = tuple(chain(row, lsv_f))

        if q:
            q.put(row)
        return row



    @staticmethod
    def _row_data(gene_id, keys):
        """
        For each row in index, zip list of keys with values in the row.
        :param gene_id: gene id
        :param keys: index field names
        :return:
        """

        index_file = Index._get_voila_index_file()

        try:
            gene_id = gene_id.encode('utf-8')
        except AttributeError:
            pass

        with h5py.File(index_file, 'r') as h:

            try:

                for row in h['index'].value:
                    if gene_id is None or gene_id == row[1]:
                        yield dict(zip(keys, row))

            except KeyError:
                raise IndexNotFound()

    @classmethod
    def psi(cls, gene_id=None):
        """
        Get PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._row_data(gene_id, psi_keys)

    @classmethod
    def delta_psi(cls, gene_id=None):
        """
        Get Delta PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._row_data(gene_id, dpsi_keys)

    @classmethod
    def heterogen(cls, gene_id=None):
        """
        Get Heterogen index data in a dictionary for each row.
        :return:
        """

        yield from cls._row_data(gene_id, het_keys)

class HDF5Index(Index):
    def __init__(self, force_create=False, voila_files=None):
        super(HDF5Index, self).__init__(force_create, voila_files)

    @staticmethod
    def _write_index(voila_file, voila_index, dtype):
        """
        Helper method to write voila index to voila file using specific numpy data type.

        :param voila_file: location and name of voila file.
        :param voila_index: Array of index data.
        :param dtype: numpy data type string.
        :return: None
        """
        voila_index = np.array(voila_index, dtype=np.dtype(dtype))

        voila_files = ViewConfig().voila_files

        with h5py.File(voila_file, 'a') as h:
            if 'index' in h:
                del h['index']
            if 'input_hash' in h:
                del h['input_hash']
            h.create_dataset('index', voila_index.shape, data=voila_index)
            hashval = Index._get_files_hash(voila_files)
            h.create_dataset("input_hash", (1,), dtype="S40", data=(hashval.encode('utf-8'),))

    @staticmethod
    def _index_in_voila(voila_file, remove_index=False):
        """
        Check if index has already been created in voila file. If remove_index has been set, then attempt to use the
        h5py api to remove the index dataset from the file.

        If the case of multiple voila files, we hash all of them and check that the hash of that group of hashes
        matches what is stored in the voila file. Here we are assuming that the 'order' of the files is not important
        (for example, if we check for a match in the first file but the previous index was stored in a different
        file, the index will be rebuilt again for the same data)
        :param voila_file:
        :param remove_index:
        :return:
        """

        with h5py.File(voila_file, 'a') as h:
            index_in_h = 'index' in h

            if remove_index and index_in_h:
                voila_log().info('Removing index from HDF5')
                del h['index']

            # voila_files = ViewConfig().voila_files

            # if 'input_hash' in h:
            #     prior = h.get('input_hash')[0].decode('utf-8')
            #     new = Index._get_files_hash(voila_files)
            #     index_in_h = (prior == new)

            return index_in_h

    def _psi(self):
        """
        Create index to PSI analysis type.
        :return: None
        """

        config = ViewConfig()
        log = voila_log()
        force_index = remove_index = self.force_index
        voila_file = self._get_voila_index_file()

        if not self._index_in_voila(voila_file, remove_index) or force_index:

            log.info('Creating index: ' + voila_file)

            with ViewSpliceGraph() as sg:
                g = sg.gene_ids2gene_names

            if config.nproc > 1:

                manager = multiprocessing.Manager()
                q = manager.Queue()

                with ViewPsis() as m:
                    lsv_ids = [(x, g, m, q) for x in m.lsv_ids()]
                p = multiprocessing.Pool(config.nproc)
                work_size = len(lsv_ids)

                voila_index = p.map_async(self._psi_pool_add_index, lsv_ids, chunksize=config.parallel_chunksize if config.parallel_chunksize > 0 else None)

                # monitor loop
                while True:
                    if voila_index.ready():
                        break
                    else:
                        size = q.qsize()
                        print("Indexing LSV IDs: %d / %d" % (size, work_size))
                        time.sleep(2)


                voila_index = voila_index.get()

            else:
                voila_index = []
                with ViewPsis() as m:
                    for x in m.lsv_ids():
                        args = (x, g, m, None)
                        voila_index.append(self._psi_pool_add_index(args))

            dtype = self._create_dtype(voila_index)
            log.info('Writing index: ' + voila_file)
            self._write_index(voila_file, voila_index, dtype)
        else:
            log.info('Using index: ' + voila_file)

    @staticmethod
    def _row_data(gene_id, keys):
        """
        For each row in index, zip list of keys with values in the row.
        :param gene_id: gene id
        :param keys: index field names
        :return:
        """

        index_file = Index._get_voila_index_file()

        try:
            gene_id = gene_id.encode('utf-8')
        except AttributeError:
            pass

        with MatrixHdf5(index_file, 'r') as m:

            try:
                for row in m.get_index():
                    if gene_id is None or gene_id == row[1]:
                        yield dict(zip(keys, row))

            except KeyError:
                raise IndexNotFound()

    @classmethod
    def psi(cls, gene_id=None):
        """
        Get PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._row_data(gene_id, psi_keys)

    @classmethod
    def delta_psi(cls, gene_id=None):
        """
        Get Delta PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._row_data(gene_id, dpsi_keys)

    @classmethod
    def heterogen(cls, gene_id=None):
        """
        Get Heterogen index data in a dictionary for each row.
        :return:
        """

        yield from cls._row_data(gene_id, het_keys)

index_constants = None
class ZarrIndex:

    cache = None

    def __init__(self):
        pass
        #super(ZarrIndex, self).__init__(force_create, voila_files)

    @classmethod
    def init_cache(cls, dpsi, het, index_file, total, force=False):
        if cls.cache is None:
            cls.cache = {}
            if index_file and os.path.exists(index_file) and not (force is True):
                voila_log().info(f'Using Cache: {index_file}')
                with open(index_file, 'rb') as f:
                    cls.cache = pickle.load(f)
            else:
                if multiprocessing.current_process().name == "MainProcess":
                    voila_log().info('Generating Caches...')
                    pbar = tqdm(total=total)
                    for gene_id, rows in cls._row_data(None, pbar, dpsi, het):
                        cls.cache[gene_id] = rows
                    pbar.close()
                    voila_log().info('Generating Caches...Done')
                    if index_file:
                        voila_log().info(f'Saving Cache: {index_file}')
                        with open(index_file, 'wb') as f:
                            pickle.dump(cls.cache, f)


    @classmethod
    def _cached_row_data(cls, gene_id, dpsi=False, het=False):
        cls.init_cache(dpsi, het, None, None)
        if gene_id:
            yield from cls.cache.get(gene_id, [])
        else:
            for rows in cls.cache.values():
                yield from rows

    @staticmethod
    def get_dpsi_precalc(cov, lsvs):
        """
        for DPSI, both the bootstrap_posterior.mean and interval_probability are very slow to be sliced and calculated
        for each lsv for each thread. To improve performance and usability, we precalculate 11 matrices by moving all
        lsv event slices into a large numpy matrix of size (lsvs, max size of lsv), zero padded. The slicing process
        is still unfortunately slow, but in total making all 11 matrices takes about a minute, and vastly speeds up
        the actual index process after it.

        Other improvements still to-do:
        1) this calculation is done in every thread before the processing starts. It should be able to be done once
        and then saved as a constant somehow for all threads to use. Shared memory?
        2) the process of moving each lsv event slice into the calc matrix is done in a for loop, which is slow. Could
        not find satisfactory answer of how to do it without a loop, need to revisit.
        """

        slices = tuple([slice(int(x), int(y)) for x, y in zip(lsvs.ec_idx_start, lsvs.ec_idx_end)])
        max_slicelen = max(s.stop - s.start for s in slices)
        calc_arr = np.zeros(shape=(len(lsvs.ec_idx_start), max_slicelen))
        means = cov.bootstrap_posterior.mean[0].as_numpy()
        for i, _slice in enumerate(slices):
            calc_arr[i, 0:_slice.stop - _slice.start] = means[_slice]
        calc_arr_abs = np.abs(calc_arr)
        calc_arr_max = np.max(calc_arr_abs, axis=1)

        confid_probs = []

        for confid in np.linspace(0, 1, 10):
            confid_prob = cov.bootstrap_posterior.interval_probability(-confid, confid)[0].as_numpy()
            confid_probs_max = np.zeros(shape=(len(lsvs.ec_idx_start), max_slicelen))
            for i, _slice in enumerate(slices):
                confid_probs_max[i, 0:_slice.stop - _slice.start] = confid_prob[_slice]
            confid_probs_max = np.max(confid_probs_max, axis=1)
            confid_probs.append(confid_probs_max)

        return calc_arr_abs, calc_arr_max, confid_probs

    @staticmethod
    def _get_index_constants(dpsi=False):
        global index_constants

        if index_constants is None:

            config = ViewConfig()

            cov = config.primary_cov_zarr

            sg = config.sg_zarr

            events = cov.get_events(sg.introns, sg.junctions)

            lsvs = sg.exon_connections.lsvs()
            if dpsi:
                dpsi_precalc = ZarrIndex.get_dpsi_precalc(cov, lsvs)
            else:
                dpsi_precalc = None

            index_constants = (events, lsvs, dpsi_precalc)

        return index_constants

    @classmethod
    def _pool_row_data(cls, args):

        het, dpsi, q, gene_id = args

        config = ViewConfig()

        cov = config.primary_cov_zarr

        sg = config.sg_zarr

        events, lsvs, dpsi_precalc = cls._get_index_constants(dpsi=dpsi)


        gene_idx = sg.genes[gene_id]
        gene_name = sg.genes.gene_name[gene_idx]
        events_slice = events.slice_for_gene(gene_idx)
        lsv_ids = sg.exon_connections.event_id(events.ref_exon_idx[events_slice], events.event_type[events_slice])
        results = []


        if config.enable_type_indexing:
            has_intron = sg.exon_connections.has_intron(events.ref_exon_idx[events_slice],
                                                        events.event_type[events_slice])
            is_source_LSV = sg.exon_connections.is_source_LSV(events.ref_exon_idx[events_slice],
                                                              events.event_type[events_slice])
            is_target_LSV = sg.exon_connections.is_target_LSV(events.ref_exon_idx[events_slice],
                                                              events.event_type[events_slice])

        for idx, e_idx in enumerate(range(events_slice.start, events_slice.stop)):
            res = {}
            lsv_id = lsv_ids[idx].encode()
            clin_denovo = ViewConfig().clin_controls.get(gene_id, False) and (
                    lsv_id in ViewConfig().clin_controls.get(gene_id, ()))
            res.update(dict(
                lsv_id=lsv_id,
                gene_id=gene_id.encode(),
                gene_name=gene_name.encode(),
                clin_denovo=clin_denovo
            ))

            # lsv_ids = sg.exon_connections.event_id(events.ref_exon_idx[events_slice], events.event_type[events_slice])

            ec_idx_s = lsvs.connections_slice_for_event(e_idx)

            # try:
            #     gene_id = gene_id.encode('utf-8')
            # except AttributeError:
            #     pass

            # here "IDX" is relative for event slice, to use the more efficient method of indexing
            # all data points at once from above
            if config.enable_type_indexing:
                first_ec_idx = lsvs.ec_idx_start[e_idx]

                if het:
                    if not cov.any_passed[first_ec_idx].all():
                        continue
                elif dpsi:
                    if not cov.event_passed[:, first_ec_idx].all():
                        continue
                else:
                    if not cov.event_passed[first_ec_idx, :].all():  # not really sure why these are reversed...
                        continue

                res.update(dict(
                    a5ss=lsvs.event_legacy_a5ss(e_idx),
                    a3ss=lsvs.event_legacy_a3ss(e_idx),
                    exon_skipping=lsvs.event_has_alt_exons(e_idx),
                    target=is_target_LSV[idx],
                    source=is_source_LSV[idx],
                    binary=lsvs.event_size[e_idx] == 2,
                    complex=lsvs.event_size[e_idx] != 2,
                    intron_retention=has_intron[idx]
                ))

            if dpsi:
                calc_arr_abs, calc_arr_max, confid_probs = dpsi_precalc

                excl_incl = calc_arr_max[e_idx]

                dpsi_thresh = calc_arr_abs[e_idx]
                dpsi_thresh = json.dumps(dpsi_thresh.tolist())

                confidence_thresh = [x[e_idx] for x in confid_probs]
                confidence_thresh = json.dumps(confidence_thresh)

                res.update(dict(
                    excl_incl=excl_incl,
                    dpsi_threshold=dpsi_thresh,
                    confidence_threshold=confidence_thresh
                ))

            if het:
                g_dpsi_thresh = 1.0

                # minimum value for every experiment and every junction in lsv
                # will be length <number of stats>
                g_stats_thresh = np.amin(cov.approximate_pvalue[:, ec_idx_s, :], axis=(0, 1,))

                # g_dpsi_thresh = g_dpsi_thresh.tolist()
                g_dpsi_thresh = json.dumps(g_dpsi_thresh)

                g_stats_thresh = list(g_stats_thresh.to_series())
                g_stats_thresh = json.dumps(g_stats_thresh)

                res.update(dict(
                    dpsi_threshold=g_dpsi_thresh,
                    stat_threshold=g_stats_thresh
                ))


            results.append(res)

        if q:
            q.put((gene_id, results))
        return (gene_id, results)


    @classmethod
    def _row_data(cls, _gene_id, pbar=None, dpsi=False, het=False):
        """
        For each row in index, zip list of keys with values in the row.
        :param gene_id: gene id
        :param keys: index field names
        :return:
        """

        #index_file = Index._get_voila_index_file()
        #print(keys)
        config = ViewConfig()
        cov = config.primary_cov_zarr
        sg = config.sg_zarr


        if _gene_id:
            gene_ids = [_gene_id]
        else:
            gene_ids = sg.genes.gene_id

        if config.nproc > 1:

            ctx = multiprocessing.get_context("forkserver")
            manager = ctx.Manager()
            q = manager.Queue()

            p = ctx.Pool(config.nproc)
            work_size = len(gene_ids)
            _gene_ids = [(het, dpsi, q, g) for g in gene_ids]

            #func_with_constants = partial(ZarrIndex._pool_row_data, het, dpsi, confid_probs, q)
            voila_index = p.map_async(cls._pool_row_data, _gene_ids,
                                      chunksize=config.parallel_chunksize if config.parallel_chunksize > 0 else None)

            # monitor loop
            while True:
                if voila_index.ready():
                    break
                else:
                    size = q.qsize()
                    if pbar:
                        pbar.update(size - pbar.n)
                    time.sleep(2)

            voila_index = voila_index.get()
            for result in voila_index:
                yield result

        else:

            for i, gene_id in enumerate(gene_ids):
                pbar.update(i - pbar.n)
                yield cls._pool_row_data((
                    het, dpsi, None, gene_id
                ))



    @classmethod
    def psi(cls, gene_id=None):
        """
        Get PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._cached_row_data(gene_id)

    @classmethod
    def delta_psi(cls, gene_id=None):
        """
        Get Delta PSI index data in a dictionary for each row.
        :param gene_id: Filter output by specific gene.
        :return: Generator
        """

        yield from cls._cached_row_data(gene_id, dpsi=True)


    @classmethod
    def heterogen(cls, gene_id=None):
        """
        Get Heterogen index data in a dictionary for each row.
        :return:
        """

        yield from cls._cached_row_data(gene_id, het=True)

def get_index(*args, **kwargs):
    config = ViewConfig()
    if config.voila_file:
        return HDF5Index(*args, **kwargs)
    else:
        return None

def get_index_class():
    config = ViewConfig()
    if config.voila_file:
        return HDF5Index
    else:
        return ZarrIndex