import os
from operator import itemgetter
from pathlib import Path
import rna_majiq as nm
from typing import Union
import zarr
import rna_voila.config

from rna_voila.constants import EXEC_DIR


class SpliceGraphZarr:
    def __init__(self):
        """
        Splice graph class to handle metadata and table information retrieval.
        :param filename: database file name
        :param delete: delete existing database file.
        """

        # try:
        #     zarr_file = zarr_file.decode("utf-8")
        # except AttributeError:
        #     zarr_file = str(zarr_file)



        self.conn = rna_voila.config.ViewConfig().sg_zarr
        self.lsvs = rna_voila.config.ViewConfig().sg_lsvs

        #self.exp_reads = nm.SpliceGraphReads.from_zarr(sgc_files)
        self.exp_reads = rna_voila.config.ViewConfig().sgc_zarr

        self._genome = None
        self._experiment_names = None
        self._file_version = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        pass
        # self.conn.commit()
        # self.conn.close()

    def get_gene_idx(self, genes: nm.Genes, gene_id: str) -> int:
        return genes[gene_id]

    def gene_regions_slice(self, regions, gene_idx: int) -> slice:
        return regions.slice_for_gene(gene_idx)

    @property
    def genome(self):
        """
        Get genome from data base.
        :return: string
        """
        if self._genome is None:
            self._genome = ''
            # query = self.conn.execute('SELECT name FROM genome')
            # genome, = query.fetchone()
            # if not genome:
            #     self._genome = ''
            # else:
            #     self._genome = genome

        return self._genome

    @genome.setter
    def genome(self, g):
        """
        Write genome to databse.
        :param g: genome value
        :return: None
        """
        pass
        # self.conn.execute('''
        #                     INSERT
        #                     INTO genome (name)
        #                     VALUES (?)
        #                     ''', (g,))

    @property
    def experiment_names(self):
        """
        Get experiment names for database.
        :return: list of strings
        """
        if self._experiment_names is None:
            self._experiment_names = self.exp_reads.prefixes

        return self._experiment_names

    @experiment_names.setter
    def experiment_names(self, names):
        """
        Write experiment names to database
        :param names: list of experiment names
        :return: None
        """
        pass
        # self.conn.executemany('''
        #                     INSERT
        #                     INTO experiment (name)
        #                     VALUES (?)
        #                     ''', tuple((n,) for n in names))

    @property
    def file_version(self):
        """
        Get file version from database.
        :return: string
        """
        return -1


    @file_version.setter
    def file_version(self, version):
        """
        Write file version from constants to database.
        :param version:
        :return: None
        """
        pass
        # self.conn.execute('''
        #                     INSERT
        #                     INTO main.file_version (value)
        #                     VALUES (?)
        #                     ''', version)

    def _iter_results(self, coordinate, fieldnames):
        """
        Query blocks of 100 entries from database and then convert each line to a dictionary.

        :param query: database query
        :param fieldnames: table fieldnames
        :return: generator of dictionaries
        """
        for d in getattr(self.conn, coordinate).gene_idx:
            yield {fn: getattr(d, fn).values for fn in fieldnames}




gene_fieldnames = ('gene_id', 'gene_name', 'strand', 'contig_idx')
junc_fieldnames = ('start', 'end', 'simplified')
junc_reads_fieldnames = ('reads', 'experiment_name')
exon_fieldnames = ('start', 'end', 'annotated_start', 'annotated_end', 'annotated')
ir_fieldnames = ('gene_id', 'start', 'end', 'has_reads', 'annotated', 'is_simplified', 'is_constitutive')
ir_reads_fieldnames = ('reads', 'experiment_name')
alt_starts_fieldnames = ('coordinate',)
alt_ends_fieldnames = ('coordinate',)


class Genes(SpliceGraphZarr):

    def gene_ids(self):
        return self.conn.genes.gene_id

    def genes(self):
        """
        Get all the genes in the database.
        :return: generator of dictionaries
        """

        for idx in getattr(self.conn, 'genes').gene_idx:

            ret = {}
            ret['id'] = self.conn.genes.gene_id[idx]
            ret['name'] = self.conn.genes.gene_name[idx]
            ret['chromosome'] = self.conn.genes.contigs.seqid[self.conn.genes.contig_idx[idx]]
            ret['strand'] = self.conn.genes.strand[idx].decode()
            yield ret

    def gene(self, gene_id):
        """
        Get gene with gene id supplied.
        :param gene_id: gene id
        :return: dictionary
        """
        idx = self.conn.genes[gene_id]

        ret = {}
        ret['id'] = gene_id
        ret['name'] = self.conn.genes.gene_name[idx]
        ret['chromosome'] = self.conn.genes.contigs.seqid[self.conn.genes.contig_idx[idx]]
        ret['strand'] = self.conn.genes.strand[idx].decode()

        return ret

    def gene_overlap(self, gene_id):
        # query = self.conn.execute("""
        #             SELECT gene_overlap.gene_id_1, g1.name, gene_overlap.gene_id_2, g2.name FROM gene_overlap
        #             INNER JOIN gene g1 on g1.id = gene_overlap.gene_id_1
		# 			INNER JOIN gene g2 on g2.id = gene_overlap.gene_id_2
        #             WHERE gene_id_1 = ?
        #             OR gene_id_2 = ?
        #         """, (gene_id, gene_id))
        # return [(g[0], g[1],) if g[0] != gene_id else (g[2], g[3],) for g in query.fetchall()]
        return []



class Exons(SpliceGraphZarr):
    def exons(self, gene_id):
        """
        Get all exons for specified gene id
        :param gene_id: gene id
        :return: list of exons
        """
        # first get gene_idx



        # gene_idx = self.conn.genes.where(self.conn.genes.gene_id == gene_id, drop=True).gene_idx.values[0]
        #
        # exons = self.conn.exons.where(self.conn.exons.gene_idx == gene_idx, drop=True)

        eslice = self.conn.exons.df.isel(exon_idx=self.conn.exons.slice_for_gene(self.conn.genes[gene_id]))

        for exon_idx in eslice.exon_idx:

            result = {}
            result['start'] = self.conn.exons.start[exon_idx]
            result['end'] = self.conn.exons.end[exon_idx]
            result['annotated_start'] = self.conn.exons.annotated_start[exon_idx]
            result['annotated_end'] = self.conn.exons.annotated_end[exon_idx]

            result['annotated'] = int(not self.conn.exons.is_denovo([exon_idx]))
            result['gene_id'] = gene_id

            result['_exon_idx'] = int(exon_idx)

            yield result




class Junctions(SpliceGraphZarr):



    def junctions(self, gene_id, omit_simplified=False, clin_denovo_conns=None):
        """
        Get list of junctions for specified gene id.
        :param gene_id: gene id
        :return: list of junction dictionaries
        """
        jslice = self.conn.junctions.slice_for_gene(self.conn.genes[gene_id])

        for junc_idx in self.conn.junctions.gj_idx[jslice]:
            if omit_simplified and self.conn.junctions.simplified[junc_idx]:
                continue

            result = {}
            result['start'] = self.conn.junctions.start[junc_idx]
            result['end'] = self.conn.junctions.end[junc_idx]
            result['gene_id'] = gene_id
            result['annotated'] = int(not self.conn.junctions.denovo[junc_idx])
            result['is_simplified'] = int(self.conn.junctions.simplified[junc_idx])

            result['has_reads'] = int(True)
            result['_junc_idx'] = junc_idx
            result['clin_denovo'] = int(True) if (clin_denovo_conns and (result['start'], result['end']) in clin_denovo_conns) else int(False)

            yield result


    def junction_reads_exp(self, junction, experiment_names):
        """
        for a junction and a set of experiment names, get a list of reads.
        :param junction: junction dictionary
        :param experiment_names: list of experiment names
        :return: list of reads dictionaries
        """

        for experiment_name in experiment_names:
            if not experiment_name:
                continue
            result = {}
            result['experiment_name'] = experiment_name
            reads = self.exp_reads.junctions_reads[junction['_junc_idx'], self.exp_reads.prefixes.index(experiment_name)].values
            #assert self.exp_reads_conns[experiment_name].junctions_hash[junc_idx]

            result['reads'] = int(reads)

            yield result


class IntronRetentions(SpliceGraphZarr):
    def intron_retentions(self, gene_id, omit_simplified=False, clin_denovo_conns=None):
        """
        Get all intron retentions for a gene id.
        :param gene_id: gene id
        :return: list of intron retention dictionaries
        """
        jslice = self.conn.introns.slice_for_gene(self.conn.genes[gene_id])

        for junc_idx in self.conn.introns.gi_idx[jslice]:
            if omit_simplified and self.conn.introns.simplified[junc_idx]:
                continue

            result = {}
            result['start'] = self.conn.introns.start[junc_idx]
            result['end'] = self.conn.introns.end[junc_idx]
            result['gene_id'] = gene_id
            result['annotated'] = int(not self.conn.introns.denovo[junc_idx])
            result['is_simplified'] = int(self.conn.introns.simplified[junc_idx])

            result['has_reads'] = int(True)

            result['_junc_idx'] = junc_idx
            result['clin_denovo'] = int(True) if (clin_denovo_conns and (result['start'], result['end']) in clin_denovo_conns) else int(False)

            yield result

    def intron_retention_reads_exp(self, intron, experiment_names):
        """
        For an intron retention dictionary and a set of experiment names, get a list of intron retention reads.
        :param ir: intron retention dictionary
        :param experiment_names: list of experiment reads.
        :return: list of intron retention reads
        """

        for experiment_name in experiment_names:
            if not experiment_name:
                continue
            result = {}
            result['experiment_name'] = experiment_name


            reads = self.exp_reads.introns_reads[intron['_junc_idx'], self.exp_reads.prefixes.index(experiment_name)].values
            result['reads'] = int(reads)

            yield result

class AltStarts(SpliceGraphZarr):
    def alt_starts(self, gene_id):
        """
        Get alternate starts for a specific gene.
        :param gene_id: gene id
        :return: list of alt starts dictionary
        """
        return []
        # query = self.conn.execute('''
        #                             SELECT coordinate
        #                             FROM alt_start
        #                             WHERE gene_id=?
        #                             ''', (gene_id,))
        # return self._iter_results(query, alt_starts_fieldnames)


class AltEnds(SpliceGraphZarr):
    def alt_ends(self, gene_id):
        """
        Get alternate ends for specific gene.
        :param gene_id: gene id
        :return: List of alt ends dictionary
        """
        return []
        # query = self.conn.execute('''
        #                             SELECT coordinate
        #                             FROM alt_end
        #                             WHERE gene_id=?
        #                             ''', (gene_id,))
        # return self._iter_results(query, alt_ends_fieldnames)
