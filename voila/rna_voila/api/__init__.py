from rna_voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
import rna_voila.api.splice_graph_sql as splice_graph_sql
import rna_voila.api.splice_graph_zarr as splice_graph_zarr


class _SpliceGraphSQL(
    splice_graph_sql.Genes,
    splice_graph_sql.Junctions,
    splice_graph_sql.Exons,
    splice_graph_sql.IntronRetentions,
    splice_graph_sql.AltStarts,
    splice_graph_sql.AltEnds):
    pass

class _SpliceGraphZarr(
    splice_graph_zarr.Genes,
    splice_graph_zarr.Junctions,
    splice_graph_zarr.Exons,
    splice_graph_zarr.IntronRetentions,
    splice_graph_zarr.AltStarts,
    splice_graph_zarr.AltEnds):
    pass


class Matrix(DeltaPsi, Psi, Heterogen):
    pass



def SpliceGraph(*args, **kwargs):
    if kwargs.get('build', False) is True:
        return _SpliceGraphSQL(*args, **kwargs)
    from rna_voila.api.view_splice_graph import get_sg_format_str
    if get_sg_format_str() == 's':
        return _SpliceGraphSQL(*args, **kwargs)
    else:
        return _SpliceGraphZarr(*args, **kwargs)