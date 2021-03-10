from rna_voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
import rna_voila.api.splice_graph_sql as splice_graph_sql
import rna_voila.api.splice_graph_netcdf as splice_graph_netcdf


class _SpliceGraphSQL(
    splice_graph_sql.Genes,
    splice_graph_sql.Junctions,
    splice_graph_sql.Exons,
    splice_graph_sql.IntronRetentions,
    splice_graph_sql.AltStarts,
    splice_graph_sql.AltEnds):
    pass

class _SpliceGraphNetCDF(
    splice_graph_netcdf.Genes,
    splice_graph_netcdf.Junctions,
    splice_graph_netcdf.Exons,
    splice_graph_netcdf.IntronRetentions,
    splice_graph_netcdf.AltStarts,
    splice_graph_netcdf.AltEnds):
    pass


class Matrix(DeltaPsi, Psi, Heterogen):
    pass



def SpliceGraph(*args, **kwargs):
    from rna_voila.api.view_splice_graph import get_sg_format_str
    if get_sg_format_str() == 's':
        return _SpliceGraphSQL(*args, **kwargs)
    else:
        return _SpliceGraphNetCDF(*args, **kwargs)