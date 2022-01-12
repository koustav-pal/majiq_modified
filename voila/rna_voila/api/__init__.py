from rna_voila.api.matrix_hdf5 import Psi, DeltaPsi, Heterogen
import rna_voila.api.splice_graph_sql as splice_graph_sql
import rna_voila.api.splice_graph_zarr as splice_graph_zarr
import rna_voila.api.view_matrix_hdf5 as view_matrix_hdf5
import rna_voila.api.view_matrix_zarr as view_matrix_zarr
from rna_voila.api.view_matrix import find_analysis_type, get_mixed_analysis_type_str

"""

What's going on in this file?

Here is where we abstract away dealing with choice of using zarr or sql for the splicegraph, and zarr or hdf5 for the 
voila files. This is mainly for an intermediate period where we would still like voila to support usage of either standard...

"""


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


# class Matrix(DeltaPsi, Psi, Heterogen):
#     pass

class _ViewMatrixHDF5(view_matrix_hdf5.ViewPsi):
    pass

class _ViewPsiHDF5(view_matrix_hdf5.ViewPsi):
    pass

class _ViewPsisHDF5(view_matrix_hdf5.ViewPsis):
    pass

class _ViewDeltaPsiHDF5(view_matrix_hdf5.ViewDeltaPsi):
    pass

class _ViewHeterogenHDF5(view_matrix_hdf5.ViewHeterogen):
    pass

class _ViewHeterogensHDF5(view_matrix_hdf5.ViewHeterogens):
    pass

class _ViewMatrixZarr(view_matrix_zarr.ViewPsi):
    pass

class _ViewPsiZarr(view_matrix_zarr.ViewPsi):
    pass

class _ViewPsisZarr(view_matrix_zarr.ViewPsis):
    pass

class _ViewDeltaPsiZarr(view_matrix_zarr.ViewDeltaPsi):
    pass

class _ViewHeterogenZarr(view_matrix_zarr.ViewHeterogen):
    pass

class _ViewHeterogensZarr(view_matrix_zarr.ViewHeterogens):
    pass



def SpliceGraph(*args, **kwargs):
    if kwargs.get('build', False) is True:
        del kwargs['build']
        return _SpliceGraphSQL(*args, **kwargs)
    from rna_voila.api.view_splice_graph import get_sg_format_str
    return _SpliceGraphSQL(*args, **kwargs) if get_sg_format_str() == 's' else _SpliceGraphZarr(*args, **kwargs)

def ViewMatrix(*args, **kwargs):
    from rna_voila.api.view_matrix import get_matrix_format_str
    return _ViewMatrixHDF5(*args, **kwargs) if get_matrix_format_str() == 'h' else _ViewMatrixZarr(*args, **kwargs)

def ViewPsi(*args, **kwargs):
    from rna_voila.api.view_matrix import get_matrix_format_str
    return _ViewPsiHDF5(*args, **kwargs) if get_matrix_format_str() == 'h' else _ViewPsiZarr(*args, **kwargs)

def ViewPsis(*args, **kwargs):
    from rna_voila.api.view_matrix import get_matrix_format_str
    return _ViewPsisHDF5(*args, **kwargs) if get_matrix_format_str() == 'h' else _ViewPsisZarr(*args, **kwargs)

def ViewDeltaPsi(*args, **kwargs):
    from rna_voila.api.view_matrix import get_matrix_format_str
    return _ViewDeltaPsiHDF5(*args, **kwargs) if get_matrix_format_str() == 'h' else _ViewDeltaPsiZarr(*args, **kwargs)

def ViewHeterogen(*args, **kwargs):
    from rna_voila.api.view_matrix import get_matrix_format_str
    return _ViewHeterogenHDF5(*args, **kwargs) if get_matrix_format_str() == 'h' else _ViewHeterogenZarr(*args, **kwargs)

def ViewHeterogens(*args, **kwargs):
    from rna_voila.api.view_matrix import get_matrix_format_str
    return _ViewHeterogensHDF5(*args, **kwargs) if get_matrix_format_str() == 'h' else _ViewHeterogensZarr(*args, **kwargs)