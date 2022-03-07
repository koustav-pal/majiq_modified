
from rna_voila.api.matrix_hdf5 import MatrixHdf5
from rna_voila.exceptions import FoundNoSpliceGraphFile, FoundMoreThanOneSpliceGraph, \
    MixedAnalysisTypeVoilaFiles, FoundMoreThanOneVoilaFile, AnalysisTypeNotFound
from rna_voila import constants

def find_analysis_type(voila_files, cov_files):
    """
    Find the analysis type from the voila files.
    :param voila_files: list of voila files.
    :return: String
    """
    analysis_type = None

    if voila_files:

        for mf in voila_files:

            with MatrixHdf5(mf) as m:

                if analysis_type is None:
                    analysis_type = m.analysis_type

                if analysis_type != m.analysis_type:
                    raise MixedAnalysisTypeVoilaFiles()

        if not analysis_type:
            raise AnalysisTypeNotFound()

        if analysis_type in (constants.ANALYSIS_DELTAPSI,) and len(voila_files) > 1:
            raise FoundMoreThanOneVoilaFile()

        return analysis_type

    else:

        for mf in cov_files:
            if mf.name.endswith('.psicov') or mf.name.endswith('.dpsicov'):

                this_file_type = None
                if mf.name.endswith('.psicov'):
                    this_file_type = constants.ANALYSIS_PSI
                elif mf.name.endswith('.dpsicov'):
                    this_file_type = constants.ANALYSIS_DELTAPSI


                if analysis_type is None:
                    analysis_type = this_file_type

                if analysis_type != this_file_type:
                    raise MixedAnalysisTypeVoilaFiles()

        if not analysis_type:
            raise AnalysisTypeNotFound()

        if analysis_type in (constants.ANALYSIS_DELTAPSI,) and len(voila_files) > 1:
            raise FoundMoreThanOneVoilaFile()

        return analysis_type



def get_mixed_analysis_type_str(voila_files, cov_files):
    types = {'psi': 0, 'delta_psi': 0, 'het': 0}

    if voila_files:
        for mf in voila_files:

            with MatrixHdf5(mf) as m:

                if m.analysis_type == constants.ANALYSIS_PSI:
                    types['psi'] += 1

                elif m.analysis_type == constants.ANALYSIS_DELTAPSI:
                    types['delta_psi'] += 1

                elif m.analysis_type == constants.ANALYSIS_HETEROGEN:
                    types['het'] += 1


    else:
        for mf in cov_files:
            if mf.endswith('.psicov'):
                types['psi'] += 1

    strsout = []
    if types['psi']:
        strsout.append("PSIx%d" % types['psi'])
    if types['delta_psi']:
        strsout.append("dPSIx%d" % types['delta_psi'])
    if types['het']:
        strsout.append("HETx%d" % types['het'])
    return ' '.join(strsout)


def get_matrix_format_str():
    from rna_voila.config import ViewConfig
    config = ViewConfig()

    if config.voila_file is not None:
        return 'h'
    elif config.cov_file is not None:
        return 'z'
    raise NotImplementedError("Invalid Matrix File Format")