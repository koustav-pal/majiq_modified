import configparser
import inspect
import sqlite3
from collections import namedtuple
from pathlib import Path
import sys, os
from rna_voila import constants
from csv import DictReader

from rna_voila.api import ViewPsi, SpliceGraph, find_analysis_type, get_mixed_analysis_type_str, view_matrix_zarr, ViewMatrix
from rna_voila.api.matrix_hdf5 import MatrixHdf5

from rna_voila.api.splice_graph_lr import SpliceGraphLR
from rna_voila.clin import parse_clin_tsv
from rna_voila.exceptions import FoundNoSpliceGraphFile, FoundMoreThanOneSpliceGraph, \
    MixedAnalysisTypeVoilaFiles, FoundMoreThanOneVoilaFile, AnalysisTypeNotFound
from rna_voila.voila_log import voila_log
from rna_voila.api.view_matrix import open_cov_wrapper

import rna_majiq as nm



_log_keys = ['logger', 'silent']
_sys_keys = ['nproc', 'debug']
_global_keys = ['analysis_type', 'memory_map_hdf5', 'groups_to_voilas', 'license', 'preserve_handles_hdf5',
                'parallel_chunksize']
_v3_keys = ['cov_file', 'cov_files', 'zarr_file', 'sgc_files', 'sg_zarr', 'sgc_zarr', 'cov_zarr', 'cov_zarr_combined', 'primary_cov_zarr', 'lsvid2lsvidx', 'lsvidx2lsvid', 'lsvtype_cache', 'module_cache', 'cov_cache']

_ViewConfig = namedtuple('ViewConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['voila_file', 'voila_files',
                                        'splice_graph_file',
                                        'force_index', 'port', 'host', 'web_server', 'index_file', 'gunicorn_worker_class',
                                        'num_web_workers', 'strict_indexing', 'enable_type_indexing', 'splice_graph_only',
                                        'enable_passcode', 'ignore_inconsistent_group_errors', 'only_index',
                                        'enable_het_comparison_chooser', 'long_read_file',  'disable_reads',
                                        'group_order_override', 'clin_controls_file', 'clin_controls',
                                        'psicov_grouping_file'])
_ViewConfig.__new__.__defaults__ = (None,) * len(_ViewConfig._fields)
_TsvConfig = namedtuple('TsvConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['file_name', 'voila_files', 'voila_file',
                                      'splice_graph_file',
                                      'non_changing_threshold', 'threshold', 'show_all',
                                      'probability_threshold', 'gene_ids', 'gene_names', 'lsv_ids',
                                      'lsv_types', 'strict_indexing', 'show_read_counts',
                                      'ignore_inconsistent_group_errors', 'non_changing_pvalue_threshold',
                                      'non_changing_within_group_iqr',
                                      'non_changing_between_group_dpsi', 'changing_pvalue_threshold',
                                      'changing_between_group_dpsi', 'show_per_sample_psi'])
_TsvConfig.__new__.__defaults__ = (None,) * len(_TsvConfig._fields)
_ClassifyConfig = namedtuple('ClassifyConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['directory', 'voila_files',
                                      'voila_file', 'splice_graph_file',
                                      'decomplexify_psi_threshold', 'decomplexify_deltapsi_threshold',
                                      'decomplexify_reads_threshold', 'gene_ids',
                                      'keep_constitutive', 'keep_no_lsvs_modules', 'only_binary',
                                      'untrimmed_exons', 'putative_multi_gene_regions',
                                        'probability_changing_threshold',
                                        'probability_non_changing_threshold', 'show_all',
                                        'non_changing_pvalue_threshold', 'non_changing_within_group_iqr',
                                        'non_changing_between_group_dpsi', 'changing_pvalue_threshold',
                                        'changing_between_group_dpsi', 'changing_between_group_dpsi_secondary',
                                        'keep_no_lsvs_junctions', 'debug_num_genes', 'overwrite', 'output_mpe',
                                        'heatmap_selection', 'enabled_outputs',
                                        'ignore_inconsistent_group_errors', 'disable_metadata',
                                        'show_read_counts', 'cassettes_constitutive_column',
                                        'non_changing_median_reads_threshold', 'permissive_event_non_changing_threshold',
                                        'include_change_cases', 'junc_gene_dist_column', 'show_per_sample_psi',
                                                                    'psicov_grouping_file'])
_ClassifyConfig.__new__.__defaults__ = (None,) * len(_ClassifyConfig._fields)
_FilterConfig = namedtuple('FilterConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['directory', 'voila_files',
                                            'voila_file', 'splice_graph_file',
                                            'gene_ids', 'overwrite',
                                            'gene_ids_file', 'lsv_ids', 'lsv_ids_file', 'voila_files_only',
                                            'splice_graph_only',
                                            'changing_threshold', 'non_changing_threshold',
                                            'probability_changing_threshold',
                                            'probability_non_changing_threshold', 'changing', 'non_changing'])
_FilterConfig.__new__.__defaults__ = (None,) * len(_FilterConfig._fields)
_SplitterConfig = namedtuple('SplitterConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['directory', 'voila_files',
                                        'voila_file', 'splice_graph_file', 'num_divisions', 'copy_only', 'overwrite'])
_SplitterConfig.__new__.__defaults__ = (None,) * len(_SplitterConfig._fields)
_RecombineConfig = namedtuple('RecombineConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['directories', 'directory'])
_RecombineConfig.__new__.__defaults__ = (None,) * len(_RecombineConfig._fields)
_LongReadsConfig = namedtuple('LongReadsConfig', _global_keys + _sys_keys + _log_keys + _v3_keys + ['voila_file', 'lr_gtf_file',
                                            'lr_tsv_file', 'splice_graph_file', 'output_file', 'gene_id',
                                            'only_update_psi'])
_LongReadsConfig.__new__.__defaults__ = (None,) * len(_LongReadsConfig._fields)

# global config variable to act as the singleton instance of the config.
this_config = None
this_group_names_to_voila_files = None
this_group_names_to_cov_files = {}
this_cov_zarr_combined = None




def find_splice_graph_file(_vs, splice_graph_only):
    """
    Function that located all splice graphs from a list of files and directories.
    :param vs: list of files and directories.
    :return: location of splice graph file
    """

    sg_files = set()
    found_sql_version = False
    zarr_files = set()
    sgc_files = set()

    def iter_dir(vs):

        nonlocal found_sql_version

        for v in vs:


            v = Path(v)

            if v.is_file():

                if v.parts[-1].endswith('.sql'):
                    sg_files.add(v)
                    found_sql_version = True

            elif v.is_dir():

                if all(x in os.listdir(v) for x in ('contigs', 'exons', 'genes', 'introns', 'junctions')):
                    zarr_files.add(v)

                if 'sg_reads' in os.listdir(v):
                    sgc_files.add(v)

                iter_dir(v.iterdir())

    iter_dir(_vs)

    if len(sg_files) == 0 and len(zarr_files) == 0:
        raise FoundNoSpliceGraphFile()

    if len(sg_files) > 1:
        raise FoundMoreThanOneSpliceGraph()

    if found_sql_version and len(zarr_files):
        raise Exception("Found mixed .sql and .zarr inputs")

    if len(zarr_files) > 1:
        raise Exception("Found multiple .zarr inputs")

    if len(zarr_files) and not len(sgc_files):
        if not splice_graph_only:
            raise Exception("Could not find any .sgc experiment data inputs")

    if found_sql_version:
        sg_file = sg_files.pop()
        return sg_file.resolve()


    zarr_file = zarr_files.pop()

    return (zarr_file.resolve(), [f.resolve() for f in sgc_files])


def find_voila_files(vs):
    """
    Find all voila files in files and directories.
    :param vs: list of files and directories.
    :return: list of voila files
    """

    voila_files = []
    voila_files_to_group_names = {}

    for v in vs:
        v = Path(v)

        if v.is_file() and v.name.endswith('.voila'):

            try:
                with MatrixHdf5(v, pre_config=True) as m:
                    voila_files.append(v)
                    voila_files_to_group_names[v] = m.group_names[0]
            except OSError:
                voila_log().warning('Error opening voila file %s , skipping this file' % str(v))
                pass

        elif v.is_dir():
            x, x2 = find_voila_files(v.iterdir())
            voila_files = [*voila_files, *x]
            voila_files_to_group_names.update(x2)

    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.
    voila_files.sort()

    return voila_files, voila_files_to_group_names

def reorder_voila_files(voila_files, group_order_override, voila_files_to_group_names):
    try:
        return sorted(voila_files, key=lambda x: group_order_override.index(voila_files_to_group_names[x]))
    except ValueError:
        voila_log().critical("Could not match group order override to provided voila files")
        raise



def _is_cov_psi(path):
    if type(path) != Path:
        path = Path(path)
    return path.is_dir() and 'psi_coverage' in os.listdir(path)

def _is_cov_dpsi(path):
    if type(path) != Path:
        path = Path(path)
    return path.is_dir() and 'deltapsi' in os.listdir(path)

def _is_cov_het(path):
    if type(path) != Path:
        path = Path(path)
    return path.is_dir() and 'heterogen' in os.listdir(path)

def cov_file_analysis_type(path):
    if _is_cov_psi(path):
        return constants.ANALYSIS_PSI
    elif _is_cov_dpsi(path):
        return constants.ANALYSIS_DELTAPSI
    elif _is_cov_het(path):
        return constants.ANALYSIS_HETEROGEN
    else:
        raise AnalysisTypeNotFound()



def find_cov_files(vs):
    """
    Find all voila files in files and directories.
    :param vs: list of files and directories.
    :return: list of voila files
    """

    cov_files = []
    #cov_files_to_group_names = {}

    for v in vs:
        v = Path(v)

        """
        Group names addition are commented below as the values are only used in multiPSI mode. These lines are quite slow,
        and multiPSI mode needs to be re-evaluated for V3 usage anyway. 
        """

        if _is_cov_psi(v):
            cov_files.append(v)
            #cov_files_to_group_names[v] = view_matrix_zarr.preconfig_group_names_psi(v)[0]

        if _is_cov_dpsi(v):
            cov_files.append(v)
            #cov_files_to_group_names[v] = view_matrix_zarr.preconfig_group_names_dpsi(v)[0]

        if _is_cov_het(v):
            cov_files.append(v)
            #cov_files_to_group_names[v] = view_matrix_zarr.preconfig_group_names_het(v)[0]

        elif v.is_dir():
            x = find_cov_files(v.iterdir())
            cov_files = [*cov_files, *x]
            #cov_files_to_group_names.update(x2)

    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.

    cov_files.sort()

    return cov_files

"""
voila_files = []
    voila_files_to_group_names = {}

    for v in vs:
        v = Path(v)

        if v.is_file() and v.name.endswith('.voila'):

            try:
                with ViewMatrix(v, pre_config=True) as m:
                    voila_files.append(v)
                    voila_files_to_group_names[v] = m.group_names[0]
            except OSError:
                voila_log().warning('Error opening voila file %s , skipping this file' % str(v))
                pass

        elif v.is_dir():
            x, x2 = find_voila_files(v.iterdir())
            voila_files = [*voila_files, *x]
            voila_files_to_group_names.update(x2)

    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.
    voila_files.sort()

"""

def write(args):

    """
    Write command line argmuments into a ini file.
    :param args: argparse object
    :return: None
    """

    voila_log().info('config file: ' + constants.CONFIG_FILE)


    attrs = inspect.getmembers(args, lambda a: not inspect.isbuiltin(a))
    attrs = (a for a in attrs if not a[0].startswith('_'))
    attrs = dict(attrs)

    if not args.func.__name__ in ('recombine', 'longReadsInputsToLongReadsVoila'):
        sg_file = find_splice_graph_file(args.files, getattr(args, 'splice_graph_only', False))
    else:
        sg_file = None

    if (hasattr(args, 'splice_graph_only') and getattr(args, 'splice_graph_only', False)) or args.func.__name__ in (
    'recombine', 'longReadsInputsToLongReadsVoila'):

        analysis_type = ''
        voila_files = []
        cov_files = []

    else:

        group_order_override = getattr(args, "group_order_override", None)

        voila_files, voila_files_to_group_names = find_voila_files(args.files)
        cov_files = find_cov_files(args.files)

        if voila_files and cov_files:
            raise Exception("Found both voila files and cov files in provided paths. Only one type should be provided.")

        global this_group_names_to_voila_files
        this_group_names_to_voila_files = {v: k for k, v in voila_files_to_group_names.items()}

        global this_group_names_to_cov_files
        this_group_names_to_cov_files = {}




        if group_order_override:
            voila_files = reorder_voila_files(voila_files, group_order_override, voila_files_to_group_names)
        if args.func.__name__ in ("Filter", "Classify", 'splitter', 'recombine'):
            analysis_type, split_paths = get_mixed_analysis_type_str(voila_files, cov_files)
            psi_cov_files = split_paths['psi']
        else:
            analysis_type = find_analysis_type(voila_files, cov_files)
            psi_cov_files = cov_files if analysis_type in (constants.ANALYSIS_PSI,) else []


    # raise multi-file error if trying to run voila in TSV mode with multiple input files
    # (currently, multiple input is only supported in View mode)
    if analysis_type in (constants.ANALYSIS_PSI, ) and \
            args.func.__name__ not in ['run_service', 'Classify', 'splitter', 'recombine'] and len(voila_files) > 1:

        raise FoundMoreThanOneVoilaFile()

    # attributes that don't need to be in the ini file
    for remove_key in ['files', 'func', 'logger']:
        try:
            del attrs[remove_key]
        except KeyError:
            pass

    config_parser = configparser.ConfigParser()
    files = 'FILES'
    settings = 'SETTINGS'
    filters = 'FILTERS'

    config_parser.add_section('META')
    config_parser.set('META', 'func', args.func.__name__)

    # Get filters from arguments, add them to the appropriate section, and remove them from arguments.
    for lsv_filter in ['lsv_types', 'lsv_ids', 'gene_ids', 'gene_names']:
        if lsv_filter in attrs and attrs[lsv_filter]:
            if not config_parser.has_section(filters):
                config_parser.add_section(filters)
            config_parser.set(filters, lsv_filter, '\n'.join(attrs[lsv_filter]))

            del attrs[lsv_filter]

    if attrs.get('enable_passcode', False):
        import random, string
        attrs['enable_passcode'] = ''.join(random.choices(string.ascii_letters + string.digits, k=20))
    else:
        attrs['enable_passcode'] = ''

    # Get settings from arguments.
    config_parser.add_section(settings)
    for key, value in attrs.items():
        if isinstance(value, int) or isinstance(value, float) or value:
            config_parser.set(settings, key, str(value))

    if args.func.__name__ == "Classify":
        # check if default = default for some args depending on the type of classify being run.
        # we store a default of none and apply the actual default if the user does not specify another value
        if not config_parser.has_option(settings, 'decomplexify_psi_threshold'):
            if config_parser.getboolean(settings, 'show_all'):
                config_parser.set(settings, 'decomplexify_psi_threshold', '0.0')

        if not config_parser.has_option(settings, 'decomplexify_deltapsi_threshold'):
            if config_parser.getboolean(settings, 'show_all'):
                config_parser.set(settings, 'decomplexify_deltapsi_threshold', '0.0')

    config_parser.set(settings, 'analysis_type', analysis_type)

    # Get files from arguments
    config_parser.add_section(files)

    if voila_files:
        config_parser.set(files, 'voila', '\n'.join(str(m) for m in voila_files))
    else:
        config_parser.set(files, 'majiq', '\n'.join(str(m) for m in cov_files))

    if type(sg_file) is tuple:
        config_parser.set(files, 'zarr_file', str(sg_file[0]))
        config_parser.set(files, 'sgc_files', '\n'.join(str(m) for m in sg_file[1]))
    else:
        config_parser.set(files, 'splice_graph', str(sg_file))

    # Write ini file.
    with open(constants.CONFIG_FILE, 'w') as configfile:
        config_parser.write(configfile)

    # dump config file to log
    with open(constants.CONFIG_FILE, 'r') as configfile:
        for line in configfile:
            voila_log().debug('CFG| ' + line[:-1])

    if config_parser.has_option(settings, 'long_read_file') and config_parser.get(settings, 'long_read_file'):
        voila_log().info(f"Parsing long reads file: {config_parser.get(settings, 'long_read_file')}")
        SpliceGraphLR(config_parser.get(settings, 'long_read_file'))

def parse_psicov_grouping_file(psicov_grouping_file):

    grouping = {}
    with open(psicov_grouping_file, 'r') as fr:
        reader = DictReader(fr, delimiter='\t')
        headers = reader.fieldnames
        if not all(x in headers for x in ('group', 'prefix',)):
            voila_log().critical("psicov-grouping-file should have at least 'group' and 'prefix' columns")
            sys.exit(1)
        for row in reader:
            if row['group'] not in grouping:
                grouping[row['group']] = []
            grouping[row['group']].append(row['prefix'])
    return grouping


def _getInputFilesSet(config_parser, view=False, cov_multiarray=False):
    files = {}
    settings = dict(config_parser['SETTINGS'])


    if 'splice_graph' in config_parser['FILES']:
        files['splice_graph_file'] = config_parser['FILES']['splice_graph']
    else:
        files['zarr_file'] = config_parser['FILES']['zarr_file']
        files['sg_zarr'] = nm.SpliceGraph.from_zarr(files['zarr_file'])
        mask = nm.SpliceGraphMask.from_arrays(files['sg_zarr'].introns, files['sg_zarr'].junctions)
        files['module_cache'] = files['sg_zarr'].modules(mask)

        files['sgc_files'] = config_parser['FILES']['sgc_files'].split('\n')
        files['sgc_zarr'] = nm.SpliceGraphReads.from_zarr(files['sgc_files'])

    files['cov_cache'] = {}
    if 'voila' in config_parser['FILES']:
        files['voila_files'] = config_parser['FILES']['voila'].split('\n')
        files['voila_file'] = config_parser['FILES']['voila'].split('\n')[0]
    else:
        files['cov_files'] = config_parser['FILES']['majiq'].split('\n')
        files['cov_file'] = config_parser['FILES']['majiq'].split('\n')[0]

        if 'index_file' not in settings and view:
            settings['index_file'] = str(Path(files['zarr_file']).parent / 'voila_index.hdf5')

        if settings.get('splice_graph_only', 'False') != 'True':

            psi_cov_files = [f for f in files['cov_files'] if f.endswith('.psicov')]
            if psi_cov_files:
                global this_cov_zarr_combined
                global this_group_names_to_cov_files
                this_cov_zarr_combined = nm.PsiCoverage.from_zarr(psi_cov_files)
                if len(psi_cov_files) > 1:
                    prefixes = this_cov_zarr_combined.prefixes
                    if files.get('psicov_grouping_file', None):
                        group_defs = parse_psicov_grouping_file(files['psicov_grouping_file'])

                        prefix2cov = {}
                        for cov_file, prefix in zip(psi_cov_files, prefixes):
                            prefix2cov[prefix] = cov_file

                        for group_name, prefixes in group_defs.items():
                            this_group_names_to_cov_files[group_name] = []
                            for prefix in prefixes:
                                if prefix not in prefix2cov:
                                    voila_log().critical(
                                        f"For group {group_name} in psicov-grouping-file, prefix {prefix} was not found in any specified psicoverage files")
                                    sys.exit(1)
                                this_group_names_to_cov_files[group_name].append(prefix2cov[prefix])

                    else:
                        # no group definition provided, set each cov file to it's group

                        for cov_file, prefix in zip(psi_cov_files, prefixes):
                            this_group_names_to_cov_files[prefix] = cov_file

            if cov_multiarray:
                # pre load psi, dpsi, and het cov files into separate keys
                _psi = list(filter(_is_cov_psi, files['cov_files']))
                _dpsi = list(filter(_is_cov_dpsi, files['cov_files']))
                _het = list(filter(_is_cov_het, files['cov_files']))

                files['cov_zarr'] = dict(
                    psi=nm.PsiCoverage.from_zarr(_psi) if _psi else None,
                    dpsi=nm.DeltaPsiDataset.from_zarr(_dpsi) if _dpsi else None,
                    het=nm.HeterogenDataset.from_zarr(_het) if _het else None,
                )
                for cov_file in files['cov_files']:
                    files['cov_zarr'][cov_file] = open_cov_wrapper(cov_file)

                files['cov_zarr_combined'] = this_cov_zarr_combined
                files['primary_cov_zarr'] = None

                if 'zarr_file' in files:
                    lsvid2lsvidx = {}
                    lsvidx2lsvid = {}
                    if _psi:
                        lsvid2lsvidx = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], files['cov_zarr']['psi'], lsvid2lsvidx)
                        lsvidx2lsvid = view_matrix_zarr.get_lsvidx2lsvid(files['sg_zarr'], files['cov_zarr']['psi'], lsvidx2lsvid)
                    if _dpsi:
                        lsvid2lsvidx = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], files['cov_zarr']['dpsi'], lsvid2lsvidx)
                        lsvidx2lsvid = view_matrix_zarr.get_lsvidx2lsvid(files['sg_zarr'], files['cov_zarr']['dpsi'], lsvidx2lsvid)
                    if _het:
                        lsvid2lsvidx = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], files['cov_zarr']['het'], lsvid2lsvidx)
                        lsvidx2lsvid = view_matrix_zarr.get_lsvidx2lsvid(files['sg_zarr'], files['cov_zarr']['het'], lsvidx2lsvid)

                    files['lsvid2lsvidx'] = lsvid2lsvidx
                    files['lsvidx2lsvid'] = lsvidx2lsvid

            else:
                if settings['analysis_type'] == constants.ANALYSIS_PSI:
                    if this_cov_zarr_combined:
                        files['cov_zarr_combined'] = this_cov_zarr_combined
                    if settings.get('psicov_grouping_file', None):
                        files['cov_zarr'] = {}

                        for group_name, cov_files in this_group_names_to_cov_files.items():
                            files['cov_zarr'][group_name] = nm.PsiCoverage.from_zarr(cov_files)
                        files['primary_cov_zarr'] = files['cov_zarr_combined']
                    else:
                        files['cov_zarr'] = {'psi': nm.PsiCoverage.from_zarr(files['cov_files'])}
                        files['primary_cov_zarr'] = files['cov_zarr']['psi']

                elif settings['analysis_type'] == constants.ANALYSIS_DELTAPSI:
                    files['cov_zarr'] = {'dpsi': nm.DeltaPsiDataset.from_zarr(files['cov_file'])}
                    files['primary_cov_zarr'] = files['cov_zarr']['dpsi']
                elif settings['analysis_type'] == constants.ANALYSIS_HETEROGEN:
                    files['cov_zarr'] = {'het': nm.HeterogenDataset.from_zarr(files['cov_file'])}
                    files['primary_cov_zarr'] = files['cov_zarr']['het']


    analysis_type_to_key = {
        constants.ANALYSIS_PSI: 'psi',
        constants.ANALYSIS_DELTAPSI: 'dpsi',
        constants.ANALYSIS_HETEROGEN: 'het'
    }

    if 'cov_files' in files and 'zarr_file' in files:
        if settings.get('splice_graph_only', 'False') != 'True':

            if not cov_multiarray:
                if settings.get('psicov_grouping_file', None):
                    files['lsvid2lsvidx'] = {}
                    for group_name, cov_arr in files['cov_zarr'].items():
                        files['lsvid2lsvidx'] = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], cov_arr, files['lsvid2lsvidx'])
                    files['lsvtype_cache'] = view_matrix_zarr.get_lsvtype_cache(files['sg_zarr'], files['primary_cov_zarr'])
                else:
                    files['lsvid2lsvidx'] = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], files['primary_cov_zarr'], {})
                    files['lsvtype_cache'] = view_matrix_zarr.get_lsvtype_cache(files['sg_zarr'], files['primary_cov_zarr'])
                if view:
                    import rna_voila.index

                    rna_voila.index.ZarrIndex.init_cache(
                        dpsi=settings['analysis_type'] == constants.ANALYSIS_DELTAPSI,
                        het=settings['analysis_type'] == constants.ANALYSIS_HETEROGEN,
                        index_file=settings['index_file'],
                        total=len(files['lsvid2lsvidx']),
                        force=settings.get('force_index', "False") == "True"
                    )
                # else:
                #     files['lsvid2lsvidx'] = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], files['cov_zarr'], {})
                #     files['lsvtype_cache'] = view_matrix_zarr.get_lsvtype_cache(files['sg_zarr'], files['cov_zarr'])

    # if settings.get('psicov_grouping_file', None):
    #     settings['psicov_grouping'] = parse_psicov_grouping_file(settings['psicov_grouping_file'])

    return files, settings


class ViewConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the HTML view.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')

            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)




            files, settings = _getInputFilesSet(config_parser, view=True)



            for int_key in ['nproc', 'port', 'num_web_workers', 'parallel_chunksize']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for bool_key in ['force_index', 'silent', 'debug', 'strict_indexing', 'enable_type_indexing',
                             'ignore_inconsistent_group_errors', 'enable_het_comparison_chooser', 'memory_map_hdf5',
                             'disable_reads', 'preserve_handles_hdf5', 'only_index']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            # singleton data store properties
            if settings.get('memory_map_hdf5', False) and not 'index_file' in settings:
                voila_log().critical('To use hdf5 memory map performance mode, you must specify --index-file as well')
                sys.exit(1)


            if settings.get('clin_controls_file', None):
                settings['clin_controls'] = parse_clin_tsv(settings['clin_controls_file'])
            else:
                settings['clin_controls'] = {}

            if files.get('voila_file', None):
                settings['groups_to_voilas'] = this_group_names_to_voila_files
            else:
                settings['groups_to_voilas'] = this_group_names_to_cov_files


            this_config = _ViewConfig(**{**files, **settings})

        return this_config


class TsvConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the TSV output.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            # files = {
            #     'voila_files': config_parser['FILES']['voila'].split('\n'),
            #     'voila_file': config_parser['FILES']['voila'].split('\n')[0],
            # }
            #
            # if 'splice_graph' in config_parser['FILES']:
            #     files['splice_graph_file'] = config_parser['FILES']['splice_graph']
            # else:
            #     files['zarr_file'] = config_parser['FILES']['zarr_file']
            #     files['sgc_files'] = config_parser['FILES']['sgc_files'].split('\n')
            #
            # settings = dict(config_parser['SETTINGS'])

            files, settings = _getInputFilesSet(config_parser, view=False)

            for int_key in ['nproc', 'parallel_chunksize']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['non_changing_threshold', 'threshold', 'probability_threshold',
                              'non_changing_pvalue_threshold',
                              'non_changing_within_group_iqr', 'non_changing_between_group_dpsi',
                              'changing_pvalue_threshold', 'changing_between_group_dpsi']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['show_all', 'silent', 'debug', 'strict_indexing', 'show_read_counts',
                             'ignore_inconsistent_group_errors', 'memory_map_hdf5', 'show_per_sample_psi',
                             'preserve_handles_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _TsvConfig(**{**files, **settings, **filters})

        return this_config

class ClassifyConfig:
    def __new__(cls, *args, **kwargs):
        """

        """

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            files, settings = _getInputFilesSet(config_parser, cov_multiarray=True)

            for int_key in ['nproc', 'keep_constitutive', 'decomplexify_reads_threshold', 'debug_num_genes',
                            'non_changing_median_reads_threshold', 'parallel_chunksize']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['decomplexify_psi_threshold', 'decomplexify_deltapsi_threshold',
                              'probability_changing_threshold',
                              'probability_non_changing_threshold', 'non_changing_pvalue_threshold',
                              'non_changing_within_group_iqr', 'non_changing_between_group_dpsi',
                              'changing_pvalue_threshold', 'changing_between_group_dpsi',
                              'changing_between_group_dpsi_secondary', 'permissive_event_non_changing_threshold']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'keep_no_lsvs_modules', 'only_binary', 'untrimmed_exons', 'overwrite',
                             'putative_multi_gene_regions', 'show_all', 'keep_no_lsvs_junctions', 'output_mpe',
                             'ignore_inconsistent_group_errors', 'disable_metadata', 'show_read_counts',
                             'cassettes_constitutive_column', 'include_change_cases', 'junc_gene_dist_column',
                             'memory_map_hdf5', 'show_per_sample_psi', 'preserve_handles_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            if settings['decomplexify_reads_threshold'] == 0:
                voila_log().warning("--decomplexify-reads-threshold 0 is not recommended and not tested!")

            # implications
            if settings['putative_multi_gene_regions']:
                settings['keep_constitutive'] = True
            if settings['keep_constitutive']:
                settings['keep_no_lsvs_modules'] = True
                settings['keep_no_lsvs_junctions'] = True



            if not settings['putative_multi_gene_regions']:
                settings['enabled_outputs'] = ['summary', 'events', 'junctions', 'heatmap']
                if settings['output_mpe']:
                    settings['enabled_outputs'].append('mpe')
                    settings['keep_constitutive'] = True
                    settings['keep_no_lsvs_modules'] = True
                    settings['keep_no_lsvs_junctions'] = True

            if not settings['show_all']:
                if 'HET' not in settings['analysis_type'] and 'dPSI' not in settings['analysis_type']:
                    voila_log().warning("Only PSI files were provided, so unable to detect changing events. "
                                        "Enabling --show-all automatically. ")
                    settings['show_all'] = True

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _ClassifyConfig(**{**files, **settings, **filters})

        return this_config

class FilterConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the TSV output.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """
        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            files = {
                'voila_files': config_parser['FILES']['voila'].split('\n'),
                'voila_file': config_parser['FILES']['voila'].split('\n')[0],
                'splice_graph_file': config_parser['FILES']['splice_graph']
            }

            settings = dict(config_parser['SETTINGS'])
            for int_key in ['nproc', 'parallel_chunksize']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['non_changing_threshold', 'changing_threshold', 'probability_changing_threshold',
                              'probability_non_changing_threshold']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'overwrite', 'voila_files_only', 'splice_graph_only', 'changing', 'non_changing',
                             'memory_map_hdf5', 'preserve_handles_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _FilterConfig(**{**files, **settings, **filters})

        return this_config

class SplitterConfig:
    def __new__(cls, *args, **kwargs):

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            files = {
                'voila_files': config_parser['FILES']['voila'].split('\n'),
                'voila_file': config_parser['FILES']['voila'].split('\n')[0],
                'splice_graph_file': config_parser['FILES']['splice_graph']
            }

            settings = dict(config_parser['SETTINGS'])

            for int_key in ['nproc', 'num_divisions', 'parallel_chunksize']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'copy_only', 'overwrite', 'memory_map_hdf5', 'preserve_handles_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _SplitterConfig(**{**files, **settings, **filters})

        return this_config


class RecombineConfig:
    def __new__(cls, *args, **kwargs):

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            settings = dict(config_parser['SETTINGS'])

            for int_key in ['nproc', 'parallel_chunksize']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'memory_map_hdf5', 'preserve_handles_hdf5']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _RecombineConfig(**{**settings, **filters})

        return this_config

class LongReadsConfig:
    def __new__(cls, *args, **kwargs):

        global this_config

        if this_config is None:
            voila_log().debug('Generating config object')
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)

            settings = dict(config_parser['SETTINGS'])

            for int_key in []:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['only_update_psi']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            if settings['only_update_psi']:
                if not settings['voila_file'] or not os.path.exists(settings['output_file']):
                    voila_log().critical("--only-update-psi requires --voila-file and --output-file to be specified")
                    voila_log().critical("--output-file should point to an existing .lr.voila file")
                    sys.exit(1)
            else:
                if any(not settings[x] for x in ('splice_graph_file', 'lr_gtf_file', 'lr_tsv_file')):
                    voila_log().critical("--splice-graph-file, --lr-gtf-file, -lr-tsv-file are required")
                    sys.exit(1)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _LongReadsConfig(**{**settings, **filters})

        return this_config

class GlobalConfig:
    def __new__(cls, *args, **kwargs):
        """
        Before the object is created, we'll parse the ini file, save the named tuple to a global variable, and use it
        as the sington object. This class is specifically for the HTML view.

        :param args: arguments
        :param kwargs: keyword arguments
        :return: named tuple config
        """

        global this_config

        if this_config is None:
            config_parser = configparser.ConfigParser()
            config_parser.read(constants.CONFIG_FILE)
            #raise Exception("GlobalConfig called before any config initialized")
            return TsvConfig()

        else:
            return this_config