import configparser
import inspect
import sqlite3
from collections import namedtuple
from pathlib import Path
import sys
from rna_voila import constants

from rna_voila.api import ViewPsi, SpliceGraph, find_analysis_type, get_mixed_analysis_type_str, view_matrix_zarr

from rna_voila.exceptions import FoundNoSpliceGraphFile, FoundMoreThanOneSpliceGraph, \
    MixedAnalysisTypeVoilaFiles, FoundMoreThanOneVoilaFile, AnalysisTypeNotFound
from rna_voila.voila_log import voila_log

import new_majiq as nm

# TODO break singleton cache and singleton config into two separate objects?

_ViewConfig = namedtuple('ViewConfig', ['voila_file', 'voila_files',
                                        'cov_file', 'cov_files', 'splice_graph_file', 'zarr_file', 'sgc_files',
                                        'analysis_type', 'nproc',
                                        'force_index', 'debug', 'silent', 'port', 'host', 'web_server', 'index_file',
                                        'num_web_workers', 'strict_indexing', 'skip_type_indexing', 'splice_graph_only',
                                        'enable_passcode', 'ignore_inconsistent_group_errors',
                                        'enable_het_comparison_chooser', 'is_multipsi_view',
                                        'sg_zarr', 'sgc_zarr', 'cov_zarr', 'lsvid2lsvidx', 'lsvtype_cache', 'module_cache'
                                    ])
_ViewConfig.__new__.__defaults__ = (None,) * len(_ViewConfig._fields)
_TsvConfig = namedtuple('TsvConfig', ['file_name', 'voila_files', 'voila_file',
                                      'cov_file', 'cov_files', 'splice_graph_file',
                                      'zarr_file', 'sgc_files',
                                      'non_changing_threshold', 'nproc', 'threshold', 'analysis_type', 'show_all',
                                      'debug', 'probability_threshold', 'silent', 'gene_ids', 'gene_names', 'lsv_ids',
                                      'lsv_types', 'strict_indexing', 'show_read_counts',
                                      'ignore_inconsistent_group_errors', 'non_changing_pvalue_threshold',
                                      'non_changing_within_group_iqr',
                                      'non_changing_between_group_dpsi', 'changing_pvalue_threshold',
                                      'changing_between_group_dpsi'])
_TsvConfig.__new__.__defaults__ = (None,) * len(_TsvConfig._fields)
_ClassifyConfig = namedtuple('ClassifyConfig', ['directory', 'voila_files', 'voila_file',
                                                'cov_file', 'cov_files', 'splice_graph_file', 'zarr_file', 'sgc_files',
                                      'nproc', 'decomplexify_psi_threshold', 'decomplexify_deltapsi_threshold',
                                      'decomplexify_reads_threshold', 'analysis_type', 'gene_ids',
                                      'debug', 'silent', 'keep_constitutive', 'keep_no_lsvs_modules', 'only_binary',
                                      'untrimmed_exons', 'putative_multi_gene_regions',
                                                'probability_changing_threshold',
                                                'probability_non_changing_threshold', 'show_all',
                                                'non_changing_pvalue_threshold', 'non_changing_within_group_iqr',
                                                'non_changing_between_group_dpsi', 'changing_pvalue_threshold',
                                                'changing_between_group_dpsi', 'changing_between_group_dpsi_secondary',
                                                'keep_no_lsvs_junctions', 'debug_num_genes', 'overwrite', 'output_mpe',
                                                'heatmap_selection', 'logger', 'enabled_outputs',
                                                'ignore_inconsistent_group_errors', 'disable_metadata',
                                                'sg_zarr', 'sgc_zarr', 'cov_zarr', 'lsvid2lsvidx', 'lsvtype_cache', 'module_cache'])
_ClassifyConfig.__new__.__defaults__ = (None,) * len(_ClassifyConfig._fields)
_FilterConfig = namedtuple('FilterConfig', ['directory', 'voila_files', 'voila_file',
                                            'cov_file', 'cov_files', 'splice_graph_file',
                                            'nproc', 'gene_ids', 'debug', 'silent', 'analysis_type', 'overwrite',
                                            'gene_ids_file', 'lsv_ids', 'lsv_ids_file', 'voila_files_only',
                                            'splice_graph_only',
                                            'changing_threshold', 'non_changing_threshold',
                                            'probability_changing_threshold',
                                            'probability_non_changing_threshold', 'changing', 'non_changing',
                                            'logger'])
_FilterConfig.__new__.__defaults__ = (None,) * len(_FilterConfig._fields)
_SplitterConfig = namedtuple('SplitterConfig', ['directory', 'voila_files', 'voila_file',
                                                'cov_file', 'cov_files', 'splice_graph_file',
                                      'nproc', 'debug', 'silent', 'num_divisions', 'copy_only', 'analysis_type',
                                                'overwrite', 'logger'])
_SplitterConfig.__new__.__defaults__ = (None,) * len(_SplitterConfig._fields)
_RecombineConfig = namedtuple('RecombineConfig', ['directories', 'directory', 'nproc', 'debug', 'silent', 'analysis_type', 'logger'])
_RecombineConfig.__new__.__defaults__ = (None,) * len(_RecombineConfig._fields)

# global config variable to act as the singleton instance of the config.
this_config = None


def find_splice_graph_file(_vs):
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
                # try:
                #     with SpliceGraph(v):
                #         sg_files.add(v)
                # except sqlite3.DatabaseError:
                #     pass

            elif v.is_dir():

                #print(v.parts)
                if v.parts[-1].endswith('.zarr'):
                    zarr_files.add(v)

                if v.parts[-1].endswith('.sgc'):
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

    for v in vs:
        v = Path(v)

        if v.is_file() and v.name.endswith('.voila'):
            voila_files.append(v)

        elif v.is_dir():
            x = find_voila_files(v.iterdir())
            voila_files = [*voila_files, *x]

    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.
    voila_files.sort()

    return voila_files


def find_cov_files(vs):
    """
    Find all voila files in files and directories.
    :param vs: list of files and directories.
    :return: list of voila files
    """

    cov_files = []

    for v in vs:
        v = Path(v)

        if v.is_dir() and v.name.endswith('.psicov'):
            cov_files.append(v)

        if v.is_dir() and v.name.endswith('.dpsicov'):
            cov_files.append(v)

        if v.is_dir() and v.name.endswith('.hetcov'):
            cov_files.append(v)

        elif v.is_dir():
            x = find_cov_files(v.iterdir())
            cov_files = [*cov_files, *x]

    # We rely on the directory of voila files to store the index for het runs, therefore it would be best to
    # have the same directory every time.
    cov_files.sort()

    return cov_files



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

    if not args.func.__name__ == 'recombine':
        sg_file = find_splice_graph_file(args.files)
    else:
        sg_file = None

    if (hasattr(args, 'splice_graph_only') and args.splice_graph_only) or args.func.__name__ == 'recombine':

        analysis_type = ''
        voila_files = []
        cov_files = []

    else:

        voila_files = find_voila_files(args.files)
        cov_files = find_cov_files(args.files)

        if voila_files and cov_files:
            raise Exception("Found both voila files and cov files in provided paths. Only one type should be provided.")



        if args.func.__name__ in ("Filter", "Classify", 'splitter', 'recombine'):
            analysis_type = get_mixed_analysis_type_str(voila_files, cov_files)
        else:
            analysis_type = find_analysis_type(voila_files, cov_files)

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

    # Get filters from arguments, add them to the appropriate section, and remove them from arguments.
    for lsv_filter in ['lsv_types', 'lsv_ids', 'gene_ids', 'gene_names']:
        if lsv_filter in attrs and attrs[lsv_filter]:
            try:
                config_parser.set(filters, lsv_filter, '\n'.join(attrs[lsv_filter]))
            except configparser.NoSectionError:
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


def _getInputFilesSet(config_parser, view=False):
    files = {}
    settings = dict(config_parser['SETTINGS'])


    if 'splice_graph' in config_parser['FILES']:
        files['splice_graph_file'] = config_parser['FILES']['splice_graph']
    else:
        files['zarr_file'] = config_parser['FILES']['zarr_file']
        files['sgc_files'] = config_parser['FILES']['sgc_files'].split('\n')
        files['sg_zarr'] = nm.SpliceGraph.from_zarr(files['zarr_file'])
        files['sgc_zarr'] = nm.SpliceGraphReads.from_zarr(files['sgc_files'])
        mask = nm.SpliceGraphMask.from_arrays(files['sg_zarr'].introns, files['sg_zarr'].junctions)
        files['module_cache'] = files['sg_zarr'].modules(mask)

    if 'voila' in config_parser['FILES']:
        files['voila_files'] = config_parser['FILES']['voila'].split('\n')
        files['voila_file'] = config_parser['FILES']['voila'].split('\n')[0]
        if view:
            settings['is_multipsi_view'] = len(files['voila_files']) > 1
    else:
        files['cov_files'] = config_parser['FILES']['majiq'].split('\n')
        files['cov_file'] = config_parser['FILES']['majiq'].split('\n')[0]
        if view:
            settings['is_multipsi_view'] = len(files['cov_files']) > 1

        if settings.get('splice_graph_only', 'True') != 'True':
            if settings['analysis_type'] == constants.ANALYSIS_PSI:
                files['cov_zarr'] = nm.PsiCoverage.from_zarr(files['cov_files'])
            elif settings['analysis_type'] == constants.ANALYSIS_DELTAPSI:
                files['cov_zarr'] = nm.DeltaPsiDataset.from_zarr(files['cov_file'])
            elif settings['analysis_type'] == constants.ANALYSIS_HETEROGEN:
                files['cov_zarr'] = nm.HeterogenDataset.from_zarr(files['cov_file'])

    if 'cov_files' in files and 'zarr_file' in files:
        if settings.get('splice_graph_only', 'True') != 'True':

            files['lsvid2lsvidx'] = view_matrix_zarr.get_lsvid2lsvidx(files['sg_zarr'], files['cov_zarr'])
            files['lsvtype_cache'] = view_matrix_zarr.get_lsvtype_cache(files['sg_zarr'], files['cov_zarr'])
            import rna_voila.index
            rna_voila.index.ZarrIndex.init_cache(
                dpsi=settings['analysis_type'] == constants.ANALYSIS_DELTAPSI,
                het=settings['analysis_type'] == constants.ANALYSIS_HETEROGEN
            )

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


            for int_key in ['nproc', 'port', 'num_web_workers']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for bool_key in ['force_index', 'silent', 'debug', 'strict_indexing', 'skip_type_indexing',
                             'ignore_inconsistent_group_errors', 'enable_het_comparison_chooser']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)


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

            files = {
                'voila_files': config_parser['FILES']['voila'].split('\n'),
                'voila_file': config_parser['FILES']['voila'].split('\n')[0],
            }

            if 'splice_graph' in config_parser['FILES']:
                files['splice_graph_file'] = config_parser['FILES']['splice_graph']
            else:
                files['zarr_file'] = config_parser['FILES']['zarr_file']
                files['sgc_files'] = config_parser['FILES']['sgc_files'].split('\n')

            settings = dict(config_parser['SETTINGS'])
            for int_key in ['nproc']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['non_changing_threshold', 'threshold', 'probability_threshold',
                              'non_changing_pvalue_threshold',
                              'non_changing_within_group_iqr', 'non_changing_between_group_dpsi',
                              'changing_pvalue_threshold', 'changing_between_group_dpsi']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['show_all', 'silent', 'debug', 'strict_indexing', 'show_read_counts',
                             'ignore_inconsistent_group_errors']:
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

            files, settings = _getInputFilesSet(config_parser)



            for int_key in ['nproc', 'keep_constitutive', 'decomplexify_reads_threshold', 'debug_num_genes']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['decomplexify_psi_threshold', 'decomplexify_deltapsi_threshold',
                              'probability_changing_threshold',
                              'probability_non_changing_threshold', 'non_changing_pvalue_threshold',
                              'non_changing_within_group_iqr', 'non_changing_between_group_dpsi',
                              'changing_pvalue_threshold', 'changing_between_group_dpsi',
                              'changing_between_group_dpsi_secondary']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'keep_no_lsvs_modules', 'only_binary', 'untrimmed_exons', 'overwrite',
                             'putative_multi_gene_regions', 'show_all', 'keep_no_lsvs_junctions', 'output_mpe',
                             'ignore_inconsistent_group_errors', 'disable_metadata'
                             ]:
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
            for int_key in ['nproc']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in ['non_changing_threshold', 'changing_threshold', 'probability_changing_threshold',
                              'probability_non_changing_threshold']:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'overwrite', 'voila_files_only', 'splice_graph_only', 'changing', 'non_changing']:
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

            for int_key in ['nproc', 'num_divisions']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug', 'copy_only', 'overwrite']:
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

            for int_key in ['nproc']:
                settings[int_key] = config_parser['SETTINGS'].getint(int_key)
            for float_key in []:
                settings[float_key] = config_parser['SETTINGS'].getfloat(float_key)
            for bool_key in ['debug']:
                settings[bool_key] = config_parser['SETTINGS'].getboolean(bool_key)

            filters = {}
            if config_parser.has_section('FILTERS'):
                for key, value in config_parser['FILTERS'].items():
                    filters[key] = config_parser['FILTERS'][key].split('\n')

            this_config = _RecombineConfig(**{**settings, **filters})

        return this_config