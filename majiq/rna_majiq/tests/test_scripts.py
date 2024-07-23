"""
test_scripts.py

Test of console scripts (mostly smoke tests, but also checking against
saved/expected output in select cases)

Depends on pytest_console_scripts for script_runner fixture

Author: Joseph K Aicher
"""

import json
from pathlib import Path

import rna_majiq as nm
import pytest
import xarray as xr
from rna_majiq.experiments import bam_experiment_name

ANNOTATED_GFF3 = "data/gff3/SRSF4.gff3"
ANNOTATED_SG = "data/sg/annotated.sg.zip"
EXPERIMENT_NAMES = sorted(
    bam_experiment_name(x) for x in Path(__file__).parent.glob("data/bam/*.bam")
)
EXPERIMENT_GROUPS = [
    ("DUM", [x for x in EXPERIMENT_NAMES if "DUM" in x]),
    ("TIZ", [x for x in EXPERIMENT_NAMES if "TIZ" in x]),
    ("JBU", [x for x in EXPERIMENT_NAMES if "JBU" in x]),
    ("AOV", [x for x in EXPERIMENT_NAMES if "AOV" in x]),
]
COMBINED_SG = "data/sg/combined.sg.zip"
FACTORS_MODEL = "data/factors/factor_model.zarr.zip"
FACTORS_TSV = "data/factors/factors.tsv"
FACTORS_TSV_CONFOUNDERS = ["confounding"]
COVERAGE_MODEL = "data/factors/coverage_model.zarr.zip"


def get_path(p) -> str:
    """Path to p relative to test script"""
    return str(Path(__file__).parent / p)


def get_bam_path(name):
    """Path to BAM files"""
    return get_path(f"data/bam/{name}.bam")


def get_sj_path(name, strandness):
    """Path to SJ files as input/expected from BAM"""
    return get_path(f"data/sj/{name}.{strandness}.sj.zip")


def get_build_path(group, simplify, min_experiments):
    """Path to splicegraph files as input/expected from build command"""
    return get_path(f"data/build/{group}.{simplify}.{min_experiments}.sg.zip")


def get_psicov_path(group):
    """Path to psicoverage files as input/expected from psi-coverage command"""
    return get_path(f"data/psicov/{group}.psicov.zip")


def assert_splicegraphs(sg1: nm.SpliceGraph, sg2: nm.SpliceGraph):
    xr.testing.assert_equal(sg1.contigs.df, sg2.contigs.df)
    xr.testing.assert_equal(sg1.genes.df, sg2.genes.df)
    xr.testing.assert_equal(sg1.exons.df, sg2.exons.df)
    xr.testing.assert_equal(sg1.introns.df, sg2.introns.df)
    xr.testing.assert_equal(sg1.junctions.df, sg2.junctions.df)


def test_gff3_command(script_runner, tmp_path):
    """Test new-majiq gff3 runs and compare to expected result"""
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        ["new-majiq", "gff3", get_path(ANNOTATED_GFF3), path_result]
    )
    assert ret.success
    path_expected = get_path(ANNOTATED_SG)
    sg_result = nm.SpliceGraph.from_zarr(path_result)
    sg_expected = nm.SpliceGraph.from_zarr(path_expected)
    assert_splicegraphs(sg_result, sg_expected)
    return


@pytest.mark.parametrize("name", EXPERIMENT_NAMES)
@pytest.mark.parametrize("strandness", ["AUTO", "NONE", "FORWARD", "REVERSE"])
def test_sj_command(script_runner, name, strandness, tmp_path):
    """Test new-majiq sj command. Check results if found in data (AUTO)"""
    path_bam = get_bam_path(name)
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        ["new-majiq", "sj", path_bam, get_path(ANNOTATED_SG), path_result]
    )
    assert ret.success
    path_expected = get_sj_path(name, strandness=strandness)
    if Path(path_expected).exists():
        sj_result = nm.SJExperiment.from_zarr(path_result)
        sj_expected = nm.SJExperiment.from_zarr(path_expected)
        xr.testing.assert_equal(sj_result.junctions._df, sj_expected.junctions._df)
        xr.testing.assert_equal(sj_result.introns._df, sj_expected.introns._df)
        xr.testing.assert_equal(
            sj_result.junctions.regions.df, sj_expected.junctions.regions.df
        )
        xr.testing.assert_equal(
            sj_result.introns.regions.df, sj_expected.introns.regions.df
        )
    return


@pytest.mark.parametrize("build_group", EXPERIMENT_GROUPS)
@pytest.mark.parametrize("simplify", [False, True])
@pytest.mark.parametrize("min_experiments", [0.5, 1])
def test_build_command(script_runner, build_group, simplify, min_experiments, tmp_path):
    """Test new-majiq build command

    Test new-majiq build command. Check results if found in data
    (when simplify=True, min-experiments=1)
    """
    group, names = build_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "build",
            get_path(ANNOTATED_SG),
            path_result,
            "--sjs",
            *paths_sj,
            "--simplify" if simplify else "--no-simplify",
            "--min-experiments",
            str(min_experiments),
        ]
    )
    assert ret.success
    path_expected = get_build_path(
        group, simplify=simplify, min_experiments=min_experiments
    )
    if Path(path_expected).exists():
        sg_result = nm.SpliceGraph.from_zarr(path_result)
        sg_expected = nm.SpliceGraph.from_zarr(path_expected)
        assert_splicegraphs(sg_result, sg_expected)
    return


@pytest.mark.parametrize("min_experiments", ["1", "999"])
def test_update_vs_combine(script_runner, tmp_path, min_experiments):
    """Test majiq-build update with groups tsv vs update on sjs and combine

    We expect that running update on build groups individually and then
    combining should be equivalent to a build with them all together
    """
    grp1 = [get_sj_path(x, strandness="AUTO") for x in EXPERIMENT_GROUPS[0][1]]
    grp2 = [get_sj_path(x, strandness="AUTO") for x in EXPERIMENT_GROUPS[1][1]]
    # separate calls to majiq-build update, then majiq-build combine
    assert script_runner.run(
        [
            "majiq-build",
            "update",
            get_path(ANNOTATED_SG),
            str(tmp_path / "grp1.sg"),
            *("--sjs", *grp1),
            "--no-simplify",
            *("--min-experiments", min_experiments),
        ]
    ).success
    assert script_runner.run(
        [
            "majiq-build",
            "update",
            get_path(ANNOTATED_SG),
            str(tmp_path / "grp2.sg"),
            *("--sjs", *grp2),
            "--no-simplify",
            *("--min-experiments", min_experiments),
        ]
    ).success
    assert script_runner.run(
        [
            "majiq-build",
            "combine",
            str(tmp_path / "combined.sg"),
            *("--keep-denovo", str(tmp_path / "grp1.sg"), str(tmp_path / "grp2.sg")),
        ]
    )
    # single call to majiq-build
    with open(tmp_path / "groups.tsv", mode="w") as handle:
        print("group\tsj", file=handle)
        for x in grp1:
            print(f"grp1\t{x}", file=handle)
        for x in grp2:
            print(f"grp2\t{x}", file=handle)
    assert script_runner.run(
        [
            "majiq-build",
            "update",
            get_path(ANNOTATED_SG),
            str(tmp_path / "multigrp_build.sg"),
            *("--groups-tsv", str(tmp_path / "groups.tsv")),
            "--no-simplify",
            *("--min-experiments", min_experiments),
        ]
    ).success
    # check equivalence of splicegraphs
    assert_splicegraphs(
        nm.SpliceGraph.from_zarr(tmp_path / "combined.sg"),
        nm.SpliceGraph.from_zarr(tmp_path / "multigrp_build.sg"),
    )
    return


def test_combine_command(script_runner, tmp_path):
    """Test new-majiq combine command

    Test new-majiq combine command. Set first experiment group as annotated.
    Compare result to expected.
    """
    make_annotated = [
        get_build_path(group, simplify=True, min_experiments=1)
        for group, _ in EXPERIMENT_GROUPS[:1]
    ]
    keep_denovo = [
        get_build_path(group, simplify=True, min_experiments=1)
        for group, _ in EXPERIMENT_GROUPS[1:]
    ]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "combine",
            path_result,
            "--make-annotated",
            *make_annotated,
            "--keep-denovo",
            *keep_denovo,
        ]
    )
    assert ret.success
    path_expected = get_path(COMBINED_SG)
    sg_result = nm.SpliceGraph.from_zarr(path_result)
    sg_expected = nm.SpliceGraph.from_zarr(path_expected)
    assert_splicegraphs(sg_result, sg_expected)
    return


@pytest.mark.parametrize("batch_group", EXPERIMENT_GROUPS)
@pytest.mark.parametrize("ignore", [None, "first", "all"])
def test_psi_coverage_command(script_runner, batch_group, ignore, tmp_path):
    """Test psi-coverage command

    Test psi-coverage command, comparing to results when not ignoring first
    experiment group events
    """
    group, names = batch_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ignore_from = list()
    if ignore:
        ignore_from = [
            "--ignore-from",
            (
                get_path(COMBINED_SG)
                if ignore == "all"
                else get_build_path(
                    EXPERIMENT_GROUPS[0][0], simplify=True, min_experiments=1
                )
            ),
        ]
    ret = script_runner.run(
        [
            "new-majiq",
            "psi-coverage",
            get_path(COMBINED_SG),
            path_result,
            *paths_sj,
            *ignore_from,
        ]
    )
    assert ret.success
    if not ignore:
        path_expected = get_psicov_path(group)
        result = nm.PsiCoverage.from_zarr(path_result)
        expected = nm.PsiCoverage.from_zarr(path_expected)
        xr.testing.assert_equal(
            result.df.drop_vars(["bootstrap_psi", "bootstrap_total", "prefix_total"]),
            expected.df.drop_vars(["bootstrap_psi", "bootstrap_total", "prefix_total"]),
        )
    return


def test_psicoverage_no_duplicates():
    psicov_paths = [get_psicov_path(group) for group, names in EXPERIMENT_GROUPS]
    psi = nm.PsiCoverage.from_zarr(psicov_paths)
    psi2 = nm.PsiCoverage.from_zarr(psicov_paths * 2)
    xr.testing.assert_equal(psi.df, psi2.df)
    assert len(set(psi2.prefixes)) == psi2.num_prefixes
    return


@pytest.mark.parametrize("groups", [EXPERIMENT_GROUPS[:1], EXPERIMENT_GROUPS[:3]])
def test_psi_group_command_allpsicov(script_runner, groups, tmp_path):
    """Test psi-group command

    Test psi-group command, comparing to results when not ignoring first
    experiment group events
    """
    psicov_paths = [get_psicov_path(group) for group, names in groups]
    path_result = str(tmp_path / "all")
    ret = script_runner.run(
        [
            "new-majiq",
            "psi-group",
            path_result,
            *psicov_paths,
            "--disable-progress",
        ]
    )
    assert ret.success
    psigroup = nm.PsiGroup.from_zarr(path_result)
    psicov = nm.PsiCoverage.from_zarr(psicov_paths)
    xr.testing.assert_equal(psigroup.raw_alpha, psicov.raw_alpha)
    xr.testing.assert_equal(psigroup.raw_beta, psicov.raw_beta)
    xr.testing.assert_equal(psigroup.approximate_alpha, psicov.approximate_alpha)
    xr.testing.assert_equal(psigroup.approximate_beta, psicov.approximate_beta)
    return


@pytest.mark.filterwarnings(
    "ignore:.*finalize object.*dead.*:pytest.PytestUnraisableExceptionWarning"
)
@pytest.mark.parametrize("groups", [EXPERIMENT_GROUPS[:]])
def test_psi_group_command_mixedpsi(script_runner, groups, tmp_path):
    """Test psi-group command

    Test psi-group command on mixed input by running on half and then with the
    next half, and then vs running on all
    """
    if len(groups) < 2:
        pytest.skip(reason="cannot test mixedpsi with only one file")
    psicov_paths = [get_psicov_path(group) for group, names in groups]
    path_all = str(tmp_path / "all")
    ret = script_runner.run(
        [
            "new-majiq",
            "psi-group",
            path_all,
            *psicov_paths,
            "--disable-progress",
        ]
    )
    assert ret.success
    path1 = str(tmp_path / "part1")
    ret = script_runner.run(
        [
            "new-majiq",
            "psi-group",
            path1,
            *psicov_paths[: len(psicov_paths) // 2],
            "--disable-progress",
        ]
    )
    assert ret.success
    path_mixed = str(tmp_path / "mixed")
    ret = script_runner.run(
        [
            "new-majiq",
            "psi-group",
            path_mixed,
            # last half of PsiCoverage files
            *psicov_paths[len(psicov_paths) // 2 :],
            # PsiGroup file
            path1,
            "--disable-progress",
        ]
    )
    assert ret.success
    group_all = nm.PsiGroup.from_zarr(path_all)
    group_mixed = nm.PsiGroup.from_zarr(path_mixed)
    xr.testing.assert_equal(
        group_all.df.sortby("prefix"), group_mixed.df.sortby("prefix")
    )
    return


@pytest.mark.parametrize("batch_group", EXPERIMENT_GROUPS)
def test_sg_coverage_command(script_runner, batch_group, tmp_path):
    """Test sg-coverage command

    Smoke test of sg-coverage command
    """
    group, names = batch_group
    paths_sj = [get_sj_path(name, strandness="AUTO") for name in names]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "sg-coverage",
            get_path(COMBINED_SG),
            path_result,
            *paths_sj,
        ]
    )
    assert ret.success
    return


@pytest.mark.parametrize("min_experiments", [None, 2])
def test_quantify_command(script_runner, min_experiments, tmp_path):
    """Smoke test for quantify command"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    flag = []
    if min_experiments:
        flag = ["--min-experiments", str(min_experiments)]
    ret = script_runner.run(
        [
            "new-majiq",
            "quantify",
            *paths_psicov,
            *flag,
            *("--splicegraph", get_path(COMBINED_SG)),
            *("--quantiles", "0.1", "0.9"),
            "--disable-progress",
        ]
    )
    assert ret.success
    return


@pytest.mark.parametrize("min_experiments", [None, 2])
def test_quantify_empty_command(script_runner, tmp_path, min_experiments):
    """Test if quantify command can run on empty psicoverage files"""
    path_psicov = tmp_path / "mock.psicov"
    nm.PsiCoverage.mock_with_psi_and_total([], []).to_zarr(path_psicov)
    flag = []
    if min_experiments:
        flag = ["--min-experiments", str(min_experiments)]
    ret = script_runner.run(
        [
            "new-majiq",
            "quantify",
            path_psicov,
            *flag,
            *("--quantiles", "0.1", "0.9"),
            "--disable-progress",
        ]
    )
    assert ret.success
    return


@pytest.mark.parametrize("behavior", ["normal", "split", "downsample"])
def test_deltapsi_command(script_runner, tmp_path, behavior):
    """Smoke test for deltapsi command"""
    # split up psicoverage files into two groups arbitrarily
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    if behavior != "split":
        split = len(paths_psicov) // 3
        if split == 0:
            pytest.skip(reason="split of Psicoverage files makes empty group")
        grp1 = ("-grp1", *paths_psicov[:split])
        grp2 = ("-grp2", *paths_psicov[split:])
    else:
        grp1 = ("-grp1", *paths_psicov)
        grp2 = tuple()
    downsample_arg = ("--downsample2",) if behavior == "downsample" else tuple()
    path_voila = str(tmp_path / "result.voila")
    path_tsv = str(tmp_path / "result.tsv")
    ret = script_runner.run(
        [
            "new-majiq",
            "deltapsi",
            *grp1,
            *grp2,
            *("-n", "grp1", "grp2"),
            *("--min-experiments", "2"),
            *("--splicegraph", get_path(COMBINED_SG)),
            *("--output-voila", path_voila),
            *("--output-tsv", path_tsv),
            *downsample_arg,
            "--disable-progress",
        ]
    )
    assert ret.success
    with open(path_tsv, mode="r") as handle:
        metadata_list = list()
        for x in handle:
            if x.startswith("# "):
                metadata_list.append(x[2:])
            else:
                break
    metadata = json.loads("".join(metadata_list))
    group_sizes = metadata["group_sizes"]
    if behavior == "split":
        assert all(v == 3 for v in group_sizes.values())
    elif behavior == "normal":
        assert group_sizes["grp1"] == 2 and group_sizes["grp2"] == 4
    elif behavior == "downsample":
        assert all(v == 2 for v in group_sizes.values())
    return


def test_deltapsi_empty_command(script_runner, tmp_path):
    """Test if quantify command can run on empty psicoverage files"""
    path_voila = str(tmp_path / "result.voila")
    path_tsv = str(tmp_path / "result.tsv")
    path_psicov1 = tmp_path / "mock1.psicov"
    path_psicov2 = tmp_path / "mock2.psicov"
    nm.PsiCoverage.mock_with_psi_and_total([], [], "mock1").to_zarr(path_psicov1)
    nm.PsiCoverage.mock_with_psi_and_total([], [], "mock2").to_zarr(path_psicov2)
    ret = script_runner.run(
        [
            "new-majiq",
            "deltapsi",
            *("-grp1", path_psicov1),
            *("-grp2", path_psicov2),
            *("-n", "grp1", "grp2"),
            *("--output-voila", path_voila),
            *("--output-tsv", path_tsv),
            "--disable-progress",
        ]
    )
    assert ret.success
    return


@pytest.mark.filterwarnings(
    "ignore:.*finalize object.*dead.*:pytest.PytestUnraisableExceptionWarning"
)
@pytest.mark.parametrize("behavior", ["normal", "split", "downsample"])
@pytest.mark.parametrize("use_psigroup", [False, True])
def test_heterogen_command(script_runner, tmp_path, behavior, use_psigroup):
    """Smoke test for heterogen command"""
    # split up psicoverage files into two groups arbitrarily
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    if behavior != "split":
        split = len(paths_psicov) // 3
        if split == 0:
            pytest.skip(reason="split of Psicoverage files makes empty group")
        grp1 = ("-grp1", *paths_psicov[:split])
        grp2 = ("-grp2", *paths_psicov[split:])
    else:
        grp1 = ("-grp1", *paths_psicov)
        grp2 = tuple()
    downsample_arg = ("--downsample2",) if behavior == "downsample" else tuple()
    path_voila = str(tmp_path / "result.voila")
    if use_psigroup:
        psigroup = str(tmp_path / "psigroup")
        ret = script_runner.run(
            [
                "new-majiq",
                "psi-group",
                psigroup,
                *grp1[1:],
                "--disable-progress",
            ]
        )
        assert ret.success
        grp1 = ("-grp1", psigroup)
    ret = script_runner.run(
        [
            "new-majiq",
            "heterogen",
            *grp1,
            *grp2,
            *("-n", "grp1", "grp2"),
            *("--min-experiments", "2"),
            *("--splicegraph", get_path(COMBINED_SG)),
            *("--output-voila", path_voila),
            *downsample_arg,
            "--disable-progress",
        ]
    )
    assert ret.success
    # check that group sizes were correct
    ds = nm.HeterogenDataset.from_zarr(path_voila)
    group_sizes = {grp_name: grp.num_prefixes for grp_name, grp in ds.groups.items()}
    if behavior == "split":
        assert all(v == 3 for v in group_sizes.values())
    elif behavior == "normal":
        assert group_sizes["grp1"] == 2 and group_sizes["grp2"] == 4
    elif behavior == "downsample":
        assert all(v == 2 for v in group_sizes.values())
    return


@pytest.mark.filterwarnings(
    "ignore:.*specified chunks separate the stored chunks.*:UserWarning",
    "ignore:.*finalize object.*dead.*:pytest.PytestUnraisableExceptionWarning",
)
def test_heterogen_empty_command(script_runner, tmp_path):
    """Test if quantify command can run on empty psicoverage files"""
    path_voila = str(tmp_path / "result.voila")
    path_tsv = str(tmp_path / "result.tsv")
    path_psicov1 = tmp_path / "mock1.psicov"
    path_psicov2 = tmp_path / "mock2.psicov"
    nm.PsiCoverage.concat(
        *(
            nm.PsiCoverage.mock_with_psi_and_total([], [], f"mock1-{i}")
            for i in range(10)
        )
    ).to_zarr(path_psicov1)
    nm.PsiCoverage.concat(
        *(
            nm.PsiCoverage.mock_with_psi_and_total([], [], f"mock2-{i}")
            for i in range(10)
        )
    ).to_zarr(path_psicov2)
    ret = script_runner.run(
        [
            "new-majiq",
            "heterogen",
            *("-grp1", path_psicov1),
            *("-grp2", path_psicov2),
            *("-n", "grp1", "grp2"),
            *("--output-voila", path_voila),
            *("--output-tsv", path_tsv),
            "--disable-progress",
        ]
    )
    assert ret.success
    return


@pytest.mark.filterwarnings(
    "ignore:.*finalize object.*dead.*:pytest.PytestUnraisableExceptionWarning"
)
@pytest.mark.parametrize("use_psigroup", [False, True])
def test_psi_controls_command(script_runner, tmp_path, use_psigroup: bool):
    """Smoke test for psi-controls command (all but first group)"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS[1:]]
    use_paths = [str(tmp_path / "psigroup")] if use_psigroup else paths_psicov
    if use_psigroup:
        ret = script_runner.run(
            [
                "new-majiq",
                "psi-group",
                use_paths[0],
                *paths_psicov,
                "--disable-progress",
            ]
        )
        assert ret.success
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "psi-controls",
            path_result,
            *use_paths,
            "--disable-progress",
        ]
    )
    assert ret.success
    return


def test_moccasin_factors_model_command(script_runner, tmp_path):
    """Smoke test for moccasin-factors-model command, check close to expected"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "moccasin-factors-model",
            path_result,
            *paths_psicov,
            "--intercept-only",
            *("--ruv-max-new-factors", "1"),
            "--disable-progress",
        ]
    )
    assert ret.success
    xr.testing.assert_allclose(
        xr.open_zarr(path_result), xr.open_zarr(get_path(FACTORS_MODEL))
    )
    return


@pytest.mark.parametrize("num_new_factors", [1, 3])
def test_moccasin_factors_infer_command(script_runner, tmp_path, num_new_factors):
    """Smoke test for moccasin-factors-infer command"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "moccasin-factors-infer",
            get_path(FACTORS_MODEL),
            path_result,
            *paths_psicov,
            "--intercept-only",
            *("--num-unknown-factors", str(num_new_factors)),
            "--disable-progress",
        ]
    )
    assert ret.success
    factors = nm.moccasin.Factors.from_zarr(path_result)
    assert factors.num_factors == 1 + num_new_factors
    return


@pytest.mark.parametrize("prefix_nchunks", [1, 3])
def test_moccasin_coverage_model_command(script_runner, tmp_path, prefix_nchunks):
    """Smoke test for moccasin-coverage-model command, check close to expected"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "moccasin-coverage-model",
            path_result,
            *paths_psicov,
            *("--factors-tsv", get_path(FACTORS_TSV)),
            *("--confounding", *FACTORS_TSV_CONFOUNDERS),
            *("--prefix-nchunks", str(prefix_nchunks)),
            "--disable-progress",
        ]
    )
    assert ret.success
    df_out = xr.open_zarr(path_result).drop_vars(["bootstrap_model"])
    df_base = xr.open_zarr(get_path(COVERAGE_MODEL)).drop_vars(
        ["bootstrap_model", "lsv_idx", "event_size"]
    )
    # lowering tolerance for differences: single precision and updated algorithm
    xr.testing.assert_allclose(df_out, df_base, rtol=1e-4, atol=1e-7)
    return


def test_moccasin_coverage_infer_command(script_runner, tmp_path):
    """Smoke test for moccasin-coverage-infer command"""
    paths_psicov = [get_psicov_path(group) for group, _ in EXPERIMENT_GROUPS]
    path_result = str(tmp_path / "result")
    ret = script_runner.run(
        [
            "new-majiq",
            "moccasin-coverage-infer",
            get_path(COVERAGE_MODEL),
            path_result,
            *paths_psicov,
            *("--factors-tsv", get_path(FACTORS_TSV)),
            *("--confounding", *FACTORS_TSV_CONFOUNDERS),
            "--disable-progress",
        ]
    )
    assert ret.success
    return


def test_psicov_concat():
    psicov_paths = [get_psicov_path(group) for group, names in EXPERIMENT_GROUPS]
    opened_together = nm.PsiCoverage.from_zarr(psicov_paths)
    opened_separate = nm.PsiCoverage.concat(
        *(nm.PsiCoverage.from_zarr(p) for p in psicov_paths)
    )
    xr.testing.assert_equal(opened_together.df, opened_separate.df)
    return


def test_sj_rename(tmp_path):
    sj = nm.SJExperiment.from_zarr(get_sj_path(EXPERIMENT_NAMES[0], strandness="AUTO"))
    sj_renamed = sj.rename_prefix(f"renamed-{sj.prefix}")
    assert sj_renamed.prefix == f"renamed-{sj.prefix}"
    roundtrip_path = str(tmp_path / "result")
    sj_renamed.to_zarr(roundtrip_path)
    sj_roundtrip = nm.SJExperiment.from_zarr(roundtrip_path)
    assert sj_renamed.prefix == sj_roundtrip.prefix
    xr.testing.assert_equal(sj_renamed.junctions._df, sj_roundtrip.junctions._df)
    xr.testing.assert_equal(sj_renamed.introns._df, sj_roundtrip.introns._df)
    xr.testing.assert_equal(
        sj_renamed.junctions.regions.df, sj_roundtrip.junctions.regions.df
    )
    xr.testing.assert_equal(
        sj_renamed.introns.regions.df, sj_roundtrip.introns.regions.df
    )
    return
