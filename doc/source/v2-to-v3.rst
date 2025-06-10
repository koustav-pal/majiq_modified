Majiq v2 to v3 migration guide
==========

In the major version bump from Majiq V2 to Majiq V3, there are significant differences in the workflow and command-line interface that long-term users may find unintuitive. We provide this page as a sort of "translation" to compare a v2 run to a similar v3 run.

v2
~~~~~
.. code-block:: bash


    # Illustrates common MAJIQ build, quantifier, and visualization commands for two groups of two samples each

    ## MAJIQ BUILD
    majiq build -c /path/to/inputs/config.ini /path/to/inputs/Homo_sapiens.GRCh38.94.gff3 -o /path/to/results/build

    ## OPTIONAL STEP: MOCCASIN adjustment for confounding factors
    python /path/to/moccasin.py /path/to/inputs/model_matrix.tsv /path/to/results/build /path/to/results/build_after_moccasin confounder_column_1 confounder_column_2

    ## MAJIQ PSI
    majiq psi -n Brain_Cerebellum /path/to/results/build/sample_bc_1.majiq /path/to/results/build/sample_bc_2.majiq -o /path/to/results/psi
    majiq psi -n Muscle_Skeletal /path/to/results/build/sample_ms_1.majiq /path/to/results/build/sample_ms_2.majiq -o /path/to/results/psi

    ## MAJIQ DELTAPSI
    majiq deltapsi -o /path/to/results/dpsi -n Brain_Cerebellum Muscle_Skeletal -grp1 /path/to/results/build/sample_bc_1.majiq /path/to/results/build/sample_bc_2.majiq -grp2 /path/to/results/build/sample_ms_1.majiq /path/to/results/build/sample_ms_2.majiq

    ## HETEROGEN
    majiq heterogen -o /path/to/results/het -n Brain_Cerebellum Muscle_Skeletal -grp1 /path/to/results/build/sample_bc_1.majiq /path/to/results/build/sample_bc_2.majiq -grp2 /path/to/results/build/sample_ms_1.majiq /path/to/results/build/sample_ms_2.majiq

    ## VOILA PSI TSV

    voila tsv -f /path/to/results/tsv/Brain_Cerebellum.psi.voila.tsv /path/to/results/build/splicegraph.sql /path/to/results/psi/Brain_Cerebellum.psi.voila
    voila tsv -f /path/to/results/tsv/Muscle_Skeletal.psi.voila.tsv /path/to/results/build/splicegraph.sql /path/to/results/psi/Muscle_Skeletal.psi.voila

    ## VOILA DPSI TSV
    voila tsv -f /path/to/results/tsv/Brain_Cerebellum_Muscle_Skeletal.dpsi.voila.tsv /path/to/results/build/splicegraph.sql /path/to/results/dpsi/Brain_Cerebellum-Muscle_Skeletal.deltapsi.voila

    ## VOILA HET TSV
    voila tsv -f /path/to/results/tsv/Brain_Cerebellum_Muscle_Skeletal.het.voila.tsv /path/to/results/build/splicegraph.sql /path/to/results/het/Brain_Cerebellum-Muscle_Skeletal.het.voila

    ## VOILA visualization
    voila view /path/to/results/build/splicegraph.sql /path/to/results/psi/Brain_Cerebellum.psi.voila

v3
~~~~~
.. code-block:: bash

    ## MAJIQ BUILD

    # convert gff3 text file to annotation databse
    majiq-v3 gff3 /path/to/inputs/Homo_sapiens.GRCh38.94.gff3 /path/to/results/annotations/sg.zarr

    # convert bam files into splice junction databases
    majiq-v3 sj /path/to/inputs/bam/sample_bc_1.bam /path/to/results/annotations/sg.zarr /path/to/results/sj/sample_bc_1.sj
    majiq-v3 sj /path/to/inputs/bam/sample_bc_2.bam /path/to/results/annotations/sg.zarr /path/to/results/sj/sample_bc_2.sj
    majiq-v3 sj /path/to/inputs/bam/sample_ms_1.bam /path/to/results/annotations/sg.zarr /path/to/results/sj/sample_ms_1.sj
    majiq-v3 sj /path/to/inputs/bam/sample_ms_2.bam /path/to/results/annotations/sg.zarr /path/to/results/sj/sample_ms_2.sj

    # the main build command (similar to v2 majiq build)
    # the config file is specified as a TSV with a 'group', 'prefix' (experiment name), and 'sj' (path to sj file) columns
    majiq-v3 build /path/to/results/annotations/sg.zarr /path/to/results/build/sg.zarr --groups-tsv config.tsv

    # optionally, you can skip making the groups TSV and use a version with inline sjs and no config file (each specification of --sjs is a different group)
    # majiq-v3 build /path/to/results/annotations/sg.zarr /path/to/results/build/sg.zarr --sjs /path/to/results/sj/sample_bc_1.sj /path/to/results/sj/sample_bc_2.sj --sjs /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj

    ## CALCULATE PSI COVERAGE FOR USE BY QUANTIFIERS
    majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Brain_Cerebellum.psicov /path/to/results/sj/sample_bc_1.sj /path/to/results/sj/sample_bc_2.sj
    majiq-v3 psi-coverage /path/to/results/build/sg.zarr /path/to/results/psi/Muscle_Skeletal.psicov /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj

    ## OPTIONAL STEP: MOCCASIN adjustment for confounding factors (similar to the MOCCASIN project, which was separate in v2)
    majiq-v3 moccasin-pipeline /path/to/results/build_after_moccasin /path/to/results/psi/Brain_Cerebellum.psicov /path/to/results/psi/Muscle_Skeletal.psicov --factors_tsv /path/to/inputs/model_matrix.tsv --confounding confounder_column_1 confounder_column_2 --overwrite

    ## MAJIQ PSI and output TSV (similar to TSV output by majiq psi in v2)
    majiq-v3 quantify /path/to/results/psi/Brain_Cerebellum.psicov --min-experiments 0.01 --splicegraph /path/to/results/build/sg.zarr --output-tsv /path/to/results/psi/Brain_Cerebellum.tsv --overwrite
    majiq-v3 quantify /path/to/results/psi/Muscle_Skeletal.psicov --min-experiments 0.01 --splicegraph /path/to/results/build/sg.zarr --output-tsv /path/to/results/psi/Muscle_Skeletal.tsv --overwrite

    ## MAJIQ DPSI and output TSV (similar to TSV output by majiq deltapsi in v2)
    majiq-v3 deltapsi --splicegraph /path/to/results/build/sg.zarr --output-voila /path/to/results/dpsi/Brain_Cerebellum-Muscle_Skeletal.dpsicov --output-tsv /path/to/results/dpsi/Brain_Cerebellum-Muscle_Skeletal.tsv -psi1 /path/to/results/psi/Brain_Cerebellum.psicov -psi2 /path/to/results/psi/Muscle_Skeletal.psicov

    ## MAJIQ HET and output TSV (similar to TSV output by majiq heterogen in v2)
    majiq-v3 heterogen --stats infoscore mannwhitneyu ttest tnom --splicegraph /path/to/results/build/sg.zarr --output-voila /path/to/results/het/Brain_Cerebellum-Muscle_Skeletal.hetcov --output-tsv /path/to/results/het/Brain_Cerebellum-Muscle_Skeletal.tsv -psi1 /path/to/results/psi/Brain_Cerebellum.psicov -psi2 /path/to/results/psi/Muscle_Skeletal.psicov

    ## VOILA visualization and subcommands
    # further downstream voila usages are very similar to v2, except that you must provide paths to the splicegraph (sg.zarr), quant files, AND splicegraph coverage SGC files
    # the quant files, similar to v2, are generated by the quantify, deltapsi, or heterogen commands above in majiq-v3
    # the sgc files are a separate step which must be provided along with each each group of sj files like below
    ## OUTPUT SPLICEGRAPH COVERAGE FOR USE BY VOILA
    majiq-v3 sg-coverage /path/to/results/build/sg.zarr /path/to/results/build/Brain_Cerebellum.sgc /path/to/results/sj/sample_bc_1.sj /path/to/results/sj/sample_bc_2.sj
    majiq-v3 sg-coverage /path/to/results/build/sg.zarr /path/to/results/build/Muscle_Skeletal.sgc /path/to/results/sj/sample_ms_1.sj /path/to/results/sj/sample_ms_2.sj

    # voila view example using single PSI ; other voila commands, like tsv, modulizer work as in v2
    # note that voila commands can use v2 splicegraph and v2 hdf5 quantification files for backwards compatibility, but not BOTH v2 and v3 inputs at the same time
    voila view /path/to/results/build/sg.zarr /path/to/results/psi/Brain_Cerebellum.psicov /path/to/results/build/Brain_Cerebellum.sgc

Summary
~~~~~

Note that the commands for majiq v3 are separated out into smaller units, and there are more required to run. This is a design decision to aid with usage of the software for clinical purposes on larger compute environments with job schedulers, or for easier parallel execution across multiple machines or scripts. Quantifier, sj, tsv, gff3, and sg-coverage are some examples of commands that can be run in parallel to each other, allowing smaller memory footprints and more efficient compute usage. We would recommend putting your majiq run pipeline into a pipeline software such as snakemake, cwl, etc for reproducability.
