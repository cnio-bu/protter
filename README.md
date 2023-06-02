# protter

Pipeline to generate peptide-spectrum matches using different
search engines against target and decoy databases. Input may
be local or remote (see [Data sources](#data-sources)).

Allows for the easy addition of datasets and search engines.

**Index**
1. [Usage](#1-usage)
    * 1.1. [Create conda enviroment](#11-create-conda-enviroment)
    * 1.2. [Removing readthrough transcripts, PAR_Y genes and duplicates](#12-removing-readthrough-transcripts-par_y-genes-and-duplicates)
    * 1.3. [Preparing a sample sheet](#13-preparing-a-sample-sheet)
        - 1.3.1. [Splitting a sample sheet](#131-splitting-a-sample-sheet)
    * 1.4. [Running the workflow](#14-running-the-workflow)
    * 1.5. [Gathering output PSM files](#15-gathering-output-psm-files)
    * 1.6. [Processing PSM files for APPRIS's Proteo](#16-processing-psm-files-for-appriss-proteo)
2. [Troubleshooting](#2-troubleshooting)
    * 2.1. [download_sample](#21-download_sample)
    * 2.2. [raw_to_mzml](#22-raw_to_mzml)
    * 2.3. [comet](#23-comet)
    * 2.4. [percolator](#24-percolator)
3. [Data sources](#3-data-sources)
4. [Dependencies](#4-dependencies)


## 1. Usage

To run this pipeline, some previous steps can be performed to prepare the input files; and some subsequent steps to process the output files.

### 1.1. Create conda enviroment

Although not required, it is recommended to use the provided `conda` environment to launch the protter pipeline, in addition to the external pre and post-processing scripts. This environment contains snakemake, as well as some of the necessary dependencies that are not contained in the pipeline's own `conda` microenvironments. (See [Dependencies](#dependencies))

You need to have a version of `conda` installed. We recommend [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge), which has Mamba installed in the base environment. To create it, just run inside the base folder of the protter:
```shell
mamba env create -f res/conf/protter-env.yaml
```

And once all the dependencies have been installed, activate the environment with:
```shell
mamba activate protter
```

### 1.2. Removing readthrough transcripts, PAR_Y genes and duplicates

Sometimes it could be useful and interesting to delete the sequences tagged as 
readthrough transcripts and the chromosome Y pseudoautosomal regions (PAR) from
the translations file. Also it is useful to use a non-redundant FASTA file.


It is possible to get this file through the following script; which needs the
transcript file and the GTF of the desired gencode version. It can be run from the
workflow directory as follows:
```shell
python scripts/clean_pc_translations.py path/to/gencode.v{version}.pc_translations.fa.gz path/to/gencode.v{version}.annotation.gtf
```

A transcript file will be obtained in the same folder as the original file and with the
extension `.NoRT.u.fa.gz`. This indicates that the readthroughs and the PAR_Y (NoRT)
have been eliminated and that the sequences that appear are unique (u). This output,
file 'gencode.v{version}.pc_translations.NoRT.u.fa.gz', is the one that should be
indicated in the protter config.yaml file.

### 1.3. Preparing a sample sheet

A TSV sample sheet is necessary to run the workflow. At a minimum, this should
contain `dataset` and `sample` columns to identify each input sample, and a
`file` column to indicate the location of the sample input file, where the file
location may be given as a URL or as a local system path. A `checksum` column is
also needed for remote input files, to verify that they have been downloaded
correctly.

A workflow run may contain several datasets, each with hundreds or thousands of
input files, so it would be tedious and error-prone to generate the sample sheet
manually. When the workflow config YAML file has been prepared, there is a
script to prepare sample sheets automatically from the config file, which can be
run from the workflow directory as follows:
```shell
python scripts/prep_sample_sheet.py config.yaml samples.tsv
```
For each dataset marked as `enabled` in the config file, this obtains the sample
metadata from PRIDE or from local input files, depending on the `source`
configured for the given dataset. Having obtained the sample metadata of all
enabled datasets, it creates an initial sample sheet. For known datasets, the
sample sheet is further enhanced, filtering irrelevant samples, modifying
existing metadata, and adding new metadata columns.

A dataset is treated as known if there is a module of the same name in the
`scripts/ds` directory, and if that module contains a function
`enhance_sample_metadata` that accepts the initial sample sheet as a pandas
`DataFrame` and returns an enhanced sample sheet `DataFrame`. (For example,
see `scripts/ds/kim_2014.py`.) Dataset modules can be added as needed for any
known dataset, and will be recognised automatically by `protter` provided that
it has the same name as the dataset and contains a working
`enhance_sample_metadata` function.

Whether the datasets in a sample sheet are known or unknown, it is
recommended to review the sample sheet before running the workflow.

#### 1.3.1. Splitting a sample sheet

It may be necessary to split a sample sheet into multiple parts, if for example,
it would take too long to run a workflow for all the samples in a single run. In
such cases, it is possible to split the sample sheet as follows. For example:
```shell
python scripts/split_sample_sheet.py --split-by dataset "samples.tsv" "samples_"
```
Depending on the configuration, this script can be used to generate several
smaller sample sheets, which can then be processed in separate workflow runs.

### 1.4. Running the workflow

The `protter` workflow can be run just like any other Snakemake workflow,
whether by specifying the number of cores and other parameters as needed…
```shell
snakemake --use-conda --cores 1
```

…or if using some HPC cluster, like CNIO's, by specifying a slurm profile:
```shell
snakemake --use-conda --profile $SMK_PROFILE_SLURM -j 10
```

### 1.5. Gathering output PSM files

After completion of a workflow run, the output
PSM files can be gathered as follows:
```shell
python scripts/gather_psm_files.py config.yaml protter_psm_output.zip
```
The input config YAML file is used to determine the datasets
for which PSM files should be gathered, and the output ZIP
archive contains all the output PSM files for those datasets.

### 1.6. Processing PSM files for APPRIS's Proteo

If the protter output is going to be used as part of APPRIS it is necessary to 
perform post-processing to obtain the final CSV file. To do this, a script has 
been prepared that takes the `protter_psm_output.zip` file and the `config.yaml` 
file as input; and can be run from the
workflow directory as follows:
```shell
python scripts/postprocessing_for_proteo.py config.yaml protter_psm_output.zip
```

This script creates a CSV file for each enabled database in `config.yaml` and saves 
them in the same path that the PSMs zipped file. The files name is `proteo_{database}.csv`.

## 2. Troubleshooting

There are several steps of the workflow during which issues may arise.

### 2.1. download_sample

When a sample download fails, it should be possible to identify the cause of
failure from the Snakemake log files. However, sometimes a download fails simply
because of adverse network conditions, in which case the workflow can be rerun
when network conditions are more favourable.

In general, if there is sufficient local storage space, and especially if the
same dataset will be processed multiple times, it is recommended to download
sample files separately and configure the sample sheet with their local file
paths.

### 2.2. raw_to_mzml

In a small number of cases, an error may occur while converting a RAW file to
mzML format.

In such cases, it may be possible to convert the RAW file to an mzML file
outside of the workflow, then use that mzML file as input for the given sample.

If no way can be found to convert the RAW file to mzML, the corresponding
sample can be dropped from the workflow, either by deleting the sample from
the sample sheet, or by marking the sample as `NA` in the `subset` column,
which effectively excludes it.

### 2.3. comet

In a small number of cases, Comet does not create an output file. Sometimes this
is because Comet has not searched any spectra. If this occurs and no other error
is found, protter will generate a placeholder PIN file so that the workflow will
continue uninterrupted.

However, it can save time and resources to simply exclude such samples from
the workflow, whether by deleting them from the sample sheet or marking them
as `NA` in the `subset` column.

### 2.4. percolator

In a small number of cases, Percolator may fail during the training phase. One
way to ameliorate this issue is to group samples by `biosample`, `experiment`,
or some other appropriate grouping.

## 3. Data sources

Input mass spectrometry data may be obtained from files stored locally,
or may be downloaded directly from the
[PRIDE database](https://www.ebi.ac.uk/pride/)
([Perez-Riverol et al. 2019](https://doi.org/10.1093/nar/gky1106))
by specifying one or more PRIDE project accessions.

## 4. Dependencies

This workflow requires Python 3.6 or greater.

For conversion of mass spectrometry data from Thermo RAW
formats to the standard mzML format, this workflow depends on
[ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
([Hulstaert et al. 2020](https://doi.org/10.1021/acs.jproteome.9b00328)).
ThermoRawFileParser is currently available under an Apache 2.0 License,
subject to further restrictions imposed by the license of the Thermo
Finnigan LLC vendor library on which it depends. These licenses can be
viewed in the [ThermoRawFileParser GitHub repository](
https://github.com/compomics/ThermoRawFileParser).

The protter workflow also depends on the following software:

- [Biopython](https://biopython.org/)
  ([Cock et al. 2009](https://doi.org/10.1093/bioinformatics/btp163))
  — [Biopython License](https://github.com/biopython/biopython/blob/master/LICENSE.rst)
- [Crux toolkit](https://crux.ms/)
  ([McIlwain et al. 2014](https://doi.org/10.1021/pr500741y))
  — Apache 2.0 License
  * [Comet](https://crux.ms/commands/comet.html)
    ([Eng et al. 2012](https://doi.org/10.1002/pmic.201200439))
  * [Percolator](https://crux.ms/commands/percolator.html)
    ([Käll et al. 2007](https://doi.org/10.1038/nmeth1113);
    [Serang et al. 2010](https://doi.org/10.1021/pr100594k))
- [curl](https://curl.se)
  — [MIT/X derivative license](http://curl.haxx.se/docs/copyright.html)
- [DecoyPYrat](https://anaconda.org/bioconda/decoypyrat)
  ([Wright & Choudhary 2016](https://doi.org/10.4172/jpb.1000404))
  — MIT License
- [pandas](https://pandas.pydata.org/)
  ([McKinney 2010](https://doi.org/10.25080/Majora-92bf1922-00a))
  — BSD 3-Clause License
- [PyYAML](https://pyyaml.org)
  — MIT License
- [ratelimiter](https://github.com/RazerM/ratelimiter)
  — Apache 2.0 License
- [requests](https://requests.readthedocs.io/en/master/)
  — Apache 2.0 License
- [rhash](https://github.com/rhash/RHash)
  — BSD Zero-Clause License
- [zlib](http://www.zlib.net)
  — [zlib License](http://www.zlib.net/zlib_license.html)

With specific datasets for which metadata is obtained from an Excel
file, packages such as [xlrd](https://pypi.org/project/xlrd/) or
[openpyxl](https://openpyxl.readthedocs.io/en/stable/) may also
be required to open XLS or XLSX files, respectively.
