# protter

Pipeline to generate peptide-spectrum matches using different search engines against target and decoy databases.

Allows for the easy addition of datasets and search engines.

Note the the `msconvert` rule for converting raw files to mzML format currently
requires access to a computer with Docker and rsync installed.

## Data sources

Input mass spectrometry data may be obtained from files stored locally,
or may be downloaded directly from the
[PRIDE database](https://www.ebi.ac.uk/pride/)
([Perez-Riverol et al. 2019](https://doi.org/10.1093/nar/gky1106))
by specifying a PRIDE project accession.

## Dependencies

For conversion of mass spectrometry data from proprietary vendor
formats to the standard mzML format, this workflow depends on
[ProteoWizard](http://www.proteowizard.org)
([Kessner et al. 2008](https://doi.org/10.1093/bioinformatics/btn323);
[Chambers et al. 2012](https://doi.org/10.1038/nbt.2377)),
in particular the ProteoWizard Docker image
[pwiz-skyline-i-agree-to-the-vendor-licenses](
https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses).
As its name implies, use of this ProteoWizard Docker image entails agreement with the
ProteoWizard Apache 2.0 License, subject to further restrictions imposed by the licenses
of the vendor libraries on which it depends. These licenses can be viewed on the
[ProteoWizard licenses webpage](http://www.proteowizard.org/licenses.html).

The protter workflow also depends on the following packages:

- [Crux toolkit](https://crux.ms/)
  ([McIlwain et al. 2014](https://doi.org/10.1021/pr500741y))
  — Apache 2.0 License
- [dateutil](https://anaconda.org/anaconda/python-dateutil)
  — BSD 3-Clause License
- [DecoyPYrat](https://anaconda.org/bioconda/decoypyrat)
  ([Wright & Choudhary 2016](https://doi.org/10.4172/jpb.1000404))
  — MIT License
- [pandas](https://pandas.pydata.org/)
  ([McKinney 2010](https://doi.org/10.25080/Majora-92bf1922-00a))
  — BSD 3-Clause License
- [PyYAML](https://pyyaml.org)
  — MIT License
- [requests](https://requests.readthedocs.io/en/master/)
  — Apache 2.0 License
- [Singularity](https://sylabs.io)
  ([Kurtzer et al. 2017](https://doi.org/10.1371/journal.pone.0177459))
  — BSD 3-Clause License
- [tstk](https://anaconda.org/tdido/tstk)
  — MIT License