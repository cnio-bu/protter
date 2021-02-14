# protter

Pipeline to generate peptide-spectrum matches using different
search engines against target and decoy databases. Input may
be local or remote (see [Data sources](#data-sources)).

Allows for the easy addition of datasets and search engines.

## Data sources

Input mass spectrometry data may be obtained from files stored locally,
or may be downloaded directly from the
[PRIDE database](https://www.ebi.ac.uk/pride/)
([Perez-Riverol et al. 2019](https://doi.org/10.1093/nar/gky1106))
by specifying a PRIDE project accession.

## Dependencies

For conversion of mass spectrometry data from Thermo RAW
formats to the standard mzML format, this workflow depends on
[ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
([Hulstaert et al. 2020](https://doi.org/10.1021/acs.jproteome.9b00328)).
ThermoRawFileParser is currently available under an Apache 2.0 License,
subject to further restrictions imposed by the license of the Thermo
Finnigan LLC vendor library on which it depends. These licenses can be
viewed in the [ThermoRawFileParser GitHub repository](
https://github.com/compomics/ThermoRawFileParser).

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
- [ratelimiter](https://github.com/RazerM/ratelimiter)
  — Apache 2.0 License
- [requests](https://requests.readthedocs.io/en/master/)
  — Apache 2.0 License
- [tstk](https://anaconda.org/tdido/tstk)
  — MIT License