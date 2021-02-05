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
