# MAD input fixture

`mad_control.in` is a copy of `control.in` from the MAD-1.5 data record by Cesare Malosso, Filippo Bigi, Paolo Pegolo, Joseph W. Abbott, Philip Loche, Mariana Rossi, Michele Ceriotti, and Arslan Mazitov, with the trailing blank line removed:

- Source: https://archive.materialscloud.org/records/18tke-tt476/files/control.in
- DOI: https://doi.org/10.24435/materialscloud:jc-9f
- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- Original download MD5: `24bf9cb988bc197477cde6b910e29f40`

The MAD-1.6 [record](https://archive.materialscloud.org/records/ffd5a-1hc37) provides a byte-identical `control.in` with the same checksum. The tests compare the generated global input flags against this fixture, including the two repeated `adjust_scf_param` entries. Species files used in tests are synthetic; no FHI-aims species defaults are redistributed.
