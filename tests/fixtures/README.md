# Real-data fixture

`SRR11060618_subset.fasta` is copied verbatim from
[Zheludev/VNom](https://github.com/Zheludev/VNom) (`test_data/SRR11060618_subset.fasta`),
the reference implementation of the Viroid Nominator tool from the Obelisks
paper (Zheludev et al., *Viroid-like colonists of human microbiomes*, Cell
2024). VNom repo is MIT-licensed; copyright Zheludev.

It is a subset of rnaSPAdes contigs from real SRA run `SRR11060618` (a
peach-tissue metatranscriptome, per VNom's own README usage example: the
sample is used there as `peach_subset`). It is included here as a genuine,
non-synthetic positive-control fixture for the circularity detector -- the
original SRA/ENA hosts (`www.ncbi.nlm.nih.gov`, `www.ebi.ac.uk`) could not be
reached directly from this sandbox, but `raw.githubusercontent.com` is not on
its network blocklist, and this exact file happened to be vendored into a
public GitHub repo.
