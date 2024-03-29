# Overview

`kmers_ld` is a program used to compute the pairwise linkage disequilibrium (LD)
between all *k*-mers from a presence/absence variation (PAV) table of these
*k*-mers. Such a PAV table can be obtained from the `filter_kmers` program
of the [*k*-mer GWAS](https://github.com/voichek/kmersGWAS) suite of tools.

# Installation

## Dependencies

`kmers_ld` requires the `khash.h` header file from the
[htslib](https://github.com/samtools/htslib) library for functioning. Make sure
you have this library installed before compiling the program.

## Software compilation

After downloading the software with `git clone https://github.com/malemay/kmers_ld`,
compiling the software is as simple as running:

	cd kmers_ld
	make

The `kmers_ld` binary can then be used with the following usage:

	kmers_ld <pav_table.txt>

with the output being written to stdout.

## Testing

The software can be tested by running `make test`. The file `pav_test.txt`
shows what a typicaly input PAV table should look like.

# Citation

If you use this software, plase cite our publication:

Lemay, M.-A., de Ronne, M., Bélanger, R., & Belzile, F. (2023). k-mer-based GWAS enhances the discovery of causal variants and candidate genes in soybean. *The Plant Genome*, 16, e20374. [doi:10.1002/tpg2.20374](https://doi.org/10.1002/tpg2.20374)


# References

Voichek, Y. and Weigel, D., 2020. Identifying genetic variants underlying
phenotypic variation in plants without complete genomes. *Nature Genetics*,
52(5), pp.534-540. <https://doi.org/10.1038/s41588-020-0612-7>
