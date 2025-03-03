gds2bgen: Format Conversion from BGEN to GDS
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Description

This package provides functions for format conversion from [bgen](http://bgenformat.org) files to [SeqArray GDS](https://github.com/zhengxwen/SeqArray) files.

**Importance**: gds2bgen has been moved to https://github.com/CoreArray/gds2bgen from https://github.com/zhengxwen/gds2bgen since Feb 2025.

## Version

v0.9.4


## Package Maintainer

Dr. Xiuwen Zheng ([zhengx@u.washington.edu](zhengx@u.washington.edu))


## Installation

Requires R (≥ v3.5.0), [gdsfmt](http://www.bioconductor.org/packages/gdsfmt) (≥ v1.20.0), [SeqArray](http://www.bioconductor.org/packages/SeqArray) (≥ v1.24.0)

* Installation from Github:
```R
library("devtools")
install_github("CoreArray/gds2bgen")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.

Or manually intall the package
```sh
git clone https://github.com/CoreArray/gds2bgen
cd gds2bgen/src
unzip bgen_v1.1.8.zip
cd bgen_v1.1.8
python ./waf configure
python ./waf
cp build/libbgen.a ..
cp build/3rd_party/zstd-1.1.0/libzstd.a ..
rm -rf build
sleep 1; touch ../libbgen.a
cd ../../..
R CMD INSTALL gds2bgen
```


## Copyright Notice

This package includes the sources of the bgen library (https://enkre.net/cgi-bin/code/bgen/dir?ci=trunk), Boost (the C++
libraries, https://www.boost.org) and Zstandard (https://zstd.net).


## Citations for GDS

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012). A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. *Bioinformatics*. [DOI: 10.1093/bioinformatics/bts606](http://dx.doi.org/10.1093/bioinformatics/bts606).

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Examples

```R
library(gds2bgen)

seqBGEN_Info()  # bgen library version
## "bgen_lib_v1.1.8"

bgen_fn <- system.file("extdata", "example.8bits.bgen", package="gds2bgen")
# or bgen_fn <- "your_bgen_file.bgen"
seqBGEN_Info(bgen_fn)

## File: gds2bgen/extdata/example.8bits.bgen
## # of samples: 500
## # of variants: 199
## Compression method: zlib
## Layout version: v1.2
## Unphased: TRUE
## # of bits: 8
## Ploidy: 2
## sample id: sample_001, sample_002, sample_003, sample_004, ...


# example.8bits.bgen ==> example.gds, using 4 cores
seqBGEN2GDS(bgen_fn, "example.gds",
    storage.option="LZMA_RA",  # compression option, e.g., ZIP_RA for zlib or LZ4_RA for LZ4
    float.type="packed8",      # 8-bit packed real numbers
    geno=FALSE,     # 2-bit integer genotypes, stored in 'genotype/data'
    dosage=TRUE,    # numeric alternative allele dosages, stored in 'annotation/format/DS'
    prob=FALSE,     # numeric genotype probabilities, stored in 'annotation/format/GP'
    parallel=4      # the number of cores
)


# show file structure
library(SeqArray)
(f <- seqOpen("example.gds"))
seqClose(f)

## File: example.gds (137.7K)
## +    [  ] *
## |--+ description   [  ] *
## |--+ sample.id   { Str8 500 LZMA_ra(7.02%), 393B } *
## |--+ variant.id   { Int32 199 LZMA_ra(33.9%), 277B } *
## |--+ position   { Int32 199 LZMA_ra(60.6%), 489B } *
## |--+ chromosome   { Str8 199 LZMA_ra(15.7%), 101B } *
## |--+ allele   { Str8 199 LZMA_ra(11.8%), 101B } *
## |--+ genotype   [  ] *
## |--+ phase   [  ]
## |--+ annotation   [  ]
## |  |--+ id   { Str8 199 LZMA_ra(18.6%), 321B } *
## |  |--+ qual   { Float32 199 LZMA_ra(11.8%), 101B } *
## |  |--+ filter   { Int32 199 LZMA_ra(11.3%), 97B } *
## |  |--+ info   [  ]
## |  \--+ format   [  ]
## |     |--+ DS   [  ] *
## |     |  \--+ data   { PackedReal8U 500x199 LZMA_ra(55.6%), 54.0K } *
## \--+ sample.annotation   [  ]
```


## Also See

[seqVCF2GDS()](https://rdrr.io/bioc/SeqArray/man/seqVCF2GDS.html) in the [SeqArray](https://bioconductor.org/packages/SeqArray) package, conversion from VCF files to GDS files.

[seqBED2GDS()](https://rdrr.io/bioc/SeqArray/man/seqBED2GDS.html) in the [SeqArray](https://bioconductor.org/packages/SeqArray) package, conversion from PLINK BED files to GDS files.
