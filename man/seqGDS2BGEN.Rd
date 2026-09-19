\name{seqGDS2BGEN}
\alias{seqGDS2BGEN}
\title{Convert GDS to BGEN}
\description{
    Reformats a SeqArray GDS file to a bgen file (layout v1.2).
}

\usage{
seqGDS2BGEN(gdsfn, out.fn, prob.source=c("auto", "GP", "DS", "GT"),
    bits=8L, compression=c("zlib", "zstd", "none"), sample.id=NULL,
    block.size=1024L, verbose=TRUE)
}

\arguments{
    \item{gdsfn}{a file name of the input SeqArray GDS file, or a
        \code{\link[SeqArray]{SeqVarGDSClass}} object; if an opened object is
        passed, the current sample and variant filter is respected}
    \item{out.fn}{the file name of the output bgen file}
    \item{prob.source}{the source of genotype probabilities written to the bgen
        file: \code{"GP"} uses \code{"annotation/format/GP"} (lossless);
        \code{"DS"} derives probabilities from the alternate allele dosages in
        \code{"annotation/format/DS"} (lossy, see Details); \code{"GT"} uses the
        integer genotypes in \code{"genotype"} as hard calls (lossy);
        \code{"auto"} (the default) selects the first available among GP, DS and
        GT, in that order}
    \item{bits}{the number of bits used to store each probability (1--32);
        a larger value gives higher precision and a larger file}
    \item{compression}{the compression method for the genotype data blocks:
        \code{"zlib"} (default), \code{"zstd"}, or \code{"none"}}
    \item{sample.id}{\code{NULL} to write the sample IDs stored in the GDS file;
        a character vector to use user-defined sample IDs; or \code{FALSE} to
        write an anonymized bgen file without a sample identifier block}
    \item{block.size}{the number of variants processed in each block}
    \item{verbose}{if \code{TRUE}, show information}
}

\details{
    Only bi-allelic, diploid variants are supported, mirroring the bi-allelic
    limitation of the import path; multiallelic sites should be split first.

    When \code{prob.source="DS"}, the alternate allele dosage \eqn{d \in [0, 2]}
    is mapped deterministically to genotype probabilities that reproduce the
    dosage: \eqn{(1-d, d, 0)} when \eqn{d \le 1}, and \eqn{(0, 2-d, d-1)} when
    \eqn{d > 1}. This is lossy because the original heterozygous/homozygous split
    cannot be recovered from a dosage alone.

    When \code{prob.source="GT"}, the integer genotype is written as a
    deterministic (0/1) probability for the called genotype.

    Missing genotypes are written as missing probabilities in the bgen file.
}

\value{
    Return the file name of the output bgen file with an absolute path.
}

\references{
    \url{http://bgenformat.org}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqBGEN2GDS}}, \code{\link{seqBGEN_Info}}
}

\examples{
bgen.fn <- system.file("extdata", "example.8bits.bgen", package="gds2bgen")

# bgen ==> gds (storing dosages)
seqBGEN2GDS(bgen.fn, "test.gds", dosage=TRUE, verbose=FALSE)

# gds ==> bgen
seqGDS2BGEN("test.gds", "test.bgen", prob.source="DS")
seqBGEN_Info("test.bgen")

# delete the temporary files
unlink(c("test.gds", "test.bgen"), force=TRUE)
}

\keyword{GDS}
\keyword{BGEN}
\keyword{utilities}
