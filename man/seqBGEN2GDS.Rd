\name{seqBGEN2GDS}
\alias{seqBGEN2GDS}
\title{Reformat BGEN Files}
\description{
    Reformats bgen files to GDS format.
}

\usage{
seqBGEN2GDS(bgen.fn, out.fn, storage.option="LZMA_RA", float.type=
    c("packed8", "packed16", "single", "double", "sp.real32", "sp.real64"),
    geno=FALSE, dosage=TRUE, prob=FALSE, ignore.chr.prefix=c("chr", "0"),
    start=1L, count=-1L, sample.id=NULL, optimize=TRUE, digest=TRUE,
    parallel=FALSE, verbose=TRUE)
}

\arguments{
    \item{bgen.fn}{the file name of a bgen file}
    \item{out.fn}{the file name of output GDS file}
    \item{storage.option}{specify the storage and compression option,
        "ZIP_RA" (\code{\link{seqStorageOption}("ZIP_RA")}); or
        "LZMA_RA" to use LZMA compression algorithm with higher compression
        ratio by default; or "LZ4_RA" to use an extremely fast compression and
        decompression algorithm. "ZIP_RA.max", "LZMA_RA.max" and "LZ4_RA.max"
        correspond to the algorithms with a maximum compression level; the
        suffix "_RA" indicates that fine-level random access is available; see
        more details at \code{\link{seqStorageOption}}}
    \item{float.type}{"packed8": use packedreal8u to store real numbers (dense
        matrix), "packed16": use packedreal16u for dense matrix; "single":
        single-precision floating-point numbers (dense matrix);
        "double": double-precision floating-point numbers (dense matrix);
        "sp.real32": single-precision sparse matrix;
        "sp.real64": double-precision sparse matrix}
    \item{geno}{if TRUE, store integer genotypes}
    \item{dosage}{if TRUE, store floating-point dosages in "annotation/format/DS",
        i.e., alternative allele dosages}
    \item{prob}{if TRUE, store floating-point probabilities in "annotation/format/GP",
        i.e., prob of heterozygous geno, prob of homozygous alternative allele}
    \item{ignore.chr.prefix}{a vector of character, indicating the prefix of
        chromosome which should be ignored, e.g., \code{"chr"}; it is not
        case-sensitive}
    \item{start}{the starting variant if importing part of VCF files}
    \item{count}{the maximum count of variant if importing part of VCF files,
        -1 indicates importing to the end}
    \item{sample.id}{specify the sample IDs}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add md5 hash codes to the GDS file
        if TRUE or a digest algorithm is specified}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing), a numeric value indicating the number of cores, or a
        cluster object for parallel processing; \code{parallel} is passed to
        the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details}
    \item{verbose}{if \code{TRUE}, show information}
}

\value{
    Return the file name of SeqArray GDS file with an absolute path.
}

\references{
    \url{https://bitbucket.org/gavinband/bgen},
    \url{http://bgenformat.org}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqBGEN_Info}}
}

\examples{
bgen.fn <- system.file("extdata", "example.8bits.bgen", package="gds2bgen")

seqBGEN2GDS(bgen.fn, "test.gds")

# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{BGEN}
\keyword{utilities}
