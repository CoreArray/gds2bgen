# ===========================================================================
#
# gds2bgen.r: format conversion between GDS and BGEN
#
# Copyright (C) 2018-2019    Xiuwen Zheng (zhengxwen@gmail.com)
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# gds2bgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with gds2bgen.
# If not, see <http://www.gnu.org/licenses/>.

#############################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")


#############################################################
# Get bgen file information
#
seqBGEN_Info <- function(bgen.fn=NULL, verbose=TRUE)
{
    # check
    stopifnot(is.null(bgen.fn) | is.character(bgen.fn))
    if (is.null(bgen.fn))
        return(.Call(SEQ_BGEN_Info, NULL))

    stopifnot(length(bgen.fn)==1L, !is.na(bgen.fn))
    rv <- .Call(SEQ_BGEN_Info, bgen.fn)
    names(rv) <- c("num.sample", "num.variant", "compression", "layout",
        "unphased", "bits", "ploidy.min", "ploidy.max", "sample.id")
    if (verbose)
    {
        .cat("File: ", normalizePath(bgen.fn))
        .cat("# of samples: ", rv$num.sample)
        .cat("# of variants: ", rv$num.variant)
        .cat("Compression method: ", rv$compression)
        .cat("Layout version: ", rv$layout)
        .cat("Unphased: ", rv$unphased)
        .cat("# of bits: ", rv$bits)
        if (rv$ploidy.min == rv$ploidy.max)
            .cat("Ploidy: ", rv$ploidy.min)
        else
            .cat("Ploidy: [", rv$ploidy.min, ", ", rv$ploidy.max, "]")
       	cat("Sample id: ")
        if (is.null(rv$sample.id))
        {
            cat("<anonymized>\n")
        } else {
            if (length(rv$sample.id) > 4L)
                cat(c(rv$sample.id[seq_len(4L)], "..."), sep=", ")
            else
                cat(rv$sample.id, sep=", ")
            cat("\n")
        }
    }
    invisible(rv)
}



#############################################################
# Format conversion from BGEN to GDS
#
seqBGEN2GDS <- function(bgen.fn, out.fn, storage.option="LZMA_RA", float.type=
    c("packed8", "packed16", "single", "double", "sp.real32", "sp.real64"),
    geno=FALSE, dosage=TRUE, prob=FALSE, ignore.chr.prefix=c("chr", "0"),
    start=1L, count=-1L, sample.id=NULL, optimize=TRUE, digest=TRUE, parallel=FALSE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(bgen.fn), length(bgen.fn)==1L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    float.type <- match.arg(float.type)
    if (is.character(storage.option))
    {
        storage.option <- seqStorageOption(storage.option)
        s1 <- switch(float.type,
            packed8  = "packedreal8u:offset=0,scale=1/127",
            packed16 = "packedreal16u:offset=0,scale=1/32767",
            single   = "float32",
            double   = "float64",
            sp.real32 = "sp.real32",
            sp.real64 = "sp.real64")
        s2 <- switch(float.type,
            packed8  = "packedreal8u:offset=0,scale=1/254",
            packed16 = "packedreal16u:offset=0,scale=1/65534",
            single   = "float32",
            double   = "float64",
            sp.real32 = "sp.real32",
            sp.real64 = "sp.real64")
		storage.option$mode <- c(
			`annotation/format/DS`=s1,
			`annotation/format/GP`=s2
		)
    }
    stopifnot(inherits(storage.option, "SeqGDSStorageClass"))
    stopifnot(is.logical(geno), length(geno)==1L)
    stopifnot(is.logical(dosage), length(dosage)==1L)
    stopifnot(is.logical(prob), length(prob)==1L)
    stopifnot(is.character(ignore.chr.prefix), length(ignore.chr.prefix)>0L)
    stopifnot(is.numeric(start), length(start)==1L)
    stopifnot(is.numeric(count), length(count)==1L)

    # get bgen info
    info <- seqBGEN_Info(bgen.fn, verbose=FALSE)
    if (verbose)
    {
        .cat(date())
        # bgen information
        cat("BGEN Import:\n")
        cat("    file:", bgen.fn)
        cat(" (", SeqArray:::.pretty_size(file.size(bgen.fn)), ")\n", sep="")
        .cat("    # of samples: ", info$num.sample)
        .cat("    # of variants: ", info$num.variant)
        .cat("    bgen compression method: ", info$compression)
        .cat("    layout version: ", info$layout)
        .cat("    unphased: ", info$unphased)
        .cat("    # of bits: ", info$bits)
        if (info$ploidy.min == info$ploidy.max)
            .cat("    ploidy: ", info$ploidy.min)
        else
            .cat("    ploidy: [", info$ploidy.min, ", ", info$ploidy.max, "]")
       	cat("    sample id: ")
        if (is.null(info$sample.id))
        {
            cat("<anonymized>\n")
        } else {
            if (length(info$sample.id) > 4L)
                cat(c(info$sample.id[seq_len(4L)], "..."), sep=", ")
            else
                cat(info$sample.id, sep=", ")
            cat("\n")
        }
        # gds information
        .cat("Output:\n    ", out.fn)
        .cat("    saving genotypes [GT]: ", geno)
        .cat("    saving dosages [annotation/format/DS]: ", dosage)
        .cat("    saving probabilities [annotation/format/GP]: ", prob)
        flush.console()
    }

    # the number of samples and variants
    nSamp <- info$num.sample
    nVariant <- info$num.variant
    if (nSamp <= 0L)
        stop("No sample in the bgen file.")
    if (!is.null(sample.id))
    {
        stopifnot(is.vector(sample.id))
        if (length(sample.id) != nSamp)
        {
            stop(sprintf(
                "'sample.id' should have the same length as the bgen file (# %d).",
                nSamp))
        }
        info$sample.id <- sample.id
        if (verbose)
        {
            cat("    user-defined sample id: ")
            if (length(info$sample.id) > 4L)
                cat(c(info$sample.id[seq_len(4L)], "..."), sep=", ")
            else
                cat(info$sample.id, sep=", ")
            cat("\n")
        }
    }

    # the number of parallel tasks
    pnum <- SeqArray:::.NumParallel(parallel)
    if (pnum > 1L)
    {
        if (start < 1L)
            stop("'start' should be a positive integer if conversion in parallel.")
        else if (start > nVariant)
            stop("'start' should not be greater than the total number of variants.")
        if (count < 0L)
            count <- nVariant - start + 1L
        if (start+count > nVariant+1L)
            stop("Invalid 'count'.")

        if (count >= pnum)
        {
            fn <- sub("^([^.]*).*", "\\1", basename(out.fn))
            psplit <- SeqArray:::.file_split(count, pnum, start)

            # need unique temporary file names
            ptmpfn <- character()
            while (length(ptmpfn) < pnum)
            {
                s <- tempfile(pattern=sprintf("%s_tmp%02d_",
                    fn, length(ptmpfn)+1L), tmpdir=dirname(out.fn))
                file.create(s)
                if (!(s %in% ptmpfn)) ptmpfn <- c(ptmpfn, s)
            }
            if (verbose)
            {
                cat("    output to path: ", file.path(dirname(out.fn), ""), "\n", sep="")
                cat(sprintf("    writing to %d files:\n", pnum))
                cat(sprintf("        %s [%s..%s]\n", basename(ptmpfn),
                    SeqArray:::.pretty(psplit[[1L]]),
                    SeqArray:::.pretty(psplit[[1L]] + psplit[[2L]] - 1L)),
                    sep="")
                flush.console()
            }

            # conversion in parallel
            seqParallel(parallel, NULL,
                FUN = function(bgen.fn, storage.option, float.type,
                    dosage, geno, prob, optim, ptmpfn, psplit, verbose)
                {
                    # the process id, starting from one
                    i <- SeqArray:::process_index
                    attr(bgen.fn, "progress") <- TRUE
                    gds2bgen::seqBGEN2GDS(bgen.fn, ptmpfn[i],
                        storage.option = storage.option,
                        float.type = float.type, dosage = dosage, geno = geno, prob = prob,
                        start = psplit[[1L]][i], count = psplit[[2L]][i],
                        optimize = optim, digest = FALSE, parallel = FALSE,
                        verbose = FALSE
                    )
                    invisible()
                }, split = "none",
                bgen.fn = bgen.fn, storage.option = storage.option,
                float.type = float.type, dosage = dosage, geno = geno, prob = prob,
                optim = optimize, ptmpfn = ptmpfn, psplit = psplit, verbose = verbose
            )

            if (verbose)
                cat("    Done (", date(), ").\n", sep="")

        } else {
            pnum <- 1L
            message("No use of parallel environment!")
        }
    }


    # create a GDS file
    gfile <- createfn.gds(out.fn)
    on.exit({ if (!is.null(gfile)) closefn.gds(gfile) })

    put.attr.gdsn(gfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(gfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(gfile, "description")
    put.attr.gdsn(n, "bgen.version", info$layout)

    # add sample.id
    if (is.null(info$sample.id))
        info$sample.id <- seq_len(nSamp)
    SeqArray:::.AddVar(storage.option, gfile, "sample.id", info$sample.id,
        closezip=TRUE)

    # add variant.id
    SeqArray:::.AddVar(storage.option, gfile, "variant.id", storage="int32")

    # add position
    SeqArray:::.AddVar(storage.option, gfile, "position", storage="int32")

    # add chromosome
    SeqArray:::.AddVar(storage.option, gfile, "chromosome", storage="string")

    # add allele
    SeqArray:::.AddVar(storage.option, gfile, "allele", storage="string")

    # add a folder for genotypes
    varGeno <- addfolder.gdsn(gfile, "genotype")
    put.attr.gdsn(varGeno, "VariableName", "GT")
    put.attr.gdsn(varGeno, "Description", "Genotype")
    if (geno)
    {
        SeqArray:::.AddVar(storage.option, varGeno, "data", storage="bit2",
            valdim=c(2L, nSamp, 0L))
        SeqArray:::.AddVar(storage.option, varGeno, "@data", storage="uint8",
            visible=FALSE)
        node <- SeqArray:::.AddVar(storage.option, varGeno, "extra.index",
            storage="int32", valdim=c(3L,0L))
        put.attr.gdsn(node, "R.colnames", c("sample.index", "variant.index", "length"))
        SeqArray:::.AddVar(storage.option, varGeno, "extra", storage="int16")
    }

    # add phase folder
    varPhase <- addfolder.gdsn(gfile, "phase")
    if (geno)
    {
        SeqArray:::.AddVar(storage.option, varPhase, "data", storage="bit1",
            valdim=c(nSamp, 0L))
        node <- SeqArray:::.AddVar(storage.option, varPhase, "extra.index",
            storage="int32", valdim=c(3L,0L))
        put.attr.gdsn(node, "R.colnames", c("sample.index", "variant.index", "length"))
        SeqArray:::.AddVar(storage.option, varPhase, "extra", storage="bit1")
    }

    # add annotation folder
    varAnnot <- addfolder.gdsn(gfile, "annotation")

    # add id
    SeqArray:::.AddVar(storage.option, varAnnot, "id", storage="string")
    # add qual
    SeqArray:::.AddVar(storage.option, varAnnot, "qual", storage="float")
    # add filter
    varFilter <- SeqArray:::.AddVar(storage.option, varAnnot, "filter",
        storage="int32")

    # VCF INFO
    varInfo <- addfolder.gdsn(varAnnot, "info")

    # add the FORMAT field
    varFormat <- addfolder.gdsn(varAnnot, "format")
    if (dosage)
    {
        DS <- addfolder.gdsn(varFormat, "DS")
        put.attr.gdsn(DS, "Number", ".")
        put.attr.gdsn(DS, "Type", "Float")
        put.attr.gdsn(DS, "Description", "Estimated alternate allele dosage")
        SeqArray:::.AddVar(storage.option, DS, "data", storage="float",
            valdim=c(nSamp, 0L))
        SeqArray:::.AddVar(storage.option, DS, "@data", storage="int32",
            visible=FALSE)
        if (verbose)
            cat("    (writing to 'annotation/format/DS')\n")
    }
    if (prob)
    {
        GP <- addfolder.gdsn(varFormat, "GP")
        put.attr.gdsn(GP, "Number", ".")
        put.attr.gdsn(GP, "Type", "Float")
        put.attr.gdsn(GP, "Description",
            "Genotype probabilities (no prob of ref homozygous geno)")
        SeqArray:::.AddVar(storage.option, GP, "data", storage="float",
            valdim=c(nSamp, 0L))
        SeqArray:::.AddVar(storage.option, GP, "@data", storage="int32",
            visible=FALSE)
        if (verbose)
            cat("    (writing to 'annotation/format/GP')\n")
    }

    if (pnum <= 1L)
    {
        if (isTRUE(attr(bgen.fn, "progress")))
        {
            progfile <- file(paste0(out.fn, ".progress"), "wt")
            on.exit({
                close(progfile)
                unlink(paste0(out.fn, ".progress"), force=TRUE)
            }, add=TRUE)
        } else {
            progfile <- NULL
        }
        # call C function
        .Call(SEQ_BGEN_Import, bgen.fn, gfile$root, ignore.chr.prefix,
            start, count, progfile, verbose)
    } else {
        ## merge all temporary files
        varnm <- c("variant.id", "position", "chromosome", "allele",
            "annotation/id", "annotation/qual", "annotation/filter")
        if (dosage)
        {
            varnm <- c(varnm, c("annotation/format/DS/data",
                "annotation/format/DS/@data"))
        }
        if (prob)
        {
            varnm <- c(varnm, c("annotation/format/GP/data",
                "annotation/format/GP/@data"))
        }
        if (geno) {
            varnm <- c(varnm, c(
                "genotype/data", "genotype/@data",
                "genotype/extra.index", "genotype/extra",
                "phase/data",
                "phase/extra.index", "phase/extra" ))
        }

        if (verbose) cat("Merging:\n")

        # open all temporary files
        for (fn in ptmpfn)
        {
            if (verbose)
                cat("    adding ", sQuote(basename(fn)))
            # open the gds file
            tmpgds <- openfn.gds(fn)
            # merge variables
            for (nm in varnm)
                append.gdsn(index.gdsn(gfile, nm), index.gdsn(tmpgds, nm))
            # close the file
            closefn.gds(tmpgds)
            if (verbose) cat(" [done]\n")
        }

        # remove temporary files
        unlink(ptmpfn, force=TRUE)
    }


    # add annotation folder
    addfolder.gdsn(gfile, "sample.annotation")

    # RLE-coded chromosome
    SeqArray:::.optim_chrom(gfile)
    # create hash
    SeqArray:::.DigestFile(gfile, digest, verbose)

    # close the GDS file
    closefn.gds(gfile)
    gfile <- NULL

    if (verbose)
    {
        cat("Done.\n")
        cat(date(), "\n", sep="")
    }
    if (optimize)
    {
        if (verbose)
        {
            cat("Optimize the access efficiency ...\n")
            flush.console()
        }
        cleanup.gds(out.fn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.fn))
}
