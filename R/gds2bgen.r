# ===========================================================================
#
# gds2bgen.r: format conversion between GDS and BGEN
#
# Copyright (C) 2018-2025    Xiuwen Zheng (zhengxwen@gmail.com)
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

tm <- function() strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")


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
    c("packed8", "packed16", "single", "double", "sp.real8u", "sp.real16u",
    "sp.real32", "sp.real64"),
    geno=FALSE, dosage=TRUE, prob=FALSE, ignore.chr.prefix=c("chr", "0"),
    start=1L, count=-1L, sample.id=NULL, optimize=TRUE, digest=TRUE, parallel=FALSE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(bgen.fn), length(bgen.fn)==1L)
    stopifnot(is.character(out.fn), length(out.fn)==1L)
    float.type <- match.arg(float.type)
    if (float.type %in% c("sp.real8u", "sp.real16u"))
    {
        if (packageVersion("gdsfmt") < "1.49.9")
        {
            stop("float.type='", float.type, "' requires gdsfmt (>= 1.49.9), ",
                "but gdsfmt v", packageVersion("gdsfmt"), " is installed.")
        }
    }
    if (is.character(storage.option))
    {
        storage.option <- seqStorageOption(storage.option)
        s1 <- switch(float.type,
            packed8  = "packedreal8u:offset=0,scale=1/127",
            packed16 = "packedreal16u:offset=0,scale=1/32767",
            single   = "float32",
            double   = "float64",
            sp.real8u  = "sp.real8u:scale=1/127",
            sp.real16u = "sp.real16u:scale=1/32767",
            sp.real32  = "sp.real32",
            sp.real64  = "sp.real64")
        s2 <- switch(float.type,
            packed8  = "packedreal8u:offset=0,scale=1/254",
            packed16 = "packedreal16u:offset=0,scale=1/65534",
            single   = "float32",
            double   = "float64",
            sp.real8u  = "sp.real8u:scale=1/254",
            sp.real16u = "sp.real16u:scale=1/65534",
            sp.real32  = "sp.real32",
            sp.real64  = "sp.real64")
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
        .cat("##< ", tm())
        # bgen information
        cat("BGEN Import:\n")
        .cat("    file (", SeqArray:::.pretty_size(file.size(bgen.fn)), "):")
        .cat("        ", bgen.fn)
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
            progfile <- file(paste0(out.fn, ".progress.txt"), "wt")
            on.exit({
                close(progfile)
                unlink(paste0(out.fn, ".progress.txt"), force=TRUE)
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
                cat("    adding", sQuote(basename(fn)))
            # open the gds file
            tmpgds <- openfn.gds(fn)
            # merge variables
            for (nm in varnm)
                append.gdsn(index.gdsn(gfile, nm), index.gdsn(tmpgds, nm))
            # close the file
            closefn.gds(tmpgds)
            if (verbose) .cat(" [", tm(), " done]")
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
        if (optimize) .cat("Done.  # ", tm()) else cat("Done.\n")
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.fn, verbose=verbose)
    }
    if (verbose) .cat("##> ", tm())

    # output
    invisible(normalizePath(out.fn))
}



#############################################################
# Format conversion from GDS to BGEN
#

seqGDS2BGEN <- function(gdsfn, out.fn, prob.source=c("auto", "GP", "DS", "GT"),
    bits=8L, compression=c("zlib", "zstd", "none"), sample.id=NULL,
    block.size=1024L, verbose=TRUE)
{
    # check
    stopifnot(is.character(out.fn), length(out.fn)==1L, !is.na(out.fn))
    prob.source <- match.arg(prob.source)
    compression <- match.arg(compression)
    stopifnot(is.numeric(bits), length(bits)==1L, bits>=1L, bits<=32L)
    bits <- as.integer(bits)
    stopifnot(is.numeric(block.size), length(block.size)==1L, block.size>=1L)
    block.size <- as.integer(block.size)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # open the GDS file
    if (is.character(gdsfn))
    {
        stopifnot(length(gdsfn)==1L, !is.na(gdsfn))
        f <- SeqArray::seqOpen(gdsfn)
        on.exit(SeqArray::seqClose(f))
    } else {
        stopifnot(inherits(gdsfn, "SeqVarGDSClass"))
        f <- gdsfn
    }

    # dimension (respecting the current sample/variant filter)
    sid <- SeqArray::seqGetData(f, "sample.id")
    nSamp <- length(sid)
    nVar <- length(SeqArray::seqGetData(f, "variant.id"))
    if (nSamp <= 0L) stop("No selected sample in the GDS file.")
    if (nVar <= 0L) stop("No selected variant in the GDS file.")

    # variant identifying data
    chr <- as.character(SeqArray::seqGetData(f, "chromosome"))
    pos <- as.integer(SeqArray::seqGetData(f, "position"))
    rsid <- as.character(SeqArray::seqGetData(f, "annotation/id"))
    rsid[is.na(rsid)] <- ""
    allele <- SeqArray::seqGetData(f, "allele")
    sp <- strsplit(allele, ",", fixed=TRUE)
    if (any(lengths(sp) != 2L))
        stop("Only bi-allelic variants are supported; please split multiallelic sites first.")
    ref <- vapply(sp, `[`, "", 1L)
    alt <- vapply(sp, `[`, "", 2L)

    # sample id written into the bgen file
    if (is.null(sample.id))
    {
        bgen.sid <- as.character(sid)
    } else if (isFALSE(sample.id)) {
        bgen.sid <- NULL   # anonymized output
    } else {
        stopifnot(is.vector(sample.id), length(sample.id)==nSamp)
        bgen.sid <- as.character(sample.id)
    }

    # resolve "auto": prefer GP (lossless), then DS, then GT
    .has <- function(p) !is.null(index.gdsn(f, p, silent=TRUE))
    if (prob.source == "auto")
    {
        if (.has("annotation/format/GP")) prob.source <- "GP"
        else if (.has("annotation/format/DS")) prob.source <- "DS"
        else if (.has("genotype/data")) prob.source <- "GT"
        else stop("None of 'GP', 'DS' or 'GT' is available in the GDS file.")
    }

    # the FORMAT/genotype field to read
    varnm <- switch(prob.source,
        GP = "annotation/format/GP",
        DS = "annotation/format/DS",
        GT = "genotype")
    chk <- switch(prob.source, GT="genotype/data", varnm)
    if (!.has(chk))
        stop(sprintf("'%s' is not available; choose another 'prob.source'.", varnm))

    if (verbose)
    {
        .cat("##< ", tm())
        cat("GDS Import:\n")
        .cat("    file: ", f$filename)
        .cat("    # of samples: ", nSamp)
        .cat("    # of variants: ", nVar)
        .cat("    probability source: ", prob.source, " (",
            switch(prob.source,
                GP="lossless", DS="dosage-derived, lossy", GT="hard-call"), ")")
        cat("BGEN Output:\n")
        .cat("    file: ", out.fn)
        .cat("    layout version: v1.2")
        .cat("    compression: ", compression)
        .cat("    # of bits: ", bits)
        .cat("    sample id: ", if (is.null(bgen.sid)) "<anonymized>" else "yes")
        flush.console()
    }

    # create the bgen writer (writes the header & sample-id blocks)
    ptr <- .Call(SEQ_BGEN_ExportBegin, out.fn, nSamp, nVar, bgen.sid, bits,
        compression, FALSE)

    # running variant offset, shared with the block callback
    env <- new.env()
    env$off <- 0L

    # convert dosage -> (3 x nSamp) genotype probabilities, lossy
    .ds2prob <- function(d, k)
    {
        na <- is.na(d)
        p1 <- p2 <- p3 <- matrix(0, nSamp, k)
        lo <- !na & (d <= 1); hi <- !na & (d > 1)
        p1[lo] <- 1 - d[lo]; p2[lo] <- d[lo]
        p2[hi] <- 2 - d[hi]; p3[hi] <- d[hi] - 1
        p1[na] <- NaN; p2[na] <- NaN; p3[na] <- NaN
        array(rbind(as.vector(p1), as.vector(p2), as.vector(p3)),
            dim=c(3L, nSamp, k))
    }
    # convert alt-allele counts -> (3 x nSamp) one-hot probabilities
    .cnt2prob <- function(cnt, k)
    {
        na <- is.na(cnt)
        p1 <- p2 <- p3 <- matrix(0, nSamp, k)
        p1[cnt==0L] <- 1; p2[cnt==1L] <- 1; p3[cnt==2L] <- 1
        p1[na] <- NaN; p2[na] <- NaN; p3[na] <- NaN
        array(rbind(as.vector(p1), as.vector(p2), as.vector(p3)),
            dim=c(3L, nSamp, k))
    }

    # flush a block of variants to the bgen file
    .flush <- function(prob, k)
    {
        idx <- env$off + seq_len(k)
        .Call(SEQ_BGEN_ExportData, ptr, chr[idx], pos[idx], rsid[idx],
            ref[idx], alt[idx], prob)
        env$off <- env$off + k
        invisible()
    }

    if (prob.source == "DS")
    {
        SeqArray::seqBlockApply(f, varnm, function(x) {
            x <- as.matrix(x)            # nSamp x k
            k <- ncol(x)
            .flush(.ds2prob(x, k), k)
        }, margin="by.variant", as.is="none", bsize=block.size)

    } else if (prob.source == "GT") {
        SeqArray::seqBlockApply(f, varnm, function(x) {
            # x: array of dim c(2, nSamp, k); 0=ref, 1=alt, NA=missing
            k <- dim(x)[3L]
            cnt <- x[1L,,] + x[2L,,]     # nSamp x k
            dim(cnt) <- c(nSamp, k)
            .flush(.cnt2prob(cnt, k), k)
        }, margin="by.variant", as.is="none", bsize=block.size)

    } else {  # GP, reconstruct full genotype probabilities, lossless
        # seqBGEN2GDS stores GP without the ref-homozygous probability and with
        # an '@data' length that does not match the stored rows, so SeqArray's
        # decoder cannot read it; read the raw nodes directly instead.
        gp_data <- read.gdsn(index.gdsn(f, "annotation/format/GP/data"))  # nSampAll x totalRows
        nSampAll <- length(SeqArray::seqGetFilter(f)$sample.sel)
        if (is.null(dim(gp_data))) dim(gp_data) <- c(nSampAll, length(gp_data)/nSampAll)
        # the raw node ignores the sample filter, so subset to the selected samples
        ssel <- which(SeqArray::seqGetFilter(f)$sample.sel)
        if (length(ssel) != nSampAll) gp_data <- gp_data[ssel, , drop=FALSE]
        gp_len <- read.gdsn(index.gdsn(f, "annotation/format/GP/@data"))  # per-variant max length
        rows_v <- ifelse(gp_len >= 2L, gp_len - 1L, gp_len)               # stored rows per variant
        col_end <- cumsum(rows_v)
        col_beg <- col_end - rows_v + 1L
        # raw indices of the selected variants, in ascending (output) order
        sel <- which(SeqArray::seqGetFilter(f)$variant.sel)
        if (length(sel) != nVar)
            stop("Internal: variant selection does not match the GP data.")

        for (j in seq_len(nVar))
        {
            r <- sel[j]
            L <- rows_v[r]
            m <- gp_data[, col_beg[r]:col_end[r], drop=FALSE]   # nSamp x L
            prob <- matrix(NaN, 3L, nSamp)
            if (L >= 2L)
            {
                het <- m[, L-1L]; altp <- m[, L]
                ok <- !is.na(het) & !is.na(altp)
                prob[1L, ok] <- 1 - het[ok] - altp[ok]
                prob[2L, ok] <- het[ok]
                prob[3L, ok] <- altp[ok]
            } else {
                v <- m[, 1L]; ok <- !is.na(v)
                prob[1L, ok] <- v[ok]
            }
            .flush(array(prob, dim=c(3L, nSamp, 1L)), 1L)
        }
    }

    # close the bgen file (verifies the variant count)
    .Call(SEQ_BGEN_ExportEnd, ptr)

    if (verbose)
        .cat("Done.  ##> ", tm())

    invisible(normalizePath(out.fn))
}
