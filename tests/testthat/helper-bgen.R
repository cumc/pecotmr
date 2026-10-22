# A minimal BGEN writer for the BGEN reader's tests (R/genotypeIo.R), of
# either layout: diploid, unphased, biallelic records whose stored
# probability values are given as integers, so the expected dosages are
# known exactly.

# Little-endian unsigned integers as raw bytes.
.bgw_u16 <- function(x) {
    as.raw(c(x %% 256, (x %/% 256) %% 256))
}

.bgw_u32 <- function(x) {
    as.raw(c(
        x %% 256,
        (x %/% 256) %% 256,
        (x %/% 65536) %% 256,
        x %/% 16777216
    ))
}

# A string as a 2- or 4-byte length and its bytes.
.bgw_str <- function(x, width = 2L) {
    bytes <- charToRaw(x)
    len <- if (width == 2L) .bgw_u16(length(bytes)) else .bgw_u32(length(bytes))
    c(len, bytes)
}

# B-bit values packed least significant bit first, padded to whole bytes.
.bgw_pack <- function(values, bits) {
    bitsOf <- unlist(lapply(values, function(v) {
        as.integer(intToBits(v))[seq_len(bits)]
    }))
    pad <- (8L - length(bitsOf) %% 8L) %% 8L
    packBits(as.raw(c(bitsOf, rep(0L, pad))), "raw")
}

.bgw_compress <- function(data, compression) {
    if (compression == 0L) {
        return(data)
    }
    if (compression == 1L) {
        return(memCompress(data, type = "gzip"))
    }
    tf <- tempfile(fileext = ".zst")
    con <- archive::archive_write(tf, "data", format = "raw", filter = "zstd")
    writeBin(data, con)
    close(con)
    readBin(tf, "raw", file.size(tf))
}

# One variant record. `values` is an nSamples x 2 integer matrix of the
# stored P(AA) and P(AB) values (out of 2^bits - 1, or 32768 in layout 1);
# `missing` flags samples with no data.
.bgw_record <- function(
    chrom,
    pos,
    rsid,
    a1,
    a2,
    values,
    missing,
    bits,
    compression,
    layout = 2L
) {
    if (layout == 1L) {
        return(.bgw_record1(
            chrom,
            pos,
            rsid,
            a1,
            a2,
            values,
            missing,
            compression
        ))
    }
    n <- nrow(values)
    ploidy <- as.raw(ifelse(missing, 2L + 128L, 2L))
    stored <- values
    stored[missing, ] <- 0L
    probs <- c(
        .bgw_u32(n),
        .bgw_u16(2),
        as.raw(2L),
        as.raw(2L),
        ploidy,
        as.raw(0L),
        as.raw(bits),
        .bgw_pack(as.vector(t(stored)), bits)
    )
    packed <- .bgw_compress(probs, compression)
    genotype <- if (compression == 0L) {
        c(.bgw_u32(length(packed)), packed)
    } else {
        c(.bgw_u32(length(packed) + 4L), .bgw_u32(length(probs)), packed)
    }
    c(
        .bgw_str(rsid),
        .bgw_str(rsid),
        .bgw_str(chrom),
        .bgw_u32(pos),
        .bgw_u16(2),
        .bgw_str(a1, 4L),
        .bgw_str(a2, 4L),
        genotype
    )
}

# A layout-1 record: the sample count, the identifying data without an
# allele count, and three 16-bit probabilities (AA, AB, BB) per sample out of
# 32768, all zero for a missing sample, zlib-compressed behind their length
# or stored as they are.
.bgw_record1 <- function(
    chrom,
    pos,
    rsid,
    a1,
    a2,
    values,
    missing,
    compression
) {
    n <- nrow(values)
    stored <- cbind(values, 32768L - values[, 1L] - values[, 2L])
    stored[missing, ] <- 0L
    probs <- writeBin(
        as.integer(t(stored)),
        raw(),
        size = 2L,
        endian = "little"
    )
    genotype <- if (compression == 0L) {
        probs
    } else {
        packed <- memCompress(probs, type = "gzip")
        c(.bgw_u32(length(packed)), packed)
    }
    c(
        .bgw_u32(n),
        .bgw_str(rsid),
        .bgw_str(rsid),
        .bgw_str(chrom),
        .bgw_u32(pos),
        .bgw_str(a1, 4L),
        .bgw_str(a2, 4L),
        genotype
    )
}

# Write a BGEN file and its .bgi index. `variants` is a data frame with
# chrom, pos, rsid, a1, a2; `values` a list of nSamples x 2 integer matrices
# (one per variant); `missing` a logical nSamples x nVariants matrix.
# `storeIds = FALSE` leaves the sample ids out of the file.
.bgw_write <- function(
    path,
    variants,
    values,
    missing,
    sampleIds,
    bits = 8L,
    compression = 1L,
    storeIds = TRUE,
    layout = 2L
) {
    n <- length(sampleIds)
    flags <- compression + layout * 4L
    flagBytes <- .bgw_u32(flags)
    if (storeIds) {
        flagBytes[[4L]] <- as.raw(128L)
    }
    header <- c(
        .bgw_u32(20L),
        .bgw_u32(nrow(variants)),
        .bgw_u32(n),
        charToRaw("bgen"),
        flagBytes
    )
    sampleBlock <- if (storeIds) {
        ids <- unlist(lapply(sampleIds, .bgw_str))
        c(.bgw_u32(length(ids) + 8L), .bgw_u32(n), ids)
    }
    offset <- length(header) + length(sampleBlock)
    records <- lapply(seq_len(nrow(variants)), function(j) {
        .bgw_record(
            variants$chrom[[j]],
            variants$pos[[j]],
            variants$rsid[[j]],
            variants$a1[[j]],
            variants$a2[[j]],
            values[[j]],
            missing[, j],
            bits,
            compression,
            layout
        )
    })
    bytes <- c(.bgw_u32(offset), header, sampleBlock, unlist(records))
    writeBin(bytes, path)
    starts <- 4L + offset + c(0L, cumsum(lengths(records)))[seq_along(records)]
    con <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, ".bgi"))
    on.exit(DBI::dbDisconnect(con))
    DBI::dbWriteTable(
        con,
        "Variant",
        data.frame(
            chromosome = variants$chrom,
            position = variants$pos,
            rsid = variants$rsid,
            number_of_alleles = rep(2L, nrow(variants)),
            allele1 = variants$a1,
            allele2 = variants$a2,
            file_start_position = starts,
            size_in_bytes = lengths(records)
        )
    )
    invisible(path)
}

# The expected A1 (second allele) dosages of `values` out of `scale`: P(AB) +
# 2 P(BB), with P(BB) one minus the stored two.
.bgw_dosage <- function(values, missing, scale) {
    out <- vapply(
        seq_along(values),
        function(j) {
            p <- values[[j]] / scale
            p[, 2] + 2 * (1 - p[, 1] - p[, 2])
        },
        numeric(nrow(missing))
    )
    out[missing] <- NA
    out
}

# A small fixture: `nVar` variants of `n` samples with random stored values
# that sum to at most the layout's scale, one missing sample, written to a
# tempfile.
.bgw_fixture <- function(
    n = 6L,
    nVar = 5L,
    bits = 8L,
    compression = 1L,
    storeIds = TRUE,
    seed = 1L,
    layout = 2L
) {
    set.seed(seed)
    scale <- if (layout == 1L) 32768 else 2^bits - 1
    values <- lapply(seq_len(nVar), function(j) {
        aa <- sample(0:scale, n, replace = TRUE)
        ab <- vapply(scale - aa, function(r) sample(0:r, 1L), numeric(1))
        cbind(aa, ab)
    })
    missing <- matrix(FALSE, n, nVar)
    missing[2L, 3L] <- TRUE
    variants <- data.frame(
        chrom = "01",
        pos = 100L * seq_len(nVar),
        rsid = c(".", sprintf("rs%d", seq_len(nVar - 1L))),
        a1 = "A",
        a2 = "G"
    )
    ids <- sprintf("s%d", seq_len(n))
    dir <- withr::local_tempdir(.local_envir = parent.frame())
    path <- file.path(dir, "fixture.bgen")
    .bgw_write(
        path,
        variants,
        values,
        missing,
        ids,
        bits,
        compression,
        storeIds,
        layout
    )
    if (!storeIds) {
        writeLines(
            c("ID_1 ID_2 missing", "0 0 0", paste(ids, ids, 0)),
            file.path(dir, "fixture.sample")
        )
    }
    list(
        path = path,
        ids = ids,
        variants = variants,
        dosage = .bgw_dosage(values, missing, scale)
    )
}
