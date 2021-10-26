library(dada2)

path <- "/mnt/research/germs/jihoon/dada2"

set.seed(1390)
CORES = parallel::detectCores()

# Forward and reverse fastq filenames have the format:
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

# Filter and trim
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = CORES)

# If the number of sequence files were different after filtering
path_edited <- "/mnt/research/germs/jihoon/dada2/filtered"

fnFs_edited <- sort(list.files(path_edited, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs_edited <- sort(list.files(path_edited, pattern = "_R2_001.fastq.gz", full.names = TRUE))
errF <- learnErrors(fnFs_edited, multithread = CORES)
errR <- learnErrors(fnRs_edited, multithread = CORES)

# Dereplicate identical reads
derepFs <- derepFastq(fnFs_edited, verbose = TRUE)
derepRs <- derepFastq(fnRs_edited, verbose = TRUE)

# Name the derep-class objects by the sample names
sample.names_edited <- unname(sapply(fnFs_edited, get.sample.name))
names(derepFs) <- sample.names_edited
names(derepRs) <- sample.names_edited

# Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = CORES)
dadaRs <- dada(derepRs, err = errR, multithread = CORES)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct Sequence Table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=CORES, verbose=TRUE)

# Assign taxanomy
unite.ref <- "/mnt/research/germs/jihoon/dada2/sh_general_release_dynamic_10.05.2021.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = CORES, tryRC = TRUE)

saveRDS(t(seqtab.nochim), paste0(path_edited, "/seq_table.RDS"))
saveRDS(taxa, paste0(path_edited, "/tax_table.RDS"))





