# Investigate if it is ok to use overfit PRS for interaction analysis
library(data.table)
library(magrittr)
library(optparse)
options(scipen = 999)

# Background --------------------------------------------------------------

option_list <- list(
    make_option(c("--bfile", "-b"), type = "character", help = "PLINK prefix"),
    make_option(c("--fam", "-f"), type = "character", help = "QCed sample list"),
    make_option(c("--snp", "-s"), type = "character", help = "QCed SNP list"),
    make_option(c("--bbs", "-w"), type = "character",
                help = "Biobank simulation software"),
    make_option(c("--prsice", "-p"), type = "character",
                help = "PRSice"),
    make_option(
        c("--out", "-o"),
        type = "character",
        help = "Output prefix",
        default = "Out"
    ),
    make_option(
        c("--quick"),
        action = "store_true",
        default = F,
        help = "Skip phenotype generation"
    )
)

argv <- parse_args(OptionParser(option_list = option_list))
if (is.na(argv$bfile) |
    is.na(argv$fam) | is.na(argv$snp) | is.na(argv$bbs) ) {
    stop("Please provide all required arguments")
}

num.causal <- 10000
heritability <- 0.5
print("Start generating phenotypes")
if (!argv$quick) {
    system2(
        argv$bbs,
        args = c(
            "--input",
            argv$bfile,
            # use the biobank genotype file
            "--out",
            argv$out,
            "--nsnp",
            num.causal,
            "--effect 2",
            # Effect sizes are simulated under normal distribution
            "--extract",
            argv$snp,
            "--keep",
            argv$fam,
            # samples to be included in the simulation
            "--herit",
            paste(heritability, collapse = ",")
        )
    )
}

dat <- fread(paste0(argv$out, ".xbeta"))
dat[, h_0.5:=scale(h_0.5)]
dat[,Environment:=rnorm(.N)]
setnames(dat, "h_0.5", "Genetics")
herit <- 0.4
env <- 0.3
inter <- 0.2
dat[,Interaction:=Genetics*Environment]
dat[,Genetics:=Genetics*sqrt(herit)]
dat[,Environment:=Environment*sqrt(env)]
dat[,Interaction:=Interaction*sqrt(inter)]
dat[,Noise:=rnorm(.N, sd=sqrt(1-var(Genetics)-var(Environment)-var(Interaction)))]
dat[,Phenotype:=Genetics+Environment+Interaction+Noise]
fwrite(dat[,c("FID","IID", "Phenotype")], paste0(argv$out, ".pheno"), sep="\t")
# Then perform plink
system2(
    "plink",
    args = c(
        "-assoc",
        "--bfile ",
        argv$bfile,
        "--extract ",
        argv$snp,
        "--out ",
        argv$out,
        "--pheno ",
        paste0(argv$out, ".pheno")
    )
)

bim <- fread(paste0(argv$bfile, ".bim"))
bim <- bim[,c("V2","V5","V6")]
setnames(bim, c("V2","V5","V6"), c("SNP", "A1", "A2"))
assoc <- fread(paste0(argv$out, ".qassoc"))
assoc <- merge(assoc, bim)
fwrite(assoc, paste0(argv$out, ".assoc"), sep = "\t")
system2(
    argv$prsice,
    args = c(
        "--base",
        paste0(argv$out, ".assoc"),
        "--target",
        argv$bfile,
        "--pheno ",
        paste0(argv$out, ".pheno"),
        "--out ", argv$out
    )
)

# Read in PRS
prs <- fread(paste0(argv$out, ".best"))
info <- merge(dat, prs)
# Real interaction
summary(lm(Phenotype~Environment+PRS+PRS*Environment))
# No real interaction
summary(lm(Phenotype~Noise+PRS+PRS*Noise))
