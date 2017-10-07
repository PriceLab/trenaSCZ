library(trena)
library(trenaViz)
library(colorspace)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("trena"))
   trena <- Trena("hg38")


PORT.RANGE <- 8000:8020

if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=TRUE)
   setGenome(tv, "hg38")
   }

if(!exists("tbl.loci"))
   load("gwas/tbl.loci.RData")

if(!exists("tbl.credSnp"))
   load("tbl.credSnp.RData")

if(!exists("tbl.snps")){
   load("gwas/tbl.15.9192.hg38.RData")
   if(!any(grepl("score", colnames(tbl.snps))))
      tbl.snps$score <- -log10(tbl.snps$pval)
   tbl.snps$end <- tbl.snps$start   # makeup of my earlier misjudgement
   }

if(!exists("tbl.wg.trena"))
   tbl.wg.trena <- readRDS("mayo.tcx.rds")

colors <- rainbow_hcl(15)
next.color <- 1

#------------------------------------------------------------------------------------------------------------------------
prepareLoci <- function()
{
   tbl.loci <- read.table("gwas/108loci.tsv", sep="\t", header=TRUE, quote="", fill=TRUE)

   parsed.locs <- lapply(tbl.loci$Position..hg19., parseChromLocString)
   tbl.locs <- do.call(rbind, lapply(parsed.locs, function(loc) data.frame(chrom=loc$chrom, start19=loc$start, end19=loc$end, stringsAsFactors=FALSE)))
   tbl.loci <- cbind(tbl.locs, tbl.loci)


   # liftover the hg19 chromLocStrings
   write.table(tbl.loci$Position..hg19., sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file="hg19-locs.txt")
   # browse to https://genome.ucsc.edu/cgi-bin/hgLiftOver
   # get implicit download link from "View Conversions"

   hg38.locs <- read.table("hglft_genome_5e7_6971d0.bed", sep="\t", as.is=TRUE, stringsAsFactors=FALSE)$V1
   length(hg38.locs)
   parsed.locs <- lapply(hg38.locs, parseChromLocString)
   tbl.hg38 <- do.call(rbind, lapply(parsed.locs, function(loc) data.frame(chrom=loc$chrom, start19=loc$start, end19=loc$end, stringsAsFactors=FALSE)))
   colnames(tbl.hg38) <- c("chrom", "start38", "end38")
      # make sure the liftover left all the locs in order
   stopifnot(nrow(tbl.hg38) == nrow(tbl.loci))
   stopifnot(all(tbl.loci$chrom == tbl.hg38$chrom))

   tbl.loci <- cbind(tbl.hg38, tbl.loci)
   dup.chrom.column <- grep("chrom", colnames(tbl.loci))[2]
   tbl.loci <- tbl.loci[, -dup.chrom.column]
   size <- with(tbl.loci, 1+end38-start38)
   tbl.loci$size <- size
   tbl.loci <- tbl.loci[, c("chrom", "start38", "end38", "size", "Rank", "P.value", "SCZ", "Protein.coding.genes")]
   save(tbl.loci, file="gwas/tbl.loci.RData")

} # prepareLoci
#------------------------------------------------------------------------------------------------------------------------
prepareSnps <- function()
{
   f <- "gwas/pgc.scz52.credible.snps.txt"
   tbl <- read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE)
   tbl$chromLoc <- sprintf("chr%s:%d-%d", tbl$CHR, tbl$indexSNP_POS, tbl$indexSNP_POS)
   write.table(tbl$chromLoc, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file="snps19.text")
   snps38 <- read.table("snps38.text", sep="\t", as.is=TRUE)$V1
   parsed.locs <- lapply(snps38, parseChromLocString)
   tbl.locs <- do.call(rbind, lapply(parsed.locs, function(loc) data.frame(chrom=loc$chrom, start38=loc$start, end38=loc$end, stringsAsFactors=FALSE)))
   tbl$pos38 <- tbl.locs$start38
   tbl$chrom <- sprintf("chr%s", tbl$CHR)

   #tbl <- tbl[, c("indexSNP", "chrom", "pos38"
   colnames(tbl) <- c("id", "chrom", "pos38", "credibleSNP",  "crediblePval", "credibleProb", "credibleCumProb")
   tbl.credSnp <- tbl
   save(tbl.credSnp, file="tbl.credSnp.RData")

} # prepareSnps
#------------------------------------------------------------------------------------------------------------------------
addFootprintsCurrentRegion <- function(target.gene, target.gene.tss)
{
   #database.host <- "bddsrds.globusgenomics.org"
   database.host <- "whovian"
   dbNames <- c("brain_hint_16", "brain_hint_20", "brain_wellington_20", "brain_wellington_16")
   brain.hint16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[1])
   brain.hint20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[2])
   brain.wellington16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[3])
   brain.wellington20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[4])

   sources <- list(brain.hint16.db.uri, brain.hint20.db.uri, brain.wellington16.db.uri, brain.wellington20.db.uri)
   names(sources) <- dbNames

   current.region <- parseChromLocString(getGenomicRegion(tv))

   printf("span: %d", 1 + current.region$end - current.region$start)
   tbls <- getRegulatoryChromosomalRegions(trena,
                                           current.region$chrom, current.region$start, current.region$end,
                                           sources, target.gene, target.gene.tss)

   track.count <- length(tbls)

   for(i in seq_len(length(tbls))){
      track.name <- names(sources)[i]
      if(track.name %in% getTrackNames(tv))
         removeTracksByName(tv, track.name)
      columns.of.interest <- c("chrom", "motifStart", "motifEnd", "motifName")
      tbl <- tbls[[i]]
      stopifnot(all(columns.of.interest %in% colnames(tbl)))
      addBedTrackFromDataFrame(tv, track.name, tbl[, columns.of.interest], color=colors[next.color]);
      next.color <- next.color + 1
      } # for i

  tbls

} # addFootprintsCurrentRegion
#----------------------------------------------------------------------------------------------------
createCollapsedGRanges <- function()
{
   database.host <- "whovian"
   dbNames <- c("brain_hint_16", "brain_hint_20", "brain_wellington_20", "brain_wellington_16")
   brain.hint16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[1])
   brain.hint20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[2])
   brain.wellington16.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[3])
   brain.wellington20.db.uri <- sprintf("postgres://%s/%s", database.host, dbNames[4])

   sources <- list(brain.hint16.db.uri, brain.hint20.db.uri, brain.wellington16.db.uri, brain.wellington20.db.uri)
   names(sources) <- dbNames

   current.region <- parseChromLocString(getGenomicRegion(tv))

   printf("span: %d", 1 + current.region$end - current.region$start)
   tbls <- getRegulatoryChromosomalRegions(trena,
                                           current.region$chrom, current.region$start, current.region$end,
                                           sources, target.gene, target.gene.tss)
     # make a GRangesList out of the 4 GRanges, converting each data.frame to a GRange object
   grl <- GRangesList(lapply(tbls, GRanges))
   gr.all <- Reduce(c, grl)
   strand(gr.all) <- "*"  # ignore strand differences
   start(gr.all) <- start(gr.all) - 10
   end(gr.all) <- end(gr.all) + 10

   tbl.fpCollapsed <- as.data.frame(reduce(gr.all))  # metadata discarded
   trackName <- "fpCollapsed"
   if(trackName %in% getTrackNames(tv))
      removeTracksByName(tv, trackName)

   addBedTrackFromDataFrame(tv, trackName, tbl.fpCollapsed, color="green")

     # find the imputed snp, collapsed fp overlaps
   colnames(tmp.bed)[1] <- "seqnames"
   #tmp.bed$seqnames <- sprintf("chr%s", tmp.bed$seqnames)

   tmp.bed$seqnames <- "chr8"
   tbl.overlaps <- as.data.frame(findOverlaps(GRanges(tmp.bed), GRanges(tbl.fpCollapsed)))
   colnames(tbl.overlaps) <- c("snp", "fp")
   trackName <- "snpInFp"
   if(trackName %in% getTrackNames(tv))
      removeTracksByName(tv, trackName)

   tbl.snpsInFP <<- tmp.bed[tbl.overlaps$snp,]
   addBedTrackFromDataFrame(tv, trackName, tmp.bed[tbl.overlaps$snp,], color="magenta")




} # createCollapsedGRanges
#----------------------------------------------------------------------------------------------------
loci.5 <- function(locusNumber, display=FALSE)
{
    #                                             5
    # chrom                                    chr8
    # start38                             142228142
    # end38                               142249172
    # start19                             143309503
    # end19                               143330533
    # Rank                                        5
    # P.value                             1.737e-15
    # Position..hg19.      chr8:143309503-143330533
    # SCZ                                         Y
    # Protein.coding.genes                  TSNARE1
    # OMIM
    # NHGRI.GWAS.catalog        SCZ pmid = 23974872
    # KO.phenotype
    # size                                    21031

   locusNumber <- 5
   locus.chrom <- tbl.loci$chrom[locusNumber]
   locus.start <- tbl.loci$start38[locusNumber]
   locus.end   <- tbl.loci$end38[locusNumber]
   locusLoc <- sprintf("%s:%d-%d", locus.chrom, locus.start, locus.end)
   showGenomicRegion(tv, locusLoc)

   tbl.lociSnps <- with(tbl.loci[locusNumber,], subset(tbl.credSnp, chrom==chrom & pos38 >= start38 & pos38 <= end38))
   with(unique(tbl.lociSnps[, c("id", "chrom", "pos38")]),
      addBedTrackFromDataFrame(tv, "CoreSNP",
                               data.frame(chrom=chrom, start=pos38, end=pos38, stringsAsFactors=FALSE), color="red"))

   imputed.snps <- grep("^rs", tbl.lociSnps$credibleSNP, value=TRUE)
     # takes a long time: several minutes on first lookup, 30 seconds on subsequent

   date()
   tbl.imputed.snps <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP150.GRCh38, imputed.snps))
   date()
   tmp.bed <- tbl.imputed.snps[, c("seqnames", "pos", "pos", "RefSNP_id")]
   colnames(tmp.bed) <- c("chrom", "start", "end", "rsid")
   tmp.bed$chrom <- as.character(tmp.bed$chrom)
   addBedTrackFromDataFrame(tv, "imputed", tmp.bed, color="purple")


   tbls.fp <- addFootprintsCurrentRegion("TSNARE1", 142354831)

   # intersecting snps, from visual inspection:
   # rs12676511:
   #   brain_wellington_20:
} # snpsInLocus
#------------------------------------------------------------------------------------------------------------------------
test.snpsInLocus <- function()
{
   snpsInLocus(5)

} # test.snpsInLocus
#------------------------------------------------------------------------------------------------------------------------
bigModel <- function()
{
   showGenomicRegion(tv, "chr8:142,199,164-142,392,333")
   target.gene <- "TSNARE1"

   tbls.fp <- addFootprintsCurrentRegion("TSNARE1", 142354831)
   tbl.fp <- do.call(rbind, tbls.fp)
   rownames(tbl.fp) <- NULL
   tbl.fp$shortMotif <- tbl.fp$motifName
   tbl.fp2 <- associateTranscriptionFactors(MotifDb, tbl.fp, source="TFClass", expand.rows=TRUE)
   rownames(tbl.fp2) <- NULL
   tbl.fp2 <- unique(tbl.fp2)

   print(load("~/github/trenaHelpers/inst/demos/aqp4/labMeeting-22jun2017/mayo.rnaSeq.cer.and.tcx.matrices.RData"))
   tbl.geneModel <- createGeneModel(trena, target.gene, solver.names, tbl.fp2, mtx.tcx)
   tbl.geneModel.2 <- createGeneModel(trena, target.gene, solver.names, tbl.fp2, mtx.cer)
   fivenum(mtx.cer)
   fivenum(mtx.tcx)


} # bigModel
#------------------------------------------------------------------------------------------------------------------------
explore.REST.motif.disruption <- function()
{
   # intersect the imputed snps with the collected footprints
   gr.fp <- with(tbl.fp2, GRanges(seqnames=chrom, IRanges(start=motifStart, end=motifEnd)))
   gr.iSnps <- with(tbl.imputed.snps, GRanges(seqnames=sprintf("chr%s", seqnames), IRanges(start=pos, end=pos)))
   tbl.ov <- as.data.frame(findOverlaps(gr.iSnps, gr.fp))
   colnames(tbl.ov) <- c("snp", "fp")
   tbl.imputed.snps[unique(tbl.ov$snp),]
   seqnames       pos strand  RefSNP_id alleles_as_ambig
      # 6         8 142245289      * rs12234969                W
      # 8         8 142232025      * rs12676511                R
      # 20        8 142232793      *  rs4976977                W
      # 24        8 142235164      *  rs4976982                R
      # 34        8 142239356      *  rs7465677                Y
      # 36        8 142241748      *  rs7822538                Y

   possibly.affected.tfs <- unique(tbl.fp2[tbl.ov$fp,"geneSymbol"])


     #  [1] "NFKB1"   "PPARG"   "ZNF263"  "NFAT5"
     #  [5] "NFATC1"  "NFATC2"  "NFATC3"  "NFATC4"
     #  [9] "STAT1"   "STAT2"   "STAT3"   "STAT4"
     # [13] "STAT5A"  "STAT5B"  "STAT6"   "MYB"
     # [17] "REST"    "TBX21"   "TBR1"    "EOMES"
     # [21] "CREM"    "CREBL2"  "CREBZF"  "CREB3"
     # [25] "CREB3L1" "CREB3L2" "CREB3L3" "CREB3L4"
     # [29] "ATF1"    "ATF3"    "JDP2"    "FOSB"

    intersect(possibly.affected.tfs, tbl.geneModel$gene) # [1] "ZNF263" "NFATC2" "STAT3"  "STAT5A" "REST"
    subset(tbl.geneModel, gene %in% possibly.affected.tfs)
     #    gene    beta.lasso pearson.coeff  rf.score   beta.ridge spearman.coeff concordance    pcaMax binding.sites
     #   STAT3 -0.0002332129    -0.3776042 3.6418365 -0.009046759     -0.3820749   0.3795418 1.0586716           124
     #  NFATC2  0.0000000000    -0.2686272 3.2686985 -0.008081293     -0.2866660   0.3570942 0.8824561            31
     #    REST  0.0000000000    -0.2723167 1.2828013 -0.028025665     -0.2709925   0.3095579 0.6467796            70
     #  STAT5A  0.0000000000    -0.2962635 0.2690173 -0.023344832     -0.2635510   0.3063804 0.5641500           124
     #  ZNF263  0.0000000000    -0.1665805 1.1343419 -0.028341127     -0.1837981   0.2882019 0.5408957           882

   tbl.snpFpTf <- cbind(tbl.imputed.snps[tbl.ov$snp,], tbl.fp2[tbl.ov$fp,])
   snps.in.fp <- unique(tbl.snpFpTf$RefSNP_id)
   tbls.snpAssessment <- list()
   pfms <- as.list(query(query(MotifDb, "jaspar2016"), "sapiens"))
   for(snp in snps.in.fp){
      tbl.snp <- assessSnp(trena, pfms,  variant=snp, shoulder=7, pwmMatchMinimum=85)
      tbls.snpAssessment[[snp]] <- tbl.snp
      } # for nsp
   tbl.snpAssessment <- do.call(rbind, tbls.snpAssessment)
   tbl.snpAssessment[grep("only", tbl.snpAssessment$assessed),]

} # explore.REST.motif.disruption
#------------------------------------------------------------------------------------------------------------------------
assess.snps <- function()
{
   dim(tbl.snpsInFP)

   targetGene <- "TSNARE1"

   jaspar2016.pfms <- query(MotifDb, "jaspar2016")
   human.pfms <- query(jaspar2016.pfms, "sapiens")
   mouse.pfms <- query(jaspar2016.pfms, "mus")
   rat.pfms   <- query(jaspar2016.pfms, "rnorvegicus")
      # jaspar.vertebrates <- query(query(MotifDb, "jaspar"), "vert"); length(jaspar.vertebrates)  # 17
   pfms <- as.list(c(human.pfms, mouse.pfms, rat.pfms))   # 619

   mm <- MotifMatcher("hg38", pfms)

   matchThreshold <- 85
   shoulder <- 12

   tbl.wg.trena.targetGene <- subset(tbl.wg.trena, target.gene==targetGene)
   tbl.wg.trena.targetGene <- tbl.wg.trena.targetGene[order(tbl.wg.trena.targetGene$pcaMax, decreasing=TRUE),]
   canonical.tfs <- tbl.wg.trena.targetGene$gene


   for(r in seq_len(nrow(tbl.snpsInFP))){
      rsid <- tbl.snpsInFP$rsid[r]
      printf("----- assessing %2d) %s", r, rsid)
      tbl.region <- with(tbl.snpsInFP[r,],
                         data.frame(chrom=seqnames, start=start-shoulder, end=end+shoulder, stringsAdFactors=FALSE))
      tbl.wt  <- findMatchesByChromosomalRegion(mm, tbl.region, matchThreshold)
      tfs.wt <- c()
      if(nrow(tbl.wt) > 0){
         tbl.wt$shortMotif <- unlist(lapply(tbl.wt$motifName,
                                         function(s) {tokens <- strsplit(s, "-")[[1]]; return(tokens[length(tokens)])}))

        tfs.wt <- c(tfs.wt, associateTranscriptionFactors(MotifDb, tbl.wt, source="TFClass", expand.rows=TRUE)$geneSymbol)
        tfs.wt <- unique(c(tfs.wt, associateTranscriptionFactors(MotifDb, tbl.wt, source="MotifDb", expand.rows=TRUE)$geneSymbol))
        if(any(is.na(tfs.wt)))
           tfs.wt <- tfs.wt[-which(is.na(tfs.wt))]
        } # nrow(tbl.wt) > 0

      tbl.mut <- findMatchesByChromosomalRegion(mm, tbl.region, matchThreshold, variants=rsid)
      tfs.mut <- c()
      if(nrow(tbl.mut) > 0){
         tbl.mut$shortMotif <- unlist(lapply(tbl.mut$motifName,
                                            function(s) {tokens <- strsplit(s, "-")[[1]]; return(tokens[length(tokens)])}))
         tfs.mut <- c(tfs.mut, associateTranscriptionFactors(MotifDb, tbl.mut, source="TFClass", expand.rows=TRUE)$geneSymbol)
         tfs.mut <- unique(c(tfs.mut, associateTranscriptionFactors(MotifDb, tbl.mut, source="MotifDb", expand.rows=TRUE)$geneSymbol))
         if(any(is.na(tfs.mut)))
            tfs.mut <- tfs.mut[-which(is.na(tfs.mut))]
         } # nrow(tbl.mut) > 0
      tfs.lost <- setdiff(tfs.wt, tfs.mut)
      tfs.gained <- setdiff(tfs.mut, tfs.wt)
      tfs.inModel.lost   <- intersect(tfs.lost, canonical.tfs)
      tfs.inModel.gained <- intersect(tfs.gained, canonical.tfs)
      printf("%s %s, tfs.wt,     %2d in model: %2d", targetGene, rsid, length(tfs.wt),     length(intersect(tfs.wt, canonical.tfs)))
      printf("%s %s, tfs.mut,    %2d in model: %2d", targetGene, rsid, length(tfs.mut),    length(intersect(tfs.mut, canonical.tfs)))
      printf("%s %s, tfs.lost,   %2d in model: %2d", targetGene, rsid, length(tfs.lost),   length(tfs.inModel.lost))
      printf("%s %s, tfs.gained, %2d in model: %2d", targetGene, rsid, length(tfs.gained), length(tfs.inModel.gained))
      if(length(tfs.inModel.lost) > 0){
         indices <- match(tfs.inModel.lost, canonical.tfs)
         for(index in indices) printf("   tf %2d lost: %s", index, canonical.tfs[index])
         }
      if(length(tfs.inModel.gained) > 0){
         indices <- match(tfs.inModel.gained, canonical.tfs)
         for(index in indices) printf("   tf %2d gained: %s", index, canonical.tfs[index])
         }
      } # for r


} # assess.snps
#------------------------------------------------------------------------------------------------------------------------
rs12676511 <- function()
{
   rsid <- "rs12676511"
   showGenomicRegion(tv, "chr8:142,232,013-142,232,042")
   tbl.thisSnp <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP150.GRCh38, rsid))
   pos <- tbl.thisSnp$pos
   tbl.thisSnpFP <- subset(tbl.fp, motifStart <= tbl.thisSnp$pos & motifEnd >= tbl.thisSnp$pos)
   jaspar2016.pfms <- query(MotifDb, "jaspar2016")
   human.pfms <- query(jaspar2016.pfms, "sapiens")
   mouse.pfms <- query(jaspar2016.pfms, "mus")
   rat.pfms   <- query(jaspar2016.pfms, "rnorvegicus")
   jaspar.vertebrates <- query(query(MotifDb, "jaspar"), "vert"); length(jaspar.vertebrates)
   pfms <- as.list(c(human.pfms, mouse.pfms, rat.pfms))
   mm <- MotifMatcher("hg38", pfms)
   shoulder <- 12
   findMatchesByChromosomalRegion(mm,
                                  tbl.regions=data.frame(chrom="chr8", start=pos-shoulder, end=pos+shoulder, stringsAsFactors=FALSE),
                                  90)


   assessSnp(trena, pfms, variant=rsid, shoulder=10, pwmMatchMinimum=70)

   genome.db.uri    <- "postgres://whovian/hg38"                  # has gtf and motifsgenes tables
   footprint.db.uri <- "postgres://whovian/brain_wellington_20"            # has hits and regions tables
   fpf <- FootprintFinder(genome.db.uri, footprint.db.uri, quiet=FALSE)
   getFootprintsInRegion(fpf, chrom="chr8", start=pos-10, end=pos+10)




} # rs12676511
#------------------------------------------------------------------------------------------------------------------------
# id <- "rs12676511"
# loc <- subset(tmp.bed, rsid==id)$start # [1] 142232025
# subset(tbls.fp[["postgres://whovian/brain_wellington_20"]], motifStart <= loc & motifEnd >= loc)
#    chrom motifStart  motifEnd motifName strand    score length distance.from.tss                                    id
# 8   chr8  142232005 142232025  MA0528.1      +  9.36735     21           -122826 TSNARE1.fp.downstream.122826.MA0528.1
# 16  chr8  142232019 142232028  MA0606.1      - 13.09680     10           -122812 TSNARE1.fp.downstream.122812.MA0606.1
# 17  chr8  142232021 142232027  MA0152.1      - 12.20220      7           -122810 TSNARE1.fp.downstream.122810.MA0152.1
# 18  chr8  142232021 142232035  MA0517.1      -  9.31818     15           -122810 TSNARE1.fp.downstream.122810.MA0517.1
#
#
# id <- "rs4976977"
# loc <- subset(tmp.bed, rsid==id)$start # 142232793
# subset(tbls.fp[["postgres://whovian/brain_hint_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_hint_20"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_20"]], motifStart <= loc & motifEnd >= loc)
# # all agree:
# # chr8  142232787 142232796  MA0100.2      - 14.0492     10           -122044 TSNARE1.fp.downstream.122044.MA0100.2
#
# id <- "rs4976982"
# loc <- subset(tmp.bed, rsid==id)$start  # 142235164
# subset(tbls.fp[["postgres://whovian/brain_hint_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_hint_20"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_20"]], motifStart <= loc & motifEnd >= loc)
# # all 4 sources agree:
# # chr8  142235154 142235174  MA0138.2      - 9.57143     21           -119677 TSNARE1.fp.downstream.119677.MA0138.2
#
# id <- "rs7465677"
# loc <- subset(tmp.bed, rsid==id)$start  # 142239356
# subset(tbls.fp[["postgres://whovian/brain_hint_16"]], motifStart <= loc & motifEnd >= loc)
# # chrom motifStart  motifEnd motifName strand   score length distance.from.tss                                    id
# #  chr8  142239352 142239361  MA0690.1      + 10.8727     10           -115479 TSNARE1.fp.downstream.115479.MA0690.1
# #  chr8  142239352 142239364  MA0800.1      + 10.9756     13           -115479 TSNARE1.fp.downstream.115479.MA0800.1
# subset(tbls.fp[["postgres://whovian/brain_hint_20"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_20"]], motifStart <= loc & motifEnd >= loc)
#
# id <- "rs7822538"
# loc <- subset(tmp.bed, rsid==id)$start  # 142241748
# subset(tbls.fp[["postgres://whovian/brain_hint_16"]], motifStart <= loc & motifEnd >= loc)
# #  chrom motifStart  motifEnd motifName strand   score length distance.from.tss                                    id
# #   chr8  142241740 142241749  MA0609.1      - 12.7705     10           -113091 TSNARE1.fp.downstream.113091.MA0609.1
# #   chr8  142241741 142241748  MA0604.1      - 13.0492      8           -113090 TSNARE1.fp.downstream.113090.MA0604.1
# #   chr8  142241742 142241749  MA0605.1      - 11.9180      8           -113089 TSNARE1.fp.downstream.113089.MA0605.1
#
# subset(tbls.fp[["postgres://whovian/brain_hint_20"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_20"]], motifStart <= loc & motifEnd >= loc)
#
#
# id <- "rs12234969"
# loc <- subset(tmp.bed, rsid==id)$start  # 142245289
# subset(tbls.fp[["postgres://whovian/brain_hint_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_hint_20"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_16"]], motifStart <= loc & motifEnd >= loc)
# subset(tbls.fp[["postgres://whovian/brain_wellington_20"]], motifStart <= loc & motifEnd >= loc)
# # all 4 agree:
# # chrom motifStart  motifEnd motifName strand   score length distance.from.tss                                    id
# #  chr8  142245279 142245291  MA0105.4      + 1.16327     13           -109552 TSNARE1.fp.downstream.109552.MA0105.4
#  chr8  142245289 142245308  MA0066.1      - 8.57778     20           -109542 TSNARE1.fp.downstream.109542.MA0066.1


# collected:
# fp.snp.motifs <- c("MA0066.1", "MA0609.1", "MA0604.1", "MA0605.1", "MA0690.1", "MA0800.1", "MA0100.2", "MA0138.2")

# motifToGene(MotifDb, fp.snp.motifs, source="TFclass")
#      motif geneSymbol pubmedID  source
# 1  MA0066.1      PPARG 23180794 TFClass
# 2  MA0100.2        MYB 23180794 TFClass
# 3  MA0138.2       REST 23180794 TFClass
# 4  MA0604.1       ATF1 23180794 TFClass
# 5  MA0604.1     CREBL2 23180794 TFClass
# 6  MA0604.1     CREBZF 23180794 TFClass
# 7  MA0604.1      CREB3 23180794 TFClass
# 8  MA0604.1    CREB3L1 23180794 TFClass
# 9  MA0604.1    CREB3L2 23180794 TFClass
# 10 MA0604.1    CREB3L3 23180794 TFClass
# 11 MA0604.1    CREB3L4 23180794 TFClass
# 12 MA0604.1       CREM 23180794 TFClass
# 13 MA0605.1       ATF3 23180794 TFClass
# 14 MA0605.1       JDP2 23180794 TFClass
# 15 MA0605.1       FOSB 23180794 TFClass
# 16 MA0609.1       CREM 23180794 TFClass
# 17 MA0609.1     CREBL2 23180794 TFClass
# 18 MA0609.1     CREBZF 23180794 TFClass
# 19 MA0609.1      CREB3 23180794 TFClass
# 20 MA0609.1    CREB3L1 23180794 TFClass
# 21 MA0609.1    CREB3L2 23180794 TFClass
# 22 MA0609.1    CREB3L3 23180794 TFClass
# 23 MA0609.1    CREB3L4 23180794 TFClass
# 24 MA0609.1       ATF1 23180794 TFClass
# 25 MA0690.1      TBX21 23180794 TFClass
# 26 MA0690.1       TBR1 23180794 TFClass
# 27 MA0800.1      EOMES 23180794 TFClass

# from MotifDb motif/tf mapping
# tf.map.motifDB <- lapply(fp.snp.motifs, function(motif) unique(mcols(query(MotifDb, motif))$geneSymbol))
# names(tf.map.motifDB) <- fp.snp.motifs
# tf.map.motifDB
# MA0066.1  [1] "PPARG"
# MA0609.1  [1] "Crem"
# MA0604.1  [1] "Atf1"
# MA0605.1  [1] "Atf3"
# MA0690.1  [1] "TBX21"
# MA0800.1  [1] "EOMES"
# MA0100.2  [1] "Myb" "brk"
# MA0138.2  [1] "REST"
