library(trena)
library(trenaViz)
library(colorspace)
library(SNPlocs.Hsapiens.dbSNP150.GRCh38)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
# an exploration of data reported in:
#    "Biological Insights From 108 Schizophrenia-Associated Genetic Loci", Ripke et al, Nature 2014
#     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4112379/
#
#  from the summary:
#     "We identify 128 independent associations spanning 108 conservatively defined loci that meet genome-wide
#      significance, 83 of which have not been previously reported."
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "TSNARE1"
targetGene.tss <- 142354831

#database.host <- "whovian"
database.host <- "bddsrds.globusgenomics.org"

#------------------------------------------------------------------------------------------------------------------------
if(!exists("trena"))
   trena <- Trena("hg38")

PORT.RANGE <- 8000:8020

if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=TRUE)
   setGenome(tv, "hg38")
   }

   # tbl.loci: reformatted, with liftover hg38 coordinatess, for
   # Supplementary Table 3: Bioinformatic summary data for 108 genome-wide significant loci
   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4112379/bin/NIHMS59304-supplement-Supplementary_Information.docx

if(!exists("tbl.loci"))
   load("gwas/tbl.loci.RData")    # 108       14

   # tbl.credSnp: 20,373 "credible SNPs" associated with 266 "index SNPs"
   # lifted-over to hg38, obtained from pgc.scz52.credible.snps.txt
   # extracted from pgc.scz2.credible.SNPs.zip
   # downloaded from http://www.med.unc.edu/pgc/files/resultfiles/pgc.scz2.credible.SNPs.zip/view

if(!exists("tbl.credSnp"))
   load("tbl.credSnp.RData")      # 20373     7

# pshannnon (10 oct 2017):  i think this is no longer used, it is an artifact of my first
# exploration of the FURIN gene
# if(!exists("tbl.snps")){
#    load("gwas/tbl.15.9192.hg38.RData")
#    if(!any(grepl("score", colnames(tbl.snps))))
#       tbl.snps$score <- -log10(tbl.snps$pval)
#    tbl.snps$end <- tbl.snps$start   # makeup of my earlier misjudgement
#    }

# cory funk's whole genome trn
if(!exists("tbl.wg.trena"))
   tbl.wg.trena <- readRDS("mayo.tcx.rds")

# our pre-computed TSNARE1, 250kb, MotifDb motifs & TF mapping
if(!exists("tbl.geneModel.mtx.cer.TSNARE1.250kb"))
  load("tbl.geneModel.mtx.cer.TSNARE1.250kb.RDATA")

colors <- rainbow_hcl(15)
next.color <- 1

#------------------------------------------------------------------------------------------------------------------------
# query trenaViz (in the browser) for the currently displayed genomic region
# then use the trena convenience method "getRegulatoryChromosomalRegions" to get and
# display brain DHS footprints reported in that region
#
# targetGene and targetGene.tss are used so that the resulting footprints can be
# annotated with respect to the (presumed) gene of interest
#
addFootprintsInCurrentRegion <- function(targetGene, targetGene.tss)
{
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
                                           sources, targetGene, targetGene.tss)

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

} # addFootprintsInCurrentRegion
#----------------------------------------------------------------------------------------------------
# collapse all footprints from multiple sources into one non-redundant, non-overlapping track
# this may be a useful visualization, as well as the basis for an exhaustive search for binding
# motifs, using trena's MotifMatcher class
#
# the "expandFootprintsBy" parameter will lead to the unification of more individual footprints
# this will sometimes be useful:  calling footprints is an inexact art; allowing for extra
# width may help in generating hypotheses - implicating snps - which would otherwise be missed.
createCollapsedFootprintTrack <- function(expandFootprintsBy=10)
{
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
                                           sources, targetGene, targetGene.tss)
     # make a GRangesList out of the 4 GRanges, converting each data.frame to a GRange object
   grl <- GRangesList(lapply(tbls, GRanges))
   gr.all <- Reduce(c, grl)
   strand(gr.all) <- "*"  # ignore strand differences
   start(gr.all) <- start(gr.all) - expandFootprintsBy
   end(gr.all) <- end(gr.all) + expandFootprintsBy

   tbl.fpCollapsed <- as.data.frame(reduce(gr.all))  # metadata discarded
   trackName <- "fpCollapsed"
   if(trackName %in% getTrackNames(tv))
      removeTracksByName(tv, trackName)

   addBedTrackFromDataFrame(tv, trackName, tbl.fpCollapsed, color="green")

     # find the imputed snp, collapsed fp overlaps
   #colnames(tmp.bed)[1] <- "seqnames"
   #tmp.bed$seqnames <- sprintf("chr%s", tmp.bed$seqnames)

   #tmp.bed$seqnames <- "chr8"
   #tbl.overlaps <- as.data.frame(findOverlaps(GRanges(tmp.bed), GRanges(tbl.fpCollapsed)))
   #colnames(tbl.overlaps) <- c("snp", "fp")
   #trackName <- "snpInFp"
   #if(trackName %in% getTrackNames(tv))
   #   removeTracksByName(tv, trackName)

   #tbl.snpsInFP <<- tmp.bed[tbl.overlaps$snp,]
   #addBedTrackFromDataFrame(tv, trackName, tmp.bed[tbl.overlaps$snp,], color="magenta")

   invisible(tbl.fpCollapsed)

} # createCollapsedFootprintTrack
#----------------------------------------------------------------------------------------------------
gwas.locus.5 <- function()
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

       # extract the location of locus 5, display in trenaViz
   locusNumber <- 5
   locus.chrom <- tbl.loci$chrom[locusNumber]
   locus.start <- tbl.loci$start38[locusNumber]
   locus.end   <- tbl.loci$end38[locusNumber]
   locusLoc <- sprintf("%s:%d-%d", locus.chrom, locus.start, locus.end)
   showGenomicRegion(tv, locusLoc)

      # get all the credible snps from that region
   tbl.lociSnps <- with(tbl.loci[locusNumber,], subset(tbl.credSnp, chrom==chrom & pos38 >= start38 & pos38 <= end38))
   with(unique(tbl.lociSnps[, c("id", "chrom", "pos38")]),
      addBedTrackFromDataFrame(tv, "indexSNP",
                               data.frame(chrom=chrom, start=pos38, end=pos38, stringsAsFactors=FALSE), color="red"))

     # now get their hg38 coordinates.  for now, ignore all imputed snps without rsid
   imputed.snps <- grep("^rs", tbl.lociSnps$credibleSNP, value=TRUE)

   printf("retrieveing snp locations from Bioconductor package 'SNPlocs.Hsapiens.dbSNP150.GRch38'")
   printf("takes a few minutes on the first call, and ~30 seconds on subsequent calls...")

   tbl.imputed.snps <- as.data.frame(snpsById(SNPlocs.Hsapiens.dbSNP150.GRCh38, imputed.snps))

      # transform tbl.imputed.snps a bit for rendering in trenaViz: rearrange & rename columns
   tbl.imputedSnp.bed <- tbl.imputed.snps[, c("seqnames", "pos", "pos", "RefSNP_id")]
   colnames(tbl.imputedSnp.bed) <- c("chrom", "start", "end", "rsid")
   tbl.imputedSnp.bed$chrom <- as.character(tbl.imputedSnp.bed$chrom)
   tbl.imputedSnp.bed$chrom <- sprintf("chr%s", tbl.imputedSnp.bed$chrom)
   addBedTrackFromDataFrame(tv, "imputed", tbl.imputedSnp.bed, color="purple")

      # now add footprints found in this region
   tbls.fp <- addFootprintsInCurrentRegion(targetGene, targetGene.tss)

      # and a track showing those footprints, no padding, collapsed into one track
   tbl.fpCollapsed <- createCollapsedFootprintTrack(expandFootprintsBy=0)

      # intersect imputed snps with footprints
   tbl.overlaps <- as.data.frame(findOverlaps(GRanges(tbl.imputedSnp.bed), GRanges(tbl.fpCollapsed)))
   colnames(tbl.overlaps) <- c("snp", "fp")

   tbl.imputedSnpsInFootprints <- tbl.imputedSnp.bed[tbl.overlaps$snp,]
   addBedTrackFromDataFrame(tv, "imputedInFp", tbl.imputedSnpsInFootprints, color="darkRed")

   assess.snps(tbl.imputedSnpsInFootprints, matchThreshold=90, shoulder=8)

} # gwas.locus.5
#------------------------------------------------------------------------------------------------------------------------
# the indexSNP and the imputed snps in locus 5 - a 3' intron of TSNARE1 - suggests that regulatory regions
# for this gene may exist over a very broad region.
# accordingly, we here build a regulatory model for TSNARE1 over a region 250kb long.
# we make the further choice of only mapping motifs-in-footprints conservatively, using the MotifDb
# mappings
bigLocalModel.TSNARE1 <- function()
{
   showGenomicRegion(tv, "chr8:142,190,000-142,440,000")   # centered on TSNARE1, 250kb

     # get and display footprints from all four sources
   tbls.fp <- addFootprintsInCurrentRegion(targetGene, targetGene.tss)
     # combine them into one table
   tbl.fp <- do.call(rbind, tbls.fp)
   tbl.fp <- unique(tbl.fp)
   rownames(tbl.fp) <- NULL

   tbl.fp.mdb <- associateMotifsFromFootprintDatabasesWithMotifDbTFs(tbl.fp)
   rownames(tbl.fp.mdb) <- NULL
   tbl.fp.mdb <- unique(tbl.fp.mdb)

   print(load("~/github/trenaHelpers/inst/demos/aqp4/labMeeting-22jun2017/mayo.rnaSeq.cer.and.tcx.matrices.RData"))
   fivenum(mtx.cer)

   solver.names = c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")

   tbl.geneModel.2.mdb <- createGeneModel(trena, targetGene, solver.names, tbl.fp.mdb, mtx.cer)
   rownames(tbl.geneModel.2.mdb) <- NULL

   tbl.geneModel.mtx.cer.TSNARE1.250kb <- tbl.geneModel.2.mdb

   save(tbl.geneModel.mtx.cer.TSNARE1.250kb, file="tbl.geneModel.mtx.cer.TSNARE1.250kb.RDATA")

} # bigModel
#------------------------------------------------------------------------------------------------------------------------
associateMotifsFromFootprintDatabasesWithMotifDbTFs <- function(tbl)
{
   motifs <- tbl$motifName
   motifs <- gsub("{", "_", motifs, fixed=TRUE)
   motifs <- gsub("}", "_", motifs, fixed=TRUE)
   motifs.uniq <- unique(motifs)
   na.motifs <- grep("^NA", motifs.uniq)
   if(length(na.motifs) > 0)
      motifs.uniq <- motifs.uniq[-na.motifs]

   tfs <- lapply(motifs.uniq, function(motif) {mcols(MotifDb[motif])$geneSymbol})
   names(tfs) <- motifs.uniq
   tfs.all <- tfs[motifs]
   stopifnot(length(tfs.all) == nrow(tbl))
   geneSymbols <- as.character(tfs.all)
   empties <- which(geneSymbols == "character(0)")
   if(length(empties) > 0)
      geneSymbols[empties] <- ""
   tbl$geneSymbol <- geneSymbols

   return(tbl)

} # associateMotifsFromFootprintDatabasesWithMotifDbTFs
#------------------------------------------------------------------------------------------------------------------------
assess.snps <- function(tbl.imputedSnpsInFootprints, matchThreshold, shoulder)
{
   jaspar2016.pfms <- query(MotifDb, "jaspar2016")
   human.pfms <- query(jaspar2016.pfms, "sapiens")
   mouse.pfms <- query(jaspar2016.pfms, "mus")
   rat.pfms   <- query(jaspar2016.pfms, "rnorvegicus")
      # jaspar.vertebrates <- query(query(MotifDb, "jaspar"), "vert"); length(jaspar.vertebrates)  # 17
   pfms <- as.list(c(human.pfms, mouse.pfms, rat.pfms))   # 619

   mm <- MotifMatcher("hg38", pfms)

      # we have already read in a pre-computed whole genome trena model
      # let's complement that with a
   tbl.wg.trena.targetGene <- subset(tbl.wg.trena, target.gene==targetGene)
   tbl.wg.trena.targetGene <- tbl.wg.trena.targetGene[order(tbl.wg.trena.targetGene$pcaMax, decreasing=TRUE),]
   wholeGenomeModel.canonical.tfs <- tbl.wg.trena.targetGene$gene

   tsnare1.250kb.localModel.canonical.tfs <- tbl.geneModel.mtx.cer.TSNARE1.250kb$gene

   for(r in seq_len(nrow(tbl.imputedSnpsInFootprints))){
      rsid <- tbl.imputedSnpsInFootprints$rsid[r]
      printf("----- assessing %2d) %s", r, rsid)
      tbl.region <- with(tbl.imputedSnpsInFootprints[r,],
                         data.frame(chrom=chrom, start=start-shoulder, end=end+shoulder, stringsAdFactors=FALSE))
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
      tfs.inModel.lost   <- unique(c(intersect(tfs.lost, wholeGenomeModel.canonical.tfs),
                                     intersect(tfs.lost, tsnare1.250kb.localModel.canonical.tfs)))
      tfs.inModel.gained <- unique(c(intersect(tfs.gained, wholeGenomeModel.canonical.tfs),
                                     intersect(tfs.gained, tsnare1.250kb.localModel.canonical.tfs)))
      printf("%s %s, tfs.wt,     %2d, in model: %2d", targetGene, rsid, length(tfs.wt),     length(intersect(tfs.wt, wholeGenomeModel.canonical.tfs)))
      printf("%s %s, tfs.mut,    %2d, in model: %2d", targetGene, rsid, length(tfs.mut),    length(intersect(tfs.mut, wholeGenomeModel.canonical.tfs)))
      printf("%s %s, tfs.lost,   %2d, in model: %2d", targetGene, rsid, length(tfs.lost),   length(tfs.inModel.lost))
      printf("%s %s, tfs.gained, %2d, in model: %2d", targetGene, rsid, length(tfs.gained), length(tfs.inModel.gained))

      if(length(tfs.inModel.lost) > 0){
         wg.indices <- match(tfs.inModel.lost, wholeGenomeModel.canonical.tfs)
         if(any(is.na(wg.indices)))
            wg.indices <- wg.indices[-which(is.na(wg.indices))]
         for(index in wg.indices)
            printf("  tf %s rank %2d lost in wg model", wholeGenomeModel.canonical.tfs[index], index)
         tsnare1.250kb.indices <- match(tfs.inModel.lost, tsnare1.250kb.localModel.canonical.tfs)
         if(any(is.na(tsnare1.250kb.indices)))
            tsnare1.250kb.indices <- tsnare1.250kb.indices[-which(is.na(tsnare1.250kb.indices))]
         for(index in tsnare1.250kb.indices)
            printf("  tf %s rank %2d lost in tsnare1 model", tsnare1.250kb.localModel.canonical.tfs[index], index)
         }
      if(length(tfs.inModel.gained) > 0){
         wg.indices <- match(tfs.inModel.gained, wholeGenomeModel.canonical.tfs)
         if(any(is.na(wg.indices)))
            wg.indices <- wg.indices[-which(is.na(wg.indices))]
         for(index in wg.indices)
            printf("  tf %s rank %2d gained in wg model", wholeGenomeModel.canonical.tfs[index], index)
         tsnare1.250kb.indices <- match(tfs.inModel.gained, tsnare1.250kb.localModel.canonical.tfs)
         if(any(is.na(tsnare1.250kb.indices)))
            tsnare1.250kb.indices <- tsnare1.250kb.indices[-which(is.na(tsnare1.250kb.indices))]
         for(index in tsnare1.250kb.indices)
            printf("  tf %s rank %2d gained in tsnare1 model", tsnare1.250kb.localModel.canonical.tfs[index],  index)

         }
      } # for r

} # assess.snps
#------------------------------------------------------------------------------------------------------------------------
