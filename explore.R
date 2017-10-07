library(trena)
library(trenaViz)
library(colorspace)
#library(annotate)
library(org.Hs.eg.db)
library(MotifDb)
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
stopifnot(packageVersion("trena")    >= "0.99.182")
stopifnot(packageVersion("trenaViz") >= "0.99.21")
stopifnot(packageVersion("MotifDb")  >= "1.19.13")

if(!exists("trena"))
   trena <- Trena("hg38")

PORT.RANGE <- 8000:8020
if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=TRUE)
   setGenome(tv, "hg38")
   }

colors <- rainbow_hcl(15)
next.color <- 1

if(!exists("motifMatcher"))
   motifMatcher <- MotifMatcher("hg38", pfms=as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")))



# tbl.snps is my version of the following table, lifted over to hg38, pval < 0.01:
#    scz2.snp.results.txt.gz (header row and 9,444,230 SNPs)
#    hg19chrc hg19 chromosome as character string (chr1-chr22, chrX)
#    snpid rs ID of SNP
#    a1 reference allele for OR (may not be minor allele)    "ref"
#    a2 alternate allele                                     "variant"
#    bp hg19 base pair position of SNP
#    info imputation quality score
#    or odds ratio in PGC GWAS data
#    se standard error of ln(OR) in PGC GWAS data
#    p p-value in PGC GWAS data
#    ngt number of samples in which SNP directly genotyped   "count"

if(!exists("tbl.snps")){
   load("gwas/tbl.15.9192.hg38.RData")
   if(!any(grepl("score", colnames(tbl.snps))))
      tbl.snps$score <- -log10(tbl.snps$pval)
   tbl.snps$end <- tbl.snps$start   # makeup of my earlier misjudgement
   }

genes.of.interest <- c("FURIN", "TSNARE1", "CNTN4", "CLCN3", "SNAP91", "SOX2", "FOXG1", "SREBF1")
#----------------------------------------------------------------------------------------------------
setupDisplay <- function(geneSymbol)
{
   #showGenomicRegion(tv, geneSymbol)
   showGenomicRegion(tv, "chr15:90,851,129-90,884,862")
   Sys.sleep(5) # wait for igv to reset the region before querying it


} # setupDisplay
#----------------------------------------------------------------------------------------------------
addSnpsAsBedTrack <- function()
{
   x <- parseChromLocString(getGenomicRegion(tv))
   tbl.bedGraph <- subset(tbl.snps, chrom==x$chrom & start >= x$start & end <= x$end)
   stopifnot(nrow(tbl.bedGraph) > 0)
   tbl.bedGraph <- tbl.bedGraph[, c("chrom", "start", "end", "id", "score")]

   trackName <- "snps"
  # if(trackName %in% getTrackNames(tv))
  #    removeTracksByName(tv, trackName)

   addBedTrackFromDataFrame(tv, trackName, tbl.bedGraph, color="red",
                                 #minValue=min(tbl.bedGraph$score), maxValue=max(tbl.bedGraph$score),
                                 displayMode="EXPANDED")

} # addSnpsAsBedTrack
#----------------------------------------------------------------------------------------------------
addSnpsAsBedGraphTrack <- function()
{
   x <- parseChromLocString(getGenomicRegion(tv))
   tbl.bedGraph <- subset(tbl.snps, chrom==x$chrom & start >= x$start & end <= x$end)
   stopifnot(nrow(tbl.bedGraph) > 0)
   tbl.bedGraph <- tbl.bedGraph[, c("chrom", "start", "end", "score", "id")]

   trackName <- "gwas"
   #if(trackName %in% getTrackNames(tv))
   #   removeTracksByName(tv, trackName)

   min.score <- 0
   max.score <- max(tbl.bedGraph$score)

   addBedGraphTrackFromDataFrame(tv, trackName, tbl.bedGraph, color="red",
                                 minValue=min.score, maxValue=max.score,
                                 displayMode="EXPANDED")

} # addSnpsAsBedGraphTrack
#----------------------------------------------------------------------------------------------------
addAtac <- function()
{
   atac.urls <- list(S1="http://trena.systemsbiology.net/scz/GSM2199916_S1_Neuronal_Nuclei.bw",
                     S2="http://trena.systemsbiology.net/scz/GSM2199917_S2_Neuronal_Nuclei.bw",
                     S3="http://trena.systemsbiology.net/scz/GSM2199918_S3_Neuronal_Nuclei.bw",
                     S4="http://trena.systemsbiology.net/scz/GSM2199919_S4_Neuronal_Nuclei.bw",
                     S5="http://trena.systemsbiology.net/scz/GSM2199921_S6_Neuronal_Nuclei.bw",
                     S6="http://trena.systemsbiology.net/scz/GSM2199922_S7_Neuronal_Nuclei.bw",
                     S7="http://trena.systemsbiology.net/scz/GSM2199923_S8_Neuronal_Nuclei.bw")

   for(name in names(atac.urls)){
      url <- atac.urls[[name]]
      addBedTrackFromHostedFile(tv, name, uri=url, color=colors[next.color])
      next.color <- next.color + 1
      Sys.sleep(5)
      }

# in javascript, from the console
# igv.browser.loadTrack({url: 'https://s3.amazonaws.com/igv.broadinstitute.org/data/hg19/encode/wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig', name: 'H3k4me3'})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199916_S1_Neuronal_Nuclei.bw", name: "atac-S1", max:10})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199917_S2_Neuronal_Nuclei.bw", name: "atac-S2", max:10})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199918_S3_Neuronal_Nuclei.bw", name: "atac-S3", max:10})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199919_S4_Neuronal_Nuclei.bw", name: "atac-S4", max:10})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199921_S6_Neuronal_Nuclei.bw", name: "atac-S6", max:10})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199922_S7_Neuronal_Nuclei.bw", name: "atac-S7", max:10})
# igv.browser.loadTrack({url: "http://trena.systemsbiology.net/scz/GSM2199923_S8_Neuronal_Nuclei.bw", name: "atac-S8", max:10})


} # addAtac
#----------------------------------------------------------------------------------------------------
addDHSregions <- function()
{
   genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"   # has gtf and motifsgenes tables
   chromLocString <- getGenomicRegion(tv)

   dhsFilter <- HumanDHSFilter(genome="hg38",
                               encodeTableName="wgEncodeRegDnaseClustered",
                               pwmMatchPercentageThreshold=95L,
                               geneInfoDatabase.uri=genome.db.uri,
                               regionsSpec=chromLocString,
                               quiet=TRUE)
    tbl.dhs <- getCandidates(dhsFilter)

} # addDHSregions
#----------------------------------------------------------------------------------------------------
addDHS <- function()
{
   sources <- list(DHS="encodeHumanDHS")

   current.region <- parseChromLocString(getGenomicRegion(tv))
   target.gene <- "FURIN"
   target.gene.tss <- 90868599

   printf("span: %d", 1 + current.region$end - current.region$start)
   tbls <- getRegulatoryChromosomalRegions(trena,
                                           current.region$chrom, current.region$start, current.region$end,
                                           sources, target.gene, target.gene.tss)

   track.count <- length(tbls)

   for(i in seq_len(length(tbls))){
      track.name <- names(sources)[i]
      columns.of.interest <- c("chrom", "motifStart", "motifEnd", "motifName")
      tbl <- tbls[[i]]
      stopifnot(all(columns.of.interest %in% colnames(tbl)))
      addBedTrackFromDataFrame(tv, track.name, tbl[, columns.of.interest], color=colors[next.color]);
      next.color <- next.color + 1
      } # for i

   tbl.dhs <- do.call(rbind, tbls)
   rownames(tbl.dhs) <- NULL
   tbl.dhs

} # addDHS
#----------------------------------------------------------------------------------------------------
addFootprints <- function()
{
   database.host <- "bddsrds.globusgenomics.org"
   dbName <- "brain_hint"
   brain.hint.db.uri <- sprintf("postgres://%s/%s", database.host, dbName)

   sources <- list(brainHint20=brain.hint.db.uri)

   current.region <- parseChromLocString(getGenomicRegion(tv))
   target.gene <- "FURIN"
   target.gene.tss <- 90868599

   printf("span: %d", 1 + current.region$end - current.region$start)
   tbls <- getRegulatoryChromosomalRegions(trena,
                                           current.region$chrom, current.region$start, current.region$end,
                                           sources, target.gene, target.gene.tss)

   track.count <- length(tbls)

   for(i in seq_len(length(tbls))){
      track.name <- names(sources)[i]
      columns.of.interest <- c("chrom", "motifStart", "motifEnd", "motifName")
      tbl <- tbls[[i]]
      stopifnot(all(columns.of.interest %in% colnames(tbl)))
      addBedTrackFromDataFrame(tv, track.name, tbl[, columns.of.interest], color=colors[next.color]);
      next.color <- next.color + 1
      } # for i

  # dhsFilter <- HumanDHSFilter(genome="hg38",
  #                             encodeTableName="wgEncodeRegDnaseClustered",
  #                             pwmMatchPercentageThreshold=95L,
  #                             geneInfoDatabase.uri=genome.db.uri,
  #                             regionsSpec=chromLocString,
  #                             quiet=TRUE)
  #  tbl.dhs <- getCandidates(dhsFilter)

  tbls

} # addFootprints
#----------------------------------------------------------------------------------------------------
makeModel <- function(trena, targetGene, targetGene.tss, tbl.motifs, mtx.rna, pcaMaxThreshold=1.0)
{
   printf("possible motifs for which we have expression: %d/%d",
          length(intersect(tbl.motifs$geneSymbol, rownames(mtx.rna))), length(unique(tbl.motifs$geneSymbol)))

  chromosome <- "chr11"

  #solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
  solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")

  tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs, mtx.rna)
  tbl.geneModel.strong <- subset(tbl.geneModel, pcaMax >= pcaMaxThreshold)

  tbl.regulatoryRegions.strong <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel.strong$gene)   # 71 rows

    # add two additional column to describe the regulatory regions, given that the TCF7 tss is 52283013
    #  motifName chrom motifStart motifEnd strand
    #   MA0516.1 chr11   52282407 52282421      -
    #
    #    1) the distance (+ or -) to the targetGene's TSS
    #    2) an id, e.g., "TCF7.fp.-.000606.MA0516.1"
    #----------------------------------------------------------------------------------------------------

  distance <- tbl.regulatoryRegions.strong$motifStart - targetGene.tss
  direction <- rep("upstream", length(distance))
  direction[which(distance < 0)] <- "downstream"
  tbl.regulatoryRegions.strong$distance.from.tss <- distance
  tbl.regulatoryRegions.strong$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.regulatoryRegions.strong$motifName)

  #save(tbl.geneModel, tbl.geneModel.strong, tbl.regulatoryRegions.strong, file="tbl.geneModel.8sep2017.9am.RData")

  tbl.tfBindingFreq <- as.data.frame(table(tbl.regulatoryRegions.strong$geneSymbol))
  tbl.tfBindingFreq <- tbl.tfBindingFreq[order(tbl.tfBindingFreq$Freq, decreasing=TRUE),]
  colnames(tbl.tfBindingFreq) <- c("gene", "binding.sites")
  tbl.geneModel.strong <- tbl.geneModel.strong[order(tbl.geneModel.strong$pcaMax, decreasing=TRUE),]
  list(model=tbl.geneModel.strong, regions=tbl.regulatoryRegions.strong)

} # makeModel
#----------------------------------------------------------------------------------------------------
run <- function()
{
   setupDisplay(genes.of.interest[2])
   #showGenomicRegion(tv, "chr15:90,866,725-90,876,696")
   #showGenomicRegion(tv, "chr15:90,880,000-90,884,330")
   #showGenomicRegion(tv, "chr15:90,868,140-90,869,746")   # just two high-scoring gwas snps straddling tss
   #showGenomicRegion(tv, "chr15:90,870,000-90,876,000")


     # furin gene model with transcript support level 1, tss: chr15:90,868,581
   showGenomicRegion(tv, "chr15:90,866,000-90,869,000")

   addSnpsAsBedTrack() # selection in browser displays  rsid labels
   addSnpsAsBedGraphTrack() # height nicely a function of score, labels lost
   #addAtac()
   tbls.fp <- addFootprints()
   tbl.h20 <- tbls.fp[[1]]
   tbl.h20$shortMotif <- tbl.h20$motifName
   tbl.h20.tf <- associateTranscriptionFactors(MotifDb, tbl.h20, source="TFClass", expand.rows=TRUE)

   tbl.dhs <- addDHS()
   tbl.dhs.tf <- associateTranscriptionFactors(MotifDb, tbl.dhs, source="MotifDb", expand.rows=TRUE)

     # get all motifs in the currently displayed regions
   tbl.regions <- with(parseChromLocString(getGenomicRegion(tv)),
                       data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE))

   tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, 85)
   tbl.motifs.tf <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="MotifDb")
   addBedTrackFromDataFrame(tv, "all.motifs.85",
                            tbl.motifs[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore")],
                            color="magenta")
   targetGene <- "FURIN"
   targetGene.tss <- 90868599
   print(load("~/github/trenaHelpers/inst/demos/aqp4/labMeeting-22jun2017/mayo.rnaSeq.cer.and.tcx.matrices.RData"))
   fivenum(mtx.cer)
   fivenum(mtx.tcx)


   model <- makeModel(trena, targetGene, targetGene.tss, tbl.motifs.tf, mtx.tcx, pcaMaxThreshold=0.5)
   colnames(model$model) <- gsub(".", "_", colnames(model$model), fixed=TRUE)
   models <- list(allMotifs_3k=model, same=model)

   g <- buildMultiModelGraph(tv, targetGene=targetGene, models) # comment out for speed
   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   setGraph(tv, g.lo, names(models))
   setStyle(tv, "style.js")
   fit(tv)

   candidate.tfs <- unique(tbl.motifs.tf$geneSymbol)

   targetGene <- "FURIN"
   #candidate.tfs <- intersect(tbl.h20.tf$geneSymbol, rownames(mtx.tcx))  # 440
   #candidate.tfs <- intersect(tbl.dhs.tf$geneSymbol, rownames(mtx.tcx))  # 78

   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")

   #tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.h20.tf, mtx.tcx)
   #tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.h20.tf, mtx.tcx)
   #tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs.tf, mtx.tcx)
   #models <- list(allMotifs=




   current.region <- parseChromLocString(getGenomicRegion(tv))
   current.snps <- subset(tbl.snps, chrom==current.region$chrom &
                                    start >= current.region$start &
                                    end <= current.region$end)$id   # may wish to also filter on count & score
   printf("snps in currently display genome region: %d", length(current.snps))

   tbls.snps <- list()
   for(snp in current.snps){
      printf("--- %s", snp)
      tbl.assay <- assessSnp(trena, as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")),
                             variant=snp, shoulder=10, pwmMatchMinimumAsPercentage=80)
      tbls.snps[[snp]] <- tbl.assay
      print(tbl.assay)
      }

   tbl.snpAssessment <- do.call(rbind, tbls.snps)
   tbl.snpAssessment <- subset(tbl.snpAssessment, assessed != "in.both")
   tbl.snpAssessment <- associateTranscriptionFactors(MotifDb, tbl.snpAssessment, source="MotifDb")
   rownames(tbl.snpAssessment) <- NULL

   tfs.added <- sort(unique(subset(tbl.snpAssessment, assessed=="mut.only")$geneSymbol))
   tfs.lost <- sort(unique(subset(tbl.snpAssessment,  assessed=="wt.only")$geneSymbol))

   strong.model.tfs <- subset(tbl.geneModel, pcaMax > 0.7)$gene
   tfs.added.also.in.model <- intersect(strong.model.tfs, tfs.added)
   tfs.lost.also.in.model  <- intersect(strong.model.tfs, tfs.lost)

   length(intersect(strong.model.tfs, tfs.added))  # 3
   length(intersect(strong.model.tfs, tfs.lost))

   all.snps <-  c("rs34688434", "rs77003790", "rs11853441", "rs11855381", "rs74869799", "rs7171165",
                  "rs7166599", "rs4932370", "rs4932371", "rs7168951", "rs7167085", "rs4932177",
                  "rs4932372", "rs7497418", "chr15_91409409_D", "rs8029440", "rs3759929", "rs4932178",
                  "chr15_91412848_I", "rs17514846", "rs8032315", "rs8027450", "chr15_91419432_I",
                  "rs2071410", "rs1573643", "rs1573644", "rs8039305", "rs6224", "rs6226", "rs6227",
                  "rs4702", "rs12906125", "rs35346340", "rs2071382", "rs11539637", "chr15_91428521_D",
                  "rs7183988", "rs7177338", "rs1894400", "rs1894401", "rs7497304", "chr15_91429196_D",
                  "rs4932373", "rs2071384")

   close.snps <- c("rs17514846", "rs8032315",  "rs8027450")


} # run
#----------------------------------------------------------------------------------------------------
rs8029440 <- function()
{
   rsid <- "rs8029440"
   subset(tbl.snps, id==rsid)
   rsid.explicit <- "chr15:90866284:A:G"
   rsid.explicit.2 <- "chr15:90866284:G:A"   # same as rsid
   assessSnp(trena, as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")), variant=rsid, shoulder=10, pwmMatch=85)
   assessSnp(trena, as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")), variant=rsid.explicit, shoulder=10, pwmMatch=85)
   assessSnp(trena, as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")), variant=rsid.explicit.2, shoulder=10, pwmMatch=85)

   # rsid adds YY1.  already 30 binding sites for yy1.  no strong independent evidence of open chromatin, tfs bound
   # Hsapiens-jaspar2016-YY1-MA0095.1    mut mut.only          0.9878049 -0.2073171   Hsapiens-jaspar2016-YY1-MA0095.1;90866281;+ chr15   90866281 90866286      + ACCATC chr15:90866284:G:A


} # rs8029440
#----------------------------------------------------------------------------------------------------
#chr15_91409409_D <- function
#{
#   rsid <- "chr15_91409409_D"   # variant: I14, treat this as maximally destructive
#   subset(tbl.snps, id==rsid)
#   tbl.regions <- data.frame(chrom="chr15", start=91409409-7, end=91409409+7, stringsAsFactors=FALSE)
#   findMatchesByChromosomalRegion(motifMatcher, tbl.regions, 85)   # another YY1
#    # 1 Hsapiens-jaspar2016-YY1-MA0095.1 chr15   91409411 91409416      -   4.176471          0.8658537 CATGGG   91409402 91409416 ATTGGTGCTCATGGG     wt
#
#} # chr15_91409409_D
##----------------------------------------------------------------------------------------------------
#rs4932178 <- function()
#{
#   rsid <- "rs4932178"
#   subset(tbl.snps, id==rsid)
#   assessSnp(trena, as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")), variant=rsid, shoulder=10, pwmMatch=85)
#   # adds second binding site for SREBF1
#   #                              motifName status assessed motifRelativeScore      delta                                       signature chrom motifStart motifEnd strand      match   variant
#   # 2  Hsapiens-jaspar2016-ZNF354C-MA0130.1     wt  wt.only          0.8888889  0.1975309 Hsapiens-jaspar2016-ZNF354C-MA0130.1;90868424;+ chr15   90868424 90868429      +     ACCCAC rs4932178
#   # 21  Hsapiens-jaspar2016-SREBF2-MA0596.1    mut mut.only          0.8720627 -0.1227154  Hsapiens-jaspar2016-SREBF2-MA0596.1;90868425;- chr15   90868425 90868434      - CTCACCCCAA rs4932178
#   # 11  Hsapiens-jaspar2016-SREBF1-MA0595.1    mut mut.only          0.9137931 -0.1250000  Hsapiens-jaspar2016-SREBF1-MA0595.1;90868425;+ chr15   90868425 90868434      + CTCACCCCAA rs4932178
#   # 1      Hsapiens-jaspar2016-SP1-MA0079.1     wt  wt.only          0.8750000  0.1071429     Hsapiens-jaspar2016-SP1-MA0079.1;90868424;- chr15   90868424 90868433      - ACCCACCCCA rs4932178
#   # 3     Hsapiens-jaspar2016-KLF5-MA0599.1     wt  wt.only          0.8862790  0.0851088    Hsapiens-jaspar2016-KLF5-MA0599.1;90868423;+ chr15   90868423 90868432      + GACCCACCCC rs4932178
#
#
#} # rs4932178
##----------------------------------------------------------------------------------------------------
#rs3759929 <- function ()
#{
#  rsid <- "rs3759929"
#  subset(tbl.snps, id==rsid)
#  assessSnp(trena, as.list(query(query(MotifDb, "hsapiens"), "jaspar2016")), variant=rsid, shoulder=10, pwmMatch=85)
#   # nothing of apparent interest here.
#   #                             motifName status assessed motifRelativeScore      delta                                     signature chrom motifStart motifEnd strand        match   variant
#   # 11 Hsapiens-jaspar2016-GATA3-MA0037.2    mut mut.only          0.8634889 -0.1365111 Hsapiens-jaspar2016-GATA3-MA0037.2;90866773;- chr15   90866773 90866780      -     TCTTCTCT rs3759929
#   # 1  Hsapiens-jaspar2016-FOXI1-MA0042.1     wt  in.both          0.9386973  0.0000000 Hsapiens-jaspar2016-FOXI1-MA0042.1;90866775;+ chr15   90866775 90866786      + TTCTGTTTGTTT rs3759929
#   # 21 Hsapiens-jaspar2016-FOXI1-MA0042.1    mut  in.both          0.8697318  0.0000000 Hsapiens-jaspar2016-FOXI1-MA0042.1;90866775;+ chr15   90866775 90866786      + TTCTCTTTGTTT rs3759929
#   # 4  Hsapiens-jaspar2016-FOXA1-MA0148.2     wt  wt.only          0.8510463  0.1041871 Hsapiens-jaspar2016-FOXA1-MA0148.2;90866778;+ chr15   90866778 90866788      +  TGTTTGTTTTT rs3759929
#   # 3  Hsapiens-jaspar2016-FOXA1-MA0148.1     wt  wt.only          0.8514336  0.1042185 Hsapiens-jaspar2016-FOXA1-MA0148.1;90866778;+ chr15   90866778 90866788      +  TGTTTGTTTTT rs3759929
#
#
#} # rs3759929
##----------------------------------------------------------------------------------------------------
#
