
txdb <- makeTxDbFromGFF("combined_f1_f.cds_add_ncbiRefSeq.gtf")
promoterAnnotationData <- preparePromoterAnnotation(txdb, species = 'Rattus_norvegicus')


starJunctionFiles <- list.files(path="/mnt/data1/wenhm/project/rn7_mRatBN7.2/2rd_analysis",pattern=".pacbioAligned.sortedByCoord.out.pacbioSJ.out.tab", full.names = TRUE)


condition <- rep(c('AngII','Control'), each=2)

result <- proActiv(files = starJunctionFiles,promoterAnnotation = promoterAnnotationData,condition = condition)
assays(result)
rowData(result)
write.table(rowData(result) , file='promoter_activity.tab', sep="\t",row.names = TRUE)

geneExpression <- metadata(result)$geneExpression
promoterCounts <- assays(result)$promoterCounts
normalizedPromoterCounts <- assays(result)$normalizedPromoterCounts
absolutePromoterActivity <- assays(result)$absolutePromoterActivity
relativePromoterActivity <- assays(result)$relativePromoterActivity

#setwd("")
write.table(geneExpression , file="geneExpression.tab",sep="\t",row.names = TRUE)
write.table(promoterCounts , file="promoterCounts.tab",sep="\t",row.names = TRUE)
write.table(normalizedPromoterCounts , file="normalizedPromoterCounts.tab",sep="\t",row.names = TRUE)
write.table(absolutePromoterActivity , file="absolutePromoterActivity.tab",sep="\t",row.names = TRUE)
write.table(relativePromoterActivity , file="relativePromoterActivity.tab",sep="\t",row.names = TRUE)

result_filter <- result[complete.cases(assays(result)$promoterCounts),]
alternativePromoters <- getAlternativePromoters(result = result, referenceCondition = "AngII")

write.table(rdata , file="promoter_activity_modify.tab",sep="|",row.names = TRUE)

## alternative promoter
raw.assay <- assays(result)$relativePromoterActivity
rdata <- rowData(result)
condition <- result$condition
nonInternalId <- which(rdata$internalPromoter != 1)
pval <- rep(NaN, nrow(rdata))
#assay <- raw.assay[nonInternalId, ]
assay <- raw.assay
num.pros <- nrow(assay)
pval[nonInternalId] <- unlist(lapply(seq_len(num.pros), function(i) 
  tryCatch(summary(lm(unlist(assay[i,]) ~ condition))$coef[2,4], 
           warning = function(cond) 1,
           error = function(cond) NaN)))

referenceCondition = "Control"
id <- which(condition == referenceCondition)
mean.cond <- rowMeans(raw.assay[,id, drop = FALSE], na.rm = TRUE) 
mean.other <- rowMeans(raw.assay[,-id, drop = FALSE], na.rm = TRUE) 
gexp <- metadata(result)$geneExpression
gexp <- as.matrix(gexp[, result$sampleName])
rownames(gexp) <- rowData(result)$geneId
mean.gexp.cond <- rowMeans(gexp[,id, drop = FALSE], na.rm=TRUE)
mean.gexp.cond <- mean.gexp.cond[match(rdata$geneId, rownames(gexp))]
mean.gexp.other <- rowMeans(gexp[,-id, drop = FALSE], na.rm=TRUE)
mean.gexp.other <- mean.gexp.other[match(rdata$geneId, rownames(gexp))]

result <- data.frame(promoterId = rdata$promoterId,
                     geneId = rdata$geneId,
                     pval = pval,
                     padj = padj,
                     abs.cond = mean.cond,
                     abs.other = mean.other,
                     gexp.cond = mean.gexp.cond,
                     gexp.other = mean.gexp.other)


# Normalize promoter read counts by DESeq2 
## normalizedPromoterCounts.star <- normalizePromoterReadCounts(promoterCounts.star)

# Calculate absolute promoter activity
## absolutePromoterActivity.star <- getAbsolutePromoterActivity(normalizedPromoterCounts.star,promoterAnnotationData)


#DESeq2

## dds <- DESeqDataSetFromMatrix(countData=promoterCounts.star[activePromoters,], colData = coldata, design = ~ condition)
## dds <- DESeq(dds)
## res.GC.Normal <- results(dds, contrast=c("condition","Cancer","Normal"))
