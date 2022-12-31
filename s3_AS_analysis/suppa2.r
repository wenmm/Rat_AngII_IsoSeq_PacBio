setwd("/mnt/data1/wenhm/project/rn7_mRatBN7.2/final_pacbio_transcript_structure/suppa2_new")
dpsi <- read.table(file="events_diffSplice.dpsi",sep="\t",header = TRUE)
colnames(dpsi) <- c("dPSI","p_value")

psi_events <- read.table(file="pacbio_events.psi",sep="\t")
colnames(psi_events) <- c("A1","A2","CTRL1","CTRL2")

event_TPM <- read.table(file="events_diffSplice_avglogtpm.tab",sep="\t",header = FALSE)
colnames(event_TPM) <- c("event","mean_TPM")

merge1 <- merge(dpsi,psi_events,by="row.names")
merge2 <- merge(merge1,event_TPM,by.x="Row.names",by.y="event")

rownames(merge2) <- merge2$Row.names
final_table <- merge2
final_table <- final_table[,-1]
final_table <- final_table[!is.nan(final_table$dPSI),]
final_table$cpval <- p.adjust(final_table$p_value, method = "bonferroni")
final_table$log10pval <- -log10(final_table$p_value)
final_table$sig <- "not sig"
final_table[final_table$p_value < 0.06,]$sig <- "sig"
# final_table[, c(2:13)] <- sapply(final_table[, c(2:13)], function(x) as.numeric(gsub(",", ".", x)))
final_table$logRNAc <-final_table$mean_TPM

final_table <- cbind(final_table,rownames(final_table))
colnames(final_table)[12] <- "Name"

final_table$mean_PSI <- (final_table$A1 + final_table$A2 + final_table$CTRL1 + final_table$CTRL2) / 4

final_table$a_diff <- abs(final_table$A1 - final_table$A2)
final_table$c_diff <- abs(final_table$CTRL1 - final_table$CTRL2)

select_final_table <- subset(final_table,a_diff < 0.1 & c_diff < 0.1 & abs(dPSI) >= 0.2)


p <- ggplot(final_table, aes(x=dPSI, y=log10pval, color=sig))
p + geom_point() + geom_vline(xintercept=c(-0.5,0.5), linetype="solid", size=1) +
  geom_vline(xintercept=c(-0.5,0.5), linetype="solid", size=1) +
  geom_hline(yintercept=1.3, size=1) +
  xlab(expression(~Delta~PSI)) + ylab("-log10(p-value)") + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="none",
    legend.text=element_blank(),
    legend.title=element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)) +
  scale_color_manual(values=c("sig" = "darkblue", "not sig" = "gray40", "nan" = "gray80")) + 
  scale_x_continuous(breaks=pretty_breaks(n=5), limits = c(-1, 1)) +
  scale_y_continuous(breaks=pretty_breaks(n=5), limits = c(0, 4.5)) +
  guides(fill = guide_legend(reverse = FALSE))

write.table(final_table,file = "events_result.txt",col.names = TRUE, row.names = FALSE)



events_result <- read.csv("/mnt/data1/wenhm/project/rn7_mRatBN7.2/final_pacbio_transcript_structure/suppa2_new/events_result1.txt", sep="")
p <- ggplot(events_result, aes(x=dPSI, y=log10pval, color=sig1))
p + geom_point() + 
  geom_vline(xintercept=c(-0.1,0.1),  size=1) +
  geom_hline(yintercept=1.221849, size=1) +
  xlab(expression(~Delta~PSI)) + ylab("-log10(p-value)") + 
  theme_classic() +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    legend.position="none",
    legend.text=element_blank(),
    legend.title=element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16)) +
  theme(axis.ticks.length = unit(5, "pt"),axis.line = element_line(size = 1), axis.ticks = element_line(colour = "black", size = 1), panel.grid.major = element_line(colour = NA, size = 1), panel.grid.minor = element_line(colour = NA, size = 1)) +
  scale_color_manual(values=c( "not sig" = "gray40", "nan" = "gray80","A3_sig" = "green","A5_sig" = "#b38f79","AF_sig"= "darkblue","AL_sig" = "red" ,"MX_sig" = "#dc5939","RI_sig" = "#b86795","SE_sig" = "#ffbf01")) + 
  scale_x_continuous(breaks=pretty_breaks(n=5), limits = c(-1, 1)) +
  scale_y_continuous(breaks=pretty_breaks(n=5), limits = c(0, 4.5)) +
  guides(fill = guide_legend(reverse = FALSE))


tmp <- read.table("events_result1_sig.gname", quote="\"", comment.char="")
de <- bitr(tmp$V1, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Rn.eg.db)
eg = de$ENTREZID
b = enrichKEGG(eg,keyType = "kegg", pvalueCutoff = 0.05,organism = "rno")
barplot(b,showCategory=20)
ego <- enrichGO(eg , 'org.Rn.eg.db', ont = 'BP')
barplot(ego)


GO = ego@result

GO_for_plot = subset(GO, qvalue< 0.05)[1:10,]
GO_for_plot$Description <- factor(GO_for_plot$Description,levels = rev(GO_for_plot$Description))
bar_up <- ggplot(GO_for_plot, aes(y = Description, x = -log(pvalue,10))) +
  geom_bar(width = 1, color = "white",size = 2,stat="identity",fill="#B7D5A3") +
  theme_classic() +
  labs(x = "-Log10(qvalue)",title = "GO Enrichment Analysis") +
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5)) +
  #scale_x_continuous(position = "top",breaks = seq(0, 30, 5)) +
  theme(axis.title.x = element_text(size = 16, face = "bold",vjust = 0.9, hjust = 0.5),axis.title.y = element_text(size = 16, face = "bold",vjust = 0.9, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 16,  face = "bold"), axis.text.y = element_text(size = 16,  face = "bold")) +
  theme(axis.ticks.length = unit(7, "pt"),axis.line = element_line(size = 1), axis.ticks = element_line(colour = "black", size = 1))


write.table(b@result, file="gname_kegg.txt")
write.table(ego@result, file="gname_GO.txt")
write.table(de, file="gene.txt",row.names = FALSE)
