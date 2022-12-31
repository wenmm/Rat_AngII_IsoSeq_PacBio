#setwd("/")
PacBio_angii_control_tr_count_percent <- read.delim("PacBio_angii_control_tr_count_percent", header=TRUE)


condition <- c("control","angii")
rdata <- as_tibble(PacBio_angii_control_tr_count_percent)

for (group in unique(condition)) {
  percent <- paste0(group, '.percent')
  class <- paste0(group, '.class')
  max_rows <- rdata %>%
    group_by(.data$gname) %>%
    slice_max(!!as.name(percent), with_ties = FALSE)
  rdata[[class]] <- ifelse(rdata[[percent]] < 0.1, 'Inactive', 'Minor')
  rdata[[class]][match(max_rows$tr, rdata$tr)] <- "Major"
  rdata[[class]][which(rdata[[percent]] < 0.1)] <- "Inactive"
}

rdata$fc <- log2(rdata$angii * 0.7942351/rdata$control)
rdata$diff <- rdata$angii.percent - rdata$control.percent
write.table(rdata,file = "PacBio_angii_control_tr_count_percent.txt",row.names = FALSE)