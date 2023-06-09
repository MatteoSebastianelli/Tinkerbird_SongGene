#We want to obtain this:
# Chr_position order    ID  Ancestry
# 1_122579053 7312 11535 NA
# 2_158780 7313 11535 HET

fileList <- list.files(pattern="_ave_anc_pos.")

fileList

#Run in loop
for (i in 1:length(fileList)) {
  d <- read.table(fileList[i], fill = TRUE)
      data<-data.frame(paste(d[,1], d[,2], sep = "_"),
                       rownames(d),
                       rep(basename(getwd()), max(seq.int(nrow(d)))),
                       d[,3])
      head(data)
      names(data)<-NULL
      d1<-ifelse(data[,4] <= 0.5,"3",
                 ifelse(data[,4] >= 1.5,"1","2")) #1 = HOMR, 2 = HET, 3 = HOMY
      head(d1)
      final<-cbind(data[-c(4)],d1)
      head(final)
      names(final)<-NULL
      write.table(final, paste("Final", print(fileList)[i], sep = "_"),
                  col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

