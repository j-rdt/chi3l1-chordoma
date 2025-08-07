counts<-read.csv("ns2023/norm-counts.csv")
row_names<-counts[,1]
counts<-counts[-1]

pca_results <- mypca(counts, 
                     center = TRUE, 
                     scale = TRUE)
#rm(t_variance_stabilised_counts)

scores <- pca_results$scores

lbls<-c("Recurrent", "Primary", 
        "Primary", "Primary", 
        "Primary", "Recurrent", "Recurrent")

scores_with_conditions<-data.frame(lbls, scores)

#scores_with_conditions<-cbind(lbls, scores_with_conditions)
# explained variance
# one % variance value per PC
explained_variance <- 
  pca_results$explained_var %>% 
  pull("exp_var")

p<-ggplot(scores_with_conditions, 
          aes(PC1, PC2, color = lbls)) +
  geom_point(size = 4) +
  geom_text(label=row_names,
            #vjust=c(-0.2, 0.1, 0.1, 0.15, 0.15), 
            hjust=c(-0.1, 1.15, -0.12, 1.15, 1.15, -0.1, -0.1))+
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  theme(legend.title = element_blank()) + 
  ggtitle("Principal Component Analysis: Primary vs Recurrent")

ggsave("chi3l1-ns-pca.jpg", p, device="jpeg", width=5, height=3, units="in", dpi=600)


de<-read.csv("ns2023/DE-results.csv")
p<-EnhancedVolcano(toptable = de,
                   x = "Log2.fold.change",
                   y = "P.value",
                   lab = str_remove(de$X, "-mRNA"),
                   xlim = c(-5, +5),
                   ylim = c(0,4),
                   pCutoff = 5e-2,
                   pointSize = 2.0,
                   labSize=2.5,
                   FCcutoff = 1.5, 
                   title = NULL, 
                   subtitle = NULL, 
                   caption=NULL,
                   legendPosition = "bottom"
)
degenes<-str_remove(de$X[de$P.value<=0.05], "-mRNA")
plot(p)
ggsave("chi3l1-ns-volcano.jpg", p, device="jpeg", width=5, height=5, units="in", dpi=600)

counts<-read.csv("ns2023/norm-counts.csv")
row_names<-counts[,1]
counts<-counts[-1]
counts<-t(counts)
colnames(counts)<-row_names

counts_scaled <- 
  counts %>% 
  t(.) %>%                              # transpose to have the genes in columns 
  scale() %>%                           # scale(x, center = TRUE, scale = TRUE) 
  t(.)

genenames<-str_remove(row.names(counts_scaled), "\\.mRNA")
row.names(counts_scaled)<-genenames
ph<-pheatmap(counts_scaled[genenames %in% degenes,], 
             cluster_rows = TRUE,                      
             cluster_cols = TRUE, 
             show_rownames = TRUE, 
             show_colnames = TRUE,
             fontsize_row = 4,
             main = "Clustering by differentially expressed genes")
save_pheatmap_jpeg <- function(x, filename, width=5, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  jpeg(filename, width=width, height=height, quality=100, units="in", res=300)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_jpeg(ph, "chi3l1-ns-heat.jpeg")



