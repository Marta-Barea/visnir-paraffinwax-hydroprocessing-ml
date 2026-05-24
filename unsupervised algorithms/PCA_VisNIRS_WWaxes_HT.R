# Load Required Libraries
library(readxl)
library(prospectr)
library(cluster)
library(factoextra)
library(data.table)
library(viridis)

figures_dir <- file.path(getwd(), "Figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

save_stacked_plots <- function(plots, filename, width = 12, height = 12, res = 600) {
  png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = res,
    bg = "white",
    type = "cairo-png"
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(length(plots), 1)))

  for (i in seq_along(plots)) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
  }

  grDevices::dev.off()
}

# Load Parallelization when available
cl <- NULL
if (requireNamespace("doParallel", quietly = TRUE)) {
  cl <- parallel::makePSOCKcluster(8)
  doParallel::registerDoParallel(cl)
}

# Load data
pw_data <- read_excel("~/Documents/Doctorado/Tesis Doctoral/Investigación Cepsa/Vis-NIR/XDS-NIR_FOSS/Estudio según Tipo de Parafina e Hidrotratamiento/NIRS_HT_PW.xlsx",
                      sheet = "HT_Class")

# First derivative and Savitzky Golay Smoothing
pw_data$Sample <- as.factor(pw_data$Sample)

sgvec <- savitzkyGolay(X = pw_data[,-1], p = 3, w = 11, m = 1)
pw_sg <- cbind.data.frame(Sample = pw_data$Sample, sgvec)

pw_mw <- as.data.frame(pw_sg)
pw_mw$Sample = as.numeric(as.factor(pw_sg$Sample))

# PCA
pw_pca <- prcomp(pw_mw[,-1], scale = FALSE)

# Visualizing PCA results
pw_pca
summary(pw_pca)

# Visualizing eigenvalues (scree plot)
fviz_eig(pw_pca,
         xlab = "Principal Components (PCs)",
         ylab = "Explained Variance (%)",
         main = "",
         addlabels = TRUE,
         ggtheme = theme_minimal(),
         barcolor = "#404788FF",
         barfill = "#404788FF",
         linecolor = "#000000")

# Score plot for PC1 and PC2
scores_pca <- cbind.data.frame(predict(pw_pca),
                               waxes = pw_data$Sample)

m_labels <- as.matrix(pw_data[,-1])
rownames(m_labels) <- pw_data$Sample

# Score plot for PC1 and PC2
scatter_plot <- ggplot(scores_pca, aes(x = PC1, y = PC2, col = waxes, label = rownames(m_labels))) +
  geom_hline(yintercept = 0, lty = "dashed", alpha = 0.3) +
  geom_vline(xintercept = 0, lty = "dashed", alpha = 0.3) +
  geom_point(alpha = 0.1, size = 6, color = "black") +  
  geom_point(alpha = 1, size = 2, shape = 16) +      
  guides(color = guide_legend(title = "Hydroprocessing Grade")) + 
  labs(x = "PC1 (55.0%)", y = "PC2 (23.5%)", title = "") + 
  theme(axis.title = element_text(size = 12), 
        legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 12),
      plot.tag = element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 0.99)) +
  theme_test() +
  scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE) +
    labs(tag = "A")

scatter_plot

# Load plot
loadings <- cbind.data.frame(pw_pca$rotation[,c(1,2)])
setDT(loadings, keep.rownames = TRUE)[]

ld <- melt(loadings, "rn")

loadings_plot <- ggplot(ld, aes(x = rn, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_test()+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, hjust = 1, angle = 90),
    axis.title = element_text(size = 8),
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99)) +
  labs(x = "Wavelength (nm)", y = "Loadings PCs", title = "") +
  scale_x_discrete(limits = loadings$rn,
                   breaks = loadings$rn[seq(1, length(loadings$rn), by = 100)]) +
  scale_fill_manual(values = c("#FDE725FF", "#440154FF")) +
  geom_hline(yintercept = c(0.05, -0.05), linetype = "dotted") +
  labs(tag = "B")
  
loadings_plot

# Combining plots
save_stacked_plots(
  plots = list(scatter_plot, loadings_plot),
  filename = file.path(figures_dir, "Fig.3.png")
)

# Stop Parallelization
if (!is.null(cl)) {
  parallel::stopCluster(cl)
}