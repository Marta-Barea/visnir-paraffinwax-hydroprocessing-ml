# Load Required Libraries
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(prospectr)
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

# Apply Multiplicative Scatter Correction (MSC)
pw_data_msc <- pw_data
pw_data_msc[,-1] <- msc(pw_data[,-1])

pw_mean <- aggregate(.~ Sample, pw_data_msc, mean)

# Original NIR Spectra
df <- reshape2::melt(pw_mean, "Sample")

spectra_plot <- ggplot(data = df, aes(x = variable, y = value, color = Sample, group = Sample)) + 
  geom_line() +
  labs(x = "Wavelength (nm)", y = "Absorbance") +
  theme_test() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, hjust = 1, angle = 90),
      axis.title = element_text(size = 8),
      plot.tag = element_text(face = "bold", size = 16),
      plot.tag.position = c(0.01, 0.99)) +
  scale_x_discrete(limits = df$variable,
                   breaks = df$variable[seq(1, length(df$variable), by = 300)]) +
  scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE) +
    labs(tag = "A")

spectra_plot

# First derivative spectra
pw_data_msc$Sample <- as.factor(pw_data_msc$Sample)

sgvec <- savitzkyGolay(X = pw_data_msc[,-1], p = 3, w = 11, m = 1)
pw_sg <- cbind.data.frame(Sample = pw_data_msc$Sample, sgvec)

pw_sg_mean <- aggregate(.~ Sample, pw_sg, mean)

df_2 <- reshape2::melt(pw_sg_mean, "Sample")

savitzky_plot <- ggplot(data = df_2, aes(x = variable, y = value, color = Sample, group = Sample)) + 
  geom_line() +
  labs(x = "Wavelength (nm)", y = "Relative Signal") +
  theme_test() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 8),
        axis.text = element_text(size = 8, hjust = 1, angle = 90),
        axis.title = element_text(size = 8),
        plot.tag = element_text(face = "bold", size = 16),
        plot.tag.position = c(0.01, 0.99)) +
  scale_x_discrete(limits = df_2$variable,
                   breaks = df_2$variable[seq(1, length(df_2$variable), by = 300)]) +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE) +
  labs(tag = "B")

savitzky_plot

# Combining plots
save_stacked_plots(
  plots = list(spectra_plot, savitzky_plot),
  filename = file.path(figures_dir, "Fig.1.png")
)

# Stop Parallelization
if (!is.null(cl)) {
  parallel::stopCluster(cl)
}