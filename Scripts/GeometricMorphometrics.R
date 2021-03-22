# Script for geometric morphometrics of ant head shapes. 
  # This script requires tps formatted landmark data for input. See the data folder for an example of what this looks like. 

# Load the packages we will use:
library(geomorph)
library(tidyverse)

# All of the steps of the analysis are wrapped up in a single function, which reads in the raw landmark data, superimposes it for standardization, and performs a PCA on the standardized coordinates. 
geometricMorphometricPCA <- function(dataFile, pointColor, figureTitle) {
  data <- geomorph::readland.tps(file = dataFile, specID = "ID")
  
  # Superimpose the raw coordinate data:
  superimposition <- geomorph::gpagen(data,  Proj = TRUE, ProcD = TRUE, curves = NULL, surfaces = NULL)
  
  # Extract out the coords values from the object returned by gpagen:
  coordinates <- geomorph.data.frame(superimposition)
  
  # Convert it from an array to a matrix:
  coordinates2 <- matrix(coordinates$coords, nrow=dim(coordinates$coords)[3], byrow=TRUE)
  # Set the rownames of that matrix:
  rownames(coordinates2) <- dimnames(coordinates$coords)[[3]]
  # Convert the matrix to a dataframe:
  coordinates3 <- data.frame(coordinates2)
  
  # Run the PCA:
  regularPCA <- prcomp(coordinates3)
  
  # Plot the PCA:
  PCAplot <- regularPCA %>%
    broom::augment(coordinates3) %>% # add original dataset back in
    ggplot(aes(x = .fittedPC1, y = .fittedPC2)) + 
    geom_point(size = 5, color = pointColor) +
    theme_half_open(12) + background_grid()
  
  plot(PCAplot)
  
  regularPCA %>%
    tidy(matrix = "rotation")
  # define arrow style for plotting
  arrow_style <- arrow(
    angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
  )
  # plot rotation matrix
  rotationMatrix <- regularPCA %>%
    tidy(matrix = "rotation") %>%
    pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
    ggplot(aes(PC1, PC2)) +
    geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
    geom_text(
      aes(label = column),
      hjust = 1, nudge_x = -0.02, 
      color = "#904C2F", size = 3
    ) +
    xlim(-1.25, .5) + ylim(-.5, 1) +
    theme_minimal_grid(12)
  
  plot(rotationMatrix)
  
  # How much variance is explained by each pc:
  varianceValues <- regularPCA %>%
    tidy(matrix = "eigenvalues")
  variancePlot <- regularPCA %>%
    tidy(matrix = "eigenvalues") %>%
    ggplot(aes(x = PC, y = percent)) +
    geom_col(fill = "#56B4E9", alpha = 0.8) +
    scale_x_continuous(breaks = 1:9) +
    scale_y_continuous(
      labels = scales::percent_format(),
      expand = expansion(mult = c(0, 0.01)), limits = c(0, 1)
    ) +
    theme_minimal_hgrid(12)
  
  plot(variancePlot)
  
  xLabel <- paste("PC1, explains ", as.character(varianceValues$percent[1] * 100), "% of variance", sep = "")
  yLabel <- paste("PC2, explains ", as.character(varianceValues$percent[2] * 100), "% of variance", sep = "")

  
  PCAplotFinal <- PCAplot + labs(x = xLabel, y = yLabel)
  plot(PCAplotFinal)
  
  
  allPlots <- ggarrange(PCAplotFinal, ggarrange(rotationMatrix, variancePlot, ncol = 1, nrow = 2), ncol = 2, nrow = 1, widths = c(2, 1))
  allPlots <- annotate_figure(allPlots, top = text_grob(figureTitle, size = 14))
  plot(allPlots)
}

# To run this function, supply the data filename; the color you want your points to be; and a text string for the figure title. 
AntHeadPCA <- geometricMorphometricPCA(dataFile = "./Data/ExampleLandmarks.txt", pointColor = "#F4B266", figureTitle = "Ant Head Shape")
plot(AntHeadPCA)
ggsave(filename = "AntHeadPCA.png", device = "png", path = "./Plots/PCAs/", width = 16, height = 9, bg = "transparent")
