# This script will run principle component analyses on linear measurements, for example head length, head width, Weber's length, etc. 

# First load the packages you will use:
library(tidyverse)
library(broom)
library(cowplot)
library(scales)
library(ggpubr)
library(rlang)

# Get all of the measurements needed:
  # Read in data. Data should be stored in a comma separated values file, placed in a subdirectory called "Data". Either name your data file "Measurements.csv", or alter the filename below to match yours:
  rawMeasurements <- read_csv(file = "./Data/CephalotesMeasurements.csv")
  
  # OPTIONAL: if you have numeric columns that are NOT your measurements, you'll need to drop them from your data tibble. For example, you might have columns that record latitude and longitude; you'll need to drop those from the tibble used in the PCAs. This could be done by:
     Measurements <- subset(rawMeasurements, select = -c(Latitude, Longitude))
  
  # You also need to drop any specimens that have missing data:
  MeasurementsNoNA <- Measurements
  MeasurementsNoNA <- na.omit(MeasurementsNoNA)
  
# Now you can produce scatterplots of your data.
  # This function will produce three scatterplots arranged togther, plotting 1) head width vs. head length; 2) Weber's length vs. pronotal width; and 3) head width/ Weber's length vs. head length/ Weber's length. You can change this function to match your own data by changing the x and y inputs, or just write your own plotting code. 
  linearScatterplots <- function(data, color, shape, figureTitle) {
    headWidthVLength <- ggplot(data) + 
      geom_point(mapping = aes(x = headwidth, 
                               y = headlength, 
                               color = {{color}}, 
                               shape = {{shape}}))
    
    webersLengthVPronotum <- ggplot(data) + 
      geom_point(mapping = aes(x = weberslength, 
                               y = pronotalwidth, 
                               color = {{color}}, 
                               shape = {{shape}}))
    
    headScaled <- ggplot(data) + 
      geom_point(mapping = aes(x = headwidth/weberslength, 
                               y = headlength/weberslength, 
                               color = {{color}}, 
                               shape = {{shape}}))
    
    allPlots <- ggarrange(headWidthVLength, webersLengthVPronotum, headScaled, ncol = 3, nrow = 1)
    allPlots <- annotate_figure(allPlots, top = text_grob(figureTitle, size = 14))
    plot(allPlots)
  }
  
  # Run the function by supplying the name of your processed measurements tibble; a column that you want to use to color the points, for example caste or species, and a string for the figure title. 
  linearScatterplots(data = MeasurementsNoNA, color = Caste, figureTitle = 'Caste Morphology in Examined Species')
  
# You will also likely want to run and plot principal components analysis (PCAs) on your data. 
  # This function runs a PCA and plots the results:
  morphologyPCA <- function(data, ellipse, color, colorVector, shape, figureTitle){
    if (ellipse == TRUE) {
      # Remove non-numeric columns and scale to variance, then perform PCA:
      pca_fit <- data %>% 
        select(where(is.numeric)) %>% 
        prcomp(scale = TRUE) 
      # Plot the PCA:
      PCAplot <- pca_fit %>%
        broom::augment(data) %>% # add original dataset back in
        ggplot(aes(.fittedPC1, .fittedPC2, color = {{color}}, shape = {{shape}})) + 
        geom_point(size = 5) +
        scale_color_manual(
          values = colorVector
        ) +
        stat_ellipse() +
        theme_half_open(12) + background_grid()
      # Plot the rotation matrix:
      pca_fit %>%
        tidy(matrix = "rotation")
      # define arrow style for plotting
      arrow_style <- arrow(
        angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
      )
      # plot rotation matrix
      rotationMatrix <- pca_fit %>%
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
      # How much variance is explained by each pc:
      pca_fit %>%
        tidy(matrix = "eigenvalues")
      variancePlot <- pca_fit %>%
        tidy(matrix = "eigenvalues") %>%
        ggplot(aes(PC, percent)) +
        geom_col(fill = "#56B4E9", alpha = 0.8) +
        scale_x_continuous(breaks = 1:9) +
        scale_y_continuous(
          labels = scales::percent_format(),
          expand = expansion(mult = c(0, 0.01))
        ) +
        theme_minimal_hgrid(12)
      allPlots <- ggarrange(PCAplot, rotationMatrix, variancePlot, ncol = 3, nrow = 1, widths = c(3.5, 2, 2))
      allPlots <- annotate_figure(allPlots, top = text_grob(figureTitle, size = 14))
      plot(allPlots)} 
    else (ellipse == FALSE) 
    {
      pca_fit <- data %>% 
        select(where(is.numeric)) %>% 
        prcomp(scale = TRUE) 
      # Plot the PCA:
      PCAplot <- pca_fit %>%
        broom::augment(data) %>% # add original dataset back in
        ggplot(aes(x = .fittedPC1, y = .fittedPC2, color = {{color}}, shape = {{shape}})) + 
        geom_point(size = 5) + xlim(-3, 3) + ylim(-2, 2) +
        scale_color_manual(
          values = colorVector
        ) +
        theme_half_open(12) + background_grid()
      
      plot(PCAplot)
      
      # Plot the rotation matrix:
      pca_fit %>%
        tidy(matrix = "rotation")
      # define arrow style for plotting
      arrow_style <- arrow(
        angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
      )
      # plot rotation matrix
      rotationMatrix <- pca_fit %>%
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
      varianceValues <- pca_fit %>%
        tidy(matrix = "eigenvalues")
      variancePlot <- pca_fit %>%
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
  }
  
  # To run this function: 
    # 1) supply the tibble containing your data; 
    # 2) decide if you want the function to draw ellipses around each group of points (this may or may not be useful to you); 
    # 3) decide which column in your tibble you want the points to be colored based upon (maybe the points should be colored by species, or by caste). Omit this if you do not want to color the points; 
    # 4) supply a vector that lists the color value corresponding to each value in the columns you're coloring on (for example, minor = "blue, major = "red). Omit this if you do not want to color the points; 
    # 5) decide which column in your tibble you want the points to be shaped based upon (maybe the points should be shaped by species, or by caste). Omit this if you do not want to use different shapes;
    # 6) provide a string for the figure title. 
  
  # For example:
    # Define the color scheme:
    colorScheme <- c(`Cephalotes varians` = "#65A484", `Cephalotes atratus` = "#a87a45", `Cephalotes placidus` = "#A48EAD", `Cephalotes minutus` = "#F4B266", `Cephalotes opacus` = "#C2DEBA", `Cephalotes hamulus` = "#D079DF")
  
    # Run the PCA:
    PCA <- morphologyPCA(data = MeasurementsNoNA, ellipse = FALSE, color = Species, colorVector = colorScheme, shape = Polymorphism, figureTitle = "PCA of Cephalotes morphology, raw measurements")
  
    # You can save this figure with ggsave:
    ggsave(filename = "PCA.png", device = "png", path = "./Plots/PCAs/", width = 8, height = 6, bg = "transparent")
  
    
  






