
plot_sample <- function(pop, nx, ny, sample_result, show.grid = TRUE, quadrat.color = "#e31a1c", pointsize = 0.4, pointalpha = 0.8, quadratalpha = 0.8) {
  
  if(class(pop)[1] == "ppp") {
    
    pop_counts <- quadratcount(pop, nx = nx, ny = ny)
    
  } else if(class(pop)[1] == "quadratcount" | class(pop)[1] == "data.frame") {
    
    # check if the matrix is of dimensions nx*ny
    if(ncol(pop) != nx | nrow(pop) != ny){
      stop("If 'pop' is a of type 'quadratcount' it has to be of same dimensions as nx and ny.")
    }
    
    pop_counts <- pop
    
  } else {
    
    stop("Object given to argument 'pop' has to be of type 'ppp' or 'quadratcount'.")
    
  }
  
  # # get number of individuals per cell
  # pop_counts <- quadratcount(pop, nx = nx, ny = ny)
  
  # extract point coordinates
  pattern_df <- data.frame(x = pop$x, y = pop$y)
  
  # rescale point coordinates inside the nx*ny window
  pattern_df$x <- pattern_df$x / (pop$window$xrange[2] / nx)
  pattern_df$y <- pattern_df$y / (pop$window$yrange[2] / ny)
  
  # extract the matrix with the counts from the object of type 'quadratcount'
  pop_counts <- as.vector(t(pop_counts))
  pop_counts <- matrix(data = pop_counts, nrow = ny, ncol = nx, byrow = TRUE)
  
  melted_counts <- reshape2::melt(pop_counts) # here the first column is the row index (y) and the second the column index (x)
  colnames(melted_counts) <- c("row_ind", "col_ind", "value")
  
  melted_counts <- cbind(melted_counts, cell_loc = paste0(melted_counts$row_ind, ",", melted_counts$col_ind)) # add a column with the cell indexes
  melted_counts <- cbind(melted_counts, included = NA) # add a column with the info of which cells are included in the sample
  melted_counts$included[melted_counts$cell_loc %in% names(sample_result)] <- 1 # quadrats included in the sample are set to 1
  melted_counts$included <- as.factor(melted_counts$included)
  
  # colors to fill the cells included in the sample
  cols <- c("1" = quadrat.color, "2" = "#fec44f") # the value "2" is only used to plot ACS samples
  
  # make the plot with of without the grid
  if (show.grid == TRUE) {
    
    myplot <- ggplot(data = melted_counts, aes(x = col_ind, y = rev(row_ind))) + # we use rev() on row index because otherwise the graph is flipped
      # geom_point(data = pattern_df, aes(x = x, y = y), size = pointsize, color = "darkred", alpha = 0.8) + # plot the point pattern
      geom_point(data = pattern_df, aes(x = x, y = y), size = pointsize, color = "black", alpha = pointalpha) + # plot the point pattern
      # in the 2 lines below we subtract 0.5 to row and column index in order to align the tiles to the same coordinates than the points
      geom_tile(aes(x = col_ind - 0.5, y = row_ind - 0.5, fill = included), alpha = quadratalpha, color = "black") + # plot the tiles and fill cells that are included in the sample
      geom_text(aes(x = col_ind - 0.5, y = rev(row_ind - 0.5), label = value), color = "black", size = 4) + # label the tiles with the number of individuals they contain
      coord_fixed() + # to keep the graph square
      scale_fill_manual(values = cols, na.value = NA) +
      theme(axis.line = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = "none",
            # panel.background = element_blank(), 
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), plot.background = element_blank(),
            panel.background = element_rect(colour = "black", fill = "NA",
                                            linetype = "solid", size = 0.3))
    
  } else {
    
    myplot <- ggplot(data = melted_counts, aes(x = col_ind, y = rev(row_ind))) + # we use rev() on row index because otherwise the graph is flipped
      # geom_point(data = pattern_df, aes(x = x, y = y), size = pointsize, color = "darkred", alpha = 0.8) + # plot the point pattern
      geom_point(data = pattern_df, aes(x = x, y = y), size = pointsize, color = "black", alpha = pointalpha) + # plot the point pattern
      # in the 2 lines below we subtract 0.5 to row and column index in order to align the tiles to the same coordinates than the points
      geom_tile(aes(x = col_ind - 0.5, y = row_ind - 0.5, fill = included), alpha = quadratalpha) + # plot the tiles and fill cells that are included in the sample
      # geom_text(aes(x = col_ind - 0.5, y = rev(row_ind - 0.5), label = value), color = "black", size = 4) + # label the tiles with the number of individuals they contain
      coord_fixed() + # to keep the graph square
      scale_fill_manual(values = cols, na.value = NA) +
      theme(axis.line = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_blank(), axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = "none",
            # panel.background = element_blank(), 
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), plot.background = element_blank(),
            panel.background = element_rect(colour = "black", fill = "NA",
                                            linetype = "solid", size = 0.3))
    
    
  }
  
  
  return(myplot)
  
}

