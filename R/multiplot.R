#' @title Multiple plot function
#'
#' @description Enables multiple ggplot objects to appear in a panel.
#' @param ... ggplot objects can be passed in
#' @param plotlist A \code{list} of ggplot objects.
#' @param cols  An \code{integer} containing the number of columns in layout
#' @param layout A \code{matrix} specifying the layout. If present, 'cols' is ignored.
#' @return Creates a combination ggplot panel.
#' @keywords internal
#' @details
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' 
#' See the R Graphics Cookbook for more information
#' @author Winston Chang 
multiplot = function(..., plotlist=NULL, cols=1, layout=NULL) {
  #require(grid)
  
  grid.newpage=pushViewport=viewport=grid.layout=NULL
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}