##' function named ggplot_dual_axis()
##' Takes 2 ggplot plots and makes a dual y-axis plot
##' function takes 2 compulsory arguments and 1 optional argument
##' arg lhs is the ggplot whose y-axis is to be displayed on the left
##' arg rhs is the ggplot whose y-axis is to be displayed on the right
##' arg 'axis.title.y.rhs' takes value "rotate" to rotate right y-axis label
##' The function does as little as possible, namely:
##'  # display the lhs plot without minor grid lines and with a
##'  transparent background to allow grid lines to show
##'  # display the rhs plot without minor grid lines and with a
##'  secondary y axis, a rotated axis label, without minor grid lines
##'  # justify the y-axis label by setting 'hjust = 0' in 'axis.text.y'
##'  # rotate the right plot 'axis.title.y' by 270 degrees, for symmetry
##'  # rotation can be turned off with 'axis.title.y.rhs' option
##'

ggplot_dual_axis <- function(lhs, rhs, axis.title.y.rhs = "rotate",
                             second.lab=TRUE, show.axis=TRUE) {
  # 1. Fix the right y-axis label justification
  rhs <- rhs + theme(axis.text.y = element_text(hjust = 0))
  # 2. Rotate the right y-axis label by 270 degrees by default
  if (missing(axis.title.y.rhs) |
      axis.title.y.rhs %in% c("rotate", "rotated")) {
    rhs <- rhs + theme(axis.title.y = element_text(angle = 270))
  }
  # 3a. Use only major grid lines for the left axis
  lhs <- lhs + theme(panel.grid.minor = element_blank())
  # 3b. Use only major grid lines for the right axis
  #     force transparency of the backgrounds to allow grid lines to show
  rhs <- rhs + theme(panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = "transparent", colour = NA),
                     plot.background = element_rect(fill = "transparent", colour = NA))
  # Process gtable objects
  # 4. Extract gtable
  library("gtable") # loads the grid package
  g1 <- ggplot_gtable(ggplot_build(lhs))
  g2 <- ggplot_gtable(ggplot_build(rhs))
  # 5. Overlap the panel of the rhs plot on that of the lhs plot
  pp <- c(subset(g1$layout, name == "panel", se = t:r))
  g <- gtable_add_grob(g1,
                       g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
  # Tweak axis position and labels
  ia <- which(g2$layout$name == "axis-l")
  ga <- g2$grobs[[ia]]
  ax <- ga$children[["axis"]]  # ga$children[[2]]
  ax$widths <- rev(ax$widths)
  ax$grobs <- rev(ax$grobs)
  if(!is.null(ax$grobs[[1]]$x)) ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
  g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
  if(show.axis) g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)

  if(second.lab){
    g <- gtable_add_grob(g, g2$grobs[[7]], pp$t, length(g$widths), pp$b)
    #Jeans fix for axis label
    g$widths[[6]] = g$widths[[5]]
  }
  # Display plot with arrangeGrob wrapper arrangeGrob(g)
  return(g)
  library("gridExtra")
  grid.newpage()
  return(arrangeGrob(g))
}
