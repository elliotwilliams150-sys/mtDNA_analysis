# Plotting themes and visualization helpers
# Functions for consistent figure styling across the project

#' Transparent background theme for publication-quality plots
#'
#' Creates a ggplot2 theme with transparent backgrounds for all elements,
#' useful for plots that will be composited into figures with custom backgrounds.
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' p <- p + theme_transparent
#' ggsave('myplot.pdf', p, width = 4, height = 3, bg = 'transparent')
theme_transparent <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    plot.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    legend.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    legend.box.background = ggplot2::element_rect(fill = 'transparent', colour = NA),
    legend.key = ggplot2::element_rect(fill = 'transparent', colour = NA)
)
