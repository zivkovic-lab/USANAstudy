suppressPackageStartupMessages(library(ggplot2))

# Function to extract legend ----------------------------------------------

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# Theme for boxplot -----------------------------------------------------

theme_boxplot <- function() {
    library(grid)
    theme_bw() %+replace%
        theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1),
                hjust = 0.5,
                margin = unit(c(0,0,2,0), "mm")
            ),
            panel.background = element_rect(colour = "black"),
            panel.grid = element_blank(),
            axis.title = element_text(
                face = "plain", 
                size = rel(1)
            ),
            axis.text = element_text(
                size = rel(1),
                face = "plain"
            ),
            axis.ticks = element_line(size = 1),
            axis.title.x = element_blank(),
            strip.background = element_rect(
                fill = "grey70"
            ),
            strip.text = element_text(
                size = rel(0.7),
                face = "bold",
                margin = unit(c(1, 6, 1, 6), "mm")
            )
        )
}

# theme for pie chart -----------------------------------------------------

abund_theme <- function() {
    library(grid)
    theme_bw() %+replace%
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        legend.background = element_blank()
    )
}

# Theme for scatter plot --------------------------------------------------

theme_point <- function() {
    library(grid)
    theme_bw() %+replace%
    theme(
        plot.title = element_text(
                face = "bold",
                size = rel(1),
                hjust = 0.5,
                margin = unit(c(0,0,2,0), "mm")
            ),
        axis.title = element_text(
                face = "plain", 
                size = rel(1)
            ),
        axis.text = element_text(
                size = rel(1),
                face = "plain"
            ),
        axis.ticks = element_line(size = 1),
        panel.background = element_rect(colour = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(
            fill = "grey70"
        ),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        legend.background = element_blank()
    )
}

theme_barplot <- function() {
    library(grid)
    theme_bw() %+replace%
        theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1),
                hjust = 0.5,
                margin = unit(c(0,0,2,0), "mm")
            ),
            panel.background = element_rect(colour = "black"),
            panel.grid = element_blank(),
            axis.title = element_text(
                face = "plain", 
                size = rel(1)
            ),
            axis.text = element_text(
                size = rel(1),
                face = "plain"
            ),
            axis.ticks = element_line(size = 1),
            axis.title.x = element_blank(),
            strip.background = element_rect(
                fill = "grey70"
            ),
            strip.text = element_text(
                size = rel(0.7),
                face = "bold",
                margin = unit(c(1, 6, 1, 6), "mm")
            ),
            legend.text = element_text(size = rel(1)),
            legend.title = element_text(size = rel(1)),
        )
}
