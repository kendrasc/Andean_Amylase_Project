library(ggplot2)

# Read in Data
PEL_tajima = read.delim("PEL_Tajimas_D.Tajima.txt")

# Figure
PEL <- ggplot(PEL_tajima, aes(x=BIN_START, y=TajimaD)) +
      geom_point( alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("maroon"), 22 )) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      xlab("BP Position") + 
      ylab("Tajima values") +
      ggtitle("PEL") +
  
# Custom the theme:
      theme_bw() +
      theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text  = element_text(size= 10, 
                           color = "darkred"),
      panel.background = element_rect(fill = "white", 
                                      color="white"),  
      legend.key = element_rect(fill = "white"), 
      plot.title=element_text( size = 15), 
      axis.title.y = element_text( size = 12),
      axis.title.x = element_text( size = 12), 
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(color = "black"))


