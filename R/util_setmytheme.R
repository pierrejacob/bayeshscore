#'@rdname setmytheme
#'@title Customize graphical settings
#'@description This function customizes the theme used by ggplot.
#'@export
setmytheme <- function(){
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               panel.margin = unit(2, "lines"),
               legend.position = "bottom")
}