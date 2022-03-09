# Script is designed to be run in RStudio
# Load ggplot2
library(ggplot2)

# Create data
# Update this path to point to TSV
df <- read.table("/Users/kchao/Desktop/sections_per_transcript.tsv", header=TRUE)
df

# Create barplot
p <- ggplot(df, aes(x=num_sections, y=num_transcripts)) +
  geom_bar(stat = "identity", fill="white")
p <- p + labs(y="Number of transcripts")

# Update this seq with maximum number of breaks -- seq(1, max breaks + 1, 1)
p <- p + scale_x_continuous("Number of missense constraint regions", seq(1, 13, 1))

# Set text size to 25 and color to white
# Also make background of plot transparent, remove major/minor grids
p <- p + theme(text = element_text(size=25, color="white"),
               panel.background = element_rect(fill = "transparent"),
               plot.background = element_rect(fill = "transparent", color = NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.background = element_rect(fill = "transparent"),
               legend.box.background = element_rect(fill = "transparent"),
               )

# Format axis ticks/text for better display
p <- p + theme(axis.title.x = element_text(colour = "white",
                                           margin = margin(t = 7)),
               axis.title.y = element_text(colour = "white",
                                           margin = margin(r = 7)),
               axis.ticks.x = element_line(colour = "white"),
               axis.ticks.y = element_line(colour = "white"),
               axis.text.x = element_text(colour = "white",
                                          margin = margin(t = 7)),
               axis.text.y = element_text(colour = "white",
                                          margin = margin(r = 7)),
               axis.line.x = element_line(colour =  "white", size = 4),
               axis.line.y = element_line(colour =  "white", size = 4),
               )

# Update this to save PNG to output path
ggsave("/Users/kchao/Desktop/sections_per_transcript.png",
       p,
       width = 10, height = 5,
       dpi = 300,
       bg = "transparent")
p
