# Script is designed to be run in RStudio
# Load libraries
library(ggplot2)
library(reshape2)

# Create data
# Update this path to point to TSV with OE bins annotated with proportion coding base pairs,
# proportion de novo missense from controls/cases and proportion ClinVar
# pathogenic/likely pathogenic missense variants in haploinsufficient genes
df <- read.table("/Users/kchao/Desktop/assessment.tsv", header=TRUE)
df

# Reorder columns
# Desired order to display:
# base pairs, controls, cases, ClinVar
df <- df[, c(1,2, 4, 5, 3)]
df

# Convert dataframe into long format
longdf <- melt(data=df,
  variable.name = "x"
)
longdf

# Recreate dataframe in long format for plotting
x_display <- longdf['x']
y_display <- longdf['value']
oe_bin <- longdf['oe_bin']
df <- data.frame(x=x_display, y=y_display, group=oe_bin)
df

p <- ggplot(df, aes(x=x, y=value, fill=oe_bin)) +
  geom_bar(stat = "identity", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(
    values = c(
      "#B90101", "#C93838", "#DA6F6F",
      "#EAA5A5", "#FADCDC"
      ))
xlabels <- c("Base pairs", "Controls", "Neurodev cases", "Severe HI variants")
p <- p + scale_x_discrete(name="",  labels=xlabels)
p <- p + scale_y_continuous(name="Proportion")
p <- p + guides(fill=guide_legend(title="Obs/Exp"))

# Set text size to 20, transparent backgrounds,
# remove major/minor grids, and add x/y axis lines
p <- p + theme(text = element_text(size=20),
               panel.background = element_rect(fill = "transparent"),
               plot.background = element_rect(fill = "transparent", color = NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.line.x = element_line(size = 0.5),
               axis.line.y = element_line(size = 0.5),
)

# Edit this line to save plot to output path
ggsave("/Users/kchao/Desktop/distributions.png",
       p,
       width = 20, height = 10,
       dpi = 300,
       bg = "transparent")
p
