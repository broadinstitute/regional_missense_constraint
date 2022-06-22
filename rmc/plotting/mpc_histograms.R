## Code used to create histograms of MPC scores
## in de novo variants from
## neurodevelopmental disorders (NDD) cases vs controls
## Designed to be run in RStudio
## (And manually adjust input/output paths)
# Load ggplot2
library(ggplot2)

# Read in input data
tsv <- read.table(
  "/Users/kchao/Desktop/RMC/mpc/ndd_mpc.tsv.bgz",
  header=T,
)
head(tsv, 5)
# This should show something like
#   case_control      mpc
# 1         case 0.274310
# 2         case 0.037434
# 3         case 0.380640
# 4         case 0.020389
# 5         case 0.432420


# Create plot
p <- ggplot(tsv, aes(x=mpc, fill=case_control)) + geom_histogram(
  # Scale counts to proportions since there are many more variants from cases than controls
  aes(y=..ncount..),
  # Make histograms slightly transparent and add black outline
  alpha=0.5, bins=75, position="identity", color="black") + scale_fill_manual(
  # Manually set controls to grey and cases to light green color (to match RMC preprint)
  values=c("#69B3A2", "#C3C6C5")
)
# Remove "chart junk" (background grid, y-axis)
p <- p + theme(
  text = element_text(size=30),
  panel.background = element_rect(fill="transparent"),
  plot.background = element_rect(fill="transparent", color=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line.x = element_line(size=0.5),
  axis.line.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.text.y = element_blank(),
  axis.text.x = element_text(size=20),
  legend.title = element_blank(),
)
p <- p + scale_x_continuous(name="MPC")
p <- p + scale_y_continuous(name="")
p

# Save plot to output file
ggsave(
  "/Users/kchao/Desktop/RMC/mpc/all_ndd_mpc_histogram.png",
  p,
  height = 10,
  dpi = 300,
  bg = "transparent"
)
