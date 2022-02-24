library(brms)
library(tidyverse) 
library(ggridges)
library(ggthemes)
library(cowplot)
library(here)

file = "aseg"
indir <- paste0("../data/",file,".RData")
load(indir)

aa <- fixef(fm, summary = FALSE)/SCALE
bb <- lapply(ranef(fm, summary = FALSE), `/`, SCALE) # Extract Group-Level (or random-effect)

EOIq <- colnames(aa)

psROI <- function(aa, bb, tm) {
  ps <- apply(bb[['ROI']][,,tm], 2, '+', aa[,tm])
  return(ps)
}

outdir <- paste0("../results/",file,"/")
if (!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}

# generates ridge plots for effects (Intercept, cond, TRAIT and STATE), and saves them.
for (ii in 1:length(EOIq)){
  posteriors <- psROI(aa, bb, EOIq[ii])
  write.table(posteriors, paste0(outdir, EOIq[ii], '_post.txt'), sep='\t',
              col.names = TRUE, row.names = FALSE)
}

rois <- c("Left BA", "Left CA1", "Left CA3", "Left CeA", "Left DG", "Left LA", "Left Subiculum",
          "Right BA", "Right CA1", "Right CA3", "Right CeA", "Right DG", "Right LA", "Right Subiculum") # Change the ROI name for the image

df <- read.table(file = paste0(outdir,'Clinical_Trajectory3_post.txt'), sep = "", header = TRUE)
fig_title = "Clinical Trajectory"
fig_name = "../results/T1_group_sub.png"

df$X <- NULL
colnames(df) <- rois
iterations <- length(df[,1])
df.long <- df %>% gather(ROI)
colnames(df.long) <- c("ROI", "Y")
df.long <- df.long %>%
  mutate(index = rep(1:length(rois), each = iterations)) %>%
  group_by(ROI) %>%
  mutate(mean = mean(Y)) %>%
  mutate(p = ((sum(Y>0)/iterations))) %>%
  mutate(p.plot = p) %>%
  mutate(p.plot = replace(p.plot, p.plot > 0.15 & p.plot < 0.85, NA))

P <- df.long %>%
  select(ROI, index, mean, p) %>%
  unique() %>%
  arrange(p)

intercept <- ggplot(df.long, aes(x = Y, y = as.numeric(reorder(index,p)), group = ROI, fill = p.plot)) +
  #coord_cartesian(xlim = c(-0.14, 0.23)) +
  geom_density_ridges(quantile_lines = TRUE,
                      quantiles = 2,
                      scale = 2.3,
                      #rel_min_height = .01,
                      color = "#404040",
                      size = .75) +
  geom_vline(xintercept = 0, alpha = .85, color = "black", size = .8) +
  scale_y_continuous(breaks = 1:length(P$ROI),
                     expand = c(0,0.1),
                     labels = P$ROI,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(P$ROI),
                                         labels = format(round(P$p, 3),nsmall = 2))) +
  #scale_x_continuous(breaks = c(-0.2,-0.1, 0, 0.1, 0.2), labels = c("-0.20","-0.10", "0", "0.10", "0.20")) +
  scale_fill_gradientn(limits = c(0,1),
                       colors = c("blue","cyan",
                                  "gray","yellow","red"),
                       values = c(0,0.15,
                                  0.150000001, 0.85,
                                  0.850000001, 1.0),
                       breaks = c(0, 0.15, 0.85, 1)) +
  scale_color_gradientn(limits = c(0,1),
                        colors = c("blue","cyan",
                                   "gray","yellow","red"),
                        values = c(0,0.15,
                                   0.150000001, 0.85,
                                   0.850000001, 1.0),
                        breaks = c(0, 0.15, 0.85, 1)) +
  guides(fill = guide_colorbar(barwidth = 1.25,
                               barheight = 8,
                               nbin = 50,
                               frame.colour = "black",
                               frame.linewidth = 1.25,
                               ticks.colour = "black")) +
  theme_stata() +
  theme(
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major.y = element_line(color = "lightgray", size = .75),
    plot.title = element_text(size = 12.5, margin = unit(c(0,0.1,.25,02),"cm"), face = "plain", hjust = 0.5),
    legend.title = element_text(size = 11, hjust = 0),
    legend.text = element_text(size = 8, angle = 0),
    legend.position = c(.85,.31),
    legend.background = element_rect(),
    legend.box.background = element_rect(colour = "black", size = .5),
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 16),
    axis.line = element_line(size = .25),
    axis.text.y = element_text(size= 9, color = "black", margin = unit(c(0,-0.05,0,0.05),"cm"), angle = 0, vjust = 0),
    axis.text.y.right = element_text(size = 9, color = "black",margin = unit(c(0,0,0,-0.05),"cm"), angle = 0),
    axis.text.x = element_text(size = 9, color = "black", margin = unit(c(0.04,0,0,0),"cm")),
    axis.ticks.x = element_line(size = .5),
    axis.ticks.length=unit(.2, "cm"),
    axis.ticks.y = element_blank()) +
  labs(
    x = NULL,
    y = NULL,
    title = fig_title,
    fill = "P+")

plot <- ggdraw(intercept) +
  draw_text("P+", x = .94, y = .92, size= 11) +
  draw_text("Remission > non-Remission", x = .31, y = 0.015, size = 11) +
  draw_text("Remission < non-Remission", x = .74, y = 0.015, size = 11)

ggsave(fig_name, plot = plot, dpi = 600, height = 120, width = 180, units = "mm")

print("Done with Figure 3, onto the next...")
