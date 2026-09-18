##Required packages
library(dplyr)
library(egg)
library(ggplot2)
library(ggpubr)
library(svglite)
library(tidyverse)

##Plot settings
# create plotting base
plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
plot_theme    <- theme_bw(base_size = 10) + theme(axis.title = element_text(size = 10),
axis.text = element_text(size=10),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.background = element_rect(fill = "transparent"),
axis.line = element_blank())
null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
axis.text.x=element_blank(),
axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
axis.text.y=element_blank(),
axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))
testPlot <- function(gg, grid = TRUE) {
if(grid)
gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
else
gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
gg_cbind <- function(..., widths = NULL) {
if(length(match.call()) - 2 != length(widths))
message("Number of widths does not match number of columns")
gg <- gtable_cbind(...)
panels = gg$layout$t[grep("panel", gg$layout$name)]
gg$widths[panels] <- unit(widths, "null")
return(gg)
}
gg_rbind <- function(..., heights = NULL, ncol = 2) {
if(length(match.call()) - 3 != length(heights))
message("Number of heights does not match number of rows")
gg <- gtable_rbind(...)
panels = gg$layout$t[grep("panel", gg$layout$name)]
gg$heights[panels] <- unit(rep(heights,each = ncol), "null")
return(gg)
}
## Extract legend
gg_legend <- function(a.gplot){
tmp <- ggplot_gtable(ggplot_build(a.gplot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
return(legend)
}
gg_blank <-
ggplot(data.frame()) +
geom_blank()

##Required tables:clin_cdkn2a_other, clin_data, anpl_prop_case, df_driver_instability_anno, fga_delta_pair_all

##Save the case order
clin_cdkn2a_other <- clin_cdkn2a_other %>% arrange(-received_rt, factor(driver_status, levels = c("R","S","P","WT")), factor(other_driver_status, levels = c("R","S","P","WT")), factor(WGD, levels = c("R","S","P","WT")), factor(Kataegis, levels = c("R","S","P","WT")), -mf_b, factor(seq, levels = c("WGS","WXS")))
case_order_seq_rt_driver_instability <- as.character(clin_cdkn2a_other$case_barcode)

##Clinical annotation landscape
clin_data <- clindata %>% mutate(case_barcode = factor(case_barcode, levels = case_order_seq_rt_driver_instability))

gg_clinical_seq_rt_driver_instability <-
clin_data %>%
gather(key = "type", value = "value", grade_change, received_rt,received_alk, is_hypermutator) %>%
mutate(type = factor(type,
levels = c("grade_change", "received_rt","received_alk",  "is_hypermutator"),
labels = c("Grade increase", "Received RT","Received Alkyl.", "Hypermutator"))) %>%
ggplot(aes(x=case_barcode)) +
geom_tile(aes(fill = factor(value), y = type)) +
scale_fill_manual(values=c("white", "#377eb8"), na.value = "gray90") +
labs(y="", fill = "Event")

testPlot(gg_clinical_seq_rt_driver_instability)

ggsave("clinical_annotation_seq_rt_driver_instability.pdf", plot = testPlot(gg_clinical_seq_rt_driver_instability), width = 25, height = 2.5, units = "in", dpi = 300)

##Aneuploidy proportion bar landscape
anpl_prop_case <- anpl_prop_case %>% mutate(case_barcode = factor(case_barcode, levels = case_order_seq_rt_driver_instability))

gg_anpl_prop_case_seq_rt_driver_instability <-
ggplot(anpl_prop_case, aes(x = case_barcode, y = aneuploidy_percent, fill=aneuploidy_type)) +
geom_bar(stat="identity") +
labs(y = "% aneuploidy") +
scale_fill_manual(values=c("Recurrent" = "#2FB3CA", "Initial" ="#CA2F66", "Shared"="#CA932F")) +
plot_grid

testPlot(gg_anpl_prop_case_seq_rt_driver_instability)

ggsave("anpl_prop_seq_rt_driver_instability.pdf", plot = testPlot(gg_anpl_prop_case_seq_rt_driver_instability), width = 25, height = 5, units = "in", dpi = 300)

##Driver annotation landscape
df_driver_instability_anno <- df_driver_instability_anno %>% mutate(case_barcode = factor(case_barcode, levels = case_order_seq_rt_driver_instability))

gg_seq_rt_driver_instability <- ggplot(df_driver_instability_anno, aes(x = case_barcode)) +
geom_tile(aes(y = event, fill = status), color = "gray80") +
scale_fill_manual(
values = c(
"A" = "white",
"R" = "#2FB3CA",
"S" = "#CA932F",
"P" = "#CA2F66"
),
na.value = "gray90"
) +
labs(y = "", fill = "status") +
theme_minimal()

testPlot(gg_seq_rt_driver_instability)

ggsave("annotation_seq_rt_driver_instability.pdf", plot = testPlot(gg_seq_rt_driver_instability), width = 25, height = 3, units = "in", dpi = 300)

##FGA delta landscape
fga_delta_pair_all <- fga_delta_pair_all %>% mutate(case_barcode = factor(case_barcode, levels = case_order_seq_rt_driver_instability))

gg_fga_delta_seq_rt_driver_instability <- ggplot(fga_delta_pair_all, aes(x= case_barcode)) +
geom_bar(aes(y = fga_norm_delta, fill = fga_group), alpha = 1, stat = "identity") +
scale_fill_manual(values = c("#8FC1DA","#D47264"))
labs(y = "Fraction of genome altered \n At recurrence")

testPlot(gg_fga_delta_seq_rt_driver_instability)

ggsave("fga_delta_seq_rt_driver_instability.pdf", plot = testPlot(gg_fga_delta_seq_rt_driver_instability), width = 25, height = 5, units = "in", dpi = 300)

##Combine all plots and legends to create full landscape in illustrator/inkscape