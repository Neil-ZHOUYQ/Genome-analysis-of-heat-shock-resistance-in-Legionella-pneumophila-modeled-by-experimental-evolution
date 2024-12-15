##########fig1
library(readxl)
fig1b<- read_excel("C:/Users/97021/Desktop/CUHK/5030/group project/heat shock/data/fig 1B RT.xlsx", sheet = "Sheet1")

library(dplyr)
fig1b_filtered <- fig1b %>% filter(log10CFU != 0)

library(ggplot2)
pdf(file="fig1b.pdf", width = 12, height = 5) 
custom_colors <- c("RT" = "red", "55" = "blue", "59" = "grey", "57" = "purple")

ggplot(fig1b_filtered, aes(x = factor(Passage), y = log10CFU, color = factor(Temperature), shape = factor(Temperature), group = Temperature)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(x = "Passage Number", y = "log10(CFU)", title = "Figure 1B") +
  scale_x_discrete(breaks = unique(fig1b_filtered$Passage)) +
  scale_color_manual(values = custom_colors)


dev.off()







##############circlar graph


library(tidyverse)
library(circlize)
library(openxlsx)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::version()
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(ComplexHeatmap)
library(reshape2)


setwd("C:/Users/97021/Desktop/CUHK/5030/group project/heat shock/data")

# read data
data <- read.xlsx("key_gene_SNP.xlsx", 1) %>%
    dplyr::select(-P1) %>%
    melt(., id.vars='gene') %>%
    na.omit %>%
    separate_rows(value, sep=", ") %>%
    separate(value, into=c('position', 'mutation'), sep=':') %>%
    mutate(position = as.numeric(position)) %>%
    mutate(type = ifelse(nchar(mutation) >4, "indel", "SNP")) %>%
    mutate(group = ifelse(grepl('HA', variable), "test", "control"))


# define data frame
dotdata <- dplyr::select(data, -mutation) %>% 
    mutate(variable = factor(variable ,levels=c('HA_1_70','HA_2_30','HA_2_70','HA_3_70','HA_5_70','HA_6_70', 'C_1_70', 'C_2_70', 'C_3_70', 'C_4_70', 'C_5_70', 'C_6_70') %>% rev)) %>% 
    mutate(y=12.5-as.numeric(variable)) %>%
    arrange(position) 
genedata <- mutate(dotdata, Name=gene) %>%group_by( Name) %>% 
    summarise(position = mean(position)) %>%
    na.omit
# name with color
coldata <- mutate(dotdata, Name=gene) %>%group_by( Name, group) %>%
    summarise(n=n()) %>%
    dcast(Name~group, value.var='n') 
coldata[is.na(coldata)] <- 0
coldata <- mutate(coldata, direction = test - control ) %>%
    mutate(direction = case_when(
        direction>0 ~ 'up',
        direction<0 ~ 'down',
        direction==0 ~ 'no change')) %>%
    mutate(col = ifelse(direction=='up', '#ef4939', ifelse(direction=='down', '#816db0', 'black')))
gene_cols <- coldata$col %>%setNames(coldata$Name)

circos.clear()
pdf(file="circos.20241214.pdf", width=9, height=9)
circos.par(
    cell.padding=c(0,0,0,0), 
    canvas.xlim =c(-1.,1.),  canvas.ylim = c(-1.2,1.2),  #circle.margin = c(0.01, 0.01, 0.01, 0.01),
    start.degree = 90-30/2, #start angle
    track.margin = c(0, 0), cell.padding= c(0, 0, 0, 0),
    gap.after=30   
)


circos.initialize(c("A"), xlim = c(0, 3.5e6))

# positions
ticks <- seq(0, 3.5e6,0.1e6)
ticklabels <- ticks/1e6
ticklabels[!ticklabels %in% seq(0, 3.5,0.5)] <- NA

circos.track('A', ylim=c(0,12),panel.fun = function(x, y) {
    circos.axis(
        labels = ifelse(is.na(ticklabels), '', paste0(sprintf('%.1f',ticklabels), ' Mbp')), 
        major.at = ticks, lwd=2, col='black', labels.cex = 1,
        labels.facing = 'clockwise',major.tick.length = mm_y(2), labels.pos.adjust = FALSE)
    
    
    for(i in c(7, 9, 11)) circos.rect(0, i, 3.5e6, i+1, col = '#BCD2EE', border = NA)
    for(i in c(6, 8, 10)) circos.rect(0, i, 3.5e6, i+1, col = '#CAE1FF', border = NA)
    for(i in c(0, 2, 4)) circos.rect(0, i, 3.5e6, i+1, col = '#FFDAB9', border = NA)
    for(i in c(1, 3, 5)) circos.rect(0, i, 3.5e6, i+1, col = '#FFEFD5', border = NA)
    circos.lines(x=CELL_META$cell.xlim, y=c(6,6), col = '#737374', lwd=2)

    
    for(i in 1:nrow(dotdata)){
        circos.points(x=dotdata$position[i], y=dotdata$y[i], 
                      pch=ifelse(dotdata$type[i]=='SNP', 19, 17), 
                      col = ifelse(dotdata$group[i]=='control', '#1e71b4', '#ee3b2b'), cex=0.6)
    }
    circos.text( x = genedata$position, 
        y= -0.5, 
        adj=c(1, 0.5),
        labels = genedata$Name, cex=0.8, 
        col = gene_cols[genedata$Name], 
        facing = 'clockwise', 
        niceFacing = TRUE)

    circos.text(CELL_META$cell.xlim[2]*360/(360-30/2), 4, 'heat', cex = 1.2, adj = c(0.5, 0.5), facing = 'downward')
    circos.text(CELL_META$cell.xlim[2]*360/(360-30/2), 2, 'adapted', cex = 1.2, adj = c(0.5, 0.5), facing = 'downward')
    circos.text(CELL_META$cell.xlim[2]*360/(360-30/2), 9, 'control', cex = 1.2, adj = c(0.5, 0.5), facing = 'downward')

}, bg.col = 'white', bg.border = NA, track.height=0.32)

#legends
lgd = Legend(labels=c('point mutation', 'in/del'),type='points',
             legend_gp = gpar(col='#ef3b2c',cex=20), pch = c(19,17), by_row = TRUE, gap = unit(1, "cm"), row_gap = unit(1, "mm"),
             labels_gp = gpar(fontsize = 10), background = "white",
             title='', border=NA)
draw(lgd,x = unit(0.9,"npc"),y=unit(0.2,"npc"),just = c("right","bottom"))

dev.off()
