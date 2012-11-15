#/usr/bin/env Rscript

# Plot concordant and discordant calls for multiple variant comparisons,
# statified by variant type
#
# Usage:
#   Rscript plot_discordance.r <summary CSV file> <target call name>

library(ggplot2)
library(stringr)
library(reshape2)
library(plyr)
library(grid)

args <- commandArgs(trailingOnly=TRUE)
in_file <- args[1]
call_target <- args[2]
out_file <- str_c(str_split(in_file, "[.]")[[1]][1], ".pdf")

d <- read.csv(in_file, header=TRUE)
d <- subset(d, call1==call_target,
            select=c(call2, type, concordant, discordant1, discordant2))

d2 <- melt(d, id=c("call2", "type"), measured=c("concordant", "discordant1", "discordant2"))


pretty_type_names <- function(val) {
  switch(val,
         total = "Total",
         snp = "SNP",
         indel = "Indel")
}
pretty_variable_names <- function(val) {
  switch(val,
         concordant = "Concordant",
         discordant1 = "Discordant (missing)",
         discordant2 = "Discordant (additional)")
}
#levels(d2$variable) <- unlist(lapply(levels(d2$variable), pretty_variable_names))
#levels(d2$type) <- unlist(lapply(levels(d2$type), pretty_type_names))
pretty_caller_names <- function(val) {
  switch(val,
         ifos = "GATK",
         ifosfb = "FreeBayes",
         ifospu = "mpileup",
         ifoscx = "cortex_var",
         ifosgh = "GATK Haplotype")
}
levels(d2$call2) <- unlist(lapply(levels(d2$call2), pretty_caller_names))

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

types = c("total", "snp", "indel")
cmps = c("concordant", "discordant1", "discordant2")

pdf(out_file, width = 9, height = 11)
pushViewport(viewport(layout = grid.layout(length(types), length(cmps))))
for (i in 1:length(cmps)) {
  for (j in 1:length(types)) {
    cmp <- cmps[i]
    cur_type <- types[j]
    d.sub <- subset(d2, type==cur_type & variable==cmp,
                    select=c(call2, value))
    p <- ggplot(d.sub, aes(x=call2, weight=value)) + geom_bar() +
      labs(title = paste(pretty_type_names(cur_type), pretty_variable_names(cmp), sep=" : ")) +
      theme(axis.title.y=element_blank(), axis.title.x=element_blank())
    if (min(d.sub$value) > 110000) {
      p <- p + coord_cartesian(ylim = c(110000, max(d.sub$value) + 1000))
    } else if (min(d.sub$value) > 9000) {
      p <- p + coord_cartesian(ylim = c(9000, max(d.sub$value + 500)))
    }
    ifelse(TRUE, #i == length(types),
           p <- p + theme(axis.text.x=element_text(angle = -90, hjust = 0)),
           p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()))
    print(p, vp = vplayout(i, j))
}}
dev.off()
