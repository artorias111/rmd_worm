#manhattan plots

library("tidyverse")
library("plotly")
library("DT")
library("ggbeeswarm")
library("knitr")
library("ggrepel")
library("genetics")
library("ggnewscale")
library("cowplot")

 #c.e. Manhattan plot


# in nextflow use sed to edit this field and make a copy of the .rmd for each trait
trait_name <- "telomere length"

#load gene list file for labelling
eleg_genes <- read.delim("/projects/b1059/projects/Shriram/data/elegans/goi/eleg_goi.tsv")

# load independent tests result
total_independent_tests <- read.table("/projects/b1059/projects/Shriram/data/elegans/nemascan_runs/Analysis_Results-20211219/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_test_cutoff <- -log10(0.05/total_independent_tests[[1]])

# load processed mapping data. 
#processed_mapping <- read.delim("/projects/b1059/projects/Shriram/worm/results/result_date/10_12/Analysis_Results-20211012/Mapping/Processed/processed_telomere_length_AGGREGATE_mapping.tsv", stringsAsFactors=FALSE) %>%
processed_mapping <- read.delim("/projects/b1059/projects/Shriram/data/elegans/nemascan_runs/Analysis_Results-20211219/Mapping/Processed/processed_length_AGGREGATE_mapping.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)


for.plot <- processed_mapping %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(algorithm = as.factor(algorithm))

BF <- processed_mapping %>% 
  dplyr::group_by(trait) %>% 
  dplyr::filter(log10p != 0) %>% 
  dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
  dplyr::ungroup() %>%
  dplyr::select(BF) %>%
  unique(.) %>%
  as.numeric()


# ntests <- data.table::fread(tests) %>%
#  as.numeric()
# EIGEN <- -log10(0.05/ntests)
BF.frame <- processed_mapping %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF, EIGEN  = independent_test_cutoff)

for.plot.ann <- for.plot %>%
  dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                log10p > BF.frame$EIGEN ~ "EIGEN",
                                TRUE ~ "NONSIG"))

sig.colors <- c("red","#EE4266")
names(sig.colors) <- c("BF","EIGEN")

man.plot <-  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::geom_point(data = for.plot.ann[which(for.plot.ann$sig != "NONSIG"),], 
                      mapping = aes(x = POS/1000000, 
                                    y = log10p,
                                    colour = sig),
                                    size=0.25)+
  ggplot2::scale_colour_manual(values = sig.colors) + 
  ggplot2::geom_point(data = for.plot[which(for.plot.ann$sig == "NONSIG"),], 
                      mapping = aes(x = POS/1000000, 
                                   y = log10p), 
                                   alpha = 0.5,
                                   size=0.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0,BF + 6.8)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  ggplot2::geom_point(data = eleg_genes, 
                      mapping=aes(x = POS/1000000, 
                                  y = log10p
                                  ),
                                  shape=25,
                                  size=1,
                                  fill="blue") +
  ggplot2::geom_hline(data = BF.frame, aes(yintercept = BF), linetype = 2) + 
  ggplot2::geom_hline(data = BF.frame, aes(yintercept = EIGEN), linetype = 3) + 
  ggplot2::labs(x = "",
                y = "" ) +
  ggplot2::theme(legend.position = "none", 
                 panel.grid = element_blank()) + 
  ggplot2::facet_grid(. ~ CHROM, scales = "free_x", space = "free") 

eleg_manplot <- man.plot
eleg_manplot

# c.b. Manhattan plot

# in nextflow use sed to edit this field and make a copy of the .rmd for each trait
trait_name <- "telomere length"

#load gene list file for labelling
brig_genes <- read.delim("/projects/b1059/projects/Shriram/data/briggsae/goi/brig_goi.tsv")

# load independent tests result
total_independent_tests <- read.table("/projects/b1059/projects/Shriram/data/briggsae/nemascan_runs/Analysis_Results-20211219/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_test_cutoff <- -log10(0.05/total_independent_tests[[1]])

# load processed mapping data. 
processed_mapping <- read.delim("/projects/b1059/projects/Shriram/data/briggsae/nemascan_runs/Analysis_Results-20211219/Mapping/Processed/processed_length_AGGREGATE_mapping.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)


for.plot <- processed_mapping %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(algorithm = as.factor(algorithm))

BF <- processed_mapping %>% 
  dplyr::group_by(trait) %>% 
  dplyr::filter(log10p != 0) %>% 
  dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
  dplyr::ungroup() %>%
  dplyr::select(BF) %>%
  unique(.) %>%
  as.numeric()


# ntests <- data.table::fread(tests) %>%
#  as.numeric()
# EIGEN <- -log10(0.05/ntests)
BF.frame <- processed_mapping %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF, EIGEN  = independent_test_cutoff)

for.plot.ann <- for.plot %>%
  dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                log10p > BF.frame$EIGEN ~ "EIGEN",
                                TRUE ~ "NONSIG"))

sig.colors <- c("red","#EE4266")
names(sig.colors) <- c("BF","EIGEN")

man.plot <-  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::geom_point(data = for.plot.ann[which(for.plot.ann$sig != "NONSIG"),], 
                      mapping = aes(x = POS/1000000, 
                                    y = log10p,
                                    colour = sig),
                                    size=0.25) +
  ggplot2::scale_colour_manual(values = sig.colors) + 
  ggplot2::geom_point(data = for.plot[which(for.plot.ann$sig == "NONSIG"),], 
                      mapping = aes(x = POS/1000000, 
                                   y = log10p), 
                                   alpha = 0.5,
                                   size=0.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0,BF + 6.7)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  ggplot2::geom_point(data = brig_genes, 
                      mapping=aes(x = POS/1000000, 
                                  y = log10p
                                  ),
                                  shape=25,
                                  size=1,#2
                                  fill="blue") +
  ggplot2::geom_hline(data = BF.frame, aes(yintercept = BF), linetype = 2) + 
  ggplot2::geom_hline(data = BF.frame, aes(yintercept = EIGEN), linetype = 3) + 
  ggplot2::labs(x = "",
                y = expression(-log[10](italic(p)))) +
                #y="") +
  ggplot2::theme(legend.position = "none", 
                 panel.grid = element_blank()) + 
  ggplot2::facet_grid(. ~ CHROM, scales = "free_x", space = "free") +
  ggplot2::theme(plot.title = element_text(face = "italic"))

brig_manplot <- man.plot

# c.t. Manhattan plot

# in nextflow use sed to edit this field and make a copy of the .rmd for each trait
trait_name <- "telomere length"

#load gene list file for labelling
trop_genes <- read.delim("/projects/b1059/projects/Shriram/data/tropicalis/goi/trop_goi.tsv")

# load independent tests result
total_independent_tests <- read.table("/projects/b1059/projects/Shriram/data/tropicalis/nemascan_runs/Analysis_Results-20220209/Genotype_Matrix/total_independent_tests.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

independent_test_cutoff <- -log10(0.05/total_independent_tests[[1]])

# load processed mapping data. 
processed_mapping <- read.delim("/projects/b1059/projects/Shriram/data/tropicalis/nemascan_runs/Analysis_Results-20220209/Mapping/Processed/processed_length_AGGREGATE_mapping.tsv", stringsAsFactors=FALSE) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F)


for.plot <- processed_mapping %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(algorithm = as.factor(algorithm))

BF <- processed_mapping %>% 
  dplyr::group_by(trait) %>% 
  dplyr::filter(log10p != 0) %>% 
  dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
  dplyr::ungroup() %>%
  dplyr::select(BF) %>%
  unique(.) %>%
  as.numeric()


# ntests <- data.table::fread(tests) %>%
#  as.numeric()
# EIGEN <- -log10(0.05/ntests)
BF.frame <- processed_mapping %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF, EIGEN  = independent_test_cutoff)

for.plot.ann <- for.plot %>%
  dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                log10p > BF.frame$EIGEN ~ "EIGEN",
                                TRUE ~ "NONSIG"))

sig.colors <- c("red","#EE4266")
names(sig.colors) <- c("BF","EIGEN")

man.plot <-  ggplot2::ggplot() + 
  ggplot2::theme_bw() + 
  ggplot2::geom_point(data = for.plot.ann[which(for.plot.ann$sig != "NONSIG"),], 
                      mapping = aes(x = POS/1000000, 
                                    y = log10p,
                                    colour = sig),
                                    size=0.25) +
  ggplot2::scale_colour_manual(values = sig.colors) + 
  ggplot2::geom_point(data = for.plot[which(for.plot.ann$sig == "NONSIG"),], 
                      mapping = aes(x = POS/1000000, 
                                   y = log10p), 
                                   alpha = 0.5,
                                   size=0.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0,5.05)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  ggplot2::geom_point(data = trop_genes, 
                      mapping=aes(x = POS/1000000, 
                                  y = log10p
                                  ),
                                  shape=25,
                                  size=1,
                                  fill="blue") +
  ggplot2::geom_hline(data = BF.frame, aes(yintercept = BF), linetype = 2) + 
  ggplot2::geom_hline(data = BF.frame, aes(yintercept = EIGEN), linetype = 3) + 
  ggplot2::labs(x = "position",
                #y = expression(-log[10](italic(p)))) +
                y="")+
  ggplot2::theme(legend.position = "none", 
                 panel.grid = element_blank()) + 
  ggplot2::facet_grid(. ~ CHROM, scales = "free_x", space = "free") +
  #ggplot2::ylim(0.0,5.0)+
  ggplot2::theme(plot.title = element_text(face = "italic"))

trop_manplot <- man.plot

combine plots: 
p <- plot_grid(eleg_manplot,brig_manplot,trop_manplot,ncol=1,align = "v",labels = c("A","B","C"))
p
ggsave("a.png",width=7.5,height=8,units="in",dpi=300)
