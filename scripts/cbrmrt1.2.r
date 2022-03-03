#candidate gene - mrt-1.2
#worm - C. briggsae

cbrmrt1 <- brig_geno %>% filter(POS>=15329043 & POS<=15332022) %>% filter(CHROM=="I")
#cbrmrt1 <- brig_geno %>% filter(POS==15329877) %>% filter(CHROM=="I")
brig_temp1 <- cbrmrt1 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)
brig_temp <- brig_temp1 %>% inner_join(cb_len_isotype) %>% filter(altref==-1) %>% dplyr::select(POS,length,isotype)
mrt_gene <- read_tsv('/projects/b1059/projects/Shriram/data/briggsae/goi/regional_gffs/cbrmrt1.tsv')

position_list <- brig_temp %>% count(POS) %>% dplyr::select(POS)

cbrmrt1_plot <- brig_temp %>% ggplot(aes(x=POS,y=length,group=POS)) + 
  geom_boxplot() + 
  geom_jitter(color="black", alpha=0.1) + 
  xlim(15329043,15332022) + 
  theme(axis.title = element_blank(), legend.position = "none",) + 
  geom_vline(xintercept = c(mrt_gene$start),linetype="dashed",alpha=0.5) 

cbrmrt1_gene_plot <- ggplot(mrt_gene, aes(xmin = start, xmax = end, y = molecule, fill = gene)) + 
  gggenes::geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) + 
  theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(), legend.position = "none") +
  geom_vline(xintercept = c(mrt_gene$start),linetype="dashed",alpha=0.5) 

#green intron, pink exon


p <- plot_grid(cbrmrt1_gene_plot,cbrmrt1_plot,ncol=1,rel_widths = c(1,1),rel_heights = c(1,6))
p <- plot_grid(p)
p

iso_list <- brig_temp1 %>% filter(altref==-1) %>% dplyr::select(isotype) %>% distinct() %>% add_column(status=1)
cbrmrt_strains <- cb_len_isotype %>% mutate(isotype=fct_reorder(isotype,length))

x <- as_tibble(setdiff(cbrmrt_strains$isotype,iso_list$isotype)) 
y <- x %>% add_column(status=-1)
y <- y %>% dplyr::rename(isotype=value)

z <- full_join(iso_list,y)
k <- inner_join(cb_len_isotype,z)

k <- k %>% mutate(status=as.character(status)) %>% mutate(status=str_replace(status,"-1","alt")) %>% mutate(status=str_replace(status,"1","ref"))

cbrmrt_bar <- k %>% mutate(isotype=fct_reorder(isotype,length)) %>%
      ggplot(aes(x=isotype,y=length,fill=as.factor(status))) + 
        labs(fill="alt") + 
        geom_bar(stat='identity') +
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
        xlab("briggsae isotypes")+
        ggtitle("cbr mrt-1")

cbrmrt_bar

#each variant in each gene? 
