#Template for plotting barplot with altref variants, and the gene model with boxplots for each position in the gene

#import the genotype matrix from nemascan run
eleg_geno <- read_tsv("/projects/b1059/projects/Shriram/data/elegans/nemascan_runs/Analysis_Results-20211219/Genotype_Matrix/Genotype_Matrix.tsv")

#include regions of the gene/variant
celpot2 <- eleg_geno %>% filter(POS==14518226) %>% filter(CHROM=="II")

#list out the positions of each variant in the candidate gene here


#convert genotype matrix from wide to long format data
eleg_temp1 <- celpot2 %>% pivot_longer(
    cols=!CHROM & !POS & !REF & !ALT,
    names_to = "isotype",
    values_to = "altref"
)

#get data for ref variants in the shortlisted genomic region (status =1 implies ref sequence)
iso_list <- eleg_temp1 %>% filter(altref==-1) %>% dplyr::select(isotype) %>% distinct() %>% add_column(status=1)
celpot2_strains <- ce_len_isotype %>% mutate(isotype=fct_reorder(isotype,length))

#get data for alt variants in the shortlisted genomic region (status=-1 implies alt sequence)
x <- as_tibble(setdiff(celpot2_strains$isotype,iso_list$isotype)) 
y <- x %>% add_column(status=-1) %>% dplyr::rename(isotype=value)

#join and clean up data, replace -1 with 'alt' and 1 with 'ref'
z <- full_join(iso_list,y)
k <- inner_join(ce_len_isotype,z)
k <- k %>% mutate(status=as.character(status)) %>% mutate(status=str_replace(status,"-1","alt")) %>% mutate(status=str_replace(status,"1","ref"))

celpot2_bar <- k %>% mutate(isotype=fct_reorder(isotype,length)) %>%
      ggplot(aes(x=isotype,y=length,fill=as.factor(status))) + 
        labs(fill="isotype") + 
        geom_bar(stat='identity') +
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
        xlab("elegans isotypes")+
ggtitle("c eleg pot-2")

celpot2_bar
ggsave("celpot2.png",dpi=300,height=4,width = 8, units = "in")