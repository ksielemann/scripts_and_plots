#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")

# install.packages("devtools")
#devtools::install_github("thackl/thacklr")
#devtools::install_github("thackl/gggenomes")



library(S4Vectors)
library(gggenomes)
library(svglite) 
library(tidyverse)
library(RColorBrewer)


# parse sequence length and some metadata from fasta file
emale_seqs <- read_seq_len("merged_species.fasta")

#extract QTL region
loci <- tribble(
  ~locus_id,~seq_id, ~start, ~end,
  "1791621-1848871", "U2Bv",1791621, 1848871,
  "1906930-2004170", "KWS2320", 1906930, 2004170,
  "1735399-1808239","EL10",1735399,1808239,
  "294447-366180","Bmar",294447,366180,
  "280937-338886","Bptu",280937,338886,
  "1813613-1909925","RefBeet",1813613,1909925)


#Annotate genes
emale_genes <- read_gff("merged_species_genes.gff") %>%
  rename(feature_id=ID) #%>%                       # we'll need this later
#mutate(gc_cont=as.numeric(gc_cont))             # per gene GC-content


#add 4 letter code as separate column
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

emale_genes <- emale_genes %>% 
  mutate(letter_code = substrRight(feature_id,4))


p2 <- gggenomes(emale_genes, emale_seqs) %>%
  focus(.loci=loci)  + 
  geom_seq() + geom_seq_label() + geom_bin_label() +
  #geom_gene_label(aes(label=Name), size=4, angle=90) +
  geom_gene(aes(fill=strand), position="strand")
p2


#Compare genome syteny
#uncomment once!!!
emale_links <- read_paf("alignment_region_species.paf")

emale_links_filtered <- filter(emale_links, (start<2154734 & start>1321396)|(start2<2154734 & start2>1321396)|(end<2154734 & end>1321396)|(end2<2154734 & end2>1321396)|
                                 (start<366180 & start>280937)|(start2<366180 & start2>280937)|(end<366180 & end>280937)|(end2<366180 & end2>280937))

emale_links_filtered <- filter(emale_links_filtered, AS > 2000)

emale_links_filtered

p4 <- gggenomes(emale_genes, emale_seqs, links=emale_links_filtered) %>%
  focus(.loci=loci)  + 
  geom_seq() + geom_bin_label() +
  #geom_feature(size=5, data=use_features(features)) +
  geom_gene(aes(fill=strand), position="strand") +
  geom_link()

#p4 <- p4 %>% flip_bins(5:6)
p4




#(Add funtional annotation/description (Swissprot BLAST)/Interpro domains)
annot <- read_tsv("merged_annotation_species.txt", col_names=NULL)
colnames(annot) <- c("Name","description")

genes_and_info_and_description <- left_join(emale_genes, annot, by=c("Name"))
genes_and_info_and_description <- distinct(genes_and_info_and_description)



p9 <- gggenomes(genes_and_info_and_description, seqs = emale_seqs, links=emale_links_filtered) %>% #rgas
  focus(.loci=loci)  + 
  geom_seq() + geom_bin_label(size=10) + 
  #geom_feat(size=15) +
  #geom_feature(size=5, data=use_features(features)) +
  geom_gene(aes(fill=description), size=6) + #position="strand"
  geom_link(color="lightgray", fill="lightgray", alpha=0.3) + #add as soon as alignment is finished!
  #geom_gene_label(aes(label=name), size = 4) +
  geom_gene_tag(aes(label=letter_code), size=5, vjust=-1) +
  scale_fill_manual(values=c("skyblue4","grey","darkseagreen3","darkcyan"))+
  theme(axis.text.x=element_text(size=14)) #+
#scale_fill_distiller(palette = "Blues", direction=1)
p9

# flip manually
p10 <- p9 %>% 
  flip(5:6)# + geom_link()
p10
#p10 <- p10 + theme(legend.position="right")
#p10 <- p10 + theme(legend.text = element_text(size=20))
ggsave("species_plot.png", dpi=700, width=40, height=20, units="cm")
ggsave("species_plot.svg", dpi=700, width=40, height=20, units="cm")
