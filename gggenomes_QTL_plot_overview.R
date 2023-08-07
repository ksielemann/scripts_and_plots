if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

# install.packages("devtools")
devtools::install_github("thackl/thacklr")
devtools::install_github("thackl/gggenomes")



library(S4Vectors)
library(gggenomes)
library(svglite) 
library(tidyverse)
library(RColorBrewer)


# parse sequence length and some metadata from fasta file
emale_seqs <- read_seq_len("merged.fasta")

#extract QTL region
loci <- tribble(
  ~locus_id,~seq_id, ~start, ~end,
  "1321396-2021946", "U2Bv",1321396, 2021946,
  "1389131-2154734", "KWS2320", 1389131, 2154734)

#plot QTL region
p1 <- gggenomes(seqs=emale_seqs) %>%
  focus(.loci=loci)  + # zoom in on all regions with links
  geom_seq_label() + geom_bin_label() +
  geom_seq()
p1


#Annotate genes
emale_genes <- read_gff("merged_genes.gff") %>%
  rename(feature_id=ID) #%>%                       # we'll need this later
  #mutate(gc_cont=as.numeric(gc_cont))             # per gene GC-content

p2 <- gggenomes(emale_genes, emale_seqs) %>%
  focus(.loci=loci)  + 
  geom_seq() + geom_seq_label() + geom_bin_label() +
  #geom_gene_label(aes(label=Name), size=4, angle=90) +
  geom_gene(aes(fill=strand), position="strand")
p2


#Compare genome syteny
emale_links <- read_paf("alignment_chr5.paf")

emale_links_filtered <- filter(emale_links, (start<2154734 & start>1321396)|(start2<2154734 & start2>1321396)|(end<2154734 & end>1321396)|(end2<2154734 & end2>1321396))

emale_links_filtered <- filter(emale_links_filtered, AS > 2000)

emale_links_filtered

p4 <- gggenomes(emale_genes, emale_seqs, links=emale_links_filtered) %>%
  focus(.loci=loci)  + 
  geom_seq() + geom_bin_label() +
  #geom_feature(size=5, data=use_features(features)) +
  geom_gene(aes(fill=strand), position="strand") +
  geom_link()

#p4 <- p4 %>% flip_bins(3:5)
p4




#(Add funtional annotation/description (Swissprot BLAST)/Interpro domains)
annot <- read_tsv("merged_annotation_refined_manually_selected.txt", col_names=NULL)
colnames(annot) <- c("Name","description")

genes_and_info_and_description <- left_join(genes_and_info, annot, by=c("Name"))

p6 <- gggenomes(genes_and_info_and_description, emale_seqs, links=emale_links_filtered) %>% #rgas
  focus(.loci=loci)  + 
  geom_seq() + geom_bin_label(size=10) + 
  #geom_feat(size=15) +
  #geom_feature(size=5, data=use_features(features)) +
  geom_gene(aes(fill=description), size=5) +#, position="strand") + #position="strand"
  #geom_gene_label(aes(label=description), size = 4) +
  geom_link(color="lightgray", fill="lightgray", alpha=0.3) +
  scale_fill_manual(values=c("skyblue4","grey","darkseagreen3","darkcyan"))+
  theme(axis.text.x=element_text(size=14)) #+
  #scale_fill_distiller(palette = "Blues", direction=1)
p6
ggsave("overview.png", dpi=700, width=40, height=20, units="cm")
ggsave("overview.svg", dpi=700, width=40, height=20, units="cm")

