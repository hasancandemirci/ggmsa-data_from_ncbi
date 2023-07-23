
############### Hasan Can Demirci #################

############## GGMSA by taking data from NCBI #####


# necessary libraries

library(ggmsa)
library(compbio4all)
library(rentrez)
library(Biostrings)
library(msa)


# creating data with accession numbers

IFT140_table <- c("NP_055529", "Human-NP_055529", "IFT140",
                  "XP_016784643","Chimp-XP_016784643","IFT140",
                  "XP_001089057","Macaca-XP_001089057","IFT140",
                  "XP_002697959","Cow-XP_002697959","IFT140",
                  "NP_598887","Mouse-NP_598887","IFT140",
                  "XP_006246116","Rat-XP_006246116","IFT140",
                  "NP_001116497","Frog-NP_001116497","IFT140",
                  "XP_695732","Zebrafish-XP_695732","IFT140",
                  "NP_995608","Fruitfly-NP_995608","IFT140",
                  "NP_506047","C. elegans-NP_506047",   "IFT140",
                  "XP_042921850",  "C. reinhardtii-XP_042921850", "IFT140",
                  "XP_001020653", "T. thermophila-XP_001020653", "IFT140",
                  "XP_001347023", "Paramecium tetraurelia-XP_001347023", "IFT140")


# convert the vector to matrix using matrix()

IFT140_table_matrix <- matrix(IFT140_table, byrow = T, nrow = 13)

# convert the matrix to a dataframe using data.frame()
# stringsAsFactors: used to determine the data change into vector or factor

IFT140_table <- as.data.frame(IFT140_table_matrix, stringsAsFactors = F)

# name columns of dataframe using names() function

colnames(IFT140_table) <- c("ncbi.protein.accession", "species", "gene.name")

# get fasta sequences

IFT140s_list <- entrez_fetch_list(db = "protein", 
                                  id = IFT140_table$ncbi.protein.accession, 
                                  rettype = "fasta")

# Go through list and remove FASTA header from each sequence
# fasta_cleaner : provide reading fasta format file as DNA, RNA, or Protein depends on your choice

IFT140s_vector <- sapply(IFT140s_list, fasta_cleaner, parse = F)

# name the vector

names(IFT140s_vector) <- IFT140_table$species

# create AAStringSet for protein alignment

IFT140s_vector_ss <- AAStringSet(IFT140s_vector)

# perform MSA (visualization)

alignment <- msa(IFT140s_vector_ss, "ClustalW")

# change class of alignment to AAMultipleAlignment

class(alignment) <- "AAMultipleAlignment"

# plot the alignment using ggmsa

ggmsa::ggmsa(alignment, start = 800, end = 840, char_width = 0.5, seq_name =T) +
  geom_seqlogo() +
  geom_msaBar()



