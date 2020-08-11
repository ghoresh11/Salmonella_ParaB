library(dplyr)
library(ggplot2)
require(scales)
library(utils)
library(RColorBrewer)
library(here)

setwd(here())  ## change your working directory here


## make it possible to try different values 
args = commandArgs(trailingOnly=TRUE)
max_val = args[1]
min_val = args[2]
outdir = paste("baps",max_val,min_val, sep = "_")
max_val = as.numeric(max_val)/100
min_val = as.numeric(min_val)/100

dir.create(file.path(".", outdir), showWarnings = FALSE)
dir.create(file.path(".", outdir, "plots"), showWarnings = FALSE)

## Step 1: Read the input files
presence_absence_tab = read.table("gene_presence_absence.Rtab", sep = "\t", comment.char = "", stringsAsFactors = F, header = F, row.names = 1)
presence_absence_csv = read.table("gene_presence_absence.csv", sep = ",", comment.char = "", stringsAsFactors = F, header = F)
md = read.table("global_metadata_baps-only.csv", sep = ",", stringsAsFactors = F, header = T, comment.char = "")
md = md[,c(1,3)]
colnames(md) = c("Isolate","Biotype")

## Step 2: Divide the Rtab file into typhi and java and find genes which are core in one but missing in the other,

## Im writing a function that gets the biotype as the input, then it does the exact set of commands written above
get_gene_frequencies <- function(biotype) {
  indexes = which(presence_absence_tab[1,] %in% md$Isolate[md$Biotype == biotype]) ## get the rows that are the specific biotype
  gene_counts = presence_absence_tab[,indexes] ## get all the relveant columns
  gene_counts = gene_counts[-1,] # remove first row with the strain names
  gene_counts_numeric <- mutate_all(gene_counts, function(x) as.numeric(as.character(x))) ## annoyingly it's not numeric, so just need to convert it
  rownames(gene_counts_numeric) = rownames(gene_counts) ## give back the gene names to each row
  colnames(gene_counts_numeric) = presence_absence_tab[1,][indexes] ## sanity check, adding the strain names
  ## sanity check for ParaB -> should return 1
  ## print(gene_counts_numeric[which(rownames(gene_counts_numeric) == "group_944"),which(colnames(gene_counts_numeric) == "11893_2#70")])
  gene_frequencies = rowSums(gene_counts_numeric) / length(indexes) ## the frequency of each gene is the sum of the rows divided by the size of the  group
  return(gene_frequencies)
}

## now it's very easy to run on both:
one_gene_frequencies = get_gene_frequencies("1")
two_gene_frequencies =  get_gene_frequencies("2")
## if you print these two you see it gives you two vectors with the freq of all genes in each group and the gene names

## output the frequency of each gene in both
write.table(data.frame(gene = names(one_gene_frequencies),
                       one = one_gene_frequencies,
                       two = two_gene_frequencies), file = "freqs_per_baps.csv",
            col.names = T, sep = ",",row.names = F, quote = F)

## calculate the difference in frequency between the two
diff_freq = one_gene_frequencies - two_gene_frequencies

## define specific to java or parab
specific_to_one = names(one_gene_frequencies)[which(one_gene_frequencies >= max_val & two_gene_frequencies <= min_val)]
specific_to_two = names(two_gene_frequencies)[which(two_gene_frequencies >= max_val & one_gene_frequencies <= min_val)]


## At this point we have two lists of genes which are specific to Jav and ParatyphiB

## Step 3: For each gff file, convert the gene name to the name in the corresponding GFF file
## If it's a Java isolate, see where all the Java specific genes are on the Java genome
## Do this for all GFFs and summarise everything in one place somehow

## initiate a summary of all the distances between every two genes for all the genomes, with every combination of two genes
init_summary_df <- function(biotype, specific_to_biotype) {
  curr_summary = data.frame(cbind(t(combn(x = specific_to_biotype, m = 2)), 
                                  data.frame(matrix(nrow = choose(length(specific_to_biotype), 2), 
                                                    ncol = length(which(md$Biotype == biotype))))), stringsAsFactors = F)
  colnames(curr_summary) = c("geneA","geneB",md$Isolate[md$Biotype == biotype])
  return(curr_summary)
}

summary_one = init_summary_df("1", specific_to_one)
summary_two = init_summary_df("2", specific_to_two)


## special function to read a GFF file
## because comments are "##" but there's a hash in the lane ID I needed to hack a solution to change all "##" to "%"
## then use "%" as the comment character. Also then I filtered to only include lines with a CDS, so the sequence at the end 
## is removed, then I apply more manipulations to make it easier to work with
read_gff<- function(file) {
  clean.lines <- sub(paste0("##"), "%", readLines(file)) ## rewrite "##" as "%"
  df = read.csv(text = paste(clean.lines, collapse = "\n"), comment.char = "%", sep= "\t", stringsAsFactors = F,header = F,  quote = "") # read the file, ignore lines with "%"
  ## I'm renaming the column names so the code is cleaer
  colnames(df) = c("contig","tool","type","start","stop","nothing","strand","nothing2","desc")
  df = df[-which(df$type != "CDS"),] ## only CDS
  ## the name of the gene in this file is the in the desc column, but it needs to be manipulated to get the exact name
  df$gene_id =  sapply(X = strsplit(x = df$desc, split = ";", fixed = T ), FUN = head, n = 1) ## here are all the gene IDs, you split the column by ";" and get the first item
  df$gene_id = sapply(X = df$gene_id, FUN=gsub, pattern = "ID=",replacement = "") ## remove the "ID="
  return(df)
}

gff_files = list.files("gffs/", full.names = T)

for (curr_gff_file in gff_files) {
  print(curr_gff_file)
  ## read the GFF -> I had to write a function because the comment character is 2 hashes not one, and you can't do that
  ## in the built in read.table function
  genome_name = gsub(x = basename(curr_gff_file), pattern = ".gff", replacement = "") ## getting the current name of the genome
  if (!genome_name %in% md$Isolate) { ## if for any reason the genome isn't in the metadata file, move on
    print(paste("Could not locate ---", genome_name, "--- in the metadata file, skipping!", sep = ""))
    next
  }
  
  curr_biotype = md$Biotype[md$Isolate == genome_name] ## getting the current biotype
  curr_df = read_gff(curr_gff_file) ## read the gff file
  
  ## if it's a Java isolate, find all the java specific genes, otherwise all the paraB specific
  if (curr_biotype == "1") {
    curr_genes = specific_to_one
    curr_summary = summary_one
  } else {
    curr_genes = specific_to_two
    curr_summary = summary_two
  }
  
  ## now I can get the IDs of these genes for this isolate in from the CSV file
  isolate_column = which(presence_absence_csv[1,] == genome_name)
  gene_rows = which(presence_absence_csv[,1] %in% curr_genes)
  geneIDs_in_curr_isolate = presence_absence_csv[gene_rows, isolate_column]
  ## if you print geneIDs_in_curr_isolate you'll see it's exactly what we need, Java specific genes that are missing in ParaB
  ## now I'll add this info to the GFF df
  curr_df$specific = rep(0,dim(curr_df)[1])
  curr_df$specific[curr_df$gene_id %in% geneIDs_in_curr_isolate] = 1
  ## it's not showing presence/absence in this isolate, it's showing which genes are in this isolate but missing in the OTHER biotype
  
  ## add the actual name of the gene to the dataframe, needed to be added to the summary
  curr_df$gene_name = presence_absence_csv$V1[match(curr_df$gene_id, presence_absence_csv[,isolate_column])]
  
  ## Now we only really care about the contigs where the specific genes are, so I will filter the dataframe to only include those contigs
  interesting_contigs = unique(curr_df$contig[which(curr_df$specific == 1)]) ## list of contigs
  
  ## filter the gff dataframe to only have these contigs (could be 1 or multiple...)
  curr_df = curr_df[curr_df$contig %in% interesting_contigs,]
 
  ## the most fun part -> plot the results for a single GFF file...
  curr_df$diff = diff_freq[match(curr_df$gene_name,names(diff_freq))]
  if (curr_biotype == "2") {
    curr_df$diff = -curr_df$diff
  }
  curr_df$specific = curr_df$diff
  other_bio = "1"
  if (curr_biotype == "1") {
    other_bio = "2"
  }
  
  
  if (length(interesting_contigs) > 0) {
    # ## this needs to be fixed in case there are multiple contigs
    curr_df$index = 1:dim(curr_df)[1]
    curr_df$index = factor(curr_df$index, curr_df$index)
    p = ggplot(curr_df, aes(x = index, y = specific)) + geom_bar(stat = "identity", color = "black", size = 0.1) +
      facet_grid(~contig, scales = "free_x") + 
      xlab("Index on contig") + ylab(paste("Frequency(", curr_biotype, ")-Frequency(", other_bio, ")",sep = "")) + 
      theme_classic(base_size = 14) +
      scale_y_continuous(expand = c(0.05,0,0.05,0), breaks = c(0,0.95,1)) + 
      ggtitle(paste(genome_name, "(",curr_biotype,")", sep ="")) + 
      scale_x_discrete(labels = curr_df$gene_name) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      geom_hline(yintercept = c(0,0.95,1), lty = 2)
    #  ggsave, not so important at the moment
    width = min(40, dim(curr_df)[1]/20)
    ggsave(p, file = file.path(outdir, "plots",paste(genome_name,".pdf",sep="")),height = 6, width = width)
  }
  ## fill out the details in the summary dataframe
  curr_column = which(colnames(curr_summary) == genome_name)
  for (i in 1:dim(curr_summary)[1]) {
    ## Check if both genes are present in this isolate's GFF file
    if ((!curr_summary$geneA[i] %in% curr_df$gene_name)  || (!curr_summary$geneB[i] %in% curr_df$gene_name)) {
      curr_summary[i, curr_column] = "One_missing"
      next
    } ## check if the two genes are on the same contig
    if (curr_df$contig[which(curr_df$gene_name == curr_summary$geneA[i])] !=
        curr_df$contig[which(curr_df$gene_name == curr_summary$geneB[i])]) {
      curr_summary[i, curr_column] = "Diff_contigs"
      next
    }
    ## on same contig - measure the difference between them in the dataframe in terms of index
    curr_summary[i,curr_column] = abs(which(curr_df$gene_name == curr_summary$geneA[i]) - which(curr_df$gene_name == curr_summary$geneB[i]))
  }
  
  ## don't forget to update the relevant results
  if (curr_biotype == "1") {
    summary_one = curr_summary 
  } else {
    summary_two = curr_summary 
  }
}

write.table(summary_one, file.path(outdir,"summary_specific_one.csv"), col.names = T, row.names = F, quote = F, sep = ",")
write.table(summary_two, file.path(outdir,"summary_specific_two.csv"), col.names = T, row.names = F, quote = F, sep = ",")
