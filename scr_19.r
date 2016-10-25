setwd("C:\\Users\\QCRI-L002\\Desktop\\Collect_all_in_all\\V_hg19///")
library("biomaRt")
ensembl37 <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL') #select the ensembl version that contains hg19 (GRch37) 
ensembl <- useDataset("hsapiens_gene_ensembl",mart= ensembl37) #select human genes
#attributes = listAttributes(ensembl)
#filters = listFilters(ensembl)
#read the 8008 genes (all genes we have in the single data)
gene_id = read.csv("gene_id.txt",header = F) #upload the values of our filter
gene_id <- as.matrix(gene_id)
#The needed attributes
attributes <- c("ensembl_gene_id", "chromosome_name","hgnc_symbol","band","start_position","end_position", "external_gene_name") #select the attributes
#In the first version, we will find the gene info using the gene name 
filter <- "hgnc_symbol" #select the filter


geneInfo = data.frame() #creat an empty dataframe
#the followinf while loop to handel the gene_id arry 500 by 500 (500 based on the advise on ensmble website beacuse whene we call the getBM function by the whole gene_id array we will miss some values)
finished = FALSE
limit = 500
k=0

while(!finished)
{
  if((length(gene_id)-(k*limit))>limit)
  {
    gene_id_tmp = gene_id[c((limit*k+1):(limit*(k+1)))]
    k=k+1
  }
  else
  {
    gene_id_tmp = gene_id[c((limit*k+1):length(gene_id))]
    finished = TRUE
  }
  geneInfo <- rbind(geneInfo,getBM(attributes = attributes, filters = filter, values = gene_id_tmp, mart = ensembl))
}
chr_names=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
#remove the umbigous chromosome names
geneInfo.f <- geneInfo[(geneInfo$chromosome_name %in% chr_names), ]
#find the missing genes
missing_genes = gene_id [!(gene_id %in% geneInfo.f$hgnc_symbol)]
#read the file contains the gene names associated with the ensemble id
gene_ids_ens <- read.table("ens_exist.txt",header = T) #upload the values of our filter
rownames(gene_ids_ens) = gene_ids_ens$gene_name
#find the gene names which are not exist in the current gene Info
missing_genes_ens = gene_ids_ens[missing_genes[missing_genes %in% gene_ids_ens$gene_name],]
rownames(missing_genes_ens) = missing_genes_ens$Ensemble_id

#Now we will find the gene Info for the missing genes using the ensemble ID
filter <- "ensembl_gene_id"
gene_id_1 = missing_genes_ens$Ensemble_id
geneInfo_1 <- getBM(attributes = attributes, filters = filter, values = gene_id_1, mart = ensembl)
#remove the umbigous chromosome names
geneInfo_1.f <- geneInfo_1[(geneInfo_1$chromosome_name %in% chr_names), ]
rownames(geneInfo_1.f) = geneInfo_1.f$ensembl_gene_id
missing_genes_ens_exist = missing_genes_ens[rownames(missing_genes_ens) %in% rownames(geneInfo_1.f),]   
geneInfo_mis = geneInfo_1.f[rownames(missing_genes_ens_exist),]
geneInfo_mis$hgnc_symbol = missing_genes_ens_exist$gene_name

rownames(geneInfo_mis) = geneInfo_mis$ensembl_gene_id
geneInfo_all = rbind(geneInfo.f, geneInfo_mis)
geneInfo_all.f = geneInfo_all[!(duplicated(geneInfo_all$hgnc_symbol)),]
rownames(geneInfo_all.f) = geneInfo_all.f$hgnc_symbol
str(geneInfo_all.f)

#sorting the geneInfo
geneInfo_all.f_num <- geneInfo_all.f[!geneInfo_all.f$chromosome_name %in% c("X", "Y"), ]
geneInfo_all.f_num$chromosome_name = as.integer(geneInfo_all.f_num$chromosome_name)
geneInfo_all.f_num_order = geneInfo_all.f_num[ order(geneInfo_all.f_num[,"chromosome_name"], geneInfo_all.f_num[,"start_position"], geneInfo_all.f_num[,"end_position"]), ]
geneInfo_all.f_x <- geneInfo_all.f[geneInfo_all.f$chromosome_name %in% c("X"), ]
geneInfo_all.f_x_order = geneInfo_all.f_x[ order(geneInfo_all.f_x[,"chromosome_name"], geneInfo_all.f_x[,"start_position"], geneInfo_all.f_x[,"end_position"]), ]
geneInfo_all.f_y <- geneInfo_all.f[geneInfo_all.f$chromosome_name %in% c("Y"), ]
geneInfo_all.f_y_order = geneInfo_all.f_y[ order(geneInfo_all.f_y[,"chromosome_name"], geneInfo_all.f_y[,"start_position"], geneInfo_all.f_y[,"end_position"]), ]
#bind the gInfo and discarding the y chromosome genes
geneInfo_all.f_sorted = rbind(geneInfo_all.f_num_order, geneInfo_all.f_x_order)
gInfo = geneInfo_all.f_sorted[,c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band")]
colnames(gInfo) = c("gene_id","chr", "starting_position","ending_position", "band")
gInfo_hg19 = gInfo
save (gInfo_hg19,file = "gInfo_hg19.RData")
