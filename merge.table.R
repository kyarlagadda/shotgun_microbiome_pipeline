#check if you have dplyr package, install it not
if(!"dplyr"%in%installed.packages()) install.packages("dplyr") 
library(dplyr)

####Explanations----

## this code is designed to merge both metaphlan outputs (i.e. perSample) and paired set (i.e paired) together or within each other
## the files to be merged should be kept in sub-folder named "PerSample" and "Paired" while metaphlan outputs should be kept in "PerSample" folder
## paired set outputs should be kept in "Paired" folder.

## PerSample tables should have their sample names names as file names such as sample.name.txt

## paired outputs should contain paired files named XXX_lineages for the lineage table and XXXX_otus for the otu table. XXX can be any names.
## you can put more than one paired files as long as the naming follows xxx_lineage.txt xxx_otus.txt; yyy_lineage.txt yyy_otus.txt ...etc

## merged files would be produced in the same directory as the source codes.

## all sample names should be unique. Otherwise, the duplicated samples will be named suffix as xx__1 xx__2 ...etc

# setting working directory, this should be where this codefile locates

setwd("set working directory here") ###change this to your own source code file path 


#setting variables
##put the source tables needed to merge in a sub-folder called "PerSample", or set your own below:
merged.files.folder <- "PerSample/" #set folder name for per sample files

merged.files.folder.2 <- "Paired/" #set folder name for paired sample files

lineage.table.name <- "merged_otu_lineage.txt" #name your lineage table here
out.counts.table.name <- "merged_otus.txt" #name your otu counts table here


#identify files to combine
#DO not change the codes below
options(stringsAsFactors = F)

files <- dir(path = merged.files.folder )%>%sort()
sample.names <- gsub(".txt","", files)

iteration <- 1

if(length(files>0)){
  
  for(file in files){
    
    print(paste0("Merging sample ", file, "  ", iteration, " of ", length(files)))
    x1 <- read.table(paste0(merged.files.folder, file), sep="\t", stringsAsFactors = F, skip = 1, 
                     comment.char = "#", header = FALSE)
    header <- readLines(paste0(merged.files.folder, file), n=2)[2]%>%strsplit(split="\t")%>%unlist
    names(x1) <- header
    
    otu.linage <- x1[,1]%>%gsub("\\|", "\\;",.)%>%gsub("[[:alpha:]]__", "",.) #adjust format from raw input
    
    if(iteration==1){
      main.table <- data.frame(Consensus.Lineage = otu.linage, counts = x1$estimated_number_of_reads_from_the_clade)#change to relative_abundance if absolute abundance is not desired
      colnames(main.table)[2] <- sample.names[iteration]
      main.table %<>% 
        group_by(Consensus.Lineage) %>%
        group_by(N=n(), add=TRUE) %>%
        summarise_all(funs(sum)) %>%
        as.data.frame()
      names(main.table)[2] <- paste0(sample.names[iteration],".N")
    }
    
    else{
      x2 <- data.frame(Consensus.Lineage = otu.linage, counts = x1$estimated_number_of_reads_from_the_clade)#change to relative_abundance if absolute abundance is not desired
      colnames(x2)[2] <- sample.names[iteration] 
      x2 %<>% 
        group_by(Consensus.Lineage) %>%
        group_by(N=n(), add=TRUE) %>%
        summarise_all(funs(sum)) %>%
        as.data.frame()
      names(x2)[2] <- paste0(sample.names[iteration],".N")
      
      main.table <- dplyr::full_join(main.table, x2)
      
    }
    iteration <- iteration + 1
  }
  
  #### counts from perSample files
  count.table.1 <- dplyr::select(main.table, -contains(".N"))
  count.table.1[is.na(count.table.1)] <- 0
  
  N.table.1 <- dplyr::select(main.table, Consensus.Lineage, contains(".N"))
}

####read in paired tables

files.2 <- dir(path = merged.files.folder.2 )%>%sort()


lineage.file <- files.2[grep("_lineages", files.2)]
otu.file <- files.2[grep("_otus", files.2)]


iteration <- 1

#Read in counts

if(length(files.2)>0){
  
  for( i in 1:length(lineage.file)){
    
    
    x1 <- read.table(paste0(merged.files.folder.2 ,otu.file[i]), sep="\t", stringsAsFactors = F, skip = 1, 
                     comment.char = "", header = TRUE)
    names(x1)[1] <- gsub("X.","", names(x1)[1], fixed = TRUE)
    
    
    x2 <- read.table(paste0(merged.files.folder.2,lineage.file[i]), sep="\t", header=T,stringsAsFactors = F)
    x2$Consensus.Lineage%<>%
      sapply(function(x)gsub(" ","",x))%>% #remove blanks
      sapply(function(x)ifelse(x=="", "Unassigned", x)) #assign empty cell to unknown
    
    
    x3 <- left_join(x1,x2) #joining two tables by OTU.ID
    
    if(iteration==1){
      main.table.2 <- x3 %>% select(-OTU.ID) %>% 
        group_by(Consensus.Lineage) %>%
        group_by(NUM_LINEAGE=n(), add=TRUE) %>%
        summarise_all(funs(sum)) %>%
        as.data.frame()
      
    }
    if(iteration>1){
      x5 <- x3 %>% select(-OTU.ID) %>% 
        group_by(Consensus.Lineage) %>%
        group_by(NUM_LINEAGE=n(), add=TRUE) %>%
        summarise_all(funs(sum)) %>%
        as.data.frame() 
      main.table.2 <- full_join(main.table.2, x5, by= "Consensus.Lineage", suffix = c("__1","__2") )
      
    }
    iteration <- iteration +1
    
  }
  main.table.2[is.na(main.table.2)] <-0 
  
  count.table.2 <- main.table.2%>%select(-contains("NUM_LINEAGE"))
  count.table.2$Consensus.Lineage <- paste0("Bacteria;", count.table.2$Consensus.Lineage) #paste bacteria
  
  N.table.2 <- main.table.2%>%select(Consensus.Lineage, contains("NUM_LINEAGE"))
  N.table.2$Consensus.Lineage <- paste0("Bacteria;", N.table.2$Consensus.Lineage) #paste bacteria
  #merge 2 count tables
  
}
if(length(files>0)&length(files.2)>0){
  main.table.3 <- full_join(count.table.1, count.table.2, by="Consensus.Lineage")
  N.table.3 <- full_join(N.table.1, N.table.2)
} else if (length(files>0)){
  main.table.3  <- count.table.1
  N.table.3 <- N.table.1
} else {
  main.table.3  <- count.table.2
  N.table.3 <- N.table.2
}



#all.equal(main.table.3$Consensus.Lineage,N.table.3$Consensus.Lineage)
main.table.3[is.na(main.table.3)] <- 0
N.table.3[is.na(N.table.3)] <- 0

N.table.3$N <- rowSums(N.table.3[,-1, drop=F])
outs <- N.table.3%>%select(Consensus.Lineage, N)%>%left_join(main.table.3, by="Consensus.Lineage")%>%mutate(OTU.ID=order(Consensus.Lineage)) 
outs <- outs %>% select(OTU.ID, everything())

#output two seprate files

out.1 <- outs%>%select(OTU.ID, Consensus.Lineage, N)
write.table(out.1, file = lineage.table.name, sep="\t", row.names = FALSE, append= FALSE, col.names = TRUE, quote = FALSE)

out.2 <- outs %>% select(-Consensus.Lineage, -N) 
cat("# QIIME-formatted OTU table\n", file = out.counts.table.name, append= FALSE) #adding special OTU table header
cat("#", file = out.counts.table.name, append = TRUE)
write.table(out.2, file = out.counts.table.name, sep="\t", row.names = FALSE, append=TRUE, col.names = TRUE, quote = FALSE)





# ##codes for share number of individuals
# Consensus.Lineage <- main.table.3$Consensus.Lineage
# 
# main.table.4 <- main.table.3%>%select(-Consensus.Lineage)
# 
# rownames(main.table.4) <- main.table.3$Consensus.Lineage
# 
# main.table.4[is.na(main.table.4)] <- 0
# 
# sample_names <- colnames(main.table.4)%>%gsub("__\\d","",.) #remove __1 and __2 to show actual sample names
# 
# 
# main.table.4.t <-as.data.frame( t(main.table.4))
# main.table.4.t$sample_names <-sample_names
# 
# main.table.5 <- main.table.4.t%>%
#   group_by(sample_names) %>%
#   summarise_all(funs(sum)) %>%
#   as.data.frame()
#   
# rownames(main.table.5) <- main.table.5$sample_names
# main.table.5 <- main.table.5%>%select(-sample_names)
# 
# main.table.6 <- t(main.table.5)%>%as.data.frame()
