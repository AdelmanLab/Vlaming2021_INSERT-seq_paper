### Function to keep only inserts for which I have data from all indicated assays.
# Input:
  # data_table = dataframe with at least the columns indicated under 'data_columns
  # data_columns= vector of character strings with names of columns where you want data to be present. Defaults to Chr-RNA, Total-RNA, sortseq
# Output:
  # Same dataframe as inputted, but without the rows where one of the data types was missing
keep_complete <- function(data_table,data_columns=c('runonRNA_average','totalRNA_average','sortseq_average')) {
  data_table[is.finite(rowSums(data_table[,data_columns])),]
}

# Function to keep genomic regions only, and re-annotate a bit
# Input:
  # data_table = dataframe with at least the columns category and window and some data columns
  # complete_only = whether to filter out rows where some information is missing, using keep_complete function
  # data_forcomplete = vector of character strings that will be passed on as data_columns parameter for the keep_complete function
# Output:
  # Dataframe with a an extra 'group' column compared to before and fewer rows, because anyting that is not a genomic region (and maybe incomplete) has been filtered out. 
keep_genomicregions <- function(data_table,
                                complete_only=F,
                                data_forcomplete=c('runonRNA_average','totalRNA_average','sortseq_average')) {
  genomicregions_table <- data_table[which(data_table$category %in% 
                                             c('mRNAs','lincRNAs','uaRNAs','eRNAs','eRNAs_reallyin_SE',
                                               'eRNAs_in_SEs','mRNA_ends',"randomized_control")),]
  genomicregions_table$window[genomicregions_table$window=="6-179"] <- "6to179"
  genomicregions_table$window[genomicregions_table$window=="160-333"] <- "160to333"
  genomicregions_table$group <- paste(genomicregions_table$category, genomicregions_table$window, sep="_")
  genomicregions_table$group[genomicregions_table$group=="randomized_control_NA"] <- "randomized_control"
  genomicregions_table$group[genomicregions_table$group=="eRNAs_in_SEs_6to179"] <- "eRNAs_6to179"
  genomicregions_table$group[genomicregions_table$group=="eRNAs_reallyin_SE_6to179"] <- "eRNAs_6to179"
  genomicregions_table$group[genomicregions_table$group==
                               "mRNA_ends_-86 to +87 around PAS (or center of PASs)"]  <- "mRNA_ends"
  genomicregions_table$group <- factor(genomicregions_table$group,  
                                       levels=c('randomized_control','mRNAs_6to179','mRNAs_160to333','lincRNAs_6to179',
                                                'lincRNAs_160to333','uaRNAs_6to179','uaRNAs_160to333','eRNAs_6to179','mRNA_ends'))         
  if(complete_only==T){ keep_complete(genomicregions_table,data_columns=data_forcomplete) 
    } else {genomicregions_table}
}


# Function to go from long table to wide table, and write to txt file
# Input:
  # table = table in long format
  # data = name of column that has the data I want to look at
  # grouping_parameter = name of column(s) that I want to separate the data by (character string, or vector of strings)
  # limits = vector of length 2, giving a lower and upper limit. Outliers beyond those will be changed to limit values. Defaults to NULL -> no limit used.
  # subset = something to add to file name to indicate what data is in 'table'
  # keep = whether to output the wide table in R, defaults to FALSE
# Output
  # A txt file with the wide table containing 1 type of data divided over columns based on variables in grouping_parameter
  # If keep ==T, the wide table that is also written to file
widetable_makeandsave <- function(table,data,grouping_parameter,limits=NULL,subset,keep=F,save=T) {
  if(!is.null(limits)) {
    table[which(table[,data]<limits[1]),data] <- limits[1]
    table[which(table[,data]>limits[2]),data] <- limits[2]
    subset <- paste0('wlimits_',subset)
    }
  if(length(grouping_parameter)==1) {
    table$number <- sapply(1:nrow(table),function(rownr){
      group <- table[rownr,grouping_parameter]
      sum(table[1:rownr,grouping_parameter]==group)
    })
    table_wide <- reshape(table[,c(data,grouping_parameter,'number')],
                          timevar=grouping_parameter,'idvar'="number", direction = "wide",sep = "___")
    rownames(table_wide) <- table_wide$number
    table_wide <- table_wide[,-1]
    colnames(table_wide) <- sapply(colnames(table_wide),
                                   function(x){ strsplit(x,split="___")[[1]][2] })
    name <- paste0('_',subset,'_by',grouping_parameter)
    if(save==T){write.table(table_wide,
                paste0('PrismTables/',data,name,'.txt'),
                sep='\t',quote=F,row.names=F)}
    table_wide
  } else{
    table[,'combi'] <- paste(table[,grouping_parameter[1]],table[,grouping_parameter[2]],sep="_")
    table$number <- sapply(1:nrow(table),function(rownr){
      group <- table[rownr,'combi']
      sum(table[1:rownr,'combi']==group)
    })
    table_wide <- reshape(table[,c(data,'combi','number')],
                          timevar='combi','idvar'="number", direction = "wide",sep = "___")
    rownames(table_wide) <- table_wide$number
    table_wide <- table_wide[,-1]
    colnames(table_wide) <- sapply(colnames(table_wide),
                                   function(x){ strsplit(x,split="___")[[1]][2] })
    name <- paste0('_',subset,'_by',grouping_parameter[1],"AND",grouping_parameter[2])
    if(save==T){write.table(table_wide,
                paste0('PrismTables/',data,name,'.txt'),
                sep='\t',quote=F,row.names=F)}
    if(keep==F) {'done'
    } else{table_wide}
  }
}


# Function to run widetable_makeandsave multiple times in one command, for different data columns.
# Input:
  # x = table in long format (defaults to 'alldata')
  # data_types = vector of column names for which you want to make a wide table, defaults to Chr-RNA, Total-RNA, Sortseq
  # grouping = name of column(s) that I want to separate the data by (character string, or vector of strings)
  # boundaries = list of same length as data_types, with each element a vector of length 2, giving a lower and upper limit per data_type. Outliers beyond those will be changed to limit values. Defaults to NULL -> no limits used.
  # subset_for_name = something to add to file name to indicate what data is in 'x'
  # print_tables = whether to output the wide tables in R, defaults to FALSE
# Output
  # One txt file per data type, as made by the widetable_makeandsave function
  # If print_tables is TRUE, a list of wide tables that are also written to files
widetable_makeandsave_multi <-  function(x=alldata,
                                         data_types=c('runonRNA_average','totalRNA_average','sortseq_average'),
                                         grouping,boundaries=NULL,subset_for_name,print_tables=F,save_table=T) {
  list_tables <- lapply(1:length(data_types),
                        function(i,y=x,dt=data_types,grouping_par=grouping,bound=boundaries,
                                 subsetname=subset_for_name, keep_t=print_tables) {
    widetable_makeandsave(table=y,dt[i],
                          grouping_parameter=grouping_par,limits=bound[[i]],
                          subset=subsetname, keep=keep_t,save=save_table)})
  if(print_tables==F) {'done'
  } else{names(list_tables) <- data_types
  list_tables}
}

widetable_makeandsave_combi <-  function(x=alldata,
                                         data_types=c('runonRNA_average','totalRNA_average','sortseq_average'),
                                         grouping,grouping_name=NULL,subset_for_name,print_tables=F,save_table=T) {
  list_tables <- lapply(unique(x[,colnames(x)==grouping]),function(i) {
                          table <- x[x[,colnames(x)==grouping]==i,data_types]
                          if(save_table==T){write.table(table,
                                                  paste0('PrismTables/',subset_for_name,'_',i,'_from',grouping_name,'.txt'),
                                                  sep='\t',quote=F,row.names=F)}
                          table
                        })
  if(print_tables==F) {'done'
  } else{names(list_tables) <- unique(x[,colnames(x)==grouping])
  list_tables}
}


get_write_fastafiles <- function(ids,name,save=T) {
  seqs <- full_library$sequence_in_read[full_library$unique_id %in% ids]
  seqs <- DNAStringSet(seqs)
  names(seqs) <- full_library$unique_id[full_library$unique_id %in% ids]
  if(save==T) {writeXStringSet(seqs,paste0("output/",name,".fa"),width=175)}
  seqs
}


averagebcsandrepl <- function(x=data_introns_bcaverages,intron_list,data,data_type,repl=T) {
  intron_mut <- paste(x["name"],x["mutstatus"],sep="_")
  rows_maintable <- which(intron_list==intron_mut)
  data_average <- median(unlist(data[rows_maintable,grep(data_type,colnames(data))]),na.rm=T)
  repl_data <- sum(!is.na(data[rows_maintable,grep(data_type,colnames(data))]))
  if(repl==T) {c(repl_data,data_average)} else {data_average}
} 

# Big function to average barcoded intron data, select which introns are spliced originally, and which mutants effectively abrogating splicing.
# Input:
# data_introns: dataframe with one row for every intron-containing insert (rows with all NAs for data removed), and columns that indclude unique_id, name, fasta_names, category, barcode, and data columns for run-on, chrRNA, totalRNA,  sortseq and splicing_efficiency, all replicates separately.
# Output: Table with per average data per intron. If intron is present with multiple barcodes, average across all of them
function_intron_tables <- function(data_introns=allintrondata_noNA,
                                   data_types=c("runonRNA"="runon_._total","totalRNA"="totalRNA_._total","sortseq"="sortseq_screen"),
                                   splicingefficiencies=list("runonRNA"=c("runonRNA_._splicingeff",0.03,0.3),"totalRNA"=c("totalRNA_._splicingeff",0.1,0.7))) {
  
  # Add intron info on the intron subtable: 1 column for type of intron, 1 for mutation status
  data_introns$introntype <- NA
  introntype_options <- c('first'='first_introns','woPAS'='introns_woPAS','wPAS'='introns_wPAS')
  for(intron_type in introntype_options) {
    data_introns$introntype[grep(intron_type,data_introns$category)] <- 
      names(introntype_options)[introntype_options==intron_type]
  }
  data_introns$mutstatus <- NA
  mutstatus_options <- c('original'='original','mutant_5ss'='5ss','mutant_3ss'='3ss',
                         'ProudfootPAS_inserted'='ProudfootPAS',
                         'ProudfootPAS_mutant5ss'='ProudfootPAS_5ss',
                         'ProudfootPAS_in_scrambled'='ProudfootPAS_in_scrambled')
  for(mutstatus in mutstatus_options) {
    data_introns$mutstatus[grep(mutstatus,data_introns$category)] <- 
      names(mutstatus_options)[mutstatus_options==mutstatus]
  }
  data_introns$group <- paste(data_introns$introntype,data_introns$mutstatus,sep="_")
  
  ### Separate no bc and bc, average barcodes, then merge into one table again
  ## First get introns w/o barcode
  data_introns_nobc <- data_introns[is.na(data_introns$barcode),]
  data_introns_nobc$name[data_introns_nobc$mutstatus=="original"] <-
    paste0(data_introns_nobc$name[data_introns_nobc$mutstatus=="original"],"_nobc")
  
  # Calculate the mean per data type
  for (i in names(data_types)) {
    data_introns_nobc[paste(i,"average",sep="_")] <-
      rowMeans(data_introns_nobc[,grep(data_types[i],colnames(data_introns_nobc))],na.rm=T)
    }
  for (i in names(splicingefficiencies)) {
    data_introns_nobc[paste("splicingeff",i,"average",sep="_")] <-
      rowMeans(data_introns_nobc[,grep(splicingefficiencies[[i]][1],colnames(data_introns_nobc))],na.rm=T)
  }
  # Set replicate number by counting the non-NAs per data type
  for (i in names(data_types)) {
    data_introns_nobc[paste('repl',i,sep="_")] <- NA
    data_introns_nobc[paste('repl',i,sep="_")] <-
      apply(data_introns_nobc[,grep(data_types[i],colnames(data_introns_nobc))], MAR=1,function(x){sum(!is.na(x))})
  }
  
  # Remove the non-averaged data and barcode columns from the table and put in the desired order using new group column
  data_introns_nobc <- data_introns_nobc[,-grep('unique_id|1|2|3|4|5|barcode',colnames(data_introns_nobc))]
  data_introns_nobc$group <- factor(data_introns_nobc$group, levels=c(
    "first_original","woPAS_original","wPAS_original",
    "wPAS_ProudfootPAS_in_scrambled"))
  data_introns_nobc <- data_introns_nobc[order(as.numeric(data_introns_nobc$group)),]
  
  
  ## Get introns with barcodes
  data_introns_withbc <- data_introns[!is.na(data_introns$barcode),]
  intron_mut_maintable <- paste(data_introns_withbc$name,
                                data_introns_withbc$mutstatus,sep="_")
  data_introns_bcaverages <- data_introns_withbc[!(duplicated(intron_mut_maintable)),]
  for (i in names(data_types)) {
    data_introns_bcaverages[,c(paste(i,"average",sep="_"),paste('repl',i,sep="_"))] <- 
      t(apply(data_introns_bcaverages,MAR=1,averagebcsandrepl,
              intron_list=intron_mut_maintable,data=data_introns_withbc,data_type=data_types[i],repl=T))
  }
  for (i in names(splicingefficiencies)) {
    data_introns_bcaverages[,paste("splicingeff",i,"average",sep="_")] <- 
      apply(data_introns_bcaverages,MAR=1,averagebcsandrepl,
              intron_list=intron_mut_maintable,data=data_introns_withbc,data_type=splicingefficiencies[[i]][1],repl=F)
  }
  # Remove the non-averaged data and barcode columns from the table
  data_introns_bcaverages <- data_introns_bcaverages[,match(colnames(data_introns_nobc),colnames(data_introns_bcaverages))]
  
  # Add '_bc' to originals (to distinguish from non-bc ones), and put in the desired order using new group column
  data_introns_bcaverages$group[grep('original',data_introns_bcaverages$group)] <- 
    paste(data_introns_bcaverages$group[grep('original',data_introns_bcaverages$group)],"bc",sep="_")
  data_introns_bcaverages$group <- factor(data_introns_bcaverages$group, levels=c(
    "woPAS_original_bc","wPAS_original_bc",
    "woPAS_mutant_5ss","wPAS_mutant_5ss",
    "woPAS_mutant_3ss","wPAS_mutant_3ss",
    "wPAS_ProudfootPAS_inserted",
    "wPAS_ProudfootPAS_mutant5ss") )
  data_introns_bcaverages <- data_introns_bcaverages[order(as.numeric(data_introns_bcaverages$group)),]
  
  ## Combine bc and no bc, make new factor of group column to sort on
  data_introns_averages <- rbind(data_introns_nobc,data_introns_bcaverages)
    data_introns_averages$group <- factor(data_introns_averages$group, levels=c(
    "first_original","woPAS_original","wPAS_original",
    "woPAS_original_bc","wPAS_original_bc",
    "woPAS_mutant_5ss","wPAS_mutant_5ss",
    "woPAS_mutant_3ss","wPAS_mutant_3ss",
    "wPAS_ProudfootPAS_inserted",
    "wPAS_ProudfootPAS_mutant5ss",
    "wPAS_ProudfootPAS_in_scrambled") )
  data_introns_averages <- data_introns_averages[order(as.numeric(data_introns_averages$group)),]
  
  # Add a unique identifier column, remove 'CDS_start' column
  data_introns_averages <- cbind("name_mutstatus"=as.character(paste(data_introns_averages$name,data_introns_averages$mutstatus,sep="_")),
                                 data_introns_averages[,-grep("CDS_start",colnames(data_introns_averages))])
  # Add 'splicing_group' columns
  for (i in names(splicingefficiencies)) {
      data_introns_averages[,paste("splicing_group",i,sep="_")] <- NA
      data_introns_averages[which(data_introns_averages[,paste("splicingeff",i,"average",sep="_")] < as.numeric(splicingefficiencies[[i]][2])),
                            paste("splicing_group",i,sep="_")] <- 
                                        paste0("<",splicingefficiencies[[i]][2])
      data_introns_averages[which(data_introns_averages[,paste("splicingeff",i,"average",sep="_")] >= as.numeric(splicingefficiencies[[i]][2]) &
                              data_introns_averages[,paste("splicingeff",i,"average",sep="_")] <= as.numeric(splicingefficiencies[[i]][3])),
                            paste("splicing_group",i,sep="_")] <- 
                                        paste0(splicingefficiencies[[i]][2],"-",splicingefficiencies[[i]][3])
      data_introns_averages[which(data_introns_averages[,paste("splicingeff",i,"average",sep="_")] > as.numeric(splicingefficiencies[[i]][3])),
                            paste("splicing_group",i,sep="_")] <- 
                                        paste0(">",splicingefficiencies[[i]][3])
  }
  
  data_introns_averages
}

# Write function that takes dataframe with category/name/fasta_names columns + indicated data columns
# data indicates which data columns should be normalized (defauls to runonRNA, totalRNA, sortseq)
# Function first calculates for each background and each data type what the controls score
# It then finds which background each insert is in and divides scores by control scores
normalize_to_scr <- function(table, data=c("runonRNA_average","totalRNA_average","sortseq_average"),
                             norm_method=c('FC','FC','shift')) {
  bgscores <- sapply(paste0('scrambled_background',1:5),function(bg,x=table,datacols=data){
    scr_rows <- grep(bg,x$category)
    apply(x[grep(bg,x$category),datacols],MAR=2,median,na.rm=T)
  })
  normalized_scores <- t(apply(table,MAR=1,function(x,cols=data,method=norm_method,bg=bgscores) {
    cell_w_bg <- x[grep('background',x)][1]
    bgloc <- str_locate(cell_w_bg, 'background')[2]
    bg_nr <- substr(cell_w_bg,bgloc+1,bgloc+1)
    # sapply(length(cols),function(i,datacols=cols,meth=method,bgs=bg,bgnr=bg_nr){
    #   if(meth[i]=="FC") {
    #       as.numeric(x[datacols[i]]) / bgs[,as.numeric(bgnr)] 
    #   } else{if(meth[i]=="shift") {
    #       as.numeric(x[datacols[i]]) - bgs[,as.numeric(bgnr)] 
    #     }}
    #   })
    sapply(1:length(cols),function(i){
      if(method[i]=="FC") {
        norm_nr <- as.numeric(x[cols[i]]) / bg[i,as.numeric(bg_nr)] 
      } else{if(method[i]=="shift") {
        norm_nr <- as.numeric(x[cols[i]]) - bg[i,as.numeric(bg_nr)] 
      }}
      norm_nr
    })
  }))
  table[,paste(data,'bgnorm',sep="_")] <- normalized_scores
  table
}  
