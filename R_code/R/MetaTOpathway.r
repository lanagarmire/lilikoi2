#' A MetaTOpathway Function
#'
#' This function allows you to convert your metabolites id such as names, kegg ids, pubchem ids.
#' into pathways. Metabolites which have not pathways will be excluded from any downstream analysis
#' make sure that you have three database files which are used for exact and fuzzy matching:
#' cmpd_db.rda, syn_nms_db.rda and Sijia_pathway.rda
#' This function was modified version of the name.match function in the below link:
#' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
#' @param q.type The type of the metabolites id such as 'name', 'kegg', 'hmdb','pubchem'
#' @keywords Match
#' @export
#' @examples Metabolite_pathway_table=MetaTOpathway('name')
#' MetaTOpathway()
#'
#'
MetaTOpathway<- function(q.type, hmdb=T, pubchem=T, chebi=F, kegg=T, metlin=F){

  # record the filter for 8 major databases
  return.cols <<- c(hmdb, pubchem, chebi, kegg, metlin);

  #if(!is.null("cmpd.db")){
    #load("/home/fadl2/lilikoi/lilikoi_Fadhl/lilikoi/lilikoi/data/cmpd_db.rda", .GlobalEnv);
	
	#load(paste0(getwd(),"/","data/","cmpd_db.rda"), .GlobalEnv);
    #readRDS(file="canada/resources/libs/compound_db.rds", .GlobalEnv);
    #how to save image from the CSV
    #cmpd.db <<- read.csv("../../libs/cmpd_name.csv", as.is=T, header=T);
    #cmpd.db <<- read.csv("canada/resources/libs/cmpd_name.csv", as.is=T, header=T);
    #com.nms <- cmpd.db$name
    #syns.vec <- cmpd.db$synonym;
    #Encoding (syns.vec) <- "UTF-8"; # otherwise, encoding warnings
    #syn.db <- list();
    #syn.db$syns.vec <- cmpd.db$synonym;
    #syn.db$syns.list <- strsplit(syns.vec, "; *");
    #save(syn.db, file="syn_nms_db.rda");
    #cmpd.db$synonym <- NULL;
    #save(cmpd.db, file="cmpd_db.rda");
  #}

  # record all the data
  # if(!exists("name.map")){
  name.map <<- list();
  # }

  # distribute job
  dataSet$q.type <<- q.type;
  if(q.type == "name"){
    y=NameMappingExact();
	sanity_check(name.map,hit.inx) 
    return(y)
  }else{
    IDMapping(q.type);
	sanity_check(name.map,hit.inx) 
  }
}

  # do some sanity check
  #' A sanity_check Function
#'
#' This function return how many metabolites matched and unmatched from Metatopathway function 
#' @param name.map is the dataframe of metabolites and its corresponding ids.
#' @param hit.inx is the status of the metabolite hits (number if the metabolites matched the database or NA otherwise.) 
#' @keywords santiy search
#' @export
#' @examples sanity_check(name.map,hit.inx)
#' sanity_check
#'
#'

  sanity_check <- function(name.map,hit.inx){
        
       todo.inx <-which(is.na(name.map$hit.inx));
       if(length(todo.inx)/length(name.map$hit.inx) > 0.5){
            print("Over half of the compound IDs could not be matched to our database. Please make 
                            sure that correct compound IDs or common compound names are used.");
       }else if (length(todo.inx) > 100){
            print("There are >100 compounds without matches. You can either proceed or if necessary, 
                                    update these compound IDs and upload again.");        
       }else{
            #print("Name matching OK, please inspect (and manual correct) the results then proceed."); 
    print( paste0( (length(name.map$hit.inx)-length(todo.inx) ), " ", "out of ",length(name.map$hit.inx) 
             , " ","Matched metabolites"," ",
              round((length(name.map$hit.inx)/(length(name.map$hit.inx)+length(todo.inx))),1)*100," " 
              ,"%")  );
           cat("\n")
    print( paste0( (length(todo.inx) ), " ", "out of ",length(name.map$hit.inx)  
             , " ","UnMatched metabolites"," ",
              round((length(todo.inx)/(length(name.map$hit.inx)+length(todo.inx))),1)*100," " 
              ,"%")  );
                    
       }
    }
	
  
  

#' A NameMappingExact Function
#'
#' This function allows you to exact matching of your metabolites names with data base
#' consists of more than 20k metabolites. it allows you also to search among metabloites synanoms.
#' make sure that you have the below two databases for the exact matching:
#' cmpd_db.rda, and syn_nms_db.rda
#' This function was modified version of the name.match function in the below link:
#' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
#' @param Nothing
#' @keywords Exact Matching
#' @export
#' @examples NameMappingExact()
#' NameMappingExact()
#'
#'
NameMappingExact<-function(){

  qvec <- dataSet$cmpd;


  # variables to record results
  hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
  match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
  match.state = vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0

  # first find exact match to the common compound names
  hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));
  match.values <- cmpd.db$name[hit.inx];
  match.state[!is.na(hit.inx)] <- 1;

  # then try to find exact match to synonyms for the remaining unmatched query names one by one
  #if(!exists("syn.db")){
    #load("/home/fadl2/lilikoi/lilikoi_Fadhl/lilikoi/lilikoi/data/syn_nms_db.rda", .GlobalEnv);
	#load(paste0(getwd(),"/","data/","syn_nms_db.rda"), .GlobalEnv);
    #load("../../libs/syn_nms_db.rda", .GlobalEnv);
    #readRDS(file="canada/resources/libs/syn_nms.rds", .GlobalEnv);
  #}
  syns.list <-  syn.db$syns.list;

  todo.inx <<-which(is.na(hit.inx));
  if(length(todo.inx) > 0){
    for(i in 1:length(syns.list)){
      syns <-  syns.list[[i]];
      hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));

      hitPos <- which(!is.na(hitInx));
      if(length(hitPos)>0){
        # record matched ones
        orig.inx<-todo.inx[hitPos];
        hit.inx[orig.inx] <- i;
        # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
        match.values[orig.inx] <- cmpd.db$name[i];    # show common name
        match.state[orig.inx] <- 1;

        # update unmatched list
        todo.inx<<-todo.inx[is.na(hitInx)];

        #name.map$hit.inx <<- hit.inx;
        #name.map$hit.values <<- match.values;
        #name.map$match.state <<- match.state;
      }
      if(length(todo.inx) == 0) break;
    }
  }

  # empty memory
  gc();

  name.map$hit.inx <<- hit.inx;
  name.map$hit.values <<- match.values;
  name.map$match.state <<- match.state;

  GetMappingResultTable();

  for (xx in 1:length(todo.inx)){
    PerformApproxMatch(todo.inx[xx])
    if(is.null(candidates)){
      next;
    }else{


      SetCandidate(todo.inx[xx],1)
    }
  }
  x=PathMapping(dataSet$map.table)
  return(x)
}

#' A IDMapping Function
#'
#' This function allows you to exact matching of your metabolites hmdb, pubchem and kegg ids with data base
#' consists of more than 20k metabolites.
#' make sure that you have the below database for the exact matching:
#' cmpd_db.rda
#' This function was modified version of the name.match function in the below link:
#' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
#' @param q.type which is one of kegg id, pubchem is, hmdb id
#' @keywords Exact Matching
#' @export
#' @examples IDMapping(q.type)
#' IDMapping()
#'
IDMapping<-function(q.type){

  qvec <- dataSet$cmpd;
  # variables to record results
  hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
  match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
  match.state = vector(mode='numeric', length=length(qvec)); # match status - 0, no match; 1, exact match; 2 approximate match. initial 0

  if(q.type == "hmdb"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
  }else if(q.type == "pubchem"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$pubchem));
  }else if(q.type == "chebi"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$chebi));
  }else if(q.type == "metlin"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$metlin));
  }else if(q.type == "kegg"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
  }else{ # hmdb + kegg
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx2 <- match(tolower(qvec), tolower(cmpd.db$kegg));
    nohmdbInx <- is.na(hit.inx);
    hit.inx[nohmdbInx]<-hit.inx2[nohmdbInx]
  }

  match.values <- cmpd.db$name[hit.inx];
  match.state[!is.na(hit.inx)] <- 1;

  name.map$hit.inx <<- hit.inx;
  name.map$hit.values <<- match.values;
  name.map$match.state <<- match.state;

  GetMappingResultTable();


  PathMapping(dataSet$map.table)
}

#' A GetMappingResultTable Function
#'
#' This function allows you to generate csv file with user metabolits matched names, kegg, pubchem
#' and hmdb ids
#' This function was modified version of the name.match function in the below link:
#' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
#' @param nonething
#' @keywords data.map
#' @export
#' @examples GetMappingResultTable()
#' GetMappingResultTable()
#'
GetMappingResultTable<-function(){

  qvec <- dataSet$cmpd;
  if(is.null(qvec)){
    return();
  }
  # style for highlighted background for unmatched names
  pre.style<-NULL;
  post.style<-NULL;

  # style for no matches
  if(dataSet$q.type == "name"){
    no.prestyle<-"<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
    no.poststyle<-"</strong>";
  }else{
    no.prestyle<-"<strong style=\"background-color:red; font-size=125%; color=\"black\">";
    no.poststyle<-"</strong>";
  }

  hit.inx<-name.map$hit.inx;
  hit.values<-name.map$hit.values;
  match.state<-name.map$match.state;

  # contruct the result table with cells wrapped in html tags
  # the unmatched will be highlighted in different background
  html.res<-matrix("", nrow=length(qvec), ncol=8);
  csv.res<-matrix("", nrow=length(qvec), ncol=8);
  colnames(csv.res)<-c("Query", "Match", "HMDB", "PubChem", "ChEBI", "KEGG", "METLIN", "Comment");

  for (i in 1:length(qvec)){
    if(match.state[i]==1){
      pre.style<-"";
      post.style="";
    }else{ # no matches
      pre.style<-no.prestyle;
      post.style<-no.poststyle;
    }
    hit <-cmpd.db[hit.inx[i], ,drop=F];
    html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                     paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$hmdb_id) || hit$hmdb_id=="" || hit$hmdb_id=="NA","-", paste("<a href=http://www.hmdb.ca/metabolites/", hit$hmdb_id, " target='_blank'>",hit$hmdb_id,"</a>", sep="")),  sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$pubchem_id) || hit$pubchem_id=="" || hit$pubchem_id=="NA", "-", paste("<a href=http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid=", hit$pubchem_id," target='_blank'>", hit$pubchem_id,"</a>", sep="")), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$chebi_id) || hit$chebi_id==""|| hit$chebi_id=="NA","-", paste("<a href=http://www.ebi.ac.uk/chebi/searchId.do?chebiId=", hit$chebi_id, " target='_blank'>",hit$chebi_id,"</a>", sep="")), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$kegg_id) || hit$kegg_id==""|| hit$kegg_id=="NA","-",paste("<a href=http://www.genome.jp/dbget-bin/www_bget?", hit$kegg_id, " target='_blank'>", hit$kegg_id,"</a>", sep="")), sep=""),
                     paste(ifelse(match.state[i]==0 || is.na(hit$metlin_id) || hit$metlin_id==""|| hit$metlin_id=="NA","-",paste("<a href=http://metlin.scripps.edu/metabo_info.php?molid=", hit$metlin_id," target='_blank'>",hit$metlin_id,"</a>", sep="")), sep=""),
                     ifelse(match.state[i]!=1,"View",""));
    csv.res[i, ]<-c(qvec[i],
                    ifelse(match.state[i]==0, "NA", hit.values[i]),
                    ifelse(match.state[i]==0, "NA", hit$hmdb_id),
                    ifelse(match.state[i]==0, "NA", hit$pubchem_id),
                    ifelse(match.state[i]==0, "NA", hit$chebi_id),
                    ifelse(match.state[i]==0, "NA", hit$kegg_id),
                    ifelse(match.state[i]==0, "NA", hit$metlin_id),
                    match.state[i]);
  }
  # return only columns user selected

  # add query and match columns at the the beginning, and 'Detail' at the end
  return.cols <- c(TRUE, TRUE, return.cols, TRUE);
  html.res <- html.res[,return.cols, drop=F];
  csv.res <- csv.res[,return.cols, drop=F];

  # store the value for report
  dataSet$map.table <<- csv.res;
  write.csv(csv.res, file="name_map.csv", row.names=F);
  #return(as.vector(html.res));
}


#' A PerformApproxMatch Function
#'
#' This function allows you to fuzzy matching of your metabolites names with data base
#' consists of more than 20k metabolites. it allows you also to match also among metabloites synanoms.
#' make sure that you have the below two databases for the exact matching:
#' cmpd_db.rda, and syn_nms_db.rda
#' This function was modified version of the name.match function in the below link:
#' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
#' @param the index of unmatched metabolites
#' @keywords fuzzy logic Matching
#' @export
#' @examples PerformApproxMatch(inx)
#' PerformApproxMatch(inx)
#'
#'

PerformApproxMatch<-function(inx){
  #print(inx)
  q <- dataSet$cmpd[inx];
  # only for none lipids
  #nonLipidInx <- cmpd.db$lipid == 0;

  #     com.nms <- cmpd.db$name[nonLipidInx]
  #     syns.vec <- syn.db$syns.vec[nonLipidInx];
  #     syns.list <- syn.db$syns.list[nonLipidInx];

  com.nms <- cmpd.db$name;
  syns.vec <- syn.db$syns.vec;
  syns.list <- syn.db$syns.list;

  matched.dist <- NULL;
  q.length <- nchar(q);
  s <- c(0);
  for (j in s) {
    new.q <- q;
    if(q.length > 32){ # note: agrep fail for exact match when length over 32 characters
      new.q<-substr(q, 1, 32);
      #new.q <- unlist(lapply(q,function(x){substr(q, 1, 32)})) ;
    }
    matched <- FALSE;
    matched.inx <- agrep(q, syns.vec, ignore.case=T, max.distance=j, useBytes=T);
    #matched.inx <- unlist(sapply(new.q, function(x) {agrep(x, syns.vec,ignore.case =T, max.distance=j, useBytes=T ) } ))
    if(length(matched.inx) > 0) {
      # record all the candidates,
      # don't use cbind, since all will be converted to character mode
      # for data.frame specify "stringsAsFactors" to prevent convert value col into factor
      candidates <- data.frame(index=vector(mode = "numeric", length=length(matched.inx)),
                               value=vector(mode = "character", length=length(matched.inx)),
                               score=vector(mode = "numeric", length=length(matched.inx)),
                               stringsAsFactors = FALSE);

      for(n in 1:length(matched.inx)){
        nm.vec<-syns.list[[matched.inx[n]]];
        # try approximate match, note: in some cases, split into element will break match using whole string
        hit3.inx <- agrep(q, nm.vec,ignore.case=T, max.distance=j, useBytes=T);
        #hit3.inx <- unlist(sapply(q, function(x) {agrep(x, nm.vec,ignore.case =T, max.distance=j, useBytes=T ) } ))
        if(length(hit3.inx)>0){
          hit3.nm <- vector(mode = "character", length=length(hit3.inx));
          hit3.score <- vector(mode = "numeric", length=length(hit3.inx));
          for(k in 1:length(hit3.inx)){
            idx <- hit3.inx[k];
            hit3.nm[k] <- nm.vec[idx];
            hit3.score[k] <- j + abs(nchar(nm.vec[idx])-nchar(q))/(10*nchar(q));
          }

          # now get the best match, the rule is that the first two character should matches
          # first check if first two character are digits or characters, otherwise will cause error
          matches2 <- c();
          if(length(grep("^[1-9a-z]{2}", q, ignore.case=T))>0){
            matches2 <- grep(paste("^", substr(q, 1, 2), sep=""), hit3.nm);
          }else if (length(grep("^[1-9a-z]", q, ignore.case=T))>0){
            matches2 <- grep(paste("^", substr(q, 1, 1), sep=""), hit3.nm);
          }

          if(length(matches2)>0){
            hit3.score[matches2] <- hit3.score[matches2] - 0.05;
          }

          best.inx<-which(hit3.score==min(hit3.score))[1];
          candidates[n,1]<-matched.inx[n];
          #    candidates[n,2]<-hit3.nm[best.inx]; # show matched syn names
          candidates[n,2]<-com.nms[matched.inx[n]] # show common names
          candidates[n,3]<-hit3.score[best.inx];
        }
      }
      rm.inx <- is.na(candidates[,2]) | candidates[,2]=="NA" | candidates[,2]=="";
      candidates<-candidates[!rm.inx, ];
      candidates<-candidates[order(candidates[,3], decreasing=F), , drop=F];
      if(nrow(candidates) > 10){
        candidates<-candidates[1:10,];
      }
      candidates <<- candidates;
      #return();
    }
  }

  candidates <<- NULL;


}


#' A SetCandidate Function
#'
#' This function allows you to fuzzy matching of your metabolites names with data base
#' consists of more than 20k metabolites. it allows you also to match also among metabloites synanoms.
#' make sure that you have the below two databases for the exact matching:
#' cmpd_db.rda, and syn_nms_db.rda
#' it uses agrap R function
#' This function was modified version of the name.match function in the below link:
#' https://github.com/cangfengzhe/Metabo/blob/master/MetaboAnalyst/website/name_match.R
#' @param the index of unmatched metabolites and the matched metabolites with the best matching score
#' @keywords data mapping
#' @export
#' @examples SetCandidate(query_inx, can_inx)
#' SetCandidate(query_inx, can_inx)
#'
#'
SetCandidate<-function(query_inx, can_inx){
  can_mat<- candidates;
  if(can_inx <= nrow(can_mat)){
    name.map$hit.inx[query_inx] <<- can_mat[can_inx,1];
    name.map$hit.values[query_inx] <<- can_mat[can_inx,2];
    name.map$match.state[query_inx] <<- 1;

    # re-generate the CSV file
    hit <-cmpd.db[name.map$hit.inx[query_inx], ,drop=F];
    csv.res <- dataSet$map.table;
    if(ncol(csv.res) > 6){ # general utilities
      csv.res[query_inx, ]<-c(csv.res[query_inx, 1],
                              name.map$hit.values[query_inx],
                              hit$hmdb_id,
                              hit$pubchem_id,
                              hit$chebi_id,
                              hit$kegg_id,
                              hit$metlin_id,
                              1);
    }else{
      csv.res[query_inx, ]<-c(csv.res[query_inx, 1],
                              name.map$hit.values[query_inx],
                              hit$hmdb_id,
                              hit$pubchem_id,
                              hit$kegg_id,
                              1);
    }
    write.csv(csv.res, file="name_map.csv", row.names=F);
    dataSet$map.table <<- csv.res;

  }else{ #no match
    name.map$hit.inx[query_inx] <<- 0;
    name.map$hit.values[query_inx] <<- "";
    name.map$match.state[query_inx] <<- 0;
  }
}

#' A PathMapping Function
#'
#' This function allows you to map user metabolites to the its corrsponding pathways
#' make sure that you have the below database for the exact matching:
#' Sijia_pathway.rda
#' @param the list of user metabolites and its corrsponding keeg, pubchem and hmdb ids
#' @keywords pathway mapping
#' @export
#' @examples PathMapping(qvec)
#' PathMapping(qvec)
#'
#'
PathMapping<-function(qvec){
  require(dplyr)
  #if(!exists("sijia_pathway")){
    #load("/home/fadl2/lilikoi/lilikoi_Fadhl/lilikoi/lilikoi/data/Sijia_pathway.rda",.GlobalEnv);
	#load(paste0(getwd(),"/","data/","Sijia_pathway.rda"), .GlobalEnv);
  #}
  #could we add the pathway for these metabolites from Sijia pathway database
  met_W_pathway=HMDBlib[match(qvec[,'HMDB'],HMDBlib$Accession),,drop=F]
  dataSet$map.table1 <<- data.frame(qvec) %>%
    mutate(pathway=met_W_pathway$Pathway_Name) %>% arrange(pathway)
  write.csv(dataSet$map.table1, file="name_map.csv", row.names=F);
  return(dataSet$map.table1)
}
