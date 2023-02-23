# ===================================================================================================================
#
#                 Functions for Stepwise or Dynamic IMCell
# 
# ===================================================================================================================


#' Convert static to dynamic network
#' 
#'
#' @param grnDF a static GRN in dataframe form
#' @param expX
#' @param sampTab meta data on cells, including pseudotime
#' @param pseudotime_column column name in sampTab with pseudotime annotation
#' @param epoch_annotation_column column name in sampTab that has pre-assigned epochs. If not NULL, num_epochs is ignored. Otherwise, epochs are computed via pseudotime.
#' @param num_epochs number of epochs, ignored if epoch_annotation_column is not NULL.
#' @param path path through the trajectory
#' @param column_annotation
#'
#' @return 
#'
#' @export
static_to_dynamic_trajectory<-function(grnDF,
                                    expX,
                                    sampTab,
                                    pseudotime_column="pseudotime",
                                    epoch_annotation_column=NULL,
                                    num_epochs=2,
                                    path=NULL,
                                    column_annotation=NULL){

  require(epoch)

  # subset data
  if (!is.null(path) & !is.null(column_annotation)){
    subset_sampTab<-sampTab[sampTab[,column_annotation] %in% path,]
    subset_expX<-expX[,rownames(subset_sampTab)]
  }else{
    cat("Taking full dataset as trajectory.")
    subset_sampTab<-sampTab
    subset_expX<-expX
  }
  
  # clean the data
  subset_expX<-subset_expX[(rowSums(subset_expX)!=0),]
  subset_expX_raw<-subset_expX

  # Normalize and transform the data if not already normalized
  if (all(subset_expX == floor(subset_expX))){
    subset_expX<-trans_rnaseq(subset_expX,1e5)
  }

  xdyn<-findDynGenes(subset_expX,subset_sampTab,path=path,group_column=column_annotation,pseudotime_column=pseudotime_column)
  
  if (is.null(epoch_annotation_column)){
    xdyn<-define_epochs(xdyn,subset_expX,method='pseudotime',num_epochs=num_epochs)
  }else{
    xdyn$cells$epoch<-subset_sampTab[,epoch_annotation_column][match(rownames(xdyn$cells),rownames(subset_sampTab))]
  }

  epoch_assignments<-assign_epochs(subset_expX,xdyn)
  dynamicGRN<-epochGRN(grnDF,epoch_assignments)

  subset_sampTab$epoch<-xdyn$cells$epoch[match(rownames(subset_sampTab),rownames(xdyn$cells))]

  return(list(dynamic_GRN=dynamicGRN,expDat=subset_expX,expRaw=subset_expX_raw,sampTab=subset_sampTab,epoch_assignments=epoch_assignments,xdyn=xdyn))

}

#' Find targets to activate and repress in dynamic trajectories at intermediate and final states
#' 
#'
#' @param dyntraj the result of running static_to_dynamic_trajectory
#' @param state_traj starting, intermediate, and final cell states
#' @param column_annotation column name in sampTab within dyntraj specifying cell states
#'
#' @return 
#' @export
#'
find_dyn_targets<-function(dyntraj, state_traj=NULL, column_annotation="celltype"){

  expX<-dyntraj$expRaw
  sampTab<-dyntraj$sampTab

  if (!is.null(state_traj)){

    targets<-vector("list",length(state_traj)-1)
    names(targets)<-state_traj[2:length(state_traj)]

    for (i in 2:length(state_traj)){
      tgs<-find_differential_nodes(expX,sampTab,state_traj[i-1],state_traj[i],column_annotation)
      targets[[state_traj[i]]]<-tgs
    }

  }else{


  }

  return(targets)

}






#' Run IMCell successively on dynamic network
#' 
#'
#' @param dynamic_grnDF a dynamic GRN in dataframe form
#' @param network_path the names and order of the subnetworks in dynamic_grnDF
#' @param expX
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited
#' @param sampTab meta data on cells, including pseudotime
#' @param kmax max number of TFs in solution set
#' @param repressor_wins See random_ICP
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param niter niter in MC
#' @param edge_probability_method method to assign edge probabilities
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param edge_probability_multiplier optional factor to multiply edge probabilities by. 1 = no multiplier
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores
#'
#' @return
#'
#' @export
IMCell_epochnets<-function(dynamic_grnDF,
                          network_path,
                          targets=NULL,
                          tfs=NULL,
                          kmax=5,
                          niter=1000,
                          min_marginal_spread=5,
                          repressor_wins=FALSE,
                          edge_probability_method="in_degree",
                          edge_probability_multiplier=1,
                          num_cores=2,
                          return_spread=FALSE,
                          temptesting=TRUE){


  epochs<-names(dynamic_grnDF)

  dynGRN<-dynamic_grnDF[network_path]

  solution_sets<-vector("list",length(names(dynGRN)))
  names(solution_sets)<-names(dynGRN)


  for (net in names(dynGRN)){
    grn<-dynGRN[[net]]
    igrn<-igraph::graph_from_data_frame(grn,directed=TRUE)

    # let tfs be TFs in GRN, not already in the solution sets
    if (is.null(tfs)){
      tfscope<-unique(grn$TF)
    }else{
      tfscope<-intersect(tfs,grn$TF)
    }
    tfscope<-setdiff(tfscope,unlist(solution_sets))

    if (!is.null(targets)){
      to_activate<-targets[[net]]$to_activate
      to_repress<-targets[[net]]$to_repress
    }else{
      to_activate<-NULL
      to_repress<-NULL
    }

    res<-IMCell(igrn,
                kmax=kmax,
                targets_activate=to_activate,
                targets_repress=to_repress,
                tfs=tfscope,
                niter=niter,
                min_marginal_spread=min_marginal_spread,
                repressor_wins=repressor_wins,
                edge_probability_method=edge_probability_method,
                edge_probability_multiplier=edge_probability_multiplier,
                num_cores=num_cores,
                return_spread=return_spread)

    solution_sets[[net]]<-res

  }

  solution_sets

}








#' Run IMCell_expweighted successively on dynamic network
#' 
#'
#' @param dynamic_grnDF a dynamic GRN in dataframe form
#' @param network_path the names and order of the subnetworks in dynamic_grnDF
#' @param expX
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited
#' @param sampTab meta data on cells, including pseudotime
#' @param epochs epoch names. If NULL inferred from network_path
#' @param kmax max number of TFs in solution set
#' @param repressor_wins See random_ICP
#' @param node_weight_method node weighitng method
#' @param min_node_weight minimum node weight
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param niter niter in MC
#' @param edge_probability_method method to assign edge probabilities
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param edge_probability_multiplier optional factor to multiply edge probabilities by. 1 = no multiplier
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores
#'
#' @return
#' @export
#'
IMCell_epochnets_expweighted<-function(dynamic_grnDF,
                          network_path,
                          targets=NULL,
                          expX=NULL,
                          sampTab=NULL,
                          epochs=NULL,
                          tfs=NULL,
                          kmax=5,
                          niter=1000,
                          min_marginal_spread=5,
                          repressor_wins=FALSE,
                          node_weight_method="no_weight",
                          min_node_weight=0.01,
                          edge_probability_method="in_degree",
                          edge_probability_multiplier=1,
                          num_cores=2,
                          return_spread=FALSE,
                          temptesting=TRUE){


  # check for expX, sampTab, and epoch assignments for node_weighting
  if (node_weight_method!="no_weight"){
    if (is.null(expX) | is.null(sampTab)){
        message("Node weighting requires expX, sampTab. Running unweighted IMCell.")
        node_weight_method<-"no_weight"
        min_node_weight<-1
      }
  }


  dynGRN<-dynamic_grnDF[network_path]

  solution_sets<-vector("list",length(names(dynGRN)))
  names(solution_sets)<-names(dynGRN)

  # get epochs labeled
  if (is.null(epochs)){
    epochs<-names(dynGRN)
    epochs<-unlist(lapply(epochs,function(t){strsplit(t,"\\..")[[1]][2]}))
    names(epochs)<-names(dynGRN)
  }

  for (net in names(dynGRN)){
    # get expX and sampTab if running node weighting
    if (node_weight_method!="no_weight"){
      secondstate<-epochs[net]
      sampTab_subset<-sampTab[sampTab$epoch==secondstate]
      expX_subset<-expX[,rownames(sampTab_subset)]
    }else{
      expX_subset<-expX
      sampTab_subset<-sampTab
      min_node_weight<-1
    }

    # get GRN
    grn<-dynGRN[[net]]
    igrn<-igraph::graph_from_data_frame(grn,directed=TRUE)

    # let tfs be TFs in GRN, not already in the solution sets
    if (is.null(tfs)){
      tfscope<-unique(grn$TF)
    }else{
      tfscope<-intersect(tfs,grn$TF)
    }
    tfscope<-setdiff(tfscope,unlist(solution_sets))

    if (!is.null(targets)){
      to_activate<-targets[[net]]$to_activate
      to_repress<-targets[[net]]$to_repress
    }else{
      to_activate<-NULL
      to_repress<-NULL
    }

    res<-IMCell_expweighted(igrn,
                            kmax=kmax,
                            expDat=expX_subset,
                            sampTab=sampTab_subset,
                            targets_activate=to_activate,
                            targets_repress=to_repress,
                            tfs=tfscope,
                            niter=niter,
                            min_marginal_spread=min_marginal_spread,
                            repressor_wins=repressor_wins,
                            node_weight_method=node_weight_method,
                            min_node_weight=min_node_weight,
                            edge_probability_method=edge_probability_method,
                            edge_probability_multiplier=edge_probability_multiplier,
                            num_cores=num_cores,
                            return_spread=return_spread)

    solution_sets[[net]]<-res

  }

  solution_sets

}





# Normalization functions from CellNet

#' weighted subtraction from mapped reades and log applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
trans_rnaseq<-function
(expRaw,
 total
 ){
    expCountDnW<-apply(expRaw, 2, downSampleW, total)
    log(1+expCountDnW)
  }

#' weighted subtraction from mapped reades
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param vector of total mapped reads per gene/transcript
#' @param total post transformation sum of read counts
#'
#' @return vector of downsampled read mapped to genes/transcripts
#'
downSampleW<-function
(vector,
total=1e5){ 

  totalSignal<-sum(vector)
  wAve<-vector/totalSignal
  resid<-sum(vector)-total #num to subtract from sample
  residW<-wAve*resid # amount to substract from each gene
  ans<-vector-residW
  ans[which(ans<0)]<-0
  ans
}



