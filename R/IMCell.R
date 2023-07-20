# ===================================================================================================================
#
#									IMCell functions + related helper functions
#	--> IMCell_CELF
#	--> IMCell
#	--> IMCell_expweighted
#	
# ===================================================================================================================



# ===================================================================================================================
#
#										Activation only functions
#
#
# ===================================================================================================================

#' Cost Effective Lazy Forward algorithm (CELF)
#' 
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities
#' @param kmax maximum number of TFs in the solution set
#' @param tfs TFlist
#' @param targets consider activation of specified targets only. If NULL activation is based on total spread
#' @param niter number of random cascades to run
#' @param edge_probability_method method to convert edge weights to probabilities
#' @param edge_probability_multiplier optional factor to multiply edge probabilities by. 1 = no multiplier
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores to use if running parallel
#'
#' @return 
#'
#' @export
#'
IMCell_CELF<-function(ig,
					kmax=5,
					tfs=NULL,
					targets=NULL,
					niter=1000,
					edge_probability_method="in_degree",
					edge_probability_multiplier=1,
					min_marginal_spread=5,
					return_spread=FALSE,
					num_cores=2){

	# set cores
	num_cores<-min(num_cores,detectCores())

	# Set weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	solution_set<-c()
	prior_spread<-0

	# --------------------- First round pass and sort: ---------------------
	# Keep track of spread of each node

	sort_TFs<-data.frame(TF=character(),spread=numeric())

	best_spread<- -0.1 # set negative to ensure replacement
	best_node<-NULL
	expected_activation<-c()
	for (node in tfs){
		print(node)
		mc_res<-mc_activation_cascade_parallel(ig,node,targets=targets,niter=niter,return_spread=return_spread,num_cores=num_cores)
		spread<-mc_res$predicted_spread

		sort_TFs<-rbind(sort_TFs,data.frame(TF=node,spread=spread))
		
		if(spread>best_spread){
			if (return_spread){
				expected_activation<-mc_res$expected_activation
			}
			best_spread<-spread
		}
	}

	sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]

	# add the first TF to the solution set, set the prior spread
	solution_set<-c(solution_set,sort_TFs$TF[1])
	prior_spread<-sort_TFs$spread[1]

	# remove any TFs that have spread < min_marginal_spread -- these aren't worth testing any further
	sort_TFs<-sort_TFs[sort_TFs$spread>=min_marginal_spread,]

	if(round(prior_spread)<1){
		warning("No predicted TFs.")
		return(list(solution_set=NA,expected_activation=NULL))
	}

	# Check if 1 TF covered all of the targets
	if (round(prior_spread)==length(targets) | nrow(sort_TFs)<=1){
		if(return_spread){
			return(list(solution_set=solution_set,expected_activation=targets))
		}else{
			return(list(solution_set=solution_set))
		}
	}

	
	# --------------------------- Rounds 2-k -------------------------------

	for (k in 2:kmax){
		print(k)
		sort_TFs<-sort_TFs[-1,]

		# remove the top spread (solution set);
		# then for each iteration of the while loop
		# re-compute marginal spread of the top TF in sort_TFs;
		# then add (replace) the marginal spread into sort_TFs and resort;
		# repeat the while loop until the marginal spread still ends up
		# at the top -- add this to the solution set.

		is_top<-FALSE
		while(!is_top){
			top_TF<-sort_TFs$TF[1]
			print(top_TF)

			mc_res<-mc_activation_cascade_parallel(ig,c(solution_set,top_TF),targets=targets,niter=niter,return_spread=return_spread,num_cores=num_cores)
			marginal_spread<-mc_res$predicted_spread - prior_spread

			sort_TFs$spread[1]<-marginal_spread
			sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]

			is_top <- sort_TFs$TF[1]==top_TF
		}

		if (round(mc_res$predicted_spread)-round(prior_spread)<min_marginal_spread){
			break
		}

		solution_set<-c(solution_set,top_TF)
		prior_spread<-mc_res$predicted_spread
		if (return_spread){
			expected_activation<-mc_res$expected_activation
		}

		if (nrow(sort_TFs)<=1){
			break
		}

	}


	if (return_spread){
		return(list(solution_set=solution_set,expected_activation=expected_activation))
	}else{
		return(list(solution_set=solution_set))
	}

}




#' run MC cascade and return expected spread, assumes edge weights are positive
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities
#' @param seed starting nodes to activate
#' @param targets consider activation of specified targets only. If NULL activation is based on total spread
#' @param niter
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores to use if running parallel
#'
#'
mc_activation_cascade_parallel<-function(ig,seed,targets=NULL,niter=1000,return_spread=FALSE,num_cores=2){

	rcs<-mclapply(1:niter,function(x){random_activation_cascade_parallel(ig=ig,seed=seed,targets=targets,return_spread=return_spread)},mc.cores=num_cores)

	total_activated<-lapply(rcs,function(x){x$activated_spread})
	total_activated<-sum(unlist(total_activated))

	predicted_num_activated<-total_activated/niter

	if(return_spread){
		activated_tracker<-unlist(lapply(rcs,function(x){x$activated}))

		if (length(activated_tracker)==0){
			expected_activation<-c()
		}else{
			activated_tracker<-as.data.frame(table(activated_tracker))
			colnames(activated_tracker)<-c("gene","Freq")
			activated_tracker<-activated_tracker[order(activated_tracker$Freq,decreasing=TRUE),]
			expected_activation<-activated_tracker$gene[1:round(predicted_num_activated)]
		}

	}


	if(return_spread){
		mc_res<-list(predicted_spread=predicted_num_activated,expected_activation=expected_activation)
	}else{
		mc_res<-list(predicted_spread=predicted_num_activated)
	}

	mc_res

}



#' Simulates a random cascade (activation only)
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities
#' @param seed starting nodes to activate
#' @param targets consider activation of specified targets only. If NULL activation is based on total spread
#' @param return_spread whether or not to return the set of activated nodes
#'
#'
random_activation_cascade_parallel<-function(ig,seed,targets=NULL,return_spread=FALSE){

	all_activated<-seed
	#failed_to_activate<-c()

	activated<-seed
	while (length(activated)>0){

		nn<-lapply(adjacent_vertices(ig,activated,mode="out"),function(x){x<-x$name;
																		  x<-x[!(x %in% activated)];
																		  x<-x[!(x %in% all_activated)]})

		# keep only the elements with targets
		nn<-nn[lapply(nn,length)>0]

		all<-Reduce(c,nn)
		if (length(all)==0){
			break
		}

		# coin flips
		flips<-runif(length(all),min=0,max=1)
		# which flips were successful (which flips are less than the edge weight)
		# e.g. if prob = .8 --> generate runif between 0 and 1 --> if falls between 0-.8, call success; if between .8-1 call fail
	
		edges<-lapply(seq_along(nn),function(i){x<-as.vector(rbind(names(nn)[i],nn[[i]]));x})
		edges<-Reduce(c,edges)
		#P<-as.vector(rbind(node,nn))
		success<-flips < E(ig,P=edges)$weight

		all_activated<-c(all_activated,all[success])
		#failed_to_activate<-c(failed_to_activate,all[!success])

		activated<-all[success]
	}

	all_activated<-unique(all_activated)		


	# compute spread

	if (!is.null(targets)){
		all_activated<-intersect(all_activated,targets)
	}
	activated_spread<-length(all_activated)


	# return
	if (return_spread){
		return(list(activated_spread=activated_spread,activated=all_activated))
	}else{
		return(list(activated_spread=activated_spread))
	}


}




# Temporary testing function that adds in marginal spread tracker, otherwise same function as IMCell_CELF

im_celf_parallel_testing<-function(ig,
					kmax=5,
					tfs=NULL,
					targets=NULL,
					niter=1000,
					edge_probability_method="in_degree",
					edge_probability_multiplier=1,
					min_marginal_spread=5,
					return_spread=FALSE,
					num_cores=2,
					temptesting=TRUE){

	if(temptesting){
		return_spread=TRUE
	}


	# set cores
	num_cores<-min(num_cores,detectCores())

	# Set weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	solution_set<-c()
	prior_spread<-0

	# --------------------- First round pass and sort: ---------------------
	# Keep track of spread of each node

	# for temptesting
	marg_spread_tracker<-c()

	sort_TFs<-data.frame(TF=character(),spread=numeric())

	best_spread<- -0.1 # set negative to ensure replacement
	best_node<-NULL
	for (node in tfs){
		print(node)
		mc_res<-mc_activation_cascade_parallel(ig,node,targets=targets,niter=niter,return_spread=return_spread,num_cores=num_cores)
		spread<-mc_res$predicted_spread

		sort_TFs<-rbind(sort_TFs,data.frame(TF=node,spread=spread))

		if(spread>best_spread){
			if (return_spread){
				expected_activation<-mc_res$expected_activation
			}
			best_spread<-spread
		}
	}

	sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]

	# add the first TF to the solution set, set the prior spread
	solution_set<-c(solution_set,sort_TFs$TF[1])
	prior_spread<-sort_TFs$spread[1]

	if (temptesting){
		marg_spread_tracker<-c(marg_spread_tracker,prior_spread)
	}


	# remove any TFs that have spread < min_marginal_spread -- these aren't worth testing any further
	sort_TFs<-sort_TFs[sort_TFs$spread>=min_marginal_spread,]

	if(round(prior_spread)<1){
		warning("No predicted TFs.")
		return(list(solution_set=NA,expected_activation=NULL))
	}

	# Check if 1 TF covered all of the targets
	if (round(prior_spread)==length(targets) | nrow(sort_TFs)<=1){
		if(return_spread){
			return(list(solution_set=solution_set,expected_activation=targets))
		}else{
			return(list(solution_set=solution_set))
		}
	}

	
	# --------------------------- Rounds 2-k -------------------------------

	for (k in 2:kmax){
		print(k)
		sort_TFs<-sort_TFs[-1,]

		# remove the top spread (solution set);
		# then for each iteration of the while loop
		# re-compute marginal spread of the top TF in sort_TFs;
		# then add (replace) the marginal spread into sort_TFs and resort;
		# repeat the while loop until the marginal spread still ends up
		# at the top -- add this to the solution set.

		is_top<-FALSE
		while(!is_top){
			top_TF<-sort_TFs$TF[1]
			print(top_TF)

			mc_res<-mc_activation_cascade_parallel(ig,c(solution_set,top_TF),targets=targets,niter=niter,return_spread=return_spread,num_cores=num_cores)
			marginal_spread<-mc_res$predicted_spread - prior_spread

			sort_TFs$spread[1]<-marginal_spread
			sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]

			is_top <- sort_TFs$TF[1]==top_TF
		}

		if (temptesting){
			marg_spread_tracker<-c(marg_spread_tracker,(mc_res$predicted_spread-prior_spread))
		}

		if (round(mc_res$predicted_spread)-round(prior_spread)<min_marginal_spread){
			break
		}

		solution_set<-c(solution_set,top_TF)
		prior_spread<-mc_res$predicted_spread
		if (return_spread){
			expected_activation<-mc_res$expected_activation
		}

		if (nrow(sort_TFs)<=1){
			break
		}

	}

	if(temptesting){
		# compute total possible spread -- assume all targets provided are in the network
		possible_act<-length(targets)
		possible_rep<-0

		return(list(solution_set=solution_set,possible_act_spread=possible_act,possible_rep_spread=possible_rep,marg_spread_tracker=marg_spread_tracker,activated_targets=expected_activation,repressed_targets=c()))
	}


	if (return_spread){
		return(list(solution_set=solution_set,expected_activation=expected_activation))
	}else{
		return(list(solution_set=solution_set))
	}

}










# ===================================================================================================================
#
#										IMCell functions
#
#
# ===================================================================================================================




#' Greedy MC Algorithm for IMCell (spread = targets activated + targets repressed)
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, 'type' corresponding to interaction type,
#' @param kmax max number of TFs in solution set
#' @param repressor_wins See random_ICP
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
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
IMCell<-function(ig,
				kmax=5,
				repressor_wins=FALSE,
				tfs=NULL,
				targets_activate=NULL,
				targets_repress=NULL,
				niter=1000,
				edge_probability_method="in_degree",
				edge_probability_multiplier=1,
				min_marginal_spread=5,
				return_spread=FALSE,
				num_cores=2){

	# set cores
	num_cores<-min(num_cores,detectCores())


	# Set weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	# If no interaction type, assume all activation.
	if (is.null(E(ig)$type)){
		E(ig)$type<-1
		message("E(ig)$type is NULL, assuming all activation.")
	}

	solution_set<-c()
	prior_spread<-0

	final_activated_nodes<-c()
	final_repressed_nodes<-c()

	for (k in 1:kmax){
		print(k)
		nodes<-setdiff(tfs,solution_set)

		best_spread<- -0.1 # set negative to ensure replacement
		best_node<-NULL

		# Run mc_ICP_parallel and add node with best spread
		for (node in nodes){
			print(node)
			mc_res<-mc_ICP_parallel(ig,c(solution_set,node),repressor_wins=repressor_wins,targets_activate=targets_activate,targets_repress=targets_repress,niter=niter,return_spread=return_spread,num_cores=num_cores)
			spread<-mc_res$activation_spread+mc_res$repression_spread

			if(spread > best_spread){
				best_spread<-spread
				best_node<-node
				
				activated_nodes<-mc_res$activated_targets
				repressed_nodes<-mc_res$repressed_targets
			}
		}

		if ((best_spread-prior_spread<min_marginal_spread) & (k>1)){
			break
		}

		cat("Marginal_spread: ",best_spread-prior_spread,'\n')

		prior_spread<-best_spread
		solution_set<-c(solution_set,best_node)
		final_activated_nodes<-activated_nodes
		final_repressed_nodes<-repressed_nodes

	}

	if (return_spread){
		return(list(solution_set=solution_set,activated_targets=final_activated_nodes,repressed_targets=final_repressed_nodes))

	}else{
		return(list(solution_set=solution_set))
	}


}






#' Cost Effective Lazy Forward algorithm (CELF) for single signed IM (either PIM or NIM) in PARALLEL
#' 
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities
#' @param kmax max number of TFs in solution set
#' @param mode either "PIM" or "NIM" -- maximize activation or repression
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param repressor_wins See random_ICP
#' @param targets consider activation (or repression is mode="NIM") of specified targets only. If NULL activation is based on total spread
#' @param niter
#' @param edge_probability_method method to assign edge probabilities
#' @param edge_probability_multiplier optional factor to multiply edge probabilities by. 1 = no multiplier
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param return_spread whether or not to return the set of influenced nodes
#'
#' @return 
#' @export
#'
im_single_PRIM_celf_parallel<-function(ig,
					kmax=5,
					mode="PIM",
					tfs=NULL,
					repressor_wins=FALSE,
					targets=NULL,
					niter=1000,
					edge_probability_method="in_degree",
					edge_probability_multiplier=1,
					min_marginal_spread=5,
					return_spread=FALSE,
					num_cores=2){

	# set cores
	num_cores<-min(num_cores,detectCores())

	# Set weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	# If no interaction type, assume all activation.
	if (is.null(E(ig)$type)){
		E(ig)$type<-1
		message("E(ig)$type is NULL, assuming all activation.")
	}

	# mode must be "PIM" or "NIM"
	if (!(mode %in% c("PIM","NIM"))){
		stop("mode must be PIM or NIM.")
	}

	solution_set<-c()
	prior_spread<-0

	final_activated_nodes<-c()
	final_repressed_nodes<-c()


	# --------------------- First round pass and sort: ---------------------
	# Keep track of spread of each node

	sort_TFs<-data.frame(TF=character(),spread=numeric())

	best_spread<- -0.1 # set negative to ensure replacement
	best_node<-NULL

	for (node in tfs){
		print(node)
		if (mode=="PIM"){
			mc_res<-mc_ICP_parallel(ig,c(solution_set,node),repressor_wins=repressor_wins,targets_activate=targets,targets_repress=NULL,niter=niter,return_spread=return_spread)
		}else{
			mc_res<-mc_ICP_parallel(ig,c(solution_set,node),repressor_wins=repressor_wins,targets_activate=NULL,targets_repress=targets,niter=niter,return_spread=return_spread)
		}

		spread<-mc_res$activation_spread + mc_res$repression_spread
		sort_TFs<-rbind(sort_TFs,data.frame(TF=node,spread=spread))

		if(spread > best_spread){
			best_spread<-spread
			best_node<-node
			
			activated_nodes<-mc_res$activated_targets			# NULL if return_spread=FALSE
			repressed_nodes<-mc_res$repressed_targets
		}

	}

	sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]
	# remove any TFs that have spread < min_marginal_spread -- these aren't worth testing any further
	sort_TFs<-sort_TFs[sort_TFs$spread>=min_marginal_spread,]

	# add the first TF to the solution set, set the prior spread
	solution_set<-c(solution_set,sort_TFs$TF[1])
	prior_spread<-sort_TFs$spread[1]


	# Check if 1 TF covered all of the targets
	if (round(prior_spread) == (length(targets))){
		if(return_spread){
			return(list(solution_set=solution_set,activated_targets=activated_nodes,repressed_targets=repressed_nodes))
		}else{
			return(list(solution_set=solution_set))
		}
	}


	# --------------------------- Rounds 2-k -------------------------------

	for (k in 2:kmax){
		print(k)
		sort_TFs<-sort_TFs[-1,]

		# remove the top spread (solution set);
		# then for each iteration of the while loop
		# re-compute marginal spread of the top TF in sort_TFs;
		# then add (replace) the marginal spread into sort_TFs and resort;
		# repeat the while loop until the marginal spread still ends up
		# at the top -- add this to the solution set.

		is_top<-FALSE
		while(!is_top){
			top_TF<-sort_TFs$TF[1]
			print(top_TF)

			if (mode=="PIM"){
				mc_res<-mc_ICP_parallel(ig,c(solution_set,top_TF),repressor_wins=repressor_wins,targets_activate=targets,targets_repress=NULL,niter=niter,return_spread=return_spread)
			}else{
				mc_res<-mc_ICP_parallel(ig,c(solution_set,top_TF),repressor_wins=repressor_wins,targets_activate=NULL,targets_repress=targets,niter=niter,return_spread=return_spread)
			}

			spread<-mc_res$activation_spread + mc_res$repression_spread
			marginal_spread<-spread - prior_spread
			sort_TFs$spread[1]<-marginal_spread
			sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]

			is_top<-sort_TFs[1]==top_TF
		}

		if ((round(spread) - round(prior_spread)) < min_marginal_spread){
			break
		}

		solution_set<-c(solution_set,top_TF)
		prior_spread<-spread

		if (return_spread){
			activated_nodes<-mc_res$activated_targets
			repressed_nodes<-mc_res$repressed_targets
		}

	}

	if (return_spread){
		return(list(solution_set=solution_set,activated_targets=activated_nodes,repressed_targets=repressed_nodes))
	}else{
		return(list(solution_set=solution_set))
	}

}






# mc_ICP_parallel only compatible with random_ICP_parallel
# targets_activate/repress = NULL dealt with in random_ICP_parallel

#' run MC cascade and return expected spread on IC-P model
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, and 'type' corresponding to interaction type
#' @param seed starting nodes to activate
#' @param repressor_wins see random_ICP for description
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
#' @param niter
#' @param return_spread whether or not to return the set of activated nodes
#'
#'
mc_ICP_parallel<-function(ig,seed,repressor_wins=FALSE,targets_activate=NULL,targets_repress=NULL,niter=1000,return_spread=FALSE,num_cores=2){


	activated_tracker<-c()
	repressed_tracker<-c()


	rcs<-mclapply(1:niter,function(x){random_ICP_parallel(ig=ig,seed=seed,repressor_wins=repressor_wins,
															targets_activate=targets_activate,targets_repress=targets_repress,
															return_spread=return_spread)},mc.cores=num_cores)

	total_activated<-lapply(rcs,function(x){x$activated_spread})
	total_activated<-sum(unlist(total_activated))
	total_repressed<-lapply(rcs,function(x){x$repressed_spread})
	total_repressed<-sum(unlist(total_repressed))

	predicted_num_activated<-total_activated/niter
	predicted_num_repressed<-total_repressed/niter

	if(return_spread){
		activated_tracker<-unlist(lapply(rcs,function(x){x$activated}))
		repressed_tracker<-unlist(lapply(rcs,function(x){x$repressed}))

		if (length(activated_tracker)==0){
			expected_activation<-c()
		}else{
			activated_tracker<-as.data.frame(table(activated_tracker))
			colnames(activated_tracker)<-c("gene","Freq")
			activated_tracker<-activated_tracker[order(activated_tracker$Freq,decreasing=TRUE),]
			expected_activation<-activated_tracker$gene[1:round(predicted_num_activated)]
		}

		if (length(repressed_tracker)==0){
			expected_repression<-c()
		}else{
			repressed_tracker<-as.data.frame(table(repressed_tracker))
			colnames(repressed_tracker)<-c("gene","Freq")
			repressed_tracker<-repressed_tracker[order(repressed_tracker$Freq,decreasing=TRUE),]
			expected_repression<-repressed_tracker$gene[1:round(predicted_num_repressed)]
		}

	}


	if(return_spread){
		mc_res<-list(activation_spread=predicted_num_activated,repression_spread=predicted_num_repressed,activated_targets=expected_activation,repressed_targets=expected_repression)
	}else{
		mc_res<-list(activation_spread=predicted_num_activated,repression_spread=predicted_num_repressed)
	}

	mc_res

}



# IC-P (independent cascade - polarity)
# Unlike random_simple_cascade, the repressed node carries propagates influence
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4111297/
# repressor_wins if FALSE, follows IC-P as described in ^^ (i.e. successful edge is chosen randomly from live edges/
#					edge flipping order is random and first successful edge influences target). If TRUE, if any repressive
#					edge is successfully flipped, then target will be repressed
# 
# almost the same as random_ICP, but with return_spread option
# and computes spread within this function (rather than in MC)

#' run random cascade on IC-P model
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, and 'type' corresponding to interaction type
#' @param seed starting nodes to activate
#' @param repressor_wins see random_ICP for description
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
#' @param return_spread whether or not to return the set of activated nodes
#'
#'
random_ICP_parallel<-function(ig,seed,repressor_wins=FALSE,targets_activate=NULL,targets_repress=NULL,return_spread=FALSE){
	# Three states a node can be in: positive (activated/on), negative (repressed/off), inactive (not influenced)
	# Three states an edge can be in: live (flipped), dead (failed to flip), unflipped

	positive_nodes<-seed
	negative_nodes<-c()

	new_positive<-seed
	new_negative<-c()
	new_influenced<-union(new_positive,new_negative)
	while(length(new_influenced)>0){
		# find not-yet-influenced neighbors
		nn<-lapply(adjacent_vertices(ig,new_influenced,mode="out"),function(x){x<-x$name;
																		  		x<-x[!(x %in% positive_nodes)];
																		  		x<-x[!(x %in% negative_nodes)];
																		  		x<-x[!(x %in% new_influenced)]})

		# keep only the elements with targets (i.e. remove non-targets)
		nn<-nn[lapply(nn,length)>0]

		all<-Reduce(c,nn)
		if (length(all)==0){
			break
		}

		# Create edge DF
		edges<-lapply(seq_along(nn),function(i){x<-as.vector(rbind(names(nn)[i],nn[[i]]));x})
		edges<-Reduce(c,edges)

		edges_df<-do.call(rbind,lapply(1:length(nn),function(i){data.frame(origin=rep(names(nn)[i],length(nn[[i]])),target=nn[[i]])}))
		edges_df$origin_type<-ifelse(edges_df$origin %in% new_positive,1,-1)
		edges_df$type<-E(ig,P=edges)$type

		edges_df$target_polarity<-edges_df$origin_type * edges_df$type

		# coin flips
		flips<-runif(length(all),min=0,max=1)
		# which flips were successful (which flips are less than the edge weight)
		# e.g. if prob = .8 --> generate runif between 0 and 1 --> if falls between 0-.8, call success; if between .8-1 call fail
		edges_df$success<-flips < E(ig,P=edges)$weight

		edges_df<-edges_df[edges_df$success,]

		if (nrow(edges_df)==0){
			break
		}

		if (!repressor_wins){
			# Traditional IC-P model = only one edge influences a give target
			# Order of edge flipping is random
			# AKA succesful edge is chosen uniformly

			# scramble the order, then limit to first live edge
			edges_df<-edges_df[sample(1:nrow(edges_df)),]
			edges_df<-edges_df[!duplicated(edges_df$target),]

			new_positive<-edges_df$target[edges_df$target_polarity==1]
			new_negative<-edges_df$target[edges_df$target_polarity==-1]

		}else{
			# Alternative combinatorial rules:
			# If any repressive edge is live, the target is repressed and assigned negative state

			new_negative<-edges_df$target[edges_df$target_polarity==-1]
			new_positive<-edges_df$target[edges_df$target_polarity==1]
			new_positive<-setdiff(new_positive,new_negative)

		}

		positive_nodes<-c(positive_nodes,new_positive)
		negative_nodes<-c(negative_nodes,new_negative)
		new_influenced<-union(new_positive,new_negative)

	}

	positive_nodes<-unique(positive_nodes)
	negative_nodes<-unique(negative_nodes)


	# compute spread

	if (!is.null(targets_activate) & is.null(targets_repress)){
		positive_nodes<-intersect(positive_nodes,targets_activate)
		negative_nodes<-c()
	}else if(is.null(targets_activate) & !is.null(targets_repress)){
		positive_nodes<-c()
		negative_nodes<-intersect(negative_nodes,targets_repress)
	}else if(!is.null(targets_activate) & !is.null(targets_repress)){
		positive_nodes<-intersect(positive_nodes,targets_activate)
		negative_nodes<-intersect(negative_nodes,targets_repress)
	}else{}

	activated_spread<-length(positive_nodes)
	repressed_spread<-length(negative_nodes)


	# return
	if (return_spread){
		return(list(activated_spread=activated_spread,repressed_spread=repressed_spread,activated=positive_nodes,repressed=negative_nodes))
	}else{
		return(list(activated_spread=activated_spread,repressed_spread=repressed_spread))
	}

}





#### SAME as IMCell, but returns other metrics useful for testing/benchmarking. Temp Function (similar to current expweighted function)
###
###
#' Greedy Algorithm for combined PRIM (spread = targets activated + targets repressed)
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, 'type' corresponding to interaction type,
#' @param kmax max number of TFs in solution set
#' @param repressor_wins See random_ICP
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
#' @param niter niter in MC
#' @param edge_probability_method method to assign edge probabilities
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores
#'
#'
im_combined_PRIM_greedy_testing<-function(ig,
								kmax=5,
								repressor_wins=FALSE,
								tfs=NULL,
								targets_activate=NULL,
								targets_repress=NULL,
								niter=1000,
								edge_probability_method="in_degree",
								edge_probability_multiplier=1,
								min_marginal_spread=5,
								return_spread=FALSE,
								num_cores=2,
								temptesting=TRUE){

	if(temptesting){
		return_spread=TRUE
	}

	# set cores
	num_cores<-min(num_cores,detectCores())


	# Set weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	# If no interaction type, assume all activation.
	if (is.null(E(ig)$type)){
		E(ig)$type<-1
		message("E(ig)$type is NULL, assuming all activation.")
	}

	solution_set<-c()
	prior_spread<-0

	final_activated_nodes<-c()
	final_repressed_nodes<-c()

	# for temptesting
	marg_spread_tracker<-c()

	for (k in 1:kmax){
		print(k)
		nodes<-setdiff(tfs,solution_set)

		best_spread<- -0.1 # set negative to ensure replacement
		best_node<-NULL

		# Run mc_ICP_parallel and add node with best spread
		for (node in nodes){
			print(node)
			mc_res<-mc_ICP_parallel(ig,c(solution_set,node),repressor_wins=repressor_wins,targets_activate=targets_activate,targets_repress=targets_repress,niter=niter,return_spread=return_spread,num_cores=num_cores)
			spread<-mc_res$activation_spread+mc_res$repression_spread

			if(spread > best_spread){
				best_spread<-spread
				best_node<-node
				
				activated_nodes<-mc_res$activated_targets
				repressed_nodes<-mc_res$repressed_targets
			}
		}

		if ((best_spread-prior_spread<min_marginal_spread) & (k>1)){
			break
		}

		cat("Marginal_spread: ",best_spread-prior_spread,'\n')

		if (temptesting){
			marg_spread_tracker<-c(marg_spread_tracker,(best_spread-prior_spread))
		}

		prior_spread<-best_spread
		solution_set<-c(solution_set,best_node)
		final_activated_nodes<-activated_nodes
		final_repressed_nodes<-repressed_nodes

	}

	if(temptesting){
		# compute total possible spread -- assume all targets provided are in the network
		possible_act<-length(targets_activate)
		possible_rep<-length(targets_repress)

		return(list(solution_set=solution_set,possible_act_spread=possible_act,possible_rep_spread=possible_rep,marg_spread_tracker=marg_spread_tracker,activated_targets=final_activated_nodes,repressed_targets=final_repressed_nodes))
	}

	if (return_spread){
		return(list(solution_set=solution_set,activated_targets=final_activated_nodes,repressed_targets=final_repressed_nodes))

	}else{
		return(list(solution_set=solution_set))
	}


}











# ===================================================================================================================
#
#										Expression-Based PRIM Functions
#
#
# ===================================================================================================================




# Based off of im_combined_PRIM_greedy_parallel... but add in expression-based node weighting

#' Greedy MC Algorithm for IMCell, with optional expression weighting
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, 'type' corresponding to interaction type
#' @param expDat expression data
#' @param sampTab sample table
#' @param source the starting cell type
#' @param target the target cell type
#' @param annotation_column the column in sampTab with cell identity annotation
#' @param node_weight_method method to weight the nodes
#' @param min_node_weight minimum assigned node weight
#' @param kmax max number of TFs in solution set
#' @param repressor_wins See random_ICP
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
#' @param niter niter in MC
#' @param edge_probability_method method to assign edge probabilities
#' @param edge_probability_multiplier optional factor to multiply edge probabilities by. 1 = no multiplier
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores
#'
#' @return 
#' @export
#'
IMCell_expweighted<-function(ig,
								expDat,
								sampTab,
								source,
								target,
								annotation_column="annotation_column",
								node_weight_method="weighted_expression_difference",
								min_node_weight=0.01,
								kmax=5,
								repressor_wins=FALSE,
								tfs=NULL,
								targets_activate=NULL,
								targets_repress=NULL,
								niter=1000,
								edge_probability_method="in_degree",
								edge_probability_multiplier=1,
								min_marginal_spread=10,
								return_spread=FALSE,
								num_cores=2,
								temptesting=FALSE){

	if(temptesting){
		return_spread=TRUE
	}

	# set cores
	num_cores<-min(num_cores,detectCores())


	# Compute Node Weights for Activation+Repression

	# Normalize/Transform data if not already done
	if (all(expDat == floor(expDat))){
		expDat<-trans_rnaseq(expDat,1e5)
  	}

  	# If weighting nodes by expression, ensure that GRN is limited to nodes within expression data
  	if (node_weight_method!="no_weight"){
  		v_remove<-setdiff(V(ig)$name,rownames(expDat))
  		ig<-delete_vertices(ig,v_remove)
  	}

	node_weights_activate<-compute_node_weights(V(ig)$name,expDat,sampTab,source,target,annotation_column=annotation_column,method=node_weight_method,min_node_weight=min_node_weight)
	node_weights_repress<-compute_node_weights(V(ig)$name,expDat,sampTab,target,source,annotation_column=annotation_column,method=node_weight_method,min_node_weight=min_node_weight)

	# Set edge weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	# If no interaction type, assume all activation.
	if (is.null(E(ig)$type)){
		E(ig)$type<-1
		message("E(ig)$type is NULL, assuming all activation.")
	}

	solution_set<-c()
	prior_spread<-0

	final_activated_nodes<-c()
	final_repressed_nodes<-c()

	# for temptesting
	marg_spread_tracker<-c()


	for (k in 1:kmax){
		print(k)
		nodes<-setdiff(tfs,solution_set)

		best_spread<- -0.1 # set negative to ensure replacement
		best_node<-NULL

		# Run mc_ICP_expweighted and add node with best spread
		for (node in nodes){
			print(node)
			mc_res<-mc_ICP_expweighted(ig,c(solution_set,node),node_weights_activate=node_weights_activate,node_weights_repress=node_weights_repress,repressor_wins=repressor_wins,targets_activate=targets_activate,targets_repress=targets_repress,niter=niter,return_spread=return_spread,num_cores=num_cores)
			spread<-mc_res$activation_spread+mc_res$repression_spread

			if(spread > best_spread){
				best_spread<-spread
				best_node<-node
				
				activated_nodes<-mc_res$activated_targets
				repressed_nodes<-mc_res$repressed_targets
			}
		}

		if ((best_spread-prior_spread<min_marginal_spread) & (k>1)){
			break
		}

		cat("Marginal_spread: ",best_spread-prior_spread,'\n')

		if (temptesting){
			marg_spread_tracker<-c(marg_spread_tracker,(best_spread-prior_spread))
		}

		prior_spread<-best_spread
		solution_set<-c(solution_set,best_node)
		final_activated_nodes<-activated_nodes
		final_repressed_nodes<-repressed_nodes

	}

	if(temptesting){
		# compute total possible spread
		possible_act<-sum(node_weights_activate[intersect(targets_activate,names(node_weights_activate))])
		possible_rep<-sum(node_weights_repress[intersect(targets_repress,names(node_weights_repress))])

		return(list(solution_set=solution_set,possible_act_spread=possible_act,possible_rep_spread=possible_rep,marg_spread_tracker=marg_spread_tracker,activated_targets=final_activated_nodes,repressed_targets=final_repressed_nodes))
	}

	if (return_spread){
		return(list(solution_set=solution_set,activated_targets=final_activated_nodes,repressed_targets=final_repressed_nodes))

	}else{
		return(list(solution_set=solution_set))
	}


}






# mc_ICP_parallel only compatible with random_ICP_parallel
# targets_activate/repress = NULL dealt with in random_ICP_parallel

# similar to random_ICP_parallel -- BUT spread is computed on weighted nodes. 

#' run MC cascade and return expected spread on IC-P model
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, and 'type' corresponding to interaction type
#' @param seed starting nodes to activate
#' @param repressor_wins see random_ICP for description
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
#' @param niter
#' @param return_spread whether or not to return the set of activated nodes
#'
#' @return 
#'
mc_ICP_expweighted<-function(ig,seed,node_weights_activate,node_weights_repress,repressor_wins=FALSE,targets_activate=NULL,targets_repress=NULL,niter=1000,return_spread=FALSE,num_cores=2){


	activated_tracker<-c()
	repressed_tracker<-c()


	rcs<-mclapply(1:niter,function(x){random_ICP_expweighted(ig=ig,seed=seed,node_weights_activate=node_weights_activate,node_weights_repress=node_weights_repress,repressor_wins=repressor_wins,
															targets_activate=targets_activate,targets_repress=targets_repress,
															return_spread=return_spread)},mc.cores=num_cores)

	total_activated<-lapply(rcs,function(x){x$activated_spread})
	total_activated<-sum(unlist(total_activated))
	total_repressed<-lapply(rcs,function(x){x$repressed_spread})
	total_repressed<-sum(unlist(total_repressed))

	predicted_num_activated<-total_activated/niter
	predicted_num_repressed<-total_repressed/niter

	if(return_spread){
		activated_tracker<-unlist(lapply(rcs,function(x){x$activated}))
		repressed_tracker<-unlist(lapply(rcs,function(x){x$repressed}))

		if (length(activated_tracker)==0){
			expected_activation<-c()
		}else{
			activated_tracker<-as.data.frame(table(activated_tracker))
			colnames(activated_tracker)<-c("gene","Freq")
			activated_tracker<-activated_tracker[order(activated_tracker$Freq,decreasing=TRUE),]
			expected_activation<-activated_tracker$gene[1:round(predicted_num_activated)]
		}

		if (length(repressed_tracker)==0){
			expected_repression<-c()
		}else{
			repressed_tracker<-as.data.frame(table(repressed_tracker))
			colnames(repressed_tracker)<-c("gene","Freq")
			repressed_tracker<-repressed_tracker[order(repressed_tracker$Freq,decreasing=TRUE),]
			expected_repression<-repressed_tracker$gene[1:round(predicted_num_repressed)]
		}

	}


	if(return_spread){
		mc_res<-list(activation_spread=predicted_num_activated,repression_spread=predicted_num_repressed,activated_targets=expected_activation,repressed_targets=expected_repression)
	}else{
		mc_res<-list(activation_spread=predicted_num_activated,repression_spread=predicted_num_repressed)
	}

	mc_res

}



# IC-P (independent cascade - polarity)
# Unlike random_simple_cascade, the repressed node carries propagates influence
# See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4111297/
# repressor_wins if FALSE, follows IC-P as described in ^^ (i.e. successful edge is chosen randomly from live edges/
#					edge flipping order is random and first successful edge influences target). If TRUE, if any repressive
#					edge is successfully flipped, then target will be repressed
# 
# almost the same as random_ICP, but with return_spread option
# and computes spread within this function (rather than in MC)

# similar to random_ICP_parallel -- BUT spread is computed on weighted nodes. 
random_ICP_expweighted<-function(ig,seed,node_weights_activate,node_weights_repress,repressor_wins=FALSE,targets_activate=NULL,targets_repress=NULL,return_spread=FALSE){
	# Three states a node can be in: positive (activated/on), negative (repressed/off), inactive (not influenced)
	# Three states an edge can be in: live (flipped), dead (failed to flip), unflipped

	positive_nodes<-seed
	negative_nodes<-c()

	new_positive<-seed
	new_negative<-c()
	new_influenced<-union(new_positive,new_negative)
	while(length(new_influenced)>0){
		# find not-yet-influenced neighbors
		nn<-lapply(adjacent_vertices(ig,new_influenced,mode="out"),function(x){x<-x$name;
																		  		x<-x[!(x %in% positive_nodes)];
																		  		x<-x[!(x %in% negative_nodes)];
																		  		x<-x[!(x %in% new_influenced)]})

		# keep only the elements with targets (i.e. remove non-targets)
		nn<-nn[lapply(nn,length)>0]

		all<-Reduce(c,nn)
		if (length(all)==0){
			break
		}

		# Create edge DF
		edges<-lapply(seq_along(nn),function(i){x<-as.vector(rbind(names(nn)[i],nn[[i]]));x})
		edges<-Reduce(c,edges)

		edges_df<-do.call(rbind,lapply(1:length(nn),function(i){data.frame(origin=rep(names(nn)[i],length(nn[[i]])),target=nn[[i]])}))
		edges_df$origin_type<-ifelse(edges_df$origin %in% new_positive,1,-1)
		edges_df$type<-E(ig,P=edges)$type

		edges_df$target_polarity<-edges_df$origin_type * edges_df$type

		# coin flips
		flips<-runif(length(all),min=0,max=1)
		# which flips were successful (which flips are less than the edge weight)
		# e.g. if prob = .8 --> generate runif between 0 and 1 --> if falls between 0-.8, call success; if between .8-1 call fail
		edges_df$success<-flips < E(ig,P=edges)$weight

		edges_df<-edges_df[edges_df$success,]

		if (nrow(edges_df)==0){
			break
		}

		if (!repressor_wins){
			# Traditional IC-P model = only one edge influences a give target
			# Order of edge flipping is random
			# AKA succesful edge is chosen uniformly

			# scramble the order, then limit to first live edge
			edges_df<-edges_df[sample(1:nrow(edges_df)),]
			edges_df<-edges_df[!duplicated(edges_df$target),]

			new_positive<-edges_df$target[edges_df$target_polarity==1]
			new_negative<-edges_df$target[edges_df$target_polarity==-1]

		}else{
			# Alternative combinatorial rules:
			# If any repressive edge is live, the target is repressed and assigned negative state

			new_negative<-edges_df$target[edges_df$target_polarity==-1]
			new_positive<-edges_df$target[edges_df$target_polarity==1]
			new_positive<-setdiff(new_positive,new_negative)

		}

		positive_nodes<-c(positive_nodes,new_positive)
		negative_nodes<-c(negative_nodes,new_negative)
		new_influenced<-union(new_positive,new_negative)

	}

	positive_nodes<-unique(positive_nodes)
	negative_nodes<-unique(negative_nodes)


	# compute spread

	if (!is.null(targets_activate) & is.null(targets_repress)){
		positive_nodes<-intersect(positive_nodes,targets_activate)
		negative_nodes<-c()

		activated_spread<-sum(node_weights_activate[intersect(positive_nodes,names(node_weights_activate))])
		repressed_spread<-0

	}else if(is.null(targets_activate) & !is.null(targets_repress)){
		positive_nodes<-c()
		negative_nodes<-intersect(negative_nodes,targets_repress)

		activated_spread<-0
		repressed_spread<-sum(node_weights_repress[intersect(negative_nodes,names(node_weights_repress))])

	}else if(!is.null(targets_activate) & !is.null(targets_repress)){
		positive_nodes<-intersect(positive_nodes,targets_activate)
		negative_nodes<-intersect(negative_nodes,targets_repress)

		activated_spread<-sum(node_weights_activate[intersect(positive_nodes,names(node_weights_activate))])
		repressed_spread<-sum(node_weights_repress[intersect(negative_nodes,names(node_weights_repress))])

	}else{
		activated_spread<-sum(node_weights_activate[intersect(positive_nodes,names(node_weights_activate))])
		repressed_spread<-sum(node_weights_repress[intersect(negative_nodes,names(node_weights_repress))])

	}


	# return
	if (return_spread){
		return(list(activated_spread=activated_spread,repressed_spread=repressed_spread,activated=positive_nodes,repressed=negative_nodes))
	}else{
		return(list(activated_spread=activated_spread,repressed_spread=repressed_spread))
	}

}




# ===================================================================================================================
#
#										Functions to rank TFs just based on spread
#
#
# ===================================================================================================================
# These functions return a ranked list. 
# (i.e. the equivalent of one loop of IM or PRIM)

#' Compute marginal spread with combined PRIM rules (spread = targets activated + targets repressed)
#'
#' @param ig an igraph object with normalized edge attribute 'weight' corresponding to probabilities, 'type' corresponding to interaction type,
#' @param kmax max number of TFs in solution set
#' @param repressor_wins See random_ICP
#' @param tfs limit search radius to specific TFs. If NULL, search radius is not limited.
#' @param targets_activate consider activation of specified targets only. If NULL activation is based on total spread
#' @param targets_repress consider repression of specified targets only. If NULL activation is based on total spread
#' @param niter niter in MC
#' @param edge_probability_method method to assign edge probabilities
#' @param return_spread whether or not to return the set of activated nodes
#' @param num_cores number of cores
#'
#' @return 
#' @export
#'
rank_by_spread<-function(ig,
						repressor_wins=FALSE,
						tfs=NULL,
						targets_activate=NULL,
						targets_repress=NULL,
						niter=1000,
						edge_probability_method="in_degree",
						edge_probability_multiplier=1,
						return_spread=FALSE,
						num_cores=2){


	# set cores
	num_cores<-min(num_cores,detectCores())


	# Set weights and normalize
	if (is.null(E(ig)$weight)){
		E(ig)<-0.5
	}
	ig<-weights_to_probability(ig,method=edge_probability_method,multiplier=edge_probability_multiplier)

	if (is.null(tfs)){
		tfs<-V(ig)[degree(ig, mode = 'out')>0]$name
	}else{
		tfs<-intersect(tfs,V(ig)[degree(ig, mode = 'out')>0]$name)
	}
	print(paste0("TFs: ",length(tfs)))

	# If no interaction type, assume all activation.
	if (is.null(E(ig)$type)){
		E(ig)$type<-1
		message("E(ig)$type is NULL, assuming all activation.")
	}


	sort_TFs<-data.frame(TF=character(),spread=numeric())
	for (node in tfs){
		print(node)
		mc_res<-mc_ICP_parallel(ig,node,repressor_wins=repressor_wins,targets_activate=targets_activate,targets_repress=targets_repress,niter=niter,return_spread=return_spread,num_cores=num_cores)
		spread<-mc_res$activation_spread+mc_res$repression_spread

		sort_TFs<-rbind(sort_TFs,data.frame(TF=node,spread=spread))

	}

	sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]
	sort_TFs$rank<-1:nrow(sort_TFs)

	if (return_spread){
		return(sort_TFs)

	}else{
		sort_TFs$spread<-NULL
		return(sort_TFs)
	}

}





