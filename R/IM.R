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
#' @param seed starting nodes to activate
#' @param targets consider activation of specified targets only. If NULL activation is based on total spread
#' @param niter
#' @param min_marginal_spread the minimum marginal spread. Greedy algorithm stops when marginal spread < minimum marginal spread. If minimum_marginal_spread=1, algorithm doesn't terminate until k=kmax or spread reaches full coverage
#' @param return_spread whether or not to return the set of activated nodes
#'
#' @return 
#'
im_celf_parallel<-function(ig,
					kmax=5,
					tfs=NULL,
					targets=NULL,
					niter=1000,
					edge_probability_method="in_degree",
					edge_probability_multiplier=1,
					min_marginal_spread=2,
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
	for (node in tfs){
		print(node)
		mc_res<-mc_activation_cascade_parallel(ig,node,targets=targets,niter=niter,return_spread=return_spread,num_cores=num_cores)
		spread<-mc_res$predicted_spread

		sort_TFs<-rbind(sort_TFs,data.frame(TF=node,spread=spread))
	}

	sort_TFs<-sort_TFs[order(sort_TFs$spread,decreasing=TRUE),]

	# add the first TF to the solution set, set the prior spread
	solution_set<-c(solution_set,sort_TFs$TF[1])
	prior_spread<-sort_TFs$spread[1]

	# remove any TFs that have spread < min_marginal_spread -- these aren't worth testing any further
	sort_TFs<-sort_TFs[sort_TFs$spread>=min_marginal_spread,]

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
#'
#' @return 
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
#'
#' @return 
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




# Temporary testing function that adds in marginal spread tracker, otherwise same function as im_celf_parallel

im_celf_parallel_testing<-function(ig,
					kmax=5,
					tfs=NULL,
					targets=NULL,
					niter=1000,
					edge_probability_method="in_degree",
					edge_probability_multiplier=1,
					min_marginal_spread=2,
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

