########################################################################
#
#	IMCell setup functions
#
########################################################################


# ===================================================
#
# 		Converting Weights to Probabilities
#
# ===================================================


#' Converts edge weights to probabilities
#'
#' @param ig a directed igraph object
#' @param method how to convert weights to probability. Either in_degree, total_in_degree, or coin_flip.
#' @param multiplier a multiplier to apply to probabilities. Ignored if method is coin_flip.
#' @param p probability if method is coin_flip, ignored otherwise.
#' @param pmax value at which to cap the probability.
#'
#' @return updated igraph object with weight now updated as probabilities
#' 
#' @export
#'
weights_to_probability<-function(ig,method="in_degree",multiplier=1,p=0.5,pmax=0.8){
	if (method=="in_degree"){
		if (sum(E(ig)$weight<0)>0){
			temp_ig<-ig
			E(temp_ig)$weight<-abs(E(temp_ig)$weight)
			E(ig)$weight<-E(ig)$weight/strength(temp_ig,mode="in")[get.edgelist(ig)[,2]]
		}else{
			E(ig)$weight<-E(ig)$weight/strength(ig,mode="in")[get.edgelist(ig)[,2]]
		}
	}else if (method=="total_in_degree"){
		numTFs<-length(V(ig)[degree(ig, mode = 'out')>0]$name)
		numgenes<-length(V(ig))

		if (sum(E(ig)$weight<0)>0){
			temp_ig<-ig
			E(temp_ig)$weight<-abs(E(temp_ig)$weight)
			total_degree<-sum(strength(temp_ig,mode="in"))
			E(ig)$weight<-(E(ig)$weight/total_degree) * numgenes
		}else{
			total_degree<-sum(strength(ig,mode="in"))
			E(ig)$weight<-(E(ig)$weight/total_degree) * numgenes
		}

	}else if (method=="out_degree"){



	}else if (method=="total_out_degree"){


	}else if (method=="coin_flip"){
		# Assign all edges the same probability, p
		E(ig)$weight<-p
	}


	# Set multiplier
	E(ig)$weight<-E(ig)$weight*multiplier

	# Cap probability at pmax
	E(ig)$weight[E(ig)$weight > pmax]<-pmax


	ig
}


# ===================================================
#
# 		Finding differential nodes to target
#
# ===================================================


# Given starting and target cell states, finds genes to activate and repress
# For now goes by DESeq2 or t-test

#' Performs DE via negative binomial method on TFs using DESeq2
#' 
#' @param expMat Genes-by-samples expression matrix
#' @param sampTab Metadata/phenodata of the samples
#' @param source Starting cell type or background for which to compute DE
#' @param target Target cell type
#' @param limit_to limit analysis to a subset of genes, e.g. TFs
#' @param pThresh adjusted p-value threshold
#' @param lfcThresh log-fold-change threshold
#' @param topX limit to_activate and to_repress to the top X genes based on LFC
#' @param min_de_targets if DESeq2 returns fewer targets than this threshold, T-test will be run instead
#'
#'
#' @return A list of genes to activate and to repress
find_differential_nodes<-function(expMat,
									sampTab,
									source,
									target,
									annotation_column,
									limit_to=NULL,
									pThresh=0.05,
									lfcThresh=1,
									topX=NULL,
									min_de_targets=5){
	
	#require(DESeq2)

	if (!is.null(limit_to)){
		expMat<-expMat[intersect(rownames(expMat),limit_to),]
	}

	sampTab$condition<-sampTab[,annotation_column]
	sampTab<-sampTab[sampTab$condition %in% c(target,source),]
  	expMat<-expMat[,rownames(sampTab)]
  	expMat<-expMat[rowSums(expMat)!=0,]

  	de<-tryCatch(
	    expr = {
	        dds<-DESeqDataSetFromMatrix(countData=expMat,colData=sampTab,design=~condition)
		  	if(sum(rowSums(expMat==0)>0)==nrow(expMat)){
			    #dds<-estimateSizeFactors(dds, type = 'iterate')
			    geoMeans <- apply(expMat, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
			    dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
		  	}
		  	dds<-DESeq(dds)

		  	diffres<-results(dds,contrast=c("condition",target,source),independentFiltering=FALSE)
		  	diffres<-as.data.frame(diffres)

		  	diffres$TF<-rownames(diffres)

		  	# Filter for significance and logfoldchange
		  	diffres<-diffres[diffres$padj<pThresh,]
		  	diffres<-diffres[abs(diffres$log2FoldChange)>lfcThresh,]


		  	# Nodes to activate/repress
		  	if(!is.null(topX)){
		  		to_activate<-diffres[diffres$log2FoldChange>0,]
		  		to_repress<-diffres[diffres$log2FoldChange<0,]

		  		to_activate<-to_activate[order(to_activate$log2FoldChange,decreasing=TRUE),]
		  		to_repress<-to_repress[order(to_repress$log2FoldChange,decreasing=FALSE),]

		  		n_act<-min(nrow(to_activate),topX)
		  		n_rep<-min(nrow(to_repress),topX)

		  		to_activate<-to_activate$TF[1:n_act]
		  		to_repress<-to_repress$TF[1:n_rep]

		  	}else{
		  		to_activate<-diffres$TF[diffres$log2FoldChange>0]
		  		to_repress<-diffres$TF[diffres$log2FoldChange<0]
		  	}

		  	list(to_activate=to_activate,to_repress=to_repress)

	    },
	    error = function(e){ 
	    	message("DESeq2 error; Running t-test...")
	        # If DEseq2 errors, normalize and run ttest instead
	        expNorm<-downsample_and_transform(expMat,observations_by_features=FALSE)

	        g1<-expNorm[,rownames(sampTab)[sampTab$condition==target]]
	        g2<-expNorm[,rownames(sampTab)[sampTab$condition==source]]

	        diffres<-data.frame(gene=character(),start_mean=double(),target_mean=double(),mean_diff=double(),pval=double())

	        for (gene in rownames(expNorm)){
				t<-t.test(g1[gene,],g2[gene,])
				ans<-data.frame(gene=gene,start_mean=t$estimate[2],target_mean=t$estimate[1],mean_diff=(t$estimate[1]-t$estimate[2]),pval=t$p.value)
				diffres<-rbind(diffres,ans)
			}

			diffres$padj<-p.adjust(diffres$pval,method="BH")
			diffres<-diffres[order(diffres$padj,decreasing=FALSE),]

			# Filter for significance and difference
			diffres<-diffres[diffres$padj<pThresh,]
			diffres<-diffres[abs(diffres$mean_diff)>lfcThresh,]

			rownames(diffres)<-diffres$gene

			if(!is.null(topX)){
		  		to_activate<-diffres[diffres$mean_diff>0,]
		  		to_repress<-diffres[diffres$mean_diff<0,]

		  		to_activate<-to_activate[order(to_activate$mean_diff,decreasing=TRUE),]
		  		to_repress<-to_repress[order(to_repress$mean_diff,decreasing=FALSE),]

		  		n_act<-min(nrow(to_activate),topX)
		  		n_rep<-min(nrow(to_repress),topX)

		  		to_activate<-to_activate$gene[1:n_act]
		  		to_repress<-to_repress$gene[1:n_rep]

		  	}else{
		  		to_activate<-diffres$gene[diffres$mean_diff>0]
		  		to_repress<-diffres$gene[diffres$mean_diff<0]
		  	}

		  	return(list(to_activate=to_activate,to_repress=to_repress))
	    }
	)


  	# Stupid hacky thing right now. If DESeq returns too few targets, run t-test instead
  	if (length(de$to_activate)<min_de_targets){
  		message("Minimal targets; Running t-test...")

        expNorm<-downsample_and_transform(expMat,observations_by_features=FALSE)

        g1<-expNorm[,rownames(sampTab)[sampTab$condition==target]]
        g2<-expNorm[,rownames(sampTab)[sampTab$condition==source]]

        diffres<-data.frame(gene=character(),start_mean=double(),target_mean=double(),mean_diff=double(),pval=double())

        for (gene in rownames(expNorm)){
			t<-t.test(g1[gene,],g2[gene,])
			ans<-data.frame(gene=gene,start_mean=t$estimate[2],target_mean=t$estimate[1],mean_diff=(t$estimate[1]-t$estimate[2]),pval=t$p.value)
			diffres<-rbind(diffres,ans)
		}

		diffres$padj<-p.adjust(diffres$pval,method="BH")
		diffres<-diffres[order(diffres$padj,decreasing=FALSE),]

		# Filter for significance and difference
		diffres<-diffres[diffres$padj<pThresh,]
		diffres<-diffres[abs(diffres$mean_diff)>lfcThresh,]

		rownames(diffres)<-diffres$gene

		if(!is.null(topX)){
	  		to_activate<-diffres[diffres$mean_diff>0,]
	  		to_repress<-diffres[diffres$mean_diff<0,]

	  		to_activate<-to_activate[order(to_activate$mean_diff,decreasing=TRUE),]
	  		to_repress<-to_repress[order(to_repress$mean_diff,decreasing=FALSE),]

	  		n_act<-min(nrow(to_activate),topX)
	  		n_rep<-min(nrow(to_repress),topX)

	  		to_activate<-to_activate$gene[1:n_act]
	  		to_repress<-to_repress$gene[1:n_rep]

	  	}else{
	  		to_activate<-diffres$gene[diffres$mean_diff>0]
	  		to_repress<-diffres$gene[diffres$mean_diff<0]
	  	}

	  	de<-list(to_activate=to_activate,to_repress=to_repress)
  	}

  	return(de)

}


# ===================================================
#
# 		Expression-based Node Weighting
#
# ===================================================


# Compute node weights based on expression
# This gives node weights with the assumption of nodes you want to activate 
# To compute separate weights for nodes to repress, swap source and target when called in the driver/main function

#' Compute node weights based on expression
#' 
#' @param nodes
#' @param expNorm Normalized and transformed expression data
#' @param sampTab Metadata/phenodata of the samples
#' @param source Starting cell type or background for which to compute DE
#' @param target Target cell type
#' @param annotation_column column in sampTab with cell identity annotation
#' @param method either target_mean_expression or weighted_expression_difference
#' @param min_node_weight the minimum node weight
#'
compute_node_weights<-function(nodes,
								expNorm,
								sampTab,
								source,
								target,
								annotation_column="cf_annotation",
								method="weighted_expression_difference",
								min_node_weight=0.01){


	tvals<-cn_make_tVals(expNorm, sampTab, dLevel=annotation_column)

	if (method=="target_mean_expression"){
		meanVect<-unlist(tvals[[target]][['mean']][nodes])
		nodeweights<-(2**meanVect)/sum(2**meanVect)
		nodeweights<-nodeweights*length(nodeweights)

		nodeweights[is.na(nodeweights)]<-min_node_weight
		nodeweights[nodeweights<0]<-min_node_weight


	}else if (method=="weighted_expression_difference"){
		# nodes score higher if they are not or lowly expressed in source but strong in target
		# similar to first half of CellNet's NIS -- zscore in source samples should be negative, and high weight in target samples

		# weights based on target cell type
		meanVect<-unlist(tvals[[target]][['mean']][nodes])
		tweights<-(2**meanVect)/sum(2**meanVect)

		# zscore expression in source cell type, in the context of target (i.e. mean and sd defined in target, where does source expression compare?)
		st_source<-sampTab[sampTab[,annotation_column]==source,]

		zmat<-matrix(0,nrow=length(nodes),ncol=nrow(st_source))
		for (i in seq(nrow(st_source))){
			sid<-rownames(st_source)[i]
			xvals<-as.vector(expNorm[nodes,sid])
			names(xvals)<-nodes

			zmat[,i]<-cn_zscoreVect(nodes,xvals,tvals,target)
		}
		rownames(zmat)<-nodes
		colnames(zmat)<-rownames(st_source)

		# multiply z-score and weight
		weighted_zmat<-zmat*as.vector(tweights)

		nodeweights<-rowMeans(weighted_zmat)*length(tweights)
		nodeweights<-nodeweights * -1

		nodeweights[is.na(nodeweights)]<-min_node_weight
		nodeweights[nodeweights<0]<-min_node_weight

	}else if (method=="no_weight"){
		nodeweights<-rep(min_node_weight,length(nodes))
		names(nodeweights)<-nodes

	}

	nodeweights

}


# ===================================================
#
# 		Functions for limiting TFscope
#
# ===================================================


#' Limit TF scope and search radius
#' 
#' @param expCounts Raw counts expression data
#' @param sampTab Metadata/phenodata of the samples
#' @param tfs list of TFs
#' @param source Starting cell type or background for which to compute DE
#' @param target Target cell type
#' @param annotation_column column in sampTab with cell identity annotation
#' @param sourceThresh Filter out TFs if they are expressed above this level in the source cell type
#' @param DE_lfcThresh Log-fold-change threshold for DE approach
#' @param DE_pThresh adjusted p-value threshold for DE approach
#' @param DE_pval_topX number of top TFs to keep by DE p-value
#' @param DE_lfc_topX number of top TFs to keep by DE LFC
#' @param TM_pThresh adjusted p-value threshold for TM approach
#' @param TM_cThresh corr threshold for TM approach
#' @param TM_topX number of top TFs to keep by TM appraoch
#' @param TM_topX_sourceThresh Filter out if considered spec gene in source type
#'
#' @return narrowed TF search radius
tfscope<-function(expCounts,sampTab,tfs,source,target,annotation_column="celltype",sourceThresh=NULL,
					DE_lfcThresh=1,DE_pThresh=0.05,DE_pval_topX=20,DE_lfc_topX=20,
					TM_pThresh=0.05,TM_cThresh=0.1,TM_topX=20,TM_topX_sourceThresh=50){


	#------------------------- RUN DE approach ------------------------------
	require(DESeq2)

	expMat<-expCounts[intersect(tfs,rownames(expCounts)),]

	st<-sampTab
    st$condition<-st[,annotation_column]
    st<-st[st$condition %in% c(target,source),]
    expMat<-expMat[,rownames(st)]
    expMat<-expMat[rowSums(expMat)!=0,]

    expMat<-expMat[,colSums(expMat)!=0]
    st<-st[colnames(expMat),]

    dds<-DESeqDataSetFromMatrix(countData=expMat,colData=st,design=~condition)
    if(sum(rowSums(expMat==0)>0)==nrow(expMat)){
        #dds<-estimateSizeFactors(dds, type = 'iterate')
        geoMeans <- apply(expMat, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
        dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
    }
    dds<-DESeq(dds)

    diffres<-results(dds,contrast=c("condition",target,source),independentFiltering=FALSE)
    diffres<-as.data.frame(diffres)
    diffres$TF<-rownames(diffres)

    # Filter for significance and logfoldchange
    diffres<-diffres[diffres$padj<DE_pThresh,]
    diffres<-diffres[abs(diffres$log2FoldChange)>DE_lfcThresh,]
    # Filter for TFs with greater expression in target, and low expression in start
    # LFC>0
    diffres<-diffres[diffres$log2FoldChange>0,]

   

    #------------------------- RUN TM approach ------------------------------

    expNorm<-trans_rnaseq(expCounts,1e5)
    expNorm<-expNorm[intersect(rownames(expNorm),tfs),]
    specgenes<-adapted_TM(expNorm, sampTab, pThresh=TM_pThresh, cThresh=TM_cThresh, annotation_column=annotation_column)
    specgenes<-lapply(specgenes,function(x){x<-x[order(x$cval,decreasing=TRUE),];x})

    specGenes_target<-specgenes[[target]]
    specGenes_source<-specgenes[[source]]

    # remove from target any TFs that are within top X of source
    toremove<-rownames(specGenes_source)[1:min(TM_topX_sourceThresh,nrow(specGenes_source))]
    specGenes_target<-specGenes_target[!(rownames(specGenes_target) %in% toremove),]




    #------------------------- If filtering by source expression ------------------------------
    if (!is.null(sourceThresh)){
    	start_exp<-expMat[,rownames(sampTab)[sampTab[annotation_column]==source]]
    	keep<-rownames(start_exp)[rowMeans(start_exp) < sourceThresh]

    	diffres<-diffres[intersect(keep,rownames(diffres)),]
    	specGenes_target<-specGenes_target[intersect(keep,rownames(specGenes_target)),]
    }


    #------------------------- Summarize results ------------------------------

    # Order DE results by pval/LFC, extract TFs to keep
    diffres<-diffres[order(diffres$padj,decreasing=FALSE),]
    tfs_de_pval<-rownames(diffres)[1:min(DE_pval_topX,nrow(diffres))]

    diffres<-diffres[order(diffres$log2FoldChange,decreasing=TRUE),]
    tfs_de_lfc<-rownames(diffres)[1:min(DE_lfc_topX,nrow(diffres))]

    # Order TM results by cval, extract TFs to keep
    specGenes_target<-specGenes_target[order(specGenes_target$cval,decreasing=TRUE),]
    tfs_tm<-rownames(specGenes_target)[1:min(TM_topX,nrow(specGenes_target))]


    tfscope<-union(tfs_de_pval,tfs_de_lfc)
    tfscope<-union(tfscope,tfs_tm)

    tfscope

}




# Some template matching functions adapted from CellNet's specGenes approach
#' Template matching, adapted
#'
#' @param expDat expression data
#' @param sampTab sample table
#' @param pThresh 
#' @param cThresh
#' @param annotation_column
#'
adapted_TM<-function(expDat,sampTab,pThresh=0.05,cThresh=0.1,annotation_column="celltype"){
    myPatternG<-TM_sampR_to_pattern(as.vector(sampTab[,annotation_column])); 
    specificSets<-apply(myPatternG, 1, TM_testPattern, expDat=expDat);
    
    # Filter
    specificSets<-lapply(specificSets,function(x){x<-x[x$holm<pThresh,];x<-x[x$cval>cThresh,];x})

}

#' Template matching, pattern
#'
#' @param sampR
#'
TM_sampR_to_pattern<-function(sampR){
    d_ids<-unique(as.vector(sampR)); 
    nnnc<-length(sampR);
    ans<-matrix(nrow=length(d_ids), ncol=nnnc);

    for(i in seq(length(d_ids))){ 
        x<-rep(0,nnnc); 
        x[which(sampR==d_ids[i])]<-1;
        ans[i,]<-x;
    }

    colnames(ans)<-as.vector(sampR); 
    rownames(ans)<-d_ids; 

    ans
} 

#' Template matching, testing
#'
#' @param pattern
#' @param expDat
#'
TM_testPattern<-function(pattern, expDat){ 
    pval<-vector(); 
    cval<-vector();
    geneids<-rownames(expDat);
    llfit<-ls.print(lsfit(pattern, t(expDat)), digits=25, print=FALSE);
    xxx<-matrix( unlist(llfit$coef), ncol=8,byrow=TRUE);
    ccorr<-xxx[,6];
    cval<- sqrt(as.numeric(llfit$summary[,2])) * sign(ccorr);
    pval<-as.numeric(xxx[,8]);

    holm<-p.adjust(pval, method='holm'); 
    data.frame(row.names=geneids, pval=pval, cval=cval,holm=holm); 
}










# ===================================================
#
# 	Helpful support functions taken from CellNet
#
# ===================================================


# tVals list of ct->mean->named vector of average gene expression, ->sd->named vector of gene standard deviation

#' From CellNet: Estimate gene expression dist in CTs
#'
#' Calculate mean and SD 
#' @param expDat training data 
#' @param sampTab, ### training sample table 
#' @param dLevel="description1", ### column to define CTs
#' @param predictSD=FALSE ### whether to predict SD based on expression level
#'
cn_make_tVals<-function### estimate gene expression dist in CTs
(expDat, ### training data 
 sampTab, ### training sample table 
 dLevel="description1", ### column to define CTs
 predictSD=FALSE ### whether to predict SD based on expression level
){
  
  if(predictSD){
    ans<-cn_make_tVals_predict(expDat, sampTab, dLevel);
  }
  else{
    # Note: returns a list of dName->gene->mean, sd, where 'dName' is a ctt or lineage 
    # make sure everything is lined up
    expDat<-expDat[,rownames(sampTab)];
    tVals<-list();
    dNames<-unique(as.vector(sampTab[,dLevel]));
    allGenes<-rownames(expDat);
    for(dName in dNames){
      #cat(dName,"\n");
      xx<-which(sampTab[,dLevel]==dName);
      sids<-rownames(sampTab[xx,]);
      xDat<-expDat[,sids];
      means<-apply(xDat, 1, mean);
      sds<-apply(xDat, 1, sd);
      tVals[[dName]][['mean']]<-as.list(means);
      tVals[[dName]][['sd']]<-as.list(sds);
    }
    ans<-tVals;
  }
  ans;
}


#' From CellNet: caculate means and SDs
#'
#' Calculate mean and SD 
#'
cn_zscoreVect<-function
### Compute the mean zscore of given genes in each sample
(genes,
 ### genes
 xvals,
 ### named vector
 tVals,
 ### tvals
 ctt
 ### ctt
 ){
  ans<-vector();
  for(gene in genes){
    ans<-append(ans, cn_zscore(xvals[gene], tVals[[ctt]][['mean']][[gene]], tVals[[ctt]][['sd']][[gene]]));
  }
  ans;
  ### zscore vector
}


# Taken from CellNet's "zscore" function
#' From CellNet: zscore
#'
#' Calculate mean and SD 
#'
cn_zscore<-function
### compute zscore
(x,
 ### numeric vector
 meanVal, 
 ### mean of distribution to compute zscore of x against 
 sdVal
 ### standard deviation of distribution to compute zscore of x agains
 ){ 
  (x-meanVal)/sdVal;
  ### zscore
}


#' weighted subtraction from mapped reades and log applied to all
#'
#' Simulate expression profile of  _total_ mapped reads
#' @param expRaw matrix of total mapped reads per gene/transcript
#' @param total numeric post transformation sum of read counts
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







