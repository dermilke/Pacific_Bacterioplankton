distance_wrapper <- function(phyloseq_obj, method = "TINA_weighted", size.thresh = 1, pseudocount = 10^-6, 
                             nblocks = 100, use.cores = 4, cor.use = "na.or.complete") {
  
  #source("Functions_Similarity_Indices.R")
  
  # For ASV x ASV Distances:
  
  #SparCC Co-Occurrences                                              -> "SparCC"
  #Phylogenetic Distances                                             -> "Phylo_Dist"
  
  # For Sample x Sample Distances:
  
  #Jaccard index, classical, unweighted                               -> "jaccard"
  #Jaccard index, classical, weighted                                 -> "jaccard.abd"
  #Jaccard index, classical, weighted, fraction-based                 -> "jaccard.abd.frac"
  #Jaccard index, classical, weighted, fraction-based, Chao's version -> "jaccard.abd.frac_chao"
  #Jaccard index, classical, weighted, alternative formula            -> "jaccard.abd.alt"
  #Jaccard index, classical, Chao' abundance correction               -> "jaccard.abd.chao"
  #Bray-Curtis, classical                                             -> "bray_curtis"
  #Morisita-Horn, classical                                           -> "morisita_horn"
  
  #Unweighted UniFrac   -> "UniFrac_unweighted"
  #Weighted UniFrac     -> "UniFrac_weighted"
  
  #unweighted TINA      -> "TINA_unweighted"
  #weighted TINA        -> "TINA_weighted"
  
  #unweighted PINA      -> "PINA_unweighted"
  #weighted PINA        -> "PINA_weighted"
  
  doMC::registerDoMC(cores = use.cores)
  
  SparCC_wrapper <- function(phyloseq_obj, size.thresh, pseudocount,
                             nblocks, use.cores, cor.use) {
    
    Count_Table <- phyloseq::otu_table(phyloseq_obj) %>%
      Matrix::Matrix()
    
    Count_Table_Filtered <- Count_Table[!(rowSums(Count_Table) < size.thresh), ] + pseudocount
    n.otu <- nrow(Count_Table_Filtered)
    
    #Preallocate blocks for parallel processing & Aitchinson's T matrix
    #=> based on https://gist.github.com/bobthecat/5024079
    nblocks <- floor(nrow(Count_Table_Filtered) / use.cores)
    size.split <- use.cores
    my.split <- list()
    length(my.split) <- nblocks
    
    my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split))
    my.split[[nblocks]] <- (size.split*(nblocks-1)+1):nrow(Count_Table_Filtered)
    
    dat.split <- parallel::mclapply(my.split, function(g) {Count_Table_Filtered[g,]}, mc.cores=use.cores);
    
    #Get combinations of splits
    my.combs <- expand.grid(1:length(my.split), 1:length(my.split)) %>%
      apply(., 1, sort) %>%
      t() %>%
      unique()
    
    #Preallocate Aitchinson's T matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
    mat.T <- bigmemory::big.matrix(nrow=nrow(Count_Table_Filtered), ncol=nrow(Count_Table_Filtered), 
                                   dimnames=c(rownames(Count_Table_Filtered), rownames(Count_Table_Filtered)), 
                                   shared=T)
    mat.T.desc <- bigmemory::describe(mat.T)
    
    #Compute Aitchinson's T matrix
    results <- foreach(i = 1:nrow(my.combs)) %dopar% {
      #Get current combination
      curr.comb <- my.combs[i, ];
      #Get current data
      g.1 <- my.split[[curr.comb[1]]];
      g.2 <- my.split[[curr.comb[2]]];
      dat.1 <- dat.split[[curr.comb[1]]];
      dat.2 <- dat.split[[curr.comb[2]]];
      #Get current part of Aitchinson's matrix
      curr.T <- apply(dat.1, 1, function(x) {apply(dat.2, 1, function(y) {var(log(x/y))})});
      #Store
      curr.mat.T <- attach.big.matrix(mat.T.desc);
      curr.mat.T[g.2, g.1] <- curr.T;
      curr.mat.T[g.1, g.2] <- t(curr.T);
      #Return
      TRUE
    }
    
    #Compute component variations ("t_i")
    var.t <- colsum(mat.T);
    
    #Estimate component variances ("omega_i") from t_i by solving a linear equation system
    mat.a <- matrix(data=1, nrow=n.otu, ncol=n.otu)
    diag(mat.a) <- n.otu-1
    omega <- sqrt(solve(a=mat.a, b=var.t))
    
    #Estimate pairwise correlations based on these values
    global.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %do% {
      (omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)
    }
    rownames(global.sparcc) <- colnames(global.sparcc) <- rownames(Count_Table_Filtered)
    
    tmp.S <- cor.par(global.sparcc, method="pearson", use=cor.use, use.cores=use.cores)
    S.sparcc <- 0.5 * (tmp.S + 1);
    
    return(S.sparcc)
    
  }
  
  Phylogeny_wrapper <- function(phyloseq_obj, use.cores, cor.use) {
    
    ot <- phyloseq::otu_table(phyloseq_obj)
    
    my.tree <- phyloseq::phy_tree(phyloseq_obj)
    my.tree$node.label <- c("Root", paste("Node", 2:my.tree$Nnode, sep = "_"))
    
    global.cophenetic_tree <- fast.cophenetic.phylo(my.tree, use.cores=use.cores)
    global.cophenetic_tree <- global.cophenetic_tree[rownames(ot), rownames(ot)]
    
    #Calculate derived OTU similarity matrix S
    tmp.S <- cor.par(global.cophenetic_tree, method="pearson", use=cor.use, use.cores=use.cores)
    S.phylo <- 0.5 * (tmp.S + 1)
    
    return(S.phylo)
    
  }
  
  if (method == "SparCC") {
    
    my.cs <- SparCC_wrapper(phyloseq_obj, size.thresh, pseudocount, nblocks, use.cores, cor.use)
    
  } else if (method == "Phylo_Dist") {
    
    my.cs <- Phylogeny_wrapper(phyloseq_obj, use.cores, cor.use)
    
  } else if (method %in% c("jaccard", "jaccard.abd", "jaccard.abd.frac",
                           "jaccard.abd.frac_chao", "jaccard.abd.alt", "jaccard.abd.chao",
                           "bray_curtis", "morisita_horn")) {
    
    my.cs <- community.similarity.par(phyloseq::otu_table(phyloseq_obj), distance=method, use.cores=use.cores)
    my.cs[my.cs < 0] <- 0
    
  } else if (method %in% c("UniFrac_unweighted", "UniFrac_weighted")) {
    
    my.cs <- fasterUniFrac(phyloseq::otu_table(phyloseq_obj), 
                           phyloseq::phy_tree(phyloseq_obj), 
                           distance=method, blocksize=nblocks, use.cores=use.cores)
    
  } else if (method %in% c("TINA_unweighted", "TINA_weighted")) {
    
    S.sparcc <- SparCC_wrapper(phyloseq_obj, size.thresh, pseudocount, nblocks, use.cores, cor.use)
    my.cs <- community.similarity.corr.par(phyloseq::otu_table(phyloseq_obj), 
                                           S=S.sparcc, distance=method, blocksize=nblocks, use.cores=use.cores)
    my.cs[my.cs < 0] <- 0
    
  } else if (method %in% c("PINA_unweighted", "PINA_weighted")) {
    
    S.phylo <- Phylogeny_wrapper(phyloseq_obj, use.cores, cor.use)
    my.cs <- community.similarity.corr.par(phyloseq::otu_table(phyloseq_obj), 
                                           S=S.phylo, distance=method, blocksize=nblocks, use.cores=use.cores)
    my.cs[my.cs < 0] <- 0
    
  }
  
  return(my.cs)
  
}