#SparCC correlation analysis
sparcc <- function(otu.table, size.thresh=1, pseudocount=10^-6, nblocks=400, use.cores=4) {

  otu.table <- as(otu.table, "matrix")
  remove.otus <- rowSums(otu.table) < size.thresh
  keep.otus <- ! remove.otus
  o.t <- otu.table[! remove.otus, ] + pseudocount
  
  otus <- rownames(o.t)
  n.otu <- length(otus)
  
  #Preallocate blocks for parallel processing & Aitchinson's T matrix
  #=> based on https://gist.github.com/bobthecat/5024079
  size.split <- floor(n.otu / nblocks)
  if (size.split < 1) {size.split <- 1}
  my.split <- list() 
  length(my.split) <- nblocks
  my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split))
  my.split[[nblocks]] <- (size.split*(nblocks-1)):n.otu
  dat.split <- parallel::mclapply(my.split, function(g) {o.t[g,]}, mc.cores=use.cores);
  #Get combinations of splits
  my.combs <- expand.grid(1:length(my.split), 1:length(my.split)) %>%
    apply(., 1, sort) %>%
    t() %>%
    unique(my.combs)
  #Preallocate Aitchinson's T matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
  mat.T <- bigmemory::big.matrix(nrow=n.otu, ncol=n.otu, dimnames=list(otus, otus), shared=T)
  mat.T.desc <- bigmemory::describe(mat.T)
  
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
    
    curr.mat.T <- bigmemory::attach.big.matrix(mat.T.desc);
    curr.mat.T[g.2, g.1] <- curr.T;
    curr.mat.T[g.1, g.2] <- t(curr.T);
    #Return
    TRUE
  }
  
  var.t <- colsum(mat.T);
  
  mat.a <- matrix(data=1, nrow=n.otu, ncol=n.otu)
  diag(mat.a) <- n.otu-1
  omega <- sqrt(solve(a=mat.a, b=var.t))
  
  mat.rho <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {
    (omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)
  }
  rownames(mat.rho) <- rownames(o.t);
  
  #Return
  mat.rho
}

#Parallel computation of correlations
cor.par <- function(mat, method="pearson", nblocks=400, use="na.or.complete", use.cores=4) {

  n.el <- ncol(mat)
  if (n.el < 2*nblocks) {return(cor(mat, method=method, use=use))}
  
  size.split <- floor(n.el / nblocks);
  if (size.split < 1) {size.split <- 1}
  
  my.split <- list()
  length(my.split) <- nblocks
  my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split))
  my.split[[nblocks]] <- (size.split*(nblocks-1)):n.el
  
  dat.split <- parallel::mclapply(my.split, function(g) {as.matrix(mat[,g])}, mc.cores=use.cores)
  #Get combinations of splits
  my.combs <- expand.grid(1:length(my.split), 1:length(my.split)) %>%
    apply(., 1, sort) %>%
    t() %>%
    unique()
  
  #Preallocate correlation matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
  mat.rho <- bigmemory::big.matrix(nrow=n.el, ncol=n.el, dimnames=list(colnames(mat), colnames(mat)), shared=T)
  mat.rho.desc <- bigmemory::describe(mat.rho)
 
  #Compute correlation matrix
  results <- foreach(i = 1:nrow(my.combs)) %dopar% {
    #Get current combination
    curr.comb <- my.combs[i, ];
    #Get current data
    g.1 <- my.split[[curr.comb[1]]];
    g.2 <- my.split[[curr.comb[2]]];
    curr.rho <- cor(x=dat.split[[curr.comb[1]]], y=dat.split[[curr.comb[2]]], method=method);
    
    curr.mat.rho <- bigmemory::attach.big.matrix(mat.rho.desc);
    curr.mat.rho[g.1, g.2] <- curr.rho;
    curr.mat.rho[g.2, g.1] <- t(curr.rho);
    #Return
    TRUE
  }
  
  as.matrix(mat.rho);
}

#Traditional indices of community similarity
community.similarity.par <- function(o.t, distance="jaccard", use.cores=4) {

  o.t <- as(o.t, "matrix")
  samples <- colnames(o.t);
  
  #Turn OTU table into list format
  ot.occ.list <- apply((o.t > 0), 2, which);
  ot.count.list <- plyr::alply(o.t, .margins=2, .fun=function(o.vec) {o.vec[o.vec > 0]}, .parallel=T)
  names(ot.occ.list) <- names(ot.count.list) <- samples;
  
  if (distance == "jaccard") {
    cs.list <- parallel::mclapply(ot.occ.list, function(a.list) {
      lapply(ot.occ.list, function(b.list) {
        1 - (length(intersect(a.list, b.list)) / length(union(a.list, b.list)))
      }) %>%
      unlist()
    }, mc.cores=use.cores)
  }
  
  #Classic, weighted Jaccard
  if (distance == "jaccard.abd") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.list <- names(a.count)
      lapply(ot.count.list, function(b.count) {
        ab.shared <- intersect(a.list, names(b.count))
        1 - (sum(a.count[ab.shared]) + sum(b.count[ab.shared])) / (a.N + sum(b.count))
      }) %>%
      unlist()
    }, mc.cores=use.cores)
  }
  
  #Classic, weighted Jaccard, based on fractions (more balancing for uneven sample sizes)
  if (distance == "jaccard.abd.frac") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.list <- names(a.count)
      a.rel <- a.count / a.N
      lapply(ot.count.list, function(b.count) {
        b.N <- sum(b.count)
        ab.shared <- intersect(a.list, names(b.count))
        b.rel <- b.count / b.N
        1 - (0.5 * (sum(a.rel[ab.shared]) + sum(b.rel[ab.shared])))
      }) %>%
      unlist()
    }, mc.cores=use.cores);
  }
  
  #Classic, weighted Jaccard, based on fractions (more balancing for uneven sample sizes), Chao's version
  if (distance == "jaccard.abd.frac_chao") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.list <- names(a.count)
      a.rel <- a.count / a.N
      lapply(ot.count.list, function(b.count) {
        b.N <- sum(b.count)
        ab.shared <- intersect(a.list, names(b.count))
        b.rel <- b.count / b.N
        U <- sum(a.rel[ab.shared])
        V <- sum(b.rel[ab.shared])
        if (U == 0 & V == 0) {return(1)} else {return(1 - ((U*V) / (U + V - (U*V))))}
      }) %>%
      unlist()
    }, mc.cores=use.cores);
  }
  
  #Classic, weighted Jaccard, alternative formulation (Bray-Curtis-like)
  if (distance == "jaccard.abd.alt") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.list <- names(a.count)
      lapply(ot.count.list, function(b.count) {
        b.list <- names(b.count);
        ab.shared <- intersect(a.list, b.list)
        ab.union <- union(a.list, b.list)
        a.freq <- b.freq <- numeric(length=length(ab.union))
        names(a.freq) <- names(b.freq) <- ab.union
        a.freq[a.list] <- a.count
        b.freq[b.list] <- b.count
        1 - (sum(pmin(a.freq, b.freq)) / sum(pmax(a.freq, b.freq)))
      }) %>%
      unlist()
    }, mc.cores=use.cores);
  }
  
  #Classic, weighted Jaccard, Chao's version
  if (distance == "jaccard.abd.chao") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.list <- names(a.count)
      lapply(ot.count.list, function(b.count) {
        b.list <- names(b.count)
        b.N <- sum(b.count)
        #Get shared taxa & union of taxa
        ab.shared <- intersect(a.list, b.list)
        ab.union <- sort(union(a.list, b.list))
        a.freq <- b.freq <- numeric(length=length(ab.union))
        names(a.freq) <- names(b.freq) <- ab.union
        a.freq[a.list] <- a.count
        b.freq[b.list] <- b.count
        #If no taxa are shared, return trivial (distance = 1), otherwise compute Chao distance
        if (length(ab.shared) == 0) {d.chao <- 1} else {
          #Which taxa observed in sample a are singletons in sample b and vice versa?
          a_shared.b_singl <- which(a.freq >= 1 & b.freq == 1);
          b_shared.a_singl <- which(a.freq == 1 & b.freq >= 1);
          #How many shared taxa are singletons / doubletons in samples a & b?
          f.a.singl <- length(b_shared.a_singl);
          f.a.doubl <- length(which(a.freq == 2 & b.freq >= 1));
          f.b.singl <- length(a_shared.b_singl);
          f.b.doubl <- length(which(a.freq >= 1 & b.freq == 2));
          if (f.a.doubl == 0) {f.a.doubl = 1}
          if (f.b.doubl == 0) {f.b.doubl = 1}
          #Get taxa abundance estimates for samples a & b
          a.adjust <- ((b.N-1)/b.N)*f.b.singl/(2*f.b.doubl)*sum(a.freq[a_shared.b_singl] / a.N);
          a.estimate <- sum(a.freq[ab.shared] / a.N) + a.adjust;
          b.adjust <- ((a.N-1)/a.N)*f.a.singl/(2*f.a.doubl)*sum(b.freq[b_shared.a_singl] / b.N);
          b.estimate <- sum(b.freq[ab.shared] / b.N) + b.adjust;
          if (a.estimate > 1) {a.estimate = 1}
          if (b.estimate > 1) {b.estimate == 1}
          d.chao <- 1 - ((a.estimate * b.estimate) / (a.estimate + b.estimate - (a.estimate * b.estimate)));
          if (d.chao < 0) {d.chao = 0}
        }
        
        d.chao
      }) %>%
      unlist()
    }, mc.cores=use.cores)
  }
  
  #Classic Bray-Curtis dissimilarity
  if (distance == "bray_curtis") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.list <- names(a.count)
      lapply(ot.count.list, function(b.count) {
        b.N <- sum(b.count)
        b.list <- names(b.count)
        ab.shared <- intersect(a.list, b.list)
        ab.union <- union(a.list, b.list)
        a.freq <- b.freq <- numeric(length=length(ab.union))
        names(a.freq) <- names(b.freq) <- ab.union
        a.freq[a.list] <- a.count
        b.freq[b.list] <- b.count
        1 - (2 * (sum(pmin(a.freq, b.freq)) / (a.N + b.N)));
      }) %>%
      unlist()
    }, mc.cores=use.cores)
  }
  
  #Classic Morisita-Horn dissimilarity
  if (distance == "morisita_horn") {
    cs.list <- parallel::mclapply(ot.count.list, function(a.count) {
      a.N <- sum(a.count)
      a.n <- length(a.count)
      a.list <- names(a.count)
      a.simpson <- sum(a.count ^ 2) / (a.N ^ 2)
      lapply(ot.count.list, function(b.count) {
        #Get current summary statistics for B
        b.N <- sum(b.count)
        b.n <- length(b.count)
        b.list <- names(b.count)
        b.simpson <- sum(b.count ^ 2) / (b.N ^ 2)
        #Get current freq vectors from intersect and union of a & b
        ab.shared <- intersect(a.list, b.list)
        ab.union <- union(a.list, b.list)
        a.freq <- b.freq <- numeric(length=length(ab.union))
        names(a.freq) <- names(b.freq) <- ab.union
        a.freq[a.list] <- a.count
        b.freq[b.list] <- b.count
        #Get summed abundance product
        ab.prod <- sum(a.freq * b.freq)
        
        1 - ((2 * ab.prod) / ((a.simpson + b.simpson) * a.N * b.N));
      }) %>%
      unlist()
    }, mc.cores=use.cores);
  }

  do.call("rbind", cs.list);
}

#Faster UniFrac
fasterUniFrac <- function(o.t, tree, distance="UniFrac_unweighted", blocksize=1000, use.cores) {
  
  sp.o.t <- Matrix::Matrix(o.t, sparse=T)
  c.samples <- colnames(sp.o.t)
  c.sample.sizes <- colSums(sp.o.t)
  n.sample <- length(c.samples)
  
  if (!all(rownames(sp.o.t) == phyloseq::taxa_names(tree)) ){
    sp.o.t <- sp.o.t[phyloseq::taxa_names(tree), ]
  }
  #Get N x 2 matrix of pairwise combinations of samples
  spn <- utils::combn(c.samples, 2, simplify=FALSE)
  
  # Create a list of descendants, starting from the first internal node (root)
  # Add the terminal edge descendants (tips). By definition, can only have one descendant
  descList <- c(as.list(1:length(tree$tip.label)), ape::prop.part(tree, check.labels = FALSE))
  
  # Convert `descList` to `edge_array` that matches the order of things in `tree$edge`
  # For each entry in the tree$edge table, sum the descendants for each sample
  # `tree$edge[i, 2]` is the node ID.
  tmp.edge.array <- parallel::mclapply(1:nrow(tree$edge), function(i) {
    colSums(sp.o.t[descList[[tree$edge[i, 2]]], ,drop=F], na.rm=T)
  }, mc.cores=use.cores)
  
  edge.array <- do.call("rbind", tmp.edge.array)
  edge.occ <- (edge.array > 0) - 0
  edge.occ.list <- apply((edge.array > 0), 2, which);
  
  z = ape::reorder.phylo(tree, order="postorder")
  #Get "tip ages", or horizontal position of edges in re-ordered phylo object
  #Keep only tips
  tip.ages = ape::node.depth.edgelength(tree)[1:length(tree$tip.label)]
  
  names(tip.ages) <- z$tip.label
  tip.ages <- tip.ages[rownames(sp.o.t)]
  
  #Unweighted UniFrac
  size.chunk <- blocksize;
  if (distance == "UniFrac_unweighted") {
    cs.list <- foreach(i=split(1:length(spn), ceiling(seq_along(1:length(spn))/size.chunk))) %dopar% {
      lapply(spn[i], function(sp) {
        #Get occurrence vectors for current samples
        a <- edge.occ.list[[sp[1]]]
        b <- edge.occ.list[[sp[2]]]
        curr.shared_edges <- intersect(a, b)
        curr.union_edges <- union(a, b)
        #Return: 1 - (shared branch length) / (total branch length)
        1 - sum(tree$edge.length[curr.shared_edges], na.rm=TRUE) / sum(tree$edge.length[curr.union_edges], na.rm=TRUE)
      }) %>%
      unlist()
    }
  }
  
  #Weighted UniFrac
  if (distance == "UniFrac_weighted") {
    cs.list <- foreach(i=split(1:length(spn), ceiling(seq_along(1:length(spn))/size.chunk))) %dopar% {
      lapply(spn[i], function(sp) {
        #Get current samples and total sample sizes
        a <- sp[1]
        b <- sp[2]
        a.T <- c.sample.sizes[a]
        b.T <- c.sample.sizes[b]
        #Calculate branchweights
        curr.bw <- abs(edge.array[, a]/a.T - edge.array[, b]/b.T)
        #Return: (unshared branch length) / (total branch length)
        sum(tree$edge.length * curr.bw, na.rm = TRUE) / sum(tip.ages * (sp.o.t[,a]/a.T + sp.o.t[,b]/b.T), na.rm = TRUE)
      }) %>%
      unlist()
    }
  }
  
  #Pass output
  mat.UF <- matrix(NA_real_, n.sample, n.sample)
  diag(mat.UF) <- 0
  rownames(mat.UF) <- colnames(mat.UF) <- c.samples
  
  # Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
  mat.idx <- do.call(rbind, spn)[, 2:1]
  #Coerce results into matrices
  mat.UF[mat.idx] <- unlist(cs.list)
  mat.UF[mat.idx[,2:1]] <- unlist(cs.list)
  mat.UF
}

#Interaction-adjusted indices of community similarity
community.similarity.corr.par <- function(o.t, S, distance="jaccard.corr.uw", blocksize=1000, use.cores=4) {

  #Get N x 2 matrix of pairwise combinations of samples
  o.t <- as(o.t, "matrix")
  
  samples <- colnames(o.t)
  ot.occ.list <- apply((o.t > 0), 2, which)
  
  ot.count.list <- plyr::alply(o.t, .margins=2, .fun=function(o.vec) {o.vec[o.vec > 0]}, .parallel=T)
  names(ot.occ.list) <- names(ot.count.list) <- samples;
  
  #Unweighted, corrected Jaccard, normalized by sample self-comparisons
  #=> this is for unweighted TINA/PINA
  if (distance %in% c("TINA_unweighted", "PINA_unweighted")) {
    #Pre-calculate sample-wise interaction term sums (normalized by number of OTUs per sample)
    smpl.colsums <- parallel::mclapply(X = ot.occ.list, FUN = function(a.list) {
      if (length(a.list) == 1) {
        tmp <- matrix(S[a.list, ], nrow = 1)
      } else {
        tmp <- S[a.list, ]
      }
      colSums(tmp)/ length(a.list)
      }, mc.cores=use.cores)
    
    #Pre-calculate sample self-comparisons
    smpl.self <- unlist(lapply(samples, function(a) {
      a.list <- ot.occ.list[[a]]
      sum(smpl.colsums[[a]][a.list]) / length(a.list)
    }))
    
    names(smpl.self) <- samples
    
    #Calculate similarity index in parallel
    cs.list <- parallel::mclapply(samples, function(a) {
      a.sums <- smpl.colsums[[a]]
      a.self <- smpl.self[a]
      lapply(samples, function(b) {
        b.list <- ot.occ.list[[b]]
        b.self <- smpl.self[b]
        1 - (sum(a.sums[b.list]) / (length(b.list) * sqrt(a.self*b.self)))
      }) %>%
      unlist()
    }, mc.cores=use.cores)
  }
  
  #Weighted, corrected Jaccard, normalized by sample self-comparisons
  #=> this is for weighted TINA/PINA
  if (distance %in% c("TINA_weighted", "PINA_weighted")) {
    #Pre-calculate relative abundances per sample
    ot.rel_count <- lapply(ot.count.list, function(a.count) {a.count / sum(a.count)})
    #Pre-calculate sample self-comparisons
    smpl.self <- parallel::mclapply(ot.rel_count, function(a.rel) {
        a.list <- names(a.rel)
        a.S <- S[a.list, a.list]
        sum(a.S * outer(a.rel, a.rel));
      }, mc.cores=use.cores) %>%
      unlist()
    #Iterate through samples in parallel
    cs.list <- parallel::mclapply(samples, function(a) {
      #Get current OTU names, relative counts and sample size
      a.rel <- ot.rel_count[[a]]
      a.list <- names(a.rel)
      #Get current interaction sub-matrix
      if (length(a.list) == 1) {
        a.S <- a.rel * matrix(S[a.list,], nrow = 1) %>%
          magrittr::set_rownames(a.list) %>%
          magrittr::set_colnames(rownames(o.t))
      } else {
        a.S <- a.rel * S[a.list,]
      }
      #Iterate through all potential partner samples
      lapply(samples, function(b) {
        b.rel <- ot.rel_count[[b]]
        b.list <- names(b.rel)
        curr.S <- a.S[,b.list]
        #Calculate sample-pair weighted interaction term
        1 - (sum(b.rel * t(curr.S)) / sqrt(smpl.self[a] * smpl.self[b]));
      }) %>%
      unlist()
    }, mc.cores=use.cores);
  }
  
  return.mat <- do.call("rbind", cs.list)
  rownames(return.mat) <- colnames(return.mat) <- samples
  return.mat
}

fast.cophenetic.phylo <- function(x, nblocks=1000, use.cores=4) {
  
  #Get tree in "postorder" order
  x = ape::reorder.phylo(x, order="postorder");
  tip.order <- order(x$tip.label)
  
  #Get "tip ages", or horizontal position of edges in re-ordered phylo object
  node.ages = ape::node.depth.edgelength(x)
  my.combs <- gRbase::combn_prim(tip.order, 2) %>%
    t()
  
  #Split combinations into processable chunks
  size.split <- floor(nrow(my.combs) / nblocks)
  if (size.split < 1) {size.split <- 1}
  
  my.split <- list()
  length(my.split) <- nblocks
  my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
  my.split[[nblocks]] <- ((size.split*(nblocks-1))+1):nrow(my.combs);
  
  #Get node paths to root for all tips
  my.node_paths <- parallel::mclapply(seq(1, ape::Ntip(x)), function (tip) {ape::nodepath(x, from=tip, to = 1)}, mc.cores=use.cores)
  
  #Loop through all pairs of tips (in parallel) and calculate pairwise distance
  #=> d.cophenetic = (d[i, root] + d[j, root]) - 2 * d[mrca.ij, root]
  split.sums <- parallel::mclapply(my.split, function(g) {node.ages[my.combs[g,1]]+node.ages[my.combs[g,2]]}, mc.cores=use.cores);
  split.mrca <- parallel::mcmapply(function(a,b) {max(intersect(my.node_paths[[a]], my.node_paths[[b]]))}, my.combs[,1], my.combs[,2], mc.cores=use.cores);
  dist.vec <- unlist(split.sums) - 2*node.ages[split.mrca];
  
  #Coerce into matrix
  mat.out <- matrix(nrow=ape::Ntip(x), ncol=ape::Ntip(x))
  diag(mat.out) <- 0
  rownames(mat.out) <- colnames(mat.out) <- x$tip.label[tip.order]
  mat.out[upper.tri(mat.out)] <- dist.vec
  mat.out[lower.tri(mat.out)] <- t(mat.out)[lower.tri(mat.out)]
  
  mat.out
}

