#updated color panel , corrected co-occurence, updated the size of the dots and stars, updated the margins 
#' Exact tests to detect mutually exclusive, co-occuring and altered genesets.
#'
#' @description Performs Pair-wise Fisher's Exact test to detect mutually exclusive or co-occuring events. Also identifies gene sets mutated significantly.
#' @details This function and plotting is inspired from genetic interaction analysis performed in the published study combining gene expression and mutation data in MDS. See reference for details.
#' @references Gerstung M, Pellagatti A, Malcovati L, et al. Combining gene mutation with gene expression data improves outcome prediction in myelodysplastic syndromes. Nature Communications. 2015;6:5901. doi:10.1038/ncomms6901.
#' @param maf an \code{\link{MAF}} object generated by \code{\link{read.maf}}
#' @param top check for interactions among top 'n' number of genes. Defaults to top 25. \code{genes}
#' @param genes List of genes among which interactions should be tested. If not provided, test will be performed between top 25 genes.
#' @param pvalue Default c(0.05, 0.01) p-value threshold. You can provide two values for upper and lower threshold.
#' @param returnAll If TRUE returns test statistics for all pair of tested genes. Default FALSE, returns for only genes below pvalue threshold.
#' @param findPathways Uses all mutually exclusive set of genes to further identify altered pathways. Default TRUE
#' @param kMax Default 3. maximum gene set size if findPathways is TRUE. This is time consuming for > 3.
#' @param fontSize cex for gene names. Default 0.8
#' @param verbose Default TRUE
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.maf)
#' somaticInteractions(maf = laml, top = 5)
#' @return list of data.tables
#' @export

somaticInteractionsUpdated = function(maf, top = 25, genes = NULL, pvalue = c(0.05, 0.01), returnAll = FALSE, findPathways = TRUE, kMax = 3, fontSize = 1.5, verbose = TRUE){
  
  if(is.null(genes)){
    genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
  }
  
  if(length(genes) < 2){
    stop("Minimum two genes required!")
  }
  
  om = createOncoMatrix(m = maf, g = genes)
  all.tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])
  
  mutMat = t(om$numericMatrix)
  missing.tsbs = all.tsbs[!all.tsbs %in% rownames(mutMat)]
  
  if(nrow(mutMat) < 2){
    stop("Minimum two genes required!")
  }
  mutMat[mutMat > 0 ] = 1
  
  if(length(missing.tsbs) > 0){
    missing.tsbs = as.data.frame(matrix(data = 0, nrow = length(missing.tsbs), ncol = ncol(mutMat)),
                                 row.names = missing.tsbs)
    colnames(missing.tsbs) = colnames(mutMat)
    mutMat = rbind(mutMat, missing.tsbs)
  }
  
  #pairwise fisher test source code borrowed from: https://www.nature.com/articles/ncomms6901
  interactions = sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") NA else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
  oddsRatio <- oddsGenes <- sapply(1:ncol(mutMat), function(i) sapply(1:ncol(mutMat), function(j) {f<- try(fisher.test(mutMat[,i], mutMat[,j]), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
  rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)
  
  if(returnAll){
    sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
  }else{
    sigPairs = which(x = 10^-abs(interactions) < max(pvalue), arr.ind = TRUE)
  }
  
  if(nrow(sigPairs) < 1){
    stop("No meaningful interactions found.")
  }
  
  sigPairsTbl = data.table::rbindlist(
    lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
      x = sigPairs[i,]
      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
      tbl = as.data.frame(table(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = "")))
      combn = data.frame(t(tbl$Freq))
      colnames(combn) = tbl$Var1
      pval = 10^-abs(interactions[x[1], x[2]])
      fest = oddsRatio[x[1], x[2]]
      d = data.table::data.table(gene1 = g1,
                                 gene2 = g2,
                                 pValue = pval, oddsRatio = fest)
      d = cbind(d, combn)
      d
    }), fill = TRUE)
  
  sigPairsTbl = sigPairsTbl[!gene1 == gene2] #Remove doagonal elements
  sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurance", no = "Mutually_Exclusive")
  sigPairsTbl$pair = apply(X = sigPairsTbl[,.(gene1, gene2)], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]
  
  #Source code borrowed from: https://www.nature.com/articles/ncomms6901
  if(nrow(interactions) >= 5){
    interactions[10^-abs(interactions) > max(pvalue)] = 0
    diag(interactions) <- 0
    m <- nrow(interactions)
    n <- ncol(interactions)
    
    interactions[interactions < -4] = -4
    interactions[interactions > 4] = 4
    # r = interactions
    # rd = hclust(dist(r))$order
    # cd = hclust(dist(t(r)))$order
    # interactions = interactions[rd, , drop = FALSE]
    # interactions = interactions[,rd, drop = FALSE]
    
    interactions[lower.tri(x = interactions)] = NA
    
    print(interactions)
    
    par(bty="n", mgp = c(2,.5,0), mar = c(2, 9, 4, 7)+.1, las=2, tcl=-.33)
    image(x=1:n, y=1:m, interactions, col=RColorBrewer::brewer.pal(9,"BrBG"),
          breaks = c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n",
          xlab="",ylab="", xlim=c(0, n+4), ylim=c(0, n+4))
    abline(h=0:n+.5, col="white", lwd=.5)
    abline(v=0:n+.5, col="white", lwd=.5)
    
   # mtext(side = 2, at = 1:m, text = colnames(interactions), las = 2, line = -0.8,cex = fontSize, font = 3)
  #  mtext(side = 3, at = 1:n, text = colnames(interactions), las = 2, line = -4.5, cex = fontSize, font = 3)
    
     mtext(side = 2, at = 1:m, text = colnames(interactions), las = 2,cex = fontSize, font = 3)
     mtext(side = 3, at = 1:n, text = colnames(interactions),line = -5, las = 2, cex = fontSize, font = 3)
    
    q <- p.adjust(10^-abs(interactions), method="BH")
    p <- p.adjust(10^-abs(interactions), method="holm")
  #  w = arrayInd(which(interactions < .05), rep(m,2))
  #  points(w, pch=".", col="white", cex=1.5)
    w = arrayInd(which(10^-abs(interactions) < min(pvalue)), rep(m,2))
    points(w, pch="*", col="black",cex=2)
    w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
    points(w, pch=".", col="black",cex=4)
    #image(y = 1:8 +6, x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"PiYG"), add=TRUE)
    image(y = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 8), x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col = RColorBrewer::brewer.pal(8,"BrBG"), add=TRUE)
    #axis(side = 4, at = seq(1,7) + 6.5,  tcl=-.15, label=seq(-3, 3), las=1, lwd=.5)
    atLims = seq(0.5*nrow(interactions), 0.9*nrow(interactions), length.out = 7)
    axis(side = 4, at = atLims,  tcl=-.15, labels =c(3:1, 0, 1:3), las=1, lwd=.5)
    mtext(side=4, at = median(atLims), "-log10 (p-value)", las=3, cex = 1.5, line = 3, font = 1)
    
    par(xpd=NA)
    text(x=n+1.8, y= max(atLims)+1.4, "Co-occurrence", pos=4, cex = 1.5, font = 1)
    text(x=n+1.8, y = min(atLims)-1.4, "Mutual Exclusion", pos=4, cex = 1.5, font = 1)
    
    points(x = n+1, y = 0.2*n, pch = "*", cex = 2)
    text(x = n+1, y = 0.2*n, paste0(" p < ", min(pvalue)), pos=4, cex = 1.5, font = 1)
    points(x = n+1, y = 0.1*n, pch = ".", cex = 6)
    text(x = n+1, y = 0.1*n, paste0("p < ", max(pvalue)), pos=4, cex = 1.5, font=1)
  }
  
  
  sig.genes.pvals = NULL
  
  if(findPathways){
    if(nrow(sigPairsTblSig[Event %in% 'Mutually_Exclusive']) > 2){
      if(verbose){
        cat("Checking for Gene sets\n")
      }
      sig.genes = unique(c(sigPairsTblSig[Event %in% 'Mutually_Exclusive', gene1], sigPairsTblSig[Event %in% 'Mutually_Exclusive', gene2]))
      sig.genes.pvals = c()
      
      for(k in 3:kMax){
        sig.genes.combn = combn(x = sig.genes, m = k)
        
        if(verbose){
          cat("------------------\n")
          cat(paste0("genes: ", length(sig.genes), "\n"))
          cat(paste0("geneset size: ", k, "\n"))
          cat(paste0(ncol(sig.genes.combn), " combinations\n"))
        }
        
        sps = lapply(seq_along(1:ncol(sig.genes.combn)), function(i){
          x = sig.genes.combn[,i]
          mm = mutMat[,x, drop = FALSE]
          grid.mat = t(expand.grid(rep(list(0:1), k)))
          #colllapse grid and get all the levels (all posiible combinations)
          lvls = names(table(apply(grid.mat, 2, paste, collapse = '')))
          mm.lvls = data.frame(table(apply(mm, 1, paste, collapse = '')))
          
          #check if for any missing combinations
          lvls.missing = lvls[!lvls %in% mm.lvls[,1]]
          
          if(length(lvls.missing) > 0){
            mm.lvls = rbind(mm.lvls, data.frame(Var1 = lvls.missing, Freq = 0)) #add missing combinations with zero values
          }
          
          #reorder
          mm.lvls = mm.lvls[order(mm.lvls$Var1),]
          if(verbose){
            #cat("Geneset: ", paste(x, collapse = ", "), "\n")
          }
          xp = cometExactTest::comet_exact_test(tbl = as.integer(as.character(mm.lvls$Freq)), mutmatplot = FALSE, pvalthresh = 0.1)
          data.table::data.table(gene_set = paste(x, collapse = ", "), pvalue = xp)
        })
        
        sig.genes.pvals = rbind(sig.genes.pvals, data.table::rbindlist(sps))
      }
      
      sig.genes.pvals = sig.genes.pvals[pvalue > 0][order(pvalue, decreasing = FALSE)]
      if(nrow(sig.genes.pvals[pvalue < 0.05]) > 0){
        if(verbose){
          cat(paste0("Signifcantly altered gene-sets: ", nrow(sig.genes.pvals[pvalue < 0.05])), "\n")
          cat("------------------\n")
        }
      }
    }
  }
  
  sigPairsTblSig[,pair := NULL]
  
  return(list(pairs = sigPairsTblSig, gene_sets = sig.genes.pvals))
}

createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
   ###### mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

#---- This is small function to sort genes according to total samples in which it is mutated.
sortByMutation = function(numMat, maf){
  
  geneOrder = getGeneSummary(x = maf)[order(MutatedSamples, decreasing = TRUE), Hugo_Symbol]
  numMat = numMat[as.character(geneOrder[geneOrder %in% rownames(numMat)]),, drop = FALSE]
  
  numMat[numMat != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tnumMat = t(numMat) #transposematrix
  numMat = t(tnumMat[do.call(order, c(as.list(as.data.frame(tnumMat)), decreasing = TRUE)), ]) #sort
  
  return(numMat)
}

#Thanks to Ryan Morin for the suggestion (https://github.com/rdmorin)
#original code has been changed with vectorized code, in-addition performs class-wise sorting.
sortByAnnotation <-function(numMat,maf, anno, annoOrder = NULL, group = TRUE, isNumeric = FALSE){
  
  if(is.numeric(anno[,1])){
    #anno.spl = split(anno, as.numeric(as.character(anno[,1]))) #sorting only first annotation
    mat_samps = colnames(numMat)[colnames(numMat) %in% rownames(anno)]
    anno = anno[mat_samps,, drop = FALSE]
    anno = anno[order(anno[,1], na.last = TRUE),, drop = FALSE]
    numMat.sorted = numMat[,rownames(anno)]
  }else{
    anno[,1] = ifelse(test = is.na(as.character(anno[,1])), yes = "NA", no = as.character(anno[,1])) #NAs are notorious; converting them to characters
    anno.spl = split(anno, as.factor(as.character(anno[,1]))) #sorting only first annotation
    
    anno.spl.sort = lapply(X = anno.spl, function(x){
      numMat[,colnames(numMat)[colnames(numMat) %in% rownames(x)], drop = FALSE]
    })
    
    if(group){
      #sort list according to number of elemnts in each classification
      anno.spl.sort = anno.spl.sort[names(sort(unlist(lapply(anno.spl.sort, ncol)), decreasing = TRUE, na.last = TRUE))]
    }
    
    if(!is.null(annoOrder)){
      annoSplOrder = names(anno.spl.sort)
      
      if(length(annoOrder[annoOrder %in% annoSplOrder]) == 0){
        message("Values in provided annotation order ", paste(annoOrder, collapse = ", ")," does not match values in clinical features. Here are the available features..")
        print(annoSplOrder)
        stop()
      }
      annoOrder = annoOrder[annoOrder %in% annoSplOrder]
      
      anno.spl.sort = anno.spl.sort[annoOrder]
      
      if(length(annoSplOrder[!annoSplOrder %in% annoOrder]) > 0){
        warning("Following levels are missing from the provided annotation order: ", paste(annoSplOrder[!annoSplOrder %in% annoOrder], collapse = ", "), immediate. = TRUE)
      }
    }
    
    numMat.sorted = c()
    for(i in 1:length(anno.spl.sort)){
      numMat.sorted  = cbind(numMat.sorted, anno.spl.sort[[i]])
    }
  }
  
  return(numMat.sorted)
}

#Sort columns while keeping gene order
sortByGeneOrder = function(m, g){
  m = m[g,, drop = FALSE]
  tsbs= colnames(m)
  tsb.anno = data.table::data.table()
  
  for(i in 1:nrow(m)){
    x = m[i,]
    tsb.anno = rbind(tsb.anno, data.table::data.table(tsbs = names(x[x!=0]), gene = g[i]))
  }
  tsb.anno = tsb.anno[!duplicated(tsbs)]
  
  if(length(tsbs[!tsbs %in% tsb.anno[,tsbs]]) > 0){
    tsb.anno = rbind(tsb.anno, data.table::data.table(tsbs = tsbs[!tsbs %in% tsb.anno[,tsbs]], gene = NA))
  }
  data.table::setDF(tsb.anno)
  #m.sorted = sortByAnnotation(numMat = m, anno = tsb.anno)
  
  anno.spl = split(tsb.anno, as.factor(as.character(tsb.anno$gene)))
  anno.spl.sort = sapply(X = anno.spl, function(x){
    m[,colnames(m)[colnames(m) %in% x$tsbs], drop = FALSE]
  })
  
  numMat.sorted = c()
  for(i in 1:length(anno.spl.sort)){
    numMat.sorted  = cbind(numMat.sorted, anno.spl.sort[[i]])
  }
  
  return(numMat.sorted)
}

#plot_dat = plotting data
#lab_dat = data to be labelled
#x_var = x
#y_var = y
#bubble_var = z (variable name for bubble size)
#bubble_size (exclusive with bubble_var)
#text_size = font size for labels
#col_var = a vector color
bubble_plot = function(plot_dat, lab_dat = NULL, x_var = NULL, y_var = NULL,
                       bubble_var = NULL, bubble_size = 1, text_var = NULL,
                       text_size = 1, col_var = NULL, return_dat = FALSE){
  
  x_col_idx = which(colnames(plot_dat) == x_var)
  y_col_idx = which(colnames(plot_dat) == y_var)
  colnames(plot_dat)[c(x_col_idx, y_col_idx)] = c("x", "y")
  
  if(!is.null(col_var)){
    col_idx = which(colnames(plot_dat) == col_var)
    colnames(plot_dat)[col_idx] = c("color_var")
  }else{
    plot_dat$color_var = grDevices::adjustcolor("black", alpha.f = "0.75")
  }
  
  if(!is.null(lab_dat)){
    x_col_idx = which(colnames(lab_dat) == x_var)
    y_col_idx = which(colnames(lab_dat) == y_var)
    colnames(lab_dat)[c(x_col_idx, y_col_idx)] = c("x", "y")
    if(is.null(text_var)){
      stop("Missing text variable")
    }else{
      colnames(lab_dat)[which(colnames(lab_dat) == text_var)] = "z_text"
    }
    if(!is.null(col_var)){
      col_idx = which(colnames(lab_dat) == col_var)
      colnames(lab_dat)[col_idx] = c("color_var")
    }else{
      lab_dat$color_var = "black"
    }
  }
  
  if(!is.null(bubble_var)){
    
    if(bubble_var %in% c(x_var, y_var)){
      if(which(bubble_var == c(x_var, y_var)) == 2){
        plot_dat$size_z = sqrt(as.numeric(plot_dat$y)/pi)
      }else if(which(bubble_var == c(x_var, y_var)) == 1){
        plot_dat$size_z = sqrt(as.numeric(plot_dat$x)/pi)
      }
    }else{
      if(length(which(colnames(plot_dat) == bubble_var)) > 0){
        colnames(plot_dat)[which(colnames(plot_dat) == bubble_var)] = "z"
        plot_dat$size_z = sqrt(as.numeric(plot_dat$z)/pi)
      }else{
        plot_dat$size_z = bubble_size
      }
    }
    
    if(!is.null(lab_dat)){
      if(bubble_var %in% c(x_var, y_var)){
        if(which(bubble_var == c(x_var, y_var)) == 2){
          lab_dat$size_z = sqrt(as.numeric(lab_dat$y)/pi) * bubble_size
        }else if(which(bubble_var == c(x_var, y_var)) == 1){
          lab_dat$size_z = sqrt(as.numeric(lab_dat$x)/pi) * bubble_size
        }
      }else{
        colnames(lab_dat)[which(colnames(lab_dat) == bubble_var)] = "z"
        lab_dat$size_z = sqrt(as.numeric(lab_dat$z)/pi) * bubble_size
      }
    }
  }else{
    plot_dat$size_z = bubble_size
    if(!is.null(lab_dat)){
      lab_dat$size_z = bubble_size
    }
  }
  
  x_lims = as.integer(seq(min(as.numeric(plot_dat$x), na.rm = TRUE),
                          max(as.numeric(plot_dat$x), na.rm = TRUE), length.out = 4))
  x_lims[4] = as.integer(ceiling(max(as.numeric(plot_dat$x), na.rm = TRUE)))
  x_lims[1] = as.integer(floor(min(as.numeric(plot_dat$x), na.rm = TRUE)))
  x_ticks = pretty(x = x_lims, na.rm = TRUE)
  x_ticks[c(1, length(x_ticks))] = x_lims[c(1, 4)]
  
  y_lims = as.integer(seq(min(as.numeric(plot_dat$y), na.rm = TRUE),
                          max(as.numeric(plot_dat$y), na.rm = TRUE), length.out = 4))
  y_lims[4] = as.integer(ceiling(max(as.numeric(plot_dat$y), na.rm = TRUE)))
  y_lims[1] = as.integer(floor(min(as.numeric(plot_dat$y), na.rm = TRUE)))
  y_ticks = pretty(x = y_lims, na.rm = TRUE)
  y_ticks[c(1, length(y_ticks))] = y_lims[c(1, 4)]
  
  if(return_dat){
    return(list(plots_data = plot_dat, label_cords = lab_dat,
                x_lims = x_lims, y_lims = y_lims))
  }
  
  # plot(x = plot_dat$x, y = plot_dat$y, cex = plot_dat$size_z,
  #      pch = 16, col = plot_dat$color_var, axes = FALSE, xlim = x_lims[c(1, 4)],
  #      ylim = y_lims[c(1, 4)], xlab = NA, ylab = NA)
  suppressWarnings(symbols(x = plot_dat$x, y = plot_dat$y, circles = plot_dat$size_z, inches = 0.1, bg = plot_dat$color_var, xlim = x_lims[c(1, 4)],
                           ylim = y_lims[c(1, 4)], xlab = NA, ylab = NA, fg = "white", axes = FALSE))
  axis(side = 1, at = x_ticks)
  axis(side = 2, at = y_ticks, las = 2)
  abline(h = y_ticks, v = x_ticks, lty = 2,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.5), lwd = 0.75)
  
  if(!is.null(lab_dat)){
    # points(x = lab_dat$x, y = lab_dat$y, cex = lab_dat$size_z,
    #        pch = 16, col = lab_dat$color_var)
    if(nrow(lab_dat) > 0 & nrow(lab_dat) < 2){
      text(x = lab_dat$x, y = lab_dat$y, labels = lab_dat$z_text, adj = 1, offset = 0.2, cex = text_size, col = lab_dat$color_var)
    }else if(nrow(lab_dat) >= 2){
      symbols(x = lab_dat$x, y = lab_dat$y, circles = lab_dat$size_z,
              bg = lab_dat$color_var, add = TRUE, fg = "white", inches = 0.1)
      
      wordcloud::textplot(x = lab_dat$x, y = lab_dat$y, words = lab_dat$z_text,
                          cex = text_size, new = FALSE, show.lines = TRUE,
                          xlim = x_lims[c(1, 4)], ylim = y_lims[c(1, 4)], font = 3, col = lab_dat$color_var)
    }
  }
}


#Get plot layout for oncoplot
plot_layout = function(clinicalFeatures = NULL, drawRowBar = TRUE,
                       drawColBar = TRUE, draw_titv = FALSE, exprsTbl = NULL){
  if(is.null(clinicalFeatures)){
    if(draw_titv){
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4, 4), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3, 4), nrow = 4, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 4, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 4, 4), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4, 4), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4, 4), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, 4, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(1,4, 1), heights = c(4, 12, 4, 4))
        }
      }
    }else{
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2), nrow = 2, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 4), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,3), nrow = 2, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,4,4), nrow = 2, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 4), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(1, 4, 1), heights = c(4, 12, 4))
        }
      }
    }
  }else{
    if(draw_titv){
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4), nrow = 4, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4, 4), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5), nrow = 5, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 1, 4, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,9), nrow = 5, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 1, 4, 4), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4, 4), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4, 4), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,9), nrow = 5, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, 1, 4, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,11,12,13,13,13), nrow = 5, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(1, 4, 1), heights = c(4, 12, 1, 4, 4))
        }
      }
    }else{
      if(!drawRowBar & !drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4), widths = c(1, 4))
        }
      }else if(!drawRowBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4), nrow = 4, ncol = 1, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 1, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(4, 12, 1, 4), widths = c(1, 4))
        }
      }else if(!drawColBar){
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4), widths = c(4, 1))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7,7), nrow = 3, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, heights = c(12, 1, 4), widths = c(1, 4, 1))
        }
      }else{
        if(is.null(exprsTbl)){
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,7), nrow = 4, ncol = 2, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(4, 1), heights = c(4, 12, 1, 4))
        }else{
          mat_lo = matrix(data = c(1,2,3,4,5,6,7,8,9,10,10,10), nrow = 4, ncol = 3, byrow = TRUE)
          lo = layout(mat = mat_lo, widths = c(1, 4, 1), heights = c(4, 12, 1, 4))
        }
      }
    }
  }
  
  lo
}

get_vcColors = function(alpha = 1){
  col = c(RColorBrewer::brewer.pal(12,name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','IGR','Missense_Mutation','Silent','Nonsense_Mutation',
                         'RNA','Splice_Site','Intron','Frame_Shift_Ins','Nonstop_Mutation','In_Frame_Del','ITD','In_Frame_Ins',
                         'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del')
  col
}

get_titvCol = function(alpha = 1){
  col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  col
}

parse_annotation_dat = function(annotationDat = NULL, clinicalFeatures = NULL){
  
  if(class(annotationDat)[1] == 'MAF'){
    annotationDat = data.table::copy(getClinicalData(x = annotationDat))
    data.table::setDF(annotationDat)
  }else if(class(annotationDat)[1] %in% c('data.frame', 'data.table')){
    data.table::setDF(annotationDat)
  }else{
    return(NULL)
  }
  
  if(length(clinicalFeatures[!clinicalFeatures %in% colnames(annotationDat)]) > 0){
    message('Following columns are missing from annotation slot of MAF. Ignoring them..')
    print(clinicalFeatures[!clinicalFeatures %in% colnames(annotationDat)])
    clinicalFeatures = clinicalFeatures[clinicalFeatures %in% colnames(annotationDat)]
    if(length(clinicalFeatures) == 0){
      message('Make sure at-least one of the values from provided clinicalFeatures are present in annotation slot of MAF. Here are available annotaions..')
      print(colnames(annotationDat))
      stop('Zero annotaions to add! You can also provide custom annotations via annotationDat argument.')
    }
  }
  annotation = data.frame(row.names = annotationDat$Tumor_Sample_Barcode ,annotationDat[,clinicalFeatures, drop = FALSE], stringsAsFactors = FALSE)
  
  return(annotation)
}

get_anno_cols = function(ann, numericAnnoCol = NULL){
  ann_cols = list()
  if(is.null(numericAnnoCol)){
    numericAnnoCol =  RColorBrewer::brewer.pal(n = 9, name = "YlOrBr")
  }else{
    numericAnnoCol =  RColorBrewer::brewer.pal(n = 9, name = numericAnnoCol)
  }
  for(i in 1:ncol(ann)){
    if(is.numeric(ann[,i])){
      x = ann[,i]
      ann_lvls_cols = colorRampPalette(numericAnnoCol)(length(x))
      names(ann_lvls_cols) = x[order(x, na.last = TRUE)]
      ann_cols[[i]] = ann_lvls_cols
    }else{
      ann_lvls = unique(as.character(ann[,i]))
      if(length(ann_lvls) <= 9){
        ann_lvls_cols = RColorBrewer::brewer.pal(n = 9, name = 'Set1')[1:length(ann_lvls)]
      }else{
        ann_lvls_cols = colors()[sample(x = 1:100, size = length(ann_lvls), replace = FALSE)]
      }
      ann_cols[[i]] = ann_lvls_cols
      names(ann_cols[[i]]) = ann_lvls
    }
  }
  
  names(ann_cols) = colnames(ann)
  
  return(ann_cols)
}
