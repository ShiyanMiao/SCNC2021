LDnClustering2 <- function (snp, map, nSNPs = 1000, w1 = 10, w2 = 100, LD_threshold1 = 0.5, 
    LD_threshold2 = 0.7, PC_threshold = 0.8, mc.cores = 1, plot.network = NULL, 
    threshold_net = 0.9) 
{
    out <- list()
    snp.id <- 1:nrow(map)
    snp.pos.org <- map[, 2]
    for (i in unique(map[, 1])) {
        cat("Finding windows \n")
        Chr <- map[, 1] == i
        snp.pos <- snp.pos.org[Chr]
        snp.id.chr <- snp.id[Chr]
        nL <- length(snp.pos)
        gds_name <- paste0(gsub(":", "", Sys.time()), ".gds")
        SNPRelate::snpgdsCreateGeno(gds_name, genmat = snp[, 
            Chr], sample.id = 1:nrow(snp), snp.pos = snp.id.chr, 
            snpfirstdim = FALSE)
        file <- SNPRelate::snpgdsOpen(gds_name)
        LDmat <- SNPRelate::snpgdsLDMat(file, method = "r", slide = 1, 
            verbose = FALSE)$LD^2
        SNPRelate::snpgdsClose(file)
        system(paste0("rm ", gds_name))
        nWindows <- as.integer(nL/nSNPs)
        if (nWindows > 1) {
            temp <- which(slideFunct(LDmat[1, ], w1) < LD_threshold1)
            hotspots <- sapply(seq(1, nL, length.out = nWindows + 
                1), function(i) {
                which.min(abs(i - temp))
            })
            hotspots <- temp[hotspots]
            hotspots[length(hotspots)] <- nL
            hotspots[1] <- 0
            Windows <- lapply(1:(length(hotspots) - 1), function(x) {
                snp.id.chr[(hotspots[x] + 1):(hotspots[x + 1])]
            })
        }
        else {
            Windows <- list(snp.id.chr)
        }
        cat(paste0("Number of windows: ", length(Windows), "; window  sizes: ", 
            paste(sapply(Windows, length), collapse = ":"), "\n"))
        out[[i]] <- parallel::mclapply(1:length(Windows), function(w) {
            cat(paste("Working on chromosome", i, ", window", 
                w, "\n"))
            Window <- Windows[[w]]
            snp.id.bin <- snp.id[which(snp.id %in% Window)]
            snp_bin <- snp[, which(snp.id %in% Window)]
            SNPRelate::snpgdsCreateGeno(paste0(w, gds_name), 
                genmat = snp_bin, sample.id = 1:nrow(snp), snp.id = snp.id.bin, 
                snpfirstdim = FALSE)
            file <- SNPRelate::snpgdsOpen(paste0(w, gds_name))
            if (w2 > length(snp.id.bin)) 
                w2 = -1
            MAT <- SNPRelate::snpgdsLDMat(file, method = "r", 
                slide = w2, verbose = FALSE)
            SNPRelate::snpgdsClose(file)
            system(paste0("rm ", paste0(w, gds_name)))
            if (w2 != -1) {
                LDmat <- data.table(remove = rep(NA, length(MAT$snp.id)))
                for (x in 1:(length(MAT$snp.id))) {
                  LDmat[(1:w2 + x)[!is.nan(MAT$LD[, x])], `:=`(as.character(x), 
                    MAT$LD[!is.nan(MAT$LD[, x]), x]^2)]
                }
                LDmat <- as.matrix(LDmat)
                LDmat <- LDmat[, -1]
            }
            else {
                for (x in 1:(length(MAT$snp.id))) {
                  LDmat <- MAT$LD^2
                }
            }
            LDmat[is.na(LDmat)] <- 0
            LDmat[upper.tri(LDmat)] <- NA
            diag(LDmat) <- NA
            colnames(LDmat) <- snp.id.bin
            rownames(LDmat) <- snp.id.bin
            g <- graph.adjacency(LDmat, mode = "lower", diag = FALSE, 
                weighted = TRUE)
            g <- delete_edges(g, which(E(g)$weight < LD_threshold1))
            d_g <- decompose.graph(g)
            PC_recursive <- function(clade) {
                lapply(clade, function(x) {
                  if (x > ntips) {
                    cl <- ape::extract.clade(tree, x)$tip.label
                  }
                  else {
                    cl <- tree$tip.label[x]
                  }
                  if (length(cl) > 1) {
                    LDmat.part <- LDmat[which(snp.id.bin %in% 
                      cl), which(snp.id.bin %in% cl)]
                    LDmat.part[upper.tri(LDmat.part)] <- t(LDmat.part)[upper.tri(LDmat.part)]
                    Median.temp <- median(LDmat.part, na.rm = TRUE)
                    if (Median.temp > LD_threshold2) {
                      snp_cl <- snp_bin[, which(snp.id.bin %in% 
                        cl)]
                      clusters.sub[[x]] <<- cl
                      PC_sc <- PC_score(snp_bin[, snp.id.bin %in% 
                        cl], PC_threshold)
                      PVE.sub[[x]] <<- PC_sc[[2]]
                      PCs.sub[[x]] <<- PC_sc[[1]]
                      MCL.sub[[x]] <<- snp_cl[, which.max(apply(LDmat.part, 
                        1, function(h) median(h, na.rm = TRUE)))]
                      Median.sub[[x]] <<- Median.temp
                      Mad.sub[[x]] <<- mad(LDmat.part, na.rm = TRUE)
                      temp <- do.call(cbind, c(list(snp_cl[, 
                        1]), lapply(2:ncol(snp_cl), function(h) {
                        if (!all(na.omit(snp_cl[, 1] == snp_cl[, 
                          h]))) {
                          if (cor.test(snp_cl[, 1], snp_cl[, 
                            h])$estimate < 0) {
                            return(ifelse(snp_cl[, h] == 2, 0, 
                              2))
                          } else {
                            return(snp_cl[, h])
                          }
                        } else {
                          return(snp_cl[, h])
                        }
                      })))
                      temp[is.na(temp)] <- -1
                      Cons.sub[[x]] <<- phyclust::find.consensus(t(temp))
                    }
                    else {
                      PC_recursive(tree$edge[, 2][tree$edge[, 
                        1] == x])
                    }
                  }
                  else {
                    clusters.sub[[x]] <<- cl
                    PVE.sub[[x]] <<- 1
                    Cons.sub[[x]] <<- MCL.sub[[x]] <<- PCs.sub[[x]] <<- as.matrix(snp_bin[, 
                      snp.id.bin %in% cl])
                    Median.sub[[x]] <<- Mad.sub[[x]] <<- NA
                  }
                })
            }
            clusters <- list()
            PVE <- list()
            PCs <- list()
            MCL <- list()
            Cons <- list()
            Median <- list()
            Mad <- list()
            for (y in 1:length(d_g)) {
                loci <- V(d_g[[y]])$name
                if (length(loci) > 1) {
                  LDmat.part <- LDmat[which(snp.id.bin %in% loci), 
                    which(snp.id.bin %in% loci)]
                  tree <- ape::as.phylo(hclust(as.dist(1 - LDmat.part), 
                    method = "complete"))
                  tree$tip.label <- loci
                  ntips <- length(tree$tip.label)
                  clusters.sub <- list()
                  PVE.sub <- list()
                  PCs.sub <- list()
                  MCL.sub <- list()
                  Cons.sub <- list()
                  Median.sub <- list()
                  Mad.sub <- list()
                  invisible(PC_recursive(tree$edge[1, 1]))
                  Keep <- !sapply(clusters.sub, is.null)
                  clusters[[y]] <- clusters.sub[Keep]
                  PVE[[y]] <- PVE.sub[Keep]
                  PCs[[y]] <- PCs.sub[Keep]
                  MCL[[y]] <- MCL.sub[Keep]
                  Cons[[y]] <- Cons.sub[Keep]
                  Median[[y]] <- Median.sub[Keep]
                  Mad[[y]] <- Mad.sub[Keep]
                }
                else {
                  clusters[[y]] <- loci
                  PVE[[y]] <- 1
                  Cons[[y]] <- MCL[[y]] <- PCs[[y]] <- as.matrix(snp_bin[, 
                    snp.id.bin %in% loci])
                  Median[[y]] <- NA
                  Mad[[y]] <- NA
                }
            }
            clusters <- flatten(clusters)
            PVE <- flatten(PVE)
            PCs <- flatten(PCs)
            MCL <- flatten(MCL)
            Cons <- flatten(Cons)
            Mad <- unlist(flatten(Mad))
            Median <- unlist(flatten(Median))
            if (!is.null(plot.network)) {
                cat("plotting network \n")
                nClust <- length(clusters)
                temp1 <- sample(rep(col_vector, ceiling(nClust/length((col_vector)))))
                Col <- unlist(lapply(1:length(clusters), function(x) {
                  ifelse(sapply(clusters, length)[x] == 1, return("black"), 
                    return(rep(temp1[1:nClust][x], sapply(clusters, 
                      length)[x])))
                }))[match(colnames(LDmat), unlist(clusters))]
                temp2 <- sample(rep(col_vector[-1], ceiling(nClust/length((col_vector)))))
                frame.col <- unlist(lapply(1:length(clusters), 
                  function(x) {
                    ifelse(sapply(clusters, length)[x] == 1, 
                      return("black"), return(rep(temp2[1:nClust][x], 
                        sapply(clusters, length)[x])))
                  }))[match(colnames(LDmat), unlist(clusters))]
                g <- graph.adjacency(LDmat, mode = "lower", diag = FALSE, 
                  weighted = TRUE)
                V(g)$color <- Col
                V(g)$frame.color <- frame.col
                E(g)$weight <- round(E(g)$weight, 5)
                g <- delete.edges(g, which(E(g)$weight <= threshold_net))
                g <- delete.vertices(g, which(degree(g) == 0))
                par(mar = c(0, 0, 2, 0))
                png(paste0(paste0(plot.network, "Chr", i, "_window", 
                  w, "_LD2:", LD_threshold2), ".png"), res = 300, 
                  width = 6, height = 6, units = "in")
                plot.igraph(g, layout = layout.fruchterman.reingold, 
                  vertex.size = 3, vertex.label.dist = NULL, 
                  edge.width = 1, vertex.label = NA)
                title(main = paste(" @", threshold_net, sep = "", 
                  "; Chr", i, "; window", w, "; LD2=", LD_threshold2))
                dev.off()
            }
            cat(paste0(length(PCs), " clusters and ", sum(sapply(PCs, 
                ncol)), " PCs extracted from of ", ncol(LDmat), 
                " original SNPs, on average explaining ", signif(mean(sapply(PVE, 
                  max)), 2) * 100, "% of the variation \n"))
            return(list(clusters = clusters, PCs = PCs, PVE = PVE, 
                MCL = MCL, Cons = Cons, Median = Median, Mad = Mad, 
                chr = i))
        }, mc.cores = mc.cores)
    }
    cat("preparing output \n")
    if (length(out) == 1) 
        out <- c(list(NULL), out)
    out <- out[sapply(out, length) > 0]
    clusters <- flatten(lapply(out, function(chromosome) {
        lapply(chromosome, function(window) {
            window[[1]]
        })
    }))
    clusters <- lapply(clusters, as.numeric)
    MCL <- do.call(cbind, flatten(lapply(out, function(chromosome) {
        lapply(chromosome, function(window) {
            window[[4]]
        })
    })))
    Cons <- do.call(cbind, flatten(lapply(out, function(chromosome) {
        lapply(chromosome, function(window) {
            window[[5]]
        })
    })))
    cluster_PCs <- do.call(cbind, flatten(lapply(out, function(chromosome) {
        lapply(chromosome, function(window) {
            window[[2]]
        })
    })))
    Window <- unlist(lapply(1:length(out), function(chromosome) {
        lapply(1:length(out[[chromosome]]), function(window) {
            rep(window, ncol(do.call(cbind, out[[chromosome]][[window]][[2]])))
        })
    }))
    PC <- unlist(sapply(out, function(chromosome) {
        sapply(chromosome, function(window) {
            sapply(sapply(window[[2]], function(bin) ifelse(is.matrix(bin), 
                ncol(bin), 1)), function(bin) 1:bin)
        })
    }))
    PVE <- unlist(lapply(out, function(chromosome) {
        lapply(chromosome, function(window) {
            window[[3]]
        })
    }))
    Grp <- rep(NA, length(PC))
    j <- 0
    for (x in 1:length(PC)) {
        if (PC[x] == 1) {
            j <- j + 1
            Grp[x] <- j
        }
        else {
            Grp[x] <- j
        }
    }
    Chr <- unlist(lapply(1:length(out), function(chromosome) {
        lapply(1:length(out[[chromosome]]), function(window) {
            rep(out[[chromosome]][[window]][[8]], ncol(do.call(cbind, 
                out[[chromosome]][[window]][[2]])))
        })
    }))
    Median <- unlist(sapply(out, function(chromosome) {
        sapply(chromosome, function(window) {
            window[[6]]
        })
    }))[Grp]
    MAD <- unlist(sapply(out, function(chromosome) {
        sapply(chromosome, function(window) {
            window[[7]]
        })
    }))[Grp]
    Pos <- sapply(clusters, function(x) mean(map[x, 2]))[Grp]
    Min <- sapply(clusters, function(x) min(map[x, 2]))[Grp]
    Max <- sapply(clusters, function(x) max(map[x, 2]))[Grp]
    Range <- abs(Min - Max)
    nSNPs <- sapply(clusters, function(x) length(x))[Grp]
    Range <- abs(Min - Max)
    cluster_summary <- data.frame(Chr, Window, Pos = Pos, Min, 
        Max, Range, nSNPs = nSNPs, Median = Median, MAD = MAD, 
        PC = PC, PVE, Grp)
    colnames(cluster_summary) <- c("Chr", "Window", "Pos", "Min", 
        "Max", "Range", "nSNPs", "Median", "MAD", "PC", "PVE", 
        "Grp")
    cat("Done")
    return(list(cluster_summary = cluster_summary, cluster_PCs = cluster_PCs, 
        clusters = clusters, MCL = MCL, Cons = Cons))
}

###################################################################################

PC_score <- function (x_A, PC_threshold) 
{
    cMean_A <- colMeans(x_A, na.rm = TRUE)
    for (i in 1:dim(x_A)[2]) x_A[, i] <- x_A[, i] - cMean_A[i]
    if (any(is.na(x_A))) {
        invisible(lapply(1:ncol(x_A), function(k) {
            x_A[is.na(x_A[, k]), k] <<- mean(x_A[, k], na.rm = TRUE)
        }))
    }
    eigen_A <- eigen(var(x_A))
    scores_A <- x_A %*% eigen_A$vectors
    percentage_A = (cumsum(eigen_A$values))/sum(eigen_A$values)
    a = percentage_A < PC_threshold
    a[which(!a)[1]] <- TRUE
    zeta = as.matrix(scores_A[, a])
    return(list(zeta, percentage_A[a]))
}


flatten <- function (x) 
{
    if (!inherits(x, "list")) 
        return(list(x))
    else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
}

#################################################################

col_vector <- c( "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
 "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
 "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99",
 "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
 "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC",
 "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
 "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
 "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
 "#CCEBC5", "#FFED6F")



