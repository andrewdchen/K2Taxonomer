#' Interactive K2 dendrogram
#'
#' Create an interactive dendrogram of the K2 Taxonomer results
#' @param K2res A list object. The output of runK2tax().
#' @param annot One of 'partition', 'pathways', 'genes, to specify annotation type.
#' @param alpha Significance level with which to filter pathways/genes.
#' @param limit Number of pathways/genes to include in output.
#' @references
#'  \insertRef{reed_2020}{K2Taxonomer}
#' @return An interactive dendrogram created by `visNetwork::visNetwork()`.
#' @export
#' @examples
#'
#' ## Read in K2 Taxonomer results
#' data(K2res)
#'
#' ## Generate interactive dendrogram
#' K2visNetwork(K2res)
#'

K2visNetwork <- function(K2res, 
                         annot = c("partition", "pathways", "genes"),
                         alpha = 0.05,
                         limit = 5,
                         recursive = FALSE) {

    ## Run checks
    .isK2(K2res)

    ## K2 algorithm
    if (length(K2results(K2res)) == 0) {
        stop("No results found. Please run K2tax() or runK2Taxonomer().\n")
    }

    ## Get results list
    dendrogram <- K2dendro(K2res)
    K2res <- K2results(K2res)
    
    ## Generate Matrix from visNetwork
    mat <- matrix(0, nrow=length(K2res), ncol=length(K2res[[1]]$obs[[1]]) +
        length(K2res[[1]]$obs[[2]]))
    colnames(mat) <- c(K2res[[1]]$obs[[1]], K2res[[1]]$obs[[2]])
    for (i in seq_len(length(K2res))) {
        mat[i, K2res[[i]]$obs[[1]]] <- 1
        mat[i, K2res[[i]]$obs[[2]]] <- 2
    }
    rownames(mat) <- names(K2res)

    ## Calculate sizes
    sizes <- apply(mat, 1, function(x) sum(x != 0))

    ## Add Labels
    annot <- match.arg(annot)
    if (annot == "partition") {
      titles <- unlist(lapply(K2res, function(x) paste("Probability:", signif(x$bootP, 2), "<br>", 
                                                       "Members(Edge:1):", length(x$obs[[1]]), "<br>", 
                                                       "Members(Edge:2):", length(x$obs[[2]]), "<br>",
                                                       "Stability(Edge:1):", signif(x$stability$clusters[[1]], 2), "<br>", 
                                                       "Stability(Edge:2):", signif(x$stability$clusters[[2]], 2))))
    } else if (annot == "pathways") {
      if (!recursive) {
        edge1_annots <- lapply(K2res,
                               function(x) x$dsse %>%
                                 dplyr::filter(edge == "1" & fdr < alpha) %>%
                                 dplyr::arrange(desc(coef)) %>%
                                 dplyr::slice(1:limit) %>%
                                 dplyr::pull(category))
        edge2_annots <- lapply(K2res,
                               function(x) x$dsse %>%
                                 dplyr::filter(edge == "2" & fdr < alpha) %>%
                                 dplyr::arrange(desc(coef)) %>%
                                 dplyr::slice(1:limit) %>%
                                 dplyr::pull(category))
        titles <- unlist(lapply(names(K2res), function(x) paste("Edge 1:", str_flatten(edge1_annots[[x]], collapse = " "), "<br>",
                                                                "Edge 2:", str_flatten(edge2_annots[[x]], collapse = " "))))
      } else {
        titles <- recursive_annotation(dendrogram = dendrogram, K2res = K2res)
        titles <- unlist(titles)
        titles <- titles[names(K2res)]
      }

    } else if (annot == "genes") {
      edge1_annots <- lapply(K2res, 
                             function(x) x$dge %>% 
                               dplyr::filter(edge == "1" & fdr < alpha) %>% 
                               dplyr::arrange(desc(coef)) %>% 
                               dplyr::slice(1:limit) %>% 
                               dplyr::pull(gene))
      edge2_annots <- lapply(K2res, 
                             function(x) x$dge %>% 
                               dplyr::filter(edge == "2" & fdr < alpha) %>% 
                               dplyr::arrange(desc(coef)) %>% 
                               dplyr::slice(1:limit) %>% 
                               dplyr::pull(gene))
      titles <- unlist(lapply(names(K2res), function(x) paste("Edge 1:", str_flatten(edge1_annots[[x]], collapse = " "), "<br>", 
                                                              "Edge 2:", str_flatten(edge2_annots[[x]], collapse = " "))))
    }

    names(titles) <- names(K2res)

    ## initialize leafe names
    len <- length(K2res)
    nalphabets <- ceiling(ncol(mat)/length(letters))
    nAlphabets <- ceiling(len/length(letters))
    alphabets <- unlist(lapply(seq_len(nalphabets), function(x) {

        vapply(letters, function(y) paste(rep(y, x), collapse=""),
            character(1))

    }))
    ALPHABETS <- unlist(lapply(seq_len(nAlphabets), function(x) vapply(LETTERS,
        function(y) paste(rep(y, x), collapse=""), character(1))))
    source <- c()
    target <- c()
    k <- 1

    ## Get edges
    for (i in seq_len(nrow(mat))) {

        source <- c(source, rep(rownames(mat)[i], 2))
        matRow <- mat[i, ]
        sub1 <- which(matRow == 1)[1]
        sub2 <- which(matRow == 2)[1]
        matSub1 <- mat[-seq_len(i), sub1]
        names(matSub1) <- rownames(mat)[-seq_len(i)]
        matSub2 <- mat[-seq_len(i), sub2]
        names(matSub2) <- rownames(mat)[-seq_len(i)]

        target1 <- names(matSub1)[which(matSub1 != 0)[1]]
        if (is.na(target1)) {
            target1 <- alphabets[k]
            sizes[target1] <- sum(matRow == 1)
            titles[target1] <- paste(colnames(mat)[matRow ==
                1], collapse="<br>")
            k <- k + 1
        }
        target2 <- names(matSub1)[which(matSub2 != 0)[1]]
        if (is.na(target2)) {
            target2 <- alphabets[k]
            sizes[target2] <- sum(matRow == 2)
            titles[target2] <- paste(colnames(mat)[matRow ==
                2], collapse="<br>")
            k <- k + 1
        }
        target <- c(target, target1, target2)
    }

    ## Set terminal nodes to 0 sizes
    sizes[names(sizes) %in% alphabets] <- 0

    ## Add Labels
    labs <- titles
    labs <- gsub("<br>", "\n", labs)
    labs[names(labs) %in% ALPHABETS] <- names(labs)[names(labs) %in%
        ALPHABETS]

    ## Add shapes
    shapes <- rep("diamond", length(sizes))
    names(shapes) <- names(sizes)
    shapes[names(shapes) %in% alphabets] <- "box"

    ## Get terminal node level
    levs <- vapply(2:nrow(mat), function(x) {
        matSub <- mat[seq_len(x), ]
        matSub <- matSub[, matSub[x, ] != 0]
        matSub <- matSub[, 1]
        sum(matSub != 0)
    }, numeric(1))
    levs <- c(1, levs)
    levs <- c(levs, rep(max(levs) + 1, length(sizes) - length(levs)))
    names(levs) <- names(sizes)

    ## Generate network
    nodes <- data.frame(id=names(levs), title=titles, level=levs,
        label=labs, shape=shapes, stringsAsFactors=FALSE)
    edges <- data.frame(from=source, to=target, stringsAsFactors=FALSE)

    ## Generate plot
    p <- visNetwork(nodes, edges) %>% visEdges(arrows="to") %>%
        visHierarchicalLayout(direction="LR")

    return(p)
}

recursive_annotation <- function(dendrogram, K2res, exclude=c(), annots = list(), alpha = 0.05, limit = 10) {
  if (is.na(attr(dendrogram, "label"))) {
    print("leaf")
  } else {
    curr_label <- attr(dendrogram, "label")
    left_label <- attr(dendrogram[[1]], "label")
    right_label <- attr(dendrogram[[2]], "label")
    
    # get annotations for Edge 1 and Edge 2
    edge1_annots <- K2res[[curr_label]][["dsse"]] %>% 
      dplyr::filter(edge == "1" & fdr < alpha) %>% 
      dplyr::arrange(desc(coef)) %>% 
      dplyr::slice(1:limit) %>% 
      dplyr::pull(category)
    edge2_annots <- K2res[[curr_label]][["dsse"]] %>% 
      dplyr::filter(edge == "2" & fdr < alpha) %>% 
      dplyr::arrange(desc(coef)) %>% 
      dplyr::slice(1:limit) %>% 
      dplyr::pull(category)
    
    # Remove excluded annotations 
    edge1_annots <- edge1_annots[!(edge1_annots %in% exclude)]
    edge2_annots <- edge2_annots[!(edge2_annots %in% exclude)]
    
    annots[[curr_label]] <- paste("Edge 1:", str_flatten(edge1_annots, collapse = " "), "<br>", 
                                  "Edge 2:", str_flatten(edge2_annots, collapse = " "))
    
    # Add current edge annotations to excluded sets
    exclude_l <- unique(c(exclude, edge2_annots))
    exclude_r <- unique(c(exclude, edge1_annots))
    
    # Keep annotating Edge 1
    if(!is.na(left_label)) {
      #print(paste0("Edge 1 node label ", left_label, " my excluded annotations are ", str_flatten(exclude_l, collapse = " "), " my annotations are ", str_flatten(edge1_annots, collapse = " ")))
      annots <- recursive_annotation(dendrogram[[1]], K2res, exclude = exclude_l, annots = annots)
    }
    
    # Keep annotating Edge 2
    if(!is.na(right_label)) {
      #print(paste0("Edge 2 node label ", right_label, " my excluded annotations are ", str_flatten(exclude_r, collapse = " "), " my annotations are ", str_flatten(edge2_annots, collapse = " ")))
      annots <- recursive_annotation(dendrogram[[2]], K2res, exclude = exclude_r, annots = annots)
    }
  }
  return(annots)
}