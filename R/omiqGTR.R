
.read_model <-
  function(path) {
    
    ## Import the omiqgt JSON file ----

    if (!file.exists(path)) {
      stop("File not found: ", path)
    }
    raw <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    if (is.null(raw$tree$nodes) || is.null(raw$tree$filterContainers)) {
      stop("Not an .omiqgt file: missing $tree$nodes / $tree$filterContainers.")
    }

    nodes <- raw$tree$nodes
    containers <- raw$tree$filterContainers # gate geometries
    node_ids <- names(nodes)

    ## Parse relationships between nodes ----

    parent_of <- # name ~ child, value ~ parent
      vapply(
        X = nodes,
        FUN = function(n) {
          p <- n$parentId
          if (is.null(p)) { "" } else { as.character(p) }
        },
      FUN.VALUE = character(1)
      )
    is_root <- # root mask (multiple roots possible if no explicit single root)
      parent_of == "" | !(parent_of %in% node_ids)
    ord_of <- # numeric flags per node for ordering within respective subtrees
      vapply(
        X = nodes,
        FUN = function(n) {
          o <- n$ord
          if (is.null(o)) { NA_real_ } else { as.numeric(o) }
        },
        FUN.VALUE = numeric(1)
      )
    {
      children <- # name ~ parent, values ~ children
        setNames(
          rep(list(character(0)), length(node_ids)),
          node_ids
        )
      for (id in node_ids) {
        if (!is_root[[id]]) {
          children[[parent_of[[id]]]] <-
            c(children[[parent_of[[id]]]], id)
        }
      }
    }
    
    ## Sort nodes ----

    ord_sort <- function(ids) { ids[order(ord_of[ids], ids)] }
    children <- lapply(children, ord_sort)
    roots    <- ord_sort(node_ids[is_root])

    res <- list(
      "meta" = raw[setdiff(names(raw), "tree")],
      "nodes" = nodes,
      "containers" = containers,
      "roots" = roots,
      "children" = children,
      "parent_of" = parent_of
    )
  }

.node_name <-
  function(model, id) {
    
    ## Grab node name based on its ID ----

    fc <- model$containers[[model$nodes[[id]]$filterContainerId]]
    nm <- if (is.null(fc)) { NULL } else { fc$name }
    if (is.null(nm) || !nzchar(nm)) { id } else { nm }
  }

.container_desc <-

  ## Make function to describe filter ----
  
  function(fc) {
    if (is.null(fc)) { return("<missing>") }
    if (identical(fc$containerType, "AtomicFilterContainer")) {
      df <- fc$defaultFilter
      sprintf("%s on (%s, %s)", df$type, df$f1, df$f2)
    } else {
      n <- length(fc$filterContainerIds)
      sprintf("%s of %d filter%s", fc$type, n, if (n == 1) "" else "s")
    }
  }

.nested_tree <-
  function(model) {

    ## Define recursive tree building function ----
    
    build <- function(id) {

      fc_id <- model$nodes[[id]]$filterContainerId
      kids <- lapply(model$children[[id]], build) # recurse
      list(
        "id" = id,
        "name" = .node_name(model, id),
        "container_id" = fc_id,
        "desc" = .container_desc(model$containers[[fc_id]]),
        "is_leaf" = length(kids) == 0L,
        "children" = kids
      )
    }

    ## Create a rooted tree ----

    list(
      "id" = "__root__",
      "name" = "(ungated)",
      "container_id" = NA_character_,
      "desc" = "root of all populations",
      "is_leaf" = FALSE,
      "children" = lapply(model$roots, build)
    )
  }

.depth_of <-
  function(model, id) {
    
    ## Count node depth (root ~ 0) ----

    d <- 1L # depth
    p <- model$parent_of[[id]]
    while (nzchar(p) && p%in%names(model$nodes)) {
      d <- d + 1L
      p <- model$parent_of[[p]]
    }
    d
  }

.node_order <-
  function(model) {

    ## Grab ordered node IDs recursively ----
    
    out <- character(0)
    rec <- function(id) {
      out[[length(out) + 1L]] <<- id
      for (c in model$children[[id]]) {
        rec(c)
      }
    }
    for (r in model$roots) {
      rec(r)
    }
    out
  }

.node_paths <-
  function(model, sep = "|") {
    
    ## Fetch full node paths recursively ----

    res <- character(0)
    rec <- function(id, prefix) {
      p <-
        paste(c(prefix, .node_name(model, id)), collapse = sep)
      res[[id]] <<- p
      for (c in model$children[[id]]) {
        rec(c, p)
      }
    }
    for (r in model$roots) {
      rec(r, character(0))
    }
    res
  }

.point_in_polygon <-
  function(x, y, vx, vy) {
    
    ## Get vector of flags: x,y in vx,vy-defined polygon? ----

    k <- length(vx)
    inside <- logical(length(x))
    j <- k
    for (i in seq_len(k)) {
      xi <- vx[i]; yi <- vy[i]; xj <- vx[j]; yj <- vy[j]
      crosses <-
        ((yi>y) != (yj > y)) &
        (x < (xj-xi)*(y-yi)/(yj-yi)+xi)
      crosses[is.na(crosses)] <- FALSE
      inside <- xor(inside, crosses)
      j <- i
    }
    inside[is.na(x) | is.na(y)] <- NA
    inside
  }

.as_gate_matrix <-
  function(mat, need, label = "") {
    
    ## Validate expression data matrix ----

    lab <- if (nzchar(label)) { paste0(" '", label, "'") } else { "" }
    cols <- colnames(mat)
    if (is.null(cols)) { 
      stop("The input must have column names (channel/marker names).")
    }
    missing <- setdiff(need, cols)
    if (length(missing)) {
      stop(
        sprintf(
            "Gate%s needs channel(s) not in the input: %s\nAvailable columns: %s",
            lab, paste(missing, collapse = ", "), paste(cols, collapse = ", ")
        )
      )
    }
      
    if (is.data.frame(mat)) {
      nonnum <-
        need[!vapply(need, function(c) is.numeric(mat[[c]]), logical(1))]
      if (length(nonnum)) {
        stop(
          sprintf(
            "Gate%s: channel column(s) are not numeric: %s",
            lab, paste(nonnum, collapse = ", ")
          )
        )
      }
        
      return(as.matrix(mat[, need, drop = FALSE]))
    }
    if (!is.numeric(mat)) {
      stop(
        sprintf(
          "Gate%s: input matrix is not numeric (storage mode '%s').",
          lab, storage.mode(mat)
        )
      )
    }
    mat
  }

.make_rectangle_fn <-
  function(fdef, label = "") {
    
    ## Extract rectangle vertices ----

    f1 <- fdef$f1; f2 <- fdef$f2
    g <- function(v, default) { # coordinate or default (fallback)
        if (is.null(v)) { default } else as.numeric(v)
    } 
    x1lo <- g(fdef$min$f1Val, -Inf); x1hi <- g(fdef$max$f1Val,  Inf)
    x2lo <- g(fdef$min$f2Val, -Inf); x2hi <- g(fdef$max$f2Val,  Inf)

    ## Ensure highs are higher than lows ----
    
    if (isTRUE(x1lo>x1hi)) { t <- x1lo; x1lo <- x1hi; x1hi <- t }
    if (isTRUE(x2lo > x2hi)) { t <- x2lo; x2lo <- x2hi; x2hi <- t }

    ## Define gating function ----

    fn <- function(mat) {
      mat <- # validated expression matrix
        .as_gate_matrix(mat, c(f1, f2), label)
      a <- mat[, f1]; b <- mat[, f2]
      (a >= x1lo) & (a <= x1hi) & (b >= x2lo) & (b <= x2hi) # membership mask
    }
    attr(fn, "gate_type") <- "RectangleGate"
    attr(fn, "channels") <- c(f1, f2)
    fn
  }

.make_polygon_fn <-
  function(fdef, label = "") {
    
    ## Extract polygon vertices ----

    f1 <- fdef$f1; f2 <- fdef$f2
    vx <- vapply(
      X = fdef$vertices,
      FUN = function(v) { as.numeric(v$f1Val) },
      FUN.VALUE = numeric(1)
    )
    vy <- vapply(
      X = fdef$vertices,
      FUN = function(v) { as.numeric(v$f2Val) },
      FUN.VALUE = numeric(1)
    )
    fn <- function(mat) {
      mat <- # validated expression matrix
        .as_gate_matrix(mat, c(f1, f2), label)
      .point_in_polygon(mat[, f1], mat[, f2], vx, vy) # membership mask
    }
    attr(fn, "gate_type") <- "PolygonGate"
    attr(fn, "channels") <- c(f1, f2)
    fn
  }

.effective_filter <- function(fc, file_id) {

  ## Fetch correct filter for specific file ----

  if (!is.null(file_id)) {
    file_id <- as.character(file_id)
    pf <- fc$perFileFilters
    if (!is.null(pf) && !is.null(pf[[file_id]])) {
      return(pf[[file_id]])
    }
  }
  fc$defaultFilter
}

.build_gates <-
  function(model, file_id = NULL) {

    ## Fetch membership function per container ID ----

    containers <- model$containers
    cache <- new.env(parent = emptyenv())
    building <- new.env(parent = emptyenv()) # "already-built" flags

    make <- function(cid) {

      if (!is.null(cache[[cid]])) {
        return(cache[[cid]])
      }
      fc <- containers[[cid]]
      if (is.null(fc)) {
        stop("Referenced filter container not found: ", cid)
      }
      if (!is.null(building[[cid]])) {
        stop("Cycle detected at container ", cid)
      }
      assign(cid, TRUE, envir = building)
      on.exit( # clear at end of function call
          rm(list = cid, envir = building),
          add = TRUE
      )
      label <- # human-readable label if exists
        if (is.null(fc$name)) { cid } else { fc$name }
      if (identical(fc$containerType, "AtomicFilterContainer")) {

        fdef <- .effective_filter(fc, file_id)
        fn <-
          switch(
            fdef$type,
            "RectangleGate" = .make_rectangle_fn(fdef, label),
            "PolygonGate"   = .make_polygon_fn(fdef, label),
            stop("Unsupported atomic gate type: ", fdef$type)
          )
      } else if (identical(fc$containerType, "CompoundFilterContainer")) {
        
        ## Note: only `AND` and `NOT` compound filters for compounds supported
        
        child_ids <- unlist(fc$filterContainerIds)
        child_fns <- lapply(child_ids, make)
        op <- fc$type
        fn <- local({ child_fns; op
          function(mat) {
            inter <- Reduce(`&`, lapply(child_fns, function(f) f(mat)))
            if (identical(op, "NOT")) !inter else inter
          } })
        attr(fn, "gate_type") <- op
        attr(fn, "channels") <-
          unique(unlist(lapply(child_fns, attr, "channels")))
      } else {
        stop("Unknown containerType: ", fc$containerType)
      }
      attr(fn, "gate_name") <- label
      attr(fn, "container_id") <- cid
      assign(cid, fn, envir = cache)
      fn
    }
    setNames(
      lapply(names(containers), make),
      names(containers)
    )
  }
  

.gate_channels <-
  function(model) {
    
    ## Fetch all channels required by model ----

    ch <- character(0)
    for (fc in model$containers) {
      if (!identical(fc$containerType, "AtomicFilterContainer")) {
        next # don't extract channels from compound gates: get them from atomics
      }
      defs <-
        c(
          list(fc$defaultFilter),
          if (!is.null(fc$perFileFilters)) { fc$perFileFilters } else { list() }
        )
      for (df in defs) {
        ch <- c(ch, df$f1, df$f2)
      }
    }
    sort(unique(ch))
  }

.as_dendrogram_node <-
  function(node, counter) {

    ## Convert subtree to `stats::dendrogram` structure ----

    kids <- node$children
    if (length(kids)==0L) { # leaf node
      idx <- counter$n <- counter$n+1L
      d <- idx
      attr(d, "members") <- 1L
      attr(d, "height") <- 0
      attr(d, "leaf") <- TRUE
      attr(d, "label") <- node$name
      attr(d, "midpoint") <- 0
      attr(d, "node_id") <- node$id
      attr(d, "gate_id") <- node$container_id
      class(d) <- "dendrogram"
      return(d)
    }
    built <- lapply(kids, .as_dendrogram_node, counter = counter)
    memb <- vapply(built, attr, integer(1), "members")
    height <- vapply(built, attr, numeric(1), "height")
    mids <- vapply(built, attr, numeric(1), "midpoint")
    M <- sum(memb)
    d <- built
    attr(d, "members") <- as.integer(M)
    attr(d, "height") <- max(height)+1
    attr(d, "midpoint") <- (mids[1]+(M-memb[length(memb)])+mids[length(mids)])/2
    attr(d, "label") <- node$name
    attr(d, "node_id") <- node$id
    attr(d, "gate_id") <- node$container_id
    class(d) <- "dendrogram"
    d
  }

#' Parse an `.omiqgt` file
#' 
#' Parses and `.omiqgt` (OMIQ gating hierarchy) file downloaded from an OMIQ 
#' workflow by opening a *Gating* task and pressing *Ctrl+Shift+D*.
#'
#' @param path Path to the `.omiqgt` file.
#' @return Object of class `c("GatingTree", "dendrogram")`: gating hierarchy as 
#'         dendrogram (the file's root populations are joined under a synthetic
#'         *"(ungated)"* root). The parsed model, per-gate membership functions,
#'         and root->node paths are carried in `attr(, "model")` and are used by
#'         `print()`, `plot()`, and `gate()`. Standard dendrogram operations
#'         (`labels()`, `order.dendrogram()`, `cut()`, `str()`) should  work too.
#' @export
parse_omiqgt <- function(path) {
  
  ## Read the gating hierarchy ----

  model <- .read_model(path)
  model$nested <- # internal nested-tree structure
    .nested_tree(model)
  model$gates <- # vectorized functions, one per gate
    .build_gates(model)
  model$order <- # ordering of nodes for printing/plotting
    .node_order(model)
  model$paths <- # unique, full paths to each node
    .node_paths(model)
  model$channels <- # all required fluorophore/isotope channels
    .gate_channels(model)

  ## Build a `stats::dendrogram` structure ----

  counter <- new.env()
  counter$n <- 0L
  dnd <- .as_dendrogram_node(model$nested, counter)
  attr(dnd, "model") <- model
  class(dnd) <- c("GatingTree", "dendrogram")
  dnd
}

.model <- function(x) {

  ## Fail if object not `GatingTree` ----

  m <- attr(x, "model")
  if (is.null(m)) {
    stop("Not a GatingTree produced by `parse_omiqgt()`.")
  }
  m
}

#' Print a parsed OMIQ gating tree object
#'
#' Outputs the gating hierarchy and other information about a `GatingTree` object
#' generated using the `omiqGTR::parse_omiqgt()` function.
#' 
#' @param x A `GatingTree` object created by `parse_omiqgt()`.
#' @param max.channels Max how many channel names to list (rest summarized).
#'                     Defaults to 12.
#' @param tree Logical; also print the indented gating hierarchy? Defaults to 
#'             TRUE.
#' @export
print.GatingTree <-
  function(x, max.channels = 12, tree = TRUE, ...) {
    m <- .model(x)
    n_atomic <- # number of atomic gates
      sum(
        vapply(
          X = m$containers,
          FUN = function(c) {
            identical(c$containerType, "AtomicFilterContainer")
          },
          FUN.VALUE = logical(1)
        )
      )
    n_compound <- # number of compound gates
      length(m$containers)-n_atomic
    depths <- # depth per node
      vapply(
        X = names(m$nodes),
        FUN = function(id) { .depth_of(m, id) },
        FUN.VALUE = integer(1)
      )

    cat("<GatingTree>\n")
    cat(
      sprintf(
        "  populations (nodes) : %d\n", length(m$nodes)
      )
    )
    cat(
      sprintf(
        "  gate definitions    : %d (%d atomic, %d compound)\n",
        length(m$containers), n_atomic, n_compound
      )
    )
    cat(
      sprintf(
        "  root populations    : %d\n", length(m$roots)
      )
    )
    cat(
      sprintf(
        "  max depth           : %d\n",
        if (length(depths)) { max(depths) } else { 0L }
      )
    )
    ch <- m$channels
    shown <- if (length(ch)>max.channels) {
      paste0(
        paste(ch[seq_len(max.channels)], collapse = ", "),
        sprintf(", ... (+%d more)", length(ch)-max.channels)
      )
    } else {
      paste(ch, collapse = ", ")
    }
      
    cat(
      sprintf(
        "  channels (%d)        : %s\n", length(ch), shown
      )
    )
    has_pf <- # are there per-file filter adjustments?
      any(vapply(m$containers, function(c) !is.null(c$perFileFilters), logical(1)))
    if (has_pf) {
      cat(
        "  note                : one or more gates have per-file overrides ",
        "(see gate(file_id=...)).\n",
        sep = ""
      )
    }
    if (isTRUE(m$meta$inverted)) {
      cat("  note                : top-level 'inverted' flag is TRUE.\n")
    }

    if (tree) {
      cat("\n")
      nested <- m$nested
      walk <-
        function(node, prefix, is_last, is_root) {
          if (is_root) {
            cat(node$name, "\n", sep = "")
          } else {
            cat(
              prefix, if (is_last) { "└─ " } else { "├─ " },
              node$name, "\n",
              sep = ""
            )
          }
          kids <- node$children
          cp <-
            if (is_root) {
              ""
            } else {
              paste0(prefix, if (is_last) "   " else "│  ")
            }
          for (i in seq_along(kids)) {
            walk(kids[[i]], cp, i == length(kids), FALSE)
          }
        }
      walk(nested, prefix = "", is_last = TRUE, is_root = TRUE)
    }
    cat(
      sprintf(
        paste0(
          "\n  use plot() for the dendrogram, ",
          "gate(x, data) to gate %d populations.\n"
        ),
        length(m$nodes)
      )
    )
    invisible(x)
  }

#' Visualize a parsed OMIQ gating tree object
#' 
#' Draws a dendrogram of gates as defined by a `GatingTree` object generated 
#' using the `omiqGTR::parse_omiqgt()` function.
#' 
#' `GatingTree` object
#' generated using the `omiqGTR::parse_omiqgt()` function.
#'
#' @param x A `GatingTree` object created by `parse_omiqgt()`.
#' @param horizontal TRUE (default) = root at left, populations read top-down.
#'                   FALSE = root at top, leaves along the bottom.
#' @param cex Base label size (leaves). Internal labels use cex*node.cex.
#'            Defaults to 0.62.
#' @param node.cex Relative size of internal-node labels. Defaults to 1.
#' @param align.leaves Align leaf labels in a column with dotted leader lines?
#'                     Defaults to TRUE.
#' @param edge.col,edge.lwd Branch color / width. Default to `"grey65"`` and 1.
#' @param leaf.col,node.col Label colours for leaves / internal nodes. Default to
#'                          `"black"` and `"#1f6feb"` (blue).
#' @param point.pch,point.cex,point.bg Node markers (set `point.pch = NA` to 
#'                                     hide). Default to 21, 0.7, and `"white"`.
#' @param main,mar Plot title and margins. Default to *"OMIQ gating hierarchy"* 
#'                 and `NULL` (default margins).
#' @param ... Extra arguments passed to `graphics::par()`.
#' @export
plot.GatingTree <-
  function(
    x,
    horizontal   = TRUE,
    cex          = 0.62,
    node.cex     = 1.0,
    align.leaves = TRUE,
    edge.col     = "grey65",
    edge.lwd     = 1,
    leaf.col     = "black",
    node.col     = "#1f6feb",
    point.pch    = 21,
    point.cex    = 0.7,
    point.bg     = "white",
    main         = "OMIQ gating hierarchy",
    mar          = NULL,
    ...
  ) {
    m <- .model(x)
    nested <- m$nested

    ## Assign depth and subtree-position per node ---

    env <- new.env()
    env$leaf <- 0L
    nodes <- list()
    assign_layout <- function(node, depth, parent) {
      kids <- node$children
      if (length(kids) == 0L) {
        env$leaf <- env$leaf + 1L
        y <- env$leaf
      } else {
        ys <-
          vapply(
            X = kids,
            FUN = function(k) { assign_layout(k, depth + 1L, node$id) },
            FUN.VALUE = numeric(1)
          )
        y <- (ys[1]+ys[length(ys)])/2 # centering over extremes
      }
      nodes[[node$id]] <<-
        list(
          "name" = node$name,
          "depth" = depth,
          "y" = y,
          "is_leaf" = length(kids)==0L,
          "parent" = parent,
          "children" =
            vapply(
              X = kids, FUN = function(k) { k$id }, FUN.VALUE = character(1))
            )
      y
    }
    assign_layout(node = nested, depth = 0L, parent = NA_character_)
    
    L <- env$leaf
    maxdepth <-
      max(
        vapply(X = nodes, FUN = function(n) { n$depth }, FUN.VALUE = numeric(1))
      )

    ## Map coordinates ----
    
    X <- function(depth, y) { if (horizontal) { depth } else { y } }
    Y <- function(depth, y) { if (horizontal) { -y } else { -depth } }

    if (is.null(mar)) {
      mar <- if (horizontal) { c(2, 1, 3, 12) } else { c(10, 2, 3, 2) }
    }
    op <- graphics::par(mar = mar, ...)
    on.exit( # reset margins at end of function call
      graphics::par(op), add = TRUE
    )

    leaf_lab_depth <- if (align.leaves) { maxdepth } else { NA }

    if (horizontal) {
      xlim <- c(-0.35, maxdepth + 0.15) # padding so there is room for labels
      ylim <- c(-(L + 0.6), -0.4)
    } else {
      xlim <- c(0.4, L + 0.6)
      ylim <- c(-(maxdepth + 0.5), 0.4)
    }
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    if (!is.null(main) && nzchar(main)) {
      graphics::title(main = main)
    }

    ## Draw edges (branches) ----
    for (id in names(nodes)) {

      nd <- nodes[[id]]
      kids <- nd$children
      if (!length(kids)) {
        next
      }
      kd <- nd$depth
      kys <-
        vapply(
          X = kids,
          FUN = function(k) { nodes[[k]]$y },
          FUN.VALUE = numeric(1)
        )
      
      graphics::segments(
        x0 = X(kd, min(kys)), y0 = Y(kd, min(kys)),
        x1 = X(kd, max(kys)), y1 = Y(kd, max(kys)),
        col = edge.col, lwd = edge.lwd
      )
      for (k in kids) {
        ky <- nodes[[k]]$y
        kdep <- nodes[[k]]$depth
        graphics::segments(
          x0 = X(kd, ky), y0 = Y(kd, ky),
          x1 = X(kdep, ky), y1 = Y(kdep, ky),
          col = edge.col, lwd = edge.lwd
        )
      }
    }

    ## Plot leader lines and node markers ----
    for (id in names(nodes)) {

      nd <- nodes[[id]]
      if (nd$is_leaf && align.leaves && nd$depth < maxdepth) {
        graphics::segments( # leader line
          x0 = X(nd$depth, nd$y), y0 = Y(nd$depth, nd$y),
          x1 = X(maxdepth, nd$y), y1 = Y(maxdepth, nd$y),
          col = edge.col, lwd = edge.lwd, lty = 3
        )
      }   
    }
    if (!is.na(point.pch)) {
      px <-
        vapply(
          X = names(nodes),
          FUN = function(id) { X(nodes[[id]]$depth, nodes[[id]]$y) },
          FUN.VALUE = numeric(1)
        )
      py <-
        vapply(
          X = names(nodes),
          FUN = function(id) { Y(nodes[[id]]$depth, nodes[[id]]$y) },
          FUN.VALUE = numeric(1)
        )
      graphics::points(
        x = px, y = py,
        pch = point.pch, cex = point.cex,
        bg = point.bg, col = edge.col
      )
    }

    ## Plot text labels per node ----
    for (id in names(nodes)) {
      nd <- nodes[[id]]
      if (nd$is_leaf) {
        ld <- if (align.leaves) { maxdepth } else { nd$depth }
        if (horizontal) {
          graphics::text(
            x = X(ld, nd$y), y = Y(ld, nd$y),
            labels = nd$name,
            pos = 4, offset = 0.3,
            cex = cex, col = leaf.col, xpd = NA
          )
        } else {
          graphics::text(
            x = X(ld, nd$y), y = Y(ld, nd$y),
            labels = nd$name,
            srt = 90, adj = c(1, 0.5),
            cex = cex, col = leaf.col, xpd = NA
          )
        } 
      } else {
        root_like <- is.na(nd$parent)
        if (horizontal) {
          pos <- if (root_like) { 4 } else { 2 }
          graphics::text(
            x = X(nd$depth, nd$y), y = Y(nd$depth, nd$y),
            labels = nd$name,
            pos = pos, offset = 0.3,
            cex = cex*node.cex, col = node.col, font = 2, xpd = NA
          )
        } else {
          graphics::text(
            x = X(nd$depth, nd$y), y = Y(nd$depth, nd$y),
            labels = nd$name,
            pos = if (root_like) { 4 } else { 3 },
            offset = 0.35, cex = cex * node.cex, col = node.col,
            font = 2, xpd = NA
          )
        }
      }
    }

    invisible(x)
  }

.flowFrame_to_matrix <-
  function(ff, need) {
    m <- flowCore::exprs(ff)
    nm <- colnames(m)
    desc <-
      tryCatch(
        as.vector(flowCore::pData(flowCore::parameters(ff))$desc),
        error = function(e) { NULL }
      )
    if (!is.null(desc)) {
      renamed <- nm
      has <- !is.na(desc)&nzchar(desc)
      renamed[has] <- desc[has]
      if (sum(need%in%renamed) > sum(need%in%nm)) {
        message("gate(): matched channels via flowFrame marker descriptions ($PnS).")
        colnames(m) <- renamed
      }
    }
    m
  }

#' Apply OMIQ gates to cytometry data
#'
#' Uses a gating tree exported from OMIQ (as an `.omiqgt` file) to gate an FCS 
#' file. If there are any adjustments to gates per individual file, the file's 
#' *OmiqID* must be specified. Otherwise, any potential adjustments will be 
#' ignored, and the default gate definitions will always be applied.
#' 
#' @param gt A `GatingTree` object created by `parse_omiqgt()` **OR** path to an
#'           `.omiqgt` file.
#' @param data Numeric matrix with named columns (expression matrix, columns
#'             named after fluorophores/channels/markers) **OR** a
#'             `flowCore::flowFrame` object **OR** path to an FCS file.
#' @param omiq_id Optional sample *OmiqID* if there are any per-file gate
#'                adjustments. Defaults to `NULL` (no adjustents).
#' @param sep Separator for the path column names. Defaults to "|".
#' @param ... Additional arguments to `flowCore::read.FCS()` if a path to FCS
#'            file was given.
#' @return A logical matrix: one row per cell, one column per path through the
#'   tree (root -> node), TRUE where the cell falls in that population. This is
#'   analogous to `flowWorkspace::gh_pop_get_indices()`.
#' @export
gate <- function(gt, data, omiq_id = NULL, sep = "|", ...) {

  is_single_string <- function(x) {
    is.atomic(x) && is.character(x)
  }
  is_ff <- function(x) {
    methods::is(x, "flowFrame") || inherits(x, "flowFrame")
  }

  ## Validate `gt` ----

  if (is_single_string(gt)) {
    stopifnot("Invalid path to .omiqgt file" = file.exists(gt))
    gt <- parse_omiqgt(gt)
  }
  m <- .model(gt)

  ## Validate `data` ----

  data_is_single_string <- is_single_string(data)

  if (data_is_single_string || is_ff(data)) {
    if (!requireNamespace("flowCore", quietly = TRUE)) {
      stop("Package 'flowCore' is required to gate a flowFrame.")
    }
    if (is_single_string(data)) {
      data <- flowCore::read.FCS(data, ...)
    }
    data <- .flowFrame_to_matrix(data, m$channels)
  }
    
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if (is.null(colnames(data))) {
    stop("`data` must have column names (channel/marker names).")
  }
  n <- nrow(data)

  ## Resolve file ID from OmiqID ----

  file_id <- NULL
  if (!is.null(omiq_id)) {
    file_id <- gsub("^F", "", omiq_id)
  }

  ## Extract gates ----

  gates <-
    if (is.null(file_id)) {
      m$gates
    } else {
      .build_gates(m, file_id = file_id)
    }

  ## Memoized memberships per node in isolation ----
  
  self_cache <- new.env(parent = emptyenv())
  self_mem <- function(id) {
    if (!is.null(self_cache[[id]])) {
      return(self_cache[[id]])
    }
    r <- gates[[m$nodes[[id]]$filterContainerId]](data)
    assign(id, r, envir = self_cache)
    r
  }

  ## Memoized cumulative memberships per node & its parents----

  cum_cache <- new.env(parent = emptyenv())
  cum_mem <- function(id) {
    if (!is.null(cum_cache[[id]])) {
      return(cum_cache[[id]])
    }
    p <- m$parent_of[[id]]
    r <- self_mem(id)
    if (nzchar(p) && p%in%names(m$nodes)) {
      r <- r & cum_mem(p)
    }
    assign(id, r, envir = cum_cache)
    r
  }

  ids <- m$order
  paths <- .node_paths(m, sep = sep)
  res <- vapply(X = ids, FUN = cum_mem, FUN.VALUE = logical(n))
  if (n == 1L) {
    res <- matrix(res, nrow = 1L)
  }
  dimnames(res) <- list(rownames(data), unname(paths[ids]))
  attr(res, "node_id") <- ids
  res
}

#' Extract all paths through a parsed OMIQ gating tree object
#'
#' Returns every full, unique gate path in a `GatingTree`. These are identical
#' to (and in the same order as) the column names of the matrix returned by
#' `omiqGTR::gate()`, provided that the same separator is used.
#'
#' @param gt A `GatingTree` object created by `parse_omiqgt()`.
#' @param sep Separator for the paths. Defaults to "|".
#' @return A character vector of full node paths.
#' @export
tree_paths <-
  function(gt, sep = "|") {

    m <- .model(gt)
    paths <- .node_paths(m, sep = sep) # named by node ID
    unname(paths[m$order]) # gate()-column order
  }

#' Restrict an OMIQ gating tree to one subtree
#'
#' Extracts the subtree whose base is the node at `path`. That node takes the
#' place of the synthetic "(ungated)" root: its own gate is dropped and its
#' children become the new top-level populations; everything outside the subtree
#' is removed.
#'
#' Any per-file gate adjustments (per-file filters) that belonged to a gate that
#' is no longer applied in the restricted tree cease to be relevant. Such
#' adjustments are dropped from the returned object, and a warning is emitted for
#' each affected node.
#'
#' @param gt A `GatingTree` object created by `parse_omiqgt()`.
#' @param sep Separator for the paths. Defaults to "|".
#' @param path Full path (as returned by `tree_paths()`) to the node that
#'             becomes the base of the extracted subtree.
#' @return A new `GatingTree` restricted to the requested subtree.
#' @export
isolate_subtree <-
  function(gt, sep = "|", path) {

    m <- .model(gt)

    ## Resolve `path` to a node ID ----

    all_paths <- .node_paths(m, sep = sep) # named by node ID
    base_id <- names(all_paths)[all_paths == path]
    if (length(base_id) == 0L) {
      stop(
        "No node matches path: ", path, "\n(see `tree_paths()` for valid paths)."
      )
    }
    if (length(base_id) > 1L) {
      stop("Path is ambiguous (", length(base_id), " nodes): ", path)
    }
    base_children <- m$children[[base_id]] # these become the new roots
    if (length(base_children) == 0L) {
      stop("Node '", path, "' is a leaf; it has no subtree to isolate.")
    }

    ## Collect nodes to keep ----

    descendants <- function(id) {
      c(id, unlist(lapply(m$children[[id]], descendants)))
    }
    keep_ids <- setdiff(descendants(base_id), base_id)

    ## Rebuild node relationships ----

    new_nodes <- m$nodes[keep_ids]
    new_children <- m$children[keep_ids]
    new_parent_of <- m$parent_of[keep_ids]
    for (r in base_children) { # detach the new roots
      new_nodes[[r]]$parentId <- ""
      new_parent_of[[r]] <- ""
    }

    ## Recursively collect used gate containers ----

    referenced <- new.env(parent = emptyenv())
    collect <- function(cid) {
      if (!is.null(referenced[[cid]])) {
        return(invisible())
      }
      assign(cid, TRUE, envir = referenced)
      fc <- m$containers[[cid]]
      if (identical(fc$containerType, "CompoundFilterContainer")) {
        for (child in unlist(fc$filterContainerIds)) {
          collect(child)
        }
      }
    }
    for (id in keep_ids) {
      collect(m$nodes[[id]]$filterContainerId)
    }
    ref_ids <- ls(referenced)

    ## Warn about vanished per-file adjustments ----

    has_pf <- function(fc) {
      !is.null(fc$perFileFilters) && length(fc$perFileFilters) > 0L
    }
    for (cid in names(m$containers)) {
      fc <- m$containers[[cid]]
      if (!has_pf(fc) || cid %in% ref_ids) {
        next
      }
      owners <- names(m$nodes)[
        vapply(
          X = m$nodes,
          FUN = function(nd) { identical(nd$filterContainerId, cid) },
          FUN.VALUE = logical(1)
        )
      ]
      owner_lbl <-
        if (length(owners)) {
          paste(
            vapply(
              X = owners,
              FUN = function(id) { .node_name(m, id) },
              FUN.VALUE = character(1)
            ),
            collapse = ", "
          )
        } else {
          cid
        }
      warning(
        sprintf(
          paste0(
            "`isolate_subtree()`: %d per-file gate adjustment(s) on node '%s' ",
            "(container %s) no longer apply in the isolated subtree."
          ),
          length(fc$perFileFilters), owner_lbl, cid
        ),
        call. = FALSE
      )
    }

    ## Assemble restricted model ----

    new_model <- list(
      "meta" = m$meta,
      "nodes" = new_nodes,
      "containers" = m$containers[ref_ids], # prune to referenced gates
      "roots" = base_children, # already ord-sorted
      "children" = new_children,
      "parent_of" = new_parent_of
    )
    new_model$nested <- .nested_tree(new_model)
    new_model$gates <- .build_gates(new_model)
    new_model$order <- .node_order(new_model)
    new_model$paths <- .node_paths(new_model, sep = sep)
    new_model$channels <- .gate_channels(new_model)

    counter <- new.env()
    counter$n <- 0L
    dnd <- .as_dendrogram_node(new_model$nested, counter)
    attr(dnd, "model") <- new_model
    class(dnd) <- c("GatingTree", "dendrogram")
    dnd
  }
