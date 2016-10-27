## Variable '.doxygen'
##   FALSE: Use pandoc to convert md -> html  (default)
##   TRUE : Use doxygen to convert md -> html
if(!exists(".doxygen")) {
    .doxygen <- FALSE
}
.doxygen_relpath <- "../html"

## Currently only works for files. We might consider setting
## Doxyfile->GENERATE_TAGFILE=tagfile.xml and use the links in this
## xml file.
doxylink <- function(x) {
    if(.doxygen) {
        ans <- paste("\\ref", x)
    } else {
        is_file <- (x != sub(".cpp$", "", x))
        ans <- ifelse(is_file,
                      paste0(
                          "[", x, "]",
                          "(", .doxygen_relpath, "/",
                          sub(".cpp$", "_8cpp-example.html", x),
                          ")" ),
                      paste0("`",x,"`") )
    }
    ans
}

## Include inline source code
include_source <- function(x, indent=0) {
    indent <- paste(rep(" ", indent), collapse="")
    type <- tolower(sub(".*\\.", "", x))
    ans <- c(paste0("```", type),
             readLines(x),
             "```")
    ans <- paste0(indent, ans)
    paste0("\n", paste(ans, collapse="\n"), "\n")
}

## Tweaks for doxygen's markdown:
doxy_markdown_tweaks <- function(file) {
    x <- readLines(file)
    ## 1: Doxygen wants code tags of the form {.cpp} (not just .cpp)
    x <- gsub("```cpp", "```{cpp}", x)
    ## 2: Doxygen does not like the 4-space rule
    x <- gsub("    ```", "  ```", x)
    ## 3: TODO: Remove empty lines before code in lists
    ## 4: Formula workaround (inline):
    x <- gsub("^\\$([^$]+?)\\$", "\\\\f$\\1\\\\f$", x)
    x <- gsub(" \\$([^$]+?)\\$", " \\\\f$\\1\\\\f$", x)
    ## Save
    writeLines(x, file)
}

## plot graph of hessian as either
## * Undirected conditinal independence graph
## * Directed graph of successive conditional distributions
plotGraph <- function(h, DAG=TRUE, group = function(x)(x-1)%%nrow,
                      nrow=1, color="black", fillcolor=NULL,
                      shape=NULL, debug=FALSE, ... ) {
    if(is.null(labels <- rownames(h)))
        if(is.null(labels <- colnames(h)))
            labels <- paste0("X", 1:nrow(h))
    rownames(h) <- colnames(h) <- NULL
    sub <- function(x){
        pattern <- as.character(0:9)
        replacement <- sapply(8320 + 0:9, intToUtf8)
        for(i in seq_along(pattern))
            x <- gsub(pattern[i], replacement[i], x)
        x
    }
    labels <- sapply(labels, sub)
    ## Draw nodes one-by-one in natural order
    dag_edges <- function(h) {
        h <- as(h, "dsCMatrix")
        h@x[] <- 0
        diag(h) <- 1
        p <- nrow(h):1
        hrev <- h[p, p]
        Lrev <- Matrix::Cholesky(hrev, super=FALSE, perm=FALSE,
                                 LDL=FALSE)
        Lrev <- as(Lrev, "sparseMatrix")
        U <- Lrev[p, p]
        U@x[] <- 1
        diag(U) <- 0
        e <- which(as.matrix(U)>0, arr=TRUE)
        e
    }
    ## Undirected case
    edges <- function(h) {
        h <- as(h, "dsCMatrix")
        h@x[] <- 1
        diag(h) <- 0
        e <- which(triu(h)>0, arr=TRUE)
        e
    }
    if(DAG)
        e <- dag_edges(h)
    else
        e <- edges(h)
    if(!is.null(group))
        spl <- split(1:nrow(h), group(1:nrow(h)))
    graph <- c(
        "digraph DAG {",
        paste0(1:nrow(h), "[label=\"",labels,"\"",
               if(!is.null(color))     paste0(" color=",     color)     else NULL,
               if(!is.null(fillcolor)) paste0(" fillcolor=", fillcolor, " style=filled") else NULL,
               if(!is.null(shape))     paste0(" shape=",     shape)     else NULL,
               "]"),
        if(!is.null(group))
            sapply(spl, function(x) paste("{rank=same",  paste(x, collapse=", ") , "}") )
        else NULL,
        paste(e[,1] , "->", e[,2]),
        "}")
    if(!DAG) {
        graph[1] <- "graph G {"
        graph <- gsub("->","--",graph)
    }
    graph <- paste(graph, collapse="\n")
    if(debug)return(graph)
    DiagrammeR::grViz(graph, ...)
}
