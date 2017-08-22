library(TMB)
formals(MakeADFun)$silent <- TRUE

## knitr global options
knitr::opts_chunk$set(fig.width=3.5, fig.height=3.5) ## half of default (both for html and pdf output)
knitr::opts_chunk$set(fig.align="center")
.latex <- identical( knitr::opts_knit$get('rmarkdown.pandoc.to') , "latex" ) ## Are we converting to latex ?
if( ! .latex ) {
    knitr::opts_chunk$set(dev="svg") ## Use svg rather than png for html output
}

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
    if(.latex) { ## Latex
        ans <- x
    }
    else if(.doxygen) { ## Doxygen
        ans <- paste("\\ref", x)
    } else { ## HTML
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

## Create new source example
create_source <- function(code, file) {
    if(!file.exists(file)) {
        writeLines(code, file)
    }
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
    ## 4: Formula workaround:
    x <- paste(x, collapse="\n ")
    x <- gsub("\\$\\$([^$]+?)\\$\\$", "\\\\f[\\1\\\\f]", x) ## centered
    x <- gsub(" \\$([^$#]+?)\\$", " \\\\f$\\1\\\\f$", x)    ## inline
    x <- strsplit(x, "\n ")[[1]]
    ## Save
    writeLines(x, file)
}

## Output graph as svg
graph2svg <- function(graph) {
    dotfile <- tempfile(fileext=".dot")
    svgfile <- tempfile(fileext=".svg")
    cat(graph, file=dotfile)
    cmd <- paste("dot", "-Tsvg", dotfile, "-o", svgfile)
    system(cmd)
    class(svgfile) <- "svgfile"
    svgfile
}

## plot graph of hessian as either
## * Undirected conditinal independence graph
## * Directed graph of successive conditional distributions
plotGraph <- function(h, DAG=TRUE, fill_in=DAG, group = function(x)(x-1)%%nrow,
                      nrow=1, color="black", fillcolor=NULL,
                      shape="circle", debug=FALSE, ... ) {
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
    fill_in_edges <- function(h) {
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
    if(fill_in)
        e <- fill_in_edges(h)
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
               " fixedsize=true",
               "]"),
        if(!is.null(group))
            sapply(spl, function(x) paste("{rank=same",  paste(x, collapse=" ") , "}") )
        else NULL,
        paste(e[,1] , "->", e[,2]),
        "}")
    if(!DAG) {
        graph[1] <- "graph G {"
        graph <- gsub("->","--",graph)
    }
    graph <- paste(graph, collapse="\n")
    if(debug)return(graph)
    graph2svg(graph)
}


############################################################################
## Computational graphs
## From https://rawgit.com/kaskr/tmb_dumpstack/master/makeDot.R
## with minor changes
dumpstack <- function(obj, type = c("value", "gradient"), ...) {
    type <- c("value"="ADdouble", "gradient"="ADGrad")[match.arg(type)]
    threads.old <- openmp()
    on.exit(openmp(threads.old))
    openmp(1)
    ans <- capture.output( out <- obj$env$f(dumpstack=TRUE, type=type) )
    attr(ans, "ninput") <- length(obj$env$par)
    attr(ans, "noutput") <- length(out)
    ans
}
makeDot <- function(dump,
                    file = as.character(substitute(dump)),
                    bold=NULL,
                    filled=NULL,
                    dotfile=paste0(file,".dot"),
                    svgfile=paste0(file,".svg"),
                    forward=TRUE,
                    rank_operators=FALSE,
                    remap_atomic=TRUE) {
    d <- dump
    d <- gsub("CExp ", "CExp", d)
    d <- d[d != ""]
    d <- sub("^o=[^ ]*[ ]*","",d)
    d <- gsub(" +"," ",d)
    d <- strsplit(d," ")
    node <- sub("v=","",sapply(d,function(x)x[1]))
    op <- sub("op=","",sapply(d,function(x)x[2]))
    ## Remove 'End' operator
    keep <- (op != "End")
    d <- d[keep]
    node <- node[keep]
    op <- op[keep]
    ## Some operators do not produce a variable
    remap <- as.character(seq_len(length(node)))
    names(remap) <- node
    node <- as.character(seq_len(length(node)))
    ## Atomic functions
    if(remap_atomic) {
        usr <- as.logical(cumsum(op=="User") %% 2)
        usr <- xor(usr, op=="User")
        newremap <- cumsum(!usr)
        remap[] <- newremap
        node[]  <- as.character(newremap)
    }
    ## Calc dependencies (and account for remapping)
    parseDep <- function(x){
        ans <- sub(".*=","",grep("v.=|.v=|v=",x,value=TRUE))
        ans <- ans[-1]
        ans <- remap[ans]
        ans
    }
    args <- lapply(d,parseDep)
    if(!forward){
        edges <- unlist(Map(paste,node,args,sep=" -> ")[sapply(args,length)>0])
    } else {
        edges <- unlist(Map(paste,args,node,sep=" -> ")[sapply(args,length)>0])
    }
    lab <- paste(op,node)
    labels <- paste0(node," [label=\"",lab,"\"]")
    input <- head(node, attr(dump, "ninput"))
    output <- tail(node, attr(dump, "noutput"))
    if(!is.null(bold))bold <- paste(eval(bold),"[style=dashed]")
    if(!is.null(filled))filled <- paste(eval(filled),"[style=filled]")
    group <- function(x,rank="same")c("{",paste0("rank=",rank),x,"}")
    spl <- split(node[op!="End"],op[op!="End"])
    graph <- c(
        "digraph graphname {",
        if(rank_operators)unlist(lapply(spl,group)) else NULL,
        edges,
        labels,
        filled,
        bold,
        "}")
    graph <- paste(graph, collapse="\n")
    graph2svg(graph)
}
splitDump <- function(dump) {
    spl <- split(dump, cumsum(dump == ""))
    spl <- spl[sapply(spl, length) > 1]
    names(spl) <- NULL
    spl <- lapply(spl, function(x){attributes(x) <- attributes(dump);x})
    spl
}
############################################################################



print.svgfile <- function(x) {
    browseURL(x)
}

knit_print.svgfile <- function(x,...,options=NULL) {
    filename <- knitr::fig_chunk(options$label, ext=options$fig.ext)
    ## print(options)
    ## print(filename)
    ## print(knitr::opts_knit$get('rmarkdown.pandoc.to'))
    if(!file.exists(options$fig.path))
        dir.create(options$fig.path, recursive=TRUE)
    if(options$dev == "pdf") {
        rsvg::rsvg_pdf(x, filename)
    } else if(options$dev == "png") {
        rsvg::rsvg_png(x, filename)
    } else if(options$dev == "svg") {
        file.copy(x, filename)
    }
    else stop("Device must be 'png', 'svg' or 'pdf'")
    knitr::include_graphics(filename)
}
