## Inject doxygen links in bookdown html version.
library(stringr)

dictionary <- function(folder = "html", what="code") {
    owd <- getwd()
    on.exit(setwd(owd))
    setwd(folder)
    ## 1. Find keywords and (html) tags
    pattern <- c('code' =     '<a class="code" [^ ]*[^;]</a>') [what]
    tag <- unlist(lapply(dir(pattern=".html"), readLines))
    tag <- unlist( str_extract_all(tag, pattern) )
    tag <- unique(tag)
    keyword <- gsub( ".*>", "", gsub("</a>", "", tag) )
    ## 2. Keep only useful ones
    discard <- keyword == "\\" | grepl("operator", keyword) | keyword=="dt"
    tag <- tag[!discard]
    keyword <- keyword[!discard]
    ## 3. Keep only unique
    keep <- !duplicated(keyword)
    tag <- tag[keep]
    keyword <- keyword[keep]
    ## 4. Modify path
    tag <- str_replace(tag, 'href="', 'href="../html/')
    ## 5. sort
    i <- rev(order(nchar(keyword)))
    list(keyword=keyword[i], tag=tag[i])
}

injectLinks <- function(text, dict) {
    keyword <- dict$keyword
    tag <- dict$tag
    ## 3. Create regexp that finds keyword in a html document
    ##    note: To ensure keyword is only substituted once regexp must
    ##    never match a tag !
    ##regexp <- paste0("([^a-zA-Z_0-9>])", keyword, "([^a-zA-Z_0-9<])")
    regexp <- paste0("\\b", keyword, "([^a-zA-Z_0-9<\\.])")
    replace <- paste0(tag, "\\1")
    names(replace) <- regexp
    text <- str_replace_all(text, replace)
    text
}

dict <- dictionary("html", "code")

injectLinksFile <- function(file) {
    doc <- readLines(file)
    doc <- paste(doc,collapse="\n")
    regexp <- "<code[^>]*>[\\s\\S]*?</code>"
    code <- str_extract_all(doc, regexp)[[1]]
    spl <- strsplit(doc, regexp, perl=TRUE)[[1]]
    code <- injectLinks(code, dict)
    doc2 <- paste0( spl[1], paste(paste0(code, spl[-1]), collapse="") )
    writeLines(doc2, file)
}

lapply(
    dir("_book", pattern=".html$", full.names=TRUE),
    injectLinksFile )
