VERSION := $(shell sed -n '/^Version: /s///p' ../TMB/DESCRIPTION)
vignettes := $(basename $(wildcard *.Rmd))
md_vignettes := $(vignettes:=.md)

## Preparations for Doxygen
$(shell cp Doxyfile Doxyfile.tmp)
$(shell echo "\n\n# Customized settings:" >> Doxyfile.tmp)
$(shell echo "PROJECT_NUMBER = v$(VERSION)" >> Doxyfile.tmp)

## Preparations for bookdown. _bookdown.yml seems to be the only way to limit the book to a subset of Rmd files ?
$(shell Rscript -e "cat(paste('rmd_files : [',paste(shQuote(dir(pattern='[0-9].*.Rmd')), collapse=','),']'), file='_bookdown.yml')")
$(shell echo "\ndownload : ['pdf']" >> _bookdown.yml)

%.md : %.Rmd
	echo ".doxygen<-TRUE; source('layout/utils.R'); knitr::knit('$(basename $<).Rmd'); doxy_markdown_tweaks('$(basename $<).md')" | R --slave

## Outcomment to get new layout
dox: $(md_vignettes)
	## echo "HTML_HEADER = layout/header.html" >> Doxyfile.tmp
	## echo "HTML_EXTRA_STYLESHEET = layout/customdoxygen.css" >> Doxyfile.tmp
	## echo "HTML_EXTRA_FILES = layout/doxy-boot.js" >> Doxyfile.tmp
	echo "USE_MATHJAX = YES" >> Doxyfile.tmp
	echo "MATHJAX_EXTENSIONS = TeX/AMSmath TeX/AMSsymbols" >> Doxyfile.tmp
	echo "MATHJAX_RELPATH = https://cdn.mathjax.org/mathjax/latest" >> Doxyfile.tmp
	echo "EXTRACT_ALL = NO" >> Doxyfile.tmp
	echo "HIDE_UNDOC_MEMBERS = YES" >> Doxyfile.tmp
	echo "HIDE_UNDOC_CLASSES = YES" >> Doxyfile.tmp
	echo "INTERNAL_DOCS = NO" >> Doxyfile.tmp
	echo "SORT_BRIEF_DOCS = YES" >> Doxyfile.tmp
	doxygen Doxyfile.tmp
	## sed -i 's/$(document).ready/$(document).skip_this_function/g' html/dynsections.js
	cp -rf figure html
	rm Doxyfile.tmp

dox-full: $(md_vignettes)
	doxygen Doxyfile.tmp
	cp -rf figure html
	rm Doxyfile.tmp

gitbook:
	echo ".doxygen<-FALSE; source('layout/utils.R'); bookdown::render_book(dir(pattern='[0-9].*.Rmd'), 'bookdown::gitbook', clean_envir=FALSE)" | R --slave

inject_links:
	R --slave < layout/inject_links.R

pdf_book:
	echo ".doxygen<-FALSE; source('layout/utils.R'); bookdown::render_book(dir(pattern='[0-9].*.Rmd'), 'bookdown::pdf_book', clean_envir=FALSE)" | R --slave

refman: dox
	cd latex; make

## For quick preview
$(vignettes) : % : %.Rmd
	echo ".doxygen<-FALSE; source('layout/utils.R'); rmarkdown::render('$(basename $<).Rmd', rmarkdown::html_document());" | R --slave
	touch $(basename $<)
	google-chrome $(basename $<).html &

## This is the version for gh-pages. We put book inside html folder and change relative links
all:
	make dox
	make gitbook
	make pdf_book
	make inject_links
	cd _book; sed -i 's|\.\./html/|\.\./|g' *.html
	rm -rf html/_book
	mv _book html
	touch html/.nojekyll
