# Build Comprehensive TMB documentation

## Overview

| Files         | Description                                              |
|---------------|----------------------------------------------------------|
|`*.Rmd`        | Human generated documentation. Auto-converted to `*.md`. |
|`Doxyfile`     | Controls the doxygen specific documentation.             |
|`mainpage.Rmd` | Controls the doxygen mainpage.                           |

- Add a new vignette by dumping a `.Rmd` file in this folder.
- Vignettes show up under 'related pages' in doxygen.

## Requirements

#### Command line tools

```shell
sudo apt-get install doxygen pandoc pandoc-citeproc graphviz librsvg2-dev
```

#### R packages

```R
install.packages(c("knitr", "bookdown", "rsvg"))
```

## Export options

| Command              | Output description                 | Output location       |
|----------------------|------------------------------------|-----------------------|
| `make dox`           | Reduced doxygen (easier to read)   | `html/index.html`     |
| `make dox-full`      | Full doxygen                       | `html/index.html`     |
| `make gitbook`       | Gitbook                            | `_book/Tutorial.html` |
| `make inject_links`  | Insert doxygen links in previous   | `_book/Tutorial.html` |
| `make pdf-book`      | Pdf book                           | `_book/_main.pdf`     |
| `make refman`        | Full reference manual              | `latex/refman.pdf`    |
