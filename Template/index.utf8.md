---
#####################
## thesis metadata ##
#####################
title: |
  Thesis title 
author: Faes, E. 
college: Master in Finance
university: Antwerp Management School
university-logo: templates/beltcrest.pdf
submitted-text: A thesis submitted for the degree of
degree: Master in Finance
degreedate: June 2021
abstract: |
  The greatest abstract all times
acknowledgements: |
  This is where you will normally thank your advisor, colleagues, family and friends, as well as funding and institutional support.  
  Thanks. Thanks
  
  We must be grateful to John Gruber for inventing the original version of Markdown, to John MacFarlane for creating Pandoc which converts Markdown to a large number of output formats, and to Yihui Xie for creating "knitr" which introduced R Markdown as a way of embedding code in Markdown documents, and "bookdown" which added tools for technical and longer-form writing.
  
  Special thanks to [Chester Ismay](http://chester.rbind.io), who created the "thesisdown" package that helped many a PhD student write their theses in R Markdown. And a very special thanks to John McManigle, whose adaption of Sam Evans' adaptation of Keith Gillow's original maths template for writing an Oxford University DPhil thesis in LaTeX provided the template that Ulrik Lyngs in turn adapted for R Markdown.
  
  Finally, profuse thanks to JJ Allaire, the founder and CEO of [RStudio](http://rstudio.com), and Hadley Wickham, the mastermind of the tidyverse without whom we'd all just given up and done data science in Python instead. Thanks for making data science easier, more accessible, and more fun for us all.
  
  \begin{flushright}
  Enjo Faes \\
  Antwerp Management School, Antwerp \\
  27 June 2021
  \end{flushright}
dedication: For Yihui Xie
abbreviations: "front-and-back-matter/abbreviations" # path to .tex file with abbreviations

#######################
## bibliography path ##
#######################
bibliography: references.bib

########################
## PDF layout options ###
#########################
abstractseparate: false  # include front page w/ abstract for examination schools?
bib-humanities: true   #set to true if you want in-text references formatted as author-year
doi-in-bibliography: false #set to true if you want DOI's to be shown in the bibliography
draft: false # add as DRAFT mark in the footer?
page-layout: nobind #'nobind' for PDF output (equal margins), 'twoside' for two-sided binding (mirror margins and blank pages), leave blank for one-sided binding (left margin > right margin)
hidelinks: true #if false, the PDF output highlights clickable links with a colored border - you will probably want to set this to true for PDF version you wish to physically print
toc-depth: 0 # depth of heading to include in table of contents
lof: false # list of figures in front matter?
lot: false # list of tables in front matter?
mini-toc: false  # mini-table of contents at start of each chapter? (this just prepares it; you must also add \minitoc after the chapter titles)
mini-lot: false  # mini-list of tables by start of each chapter?
mini-lof: false  # mini-list of figures by start of each chapter?

params:
  corrections: true # set false to stop applying blue background to blocks of corrections

#####################
## output options  ##
#####################
output:
  bookdown::pdf_book:
    template: templates/template.tex
    keep_tex: true
    citation_package: biblatex
    #pandoc_args: ["--lua-filter=scripts_and_filters/correction_filter.lua"] #remove filter to stop applying blue background to inline corrections
  bookdown::gitbook:
    css: templates/style.css
    config:
      sharing:
        facebook: false
        twitter: yes
        all: false
  bookdown::word_document2:
    toc: true   
link-citations: true
documentclass: book

# csl: /Users/Benjamin/Zotero/styles/ecology.csl

---




<!--
Include the create_chunk_options chunk above at the top of your index.Rmd file
This will include code to create additional chunk options (e.g. for adding author references to savequotes)
and to make sure lines in code soft wrap
If you need to create your own additional chunk options, edit the file scripts/create_chunk_options.R
-->

<!-- This chunk includes the front page content in HTML output -->

