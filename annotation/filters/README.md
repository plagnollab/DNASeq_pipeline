# Filters

Filters modified to work with bash pipe command.
File is read from stdin and output written to stdout, so can be piped to another filter.

Here is an example script of how to pipe filter:

```
chrCode=22

function filter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.025 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/GO-filter.R > af-csq-go-filtered-VEP_$chrCode.csv
    wc -l af-csq-go-filtered-VEP_$chrCode.csv
}

function tonyfilter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.05 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R --cadd.thresh 10 > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/expression-filter.R > af-csq-expression-filtered-VEP_$chrCode.csv
    wc -l af-csq-expression-filtered-VEP_$chrCode.csv
}


function stringent_filter() {
    DIR=~/bin/DNASeq_pipeline/annotation/filters/
    cat VEP_$chrCode.csv | Rscript $DIR/af-filter.R --ajcontrols.thresh 0.05 --uclex.thresh 0.05 --exac.thresh 0.025 > af-filtered-VEP_$chrCode.csv
    wc -l  af-filtered-VEP_$chrCode.csv
    cat af-filtered-VEP_$chrCode.csv |  Rscript $DIR/csq-filter.R > af-csq-filtered-VEP_$chrCode.csv
    wc -l af-csq-filtered-VEP_$chrCode.csv
    cat af-csq-filtered-VEP_$chrCode.csv | Rscript $DIR/GO-filter.R > af-csq-go-filtered-VEP_$chrCode.csv
    wc -l af-csq-go-filtered-VEP_$chrCode.csv
}

filter
#tonyfilter
#stringent_filter
```
