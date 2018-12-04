#!/bin/sh

pdflatex -draftmode -interaction=nonstopmode -shell-escape these.tex
bibtex these.aux
pdflatex -draftmode -interaction=nonstopmode -shell-escape these.tex
pdflatex -synctex=1 -interaction=nonstopmode -shell-escape these.tex

#Clean folder
rm ./*.aux
rm ./*.mtc*
rm ./*.toc
rm ./*.maf
rm ./*.out
rm ./*.blg
rm ./*.bbl
rm ./*.gz
rm ./*.log