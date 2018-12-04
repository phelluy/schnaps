#!/bin/sh

pdflatex -draftmode -interaction=nonstopmode -shell-escape beamer.tex
# bibtex beamer.aux
# pdflatex -draftmode -interaction=nonstopmode -shell-escape beamer.tex
pdflatex -synctex=1 -interaction=nonstopmode -shell-escape beamer.tex

#Clean folder
rm ./*.aux
rm ./*.snm
rm ./*.nav
rm ./*.mtc*
rm ./*.toc
rm ./*.maf
rm ./*.out
rm ./*.blg
rm ./*.bbl
rm ./*.gz
rm ./*.log