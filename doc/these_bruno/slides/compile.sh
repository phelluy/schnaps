#!/bin/sh

pdflatex -draftmode -interaction=nonstopmode presentation_1.tex
# bibtex presentation_1.aux
# makeindex presentation_1.idx # if needed ?
# pdflatex -draftmode -interaction=nonstopmode presentation_1.tex
pdflatex -synctex=1 -interaction=nonstopmode presentation_1.tex
