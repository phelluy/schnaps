#!/bin/sh

# Compile starpu article
# ----------------------
cd ../schnaps_opencl/
pdflatex -draftmode -interaction=nonstopmode schnaps_opencl2017.tex
bibtex schnaps_opencl2017.aux
pdflatex -draftmode -interaction=nonstopmode schnaps_opencl2017.tex
pdflatex -synctex=1 -interaction=nonstopmode schnaps_opencl2017.tex

# Compile thesis
# --------------
cd ../these_bruno/
if [ ! -d figures ]; then
  mkdir figures
fi

# Title page
pdflatex -interaction=nonstopmode 1ere.tex

# Final page
pdflatex -interaction=nonstopmode 4eme.tex

# Enable tikz externalize
sed -i 's|%\\usetikzlibrary{external}|\\usetikzlibrary{external}|g' header.tex
sed -i 's|%\\tikzexternalize\[|\\tikzexternalize\[|g' header.tex

pdflatex -draftmode -interaction=nonstopmode -shell-escape these-bwr.tex
bibtex these-bwr.aux
makeindex these-bwr.idx # if needed ?
pdflatex -draftmode -interaction=nonstopmode -shell-escape these-bwr.tex
pdflatex -synctex=1 -interaction=nonstopmode -shell-escape these-bwr.tex

# Disable tikz externalize
sed -i 's|\\usetikzlibrary{external}|%\\usetikzlibrary{external}|g' header.tex
sed -i 's|\\tikzexternalize\[|%\\tikzexternalize\[|g' header.tex

