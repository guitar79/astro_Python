#!/bin/sh
#build book html
jupyter-book clean ./ --all
#jupyter-book build ./

# build book pdf
#https://jupyterbook.org/en/stable/advanced/pdf.html
# jupyter-book build ./ --builder html
jupyter-book build ./ --builder pdfhtml
#jupyter-book build ./ --builder pdflatex 
#jupyter-book build ./ --builder latex