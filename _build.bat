::# build book html
jupyter-book clean ./ --all
:: jupyter-book build ./

::# build book pdf
::#https://jupyterbook.org/en/stable/advanced/pdf.html
jupyter-book build ./ --builder pdfhtml
::jupyter-book build ./ --builder pdflatex :: 한글 깨짐
::jupyter-book build ./ --builder latex