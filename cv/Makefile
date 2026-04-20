## Build the CV.
##
## Requires XeLaTeX (TeX Live, MacTeX, etc.) and the EB Garamond + Songti SC
## fonts. Two passes are needed for the "page X / Y" footer to settle.

TEX  := xelatex
SRC  := cv.tex
PDF  := cv.pdf
AUX  := cv.aux cv.log cv.out cv.toc

.PHONY: all clean realclean

all: $(PDF)

$(PDF): $(SRC)
	$(TEX) -interaction=nonstopmode -halt-on-error $(SRC)
	$(TEX) -interaction=nonstopmode -halt-on-error $(SRC)

clean:
	rm -f $(AUX)

realclean: clean
	rm -f $(PDF)
