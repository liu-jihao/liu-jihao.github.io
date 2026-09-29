## Local source -> CV build. No network or publishing steps.
## Requires Python 3 and XeLaTeX with EB Garamond, Fandol, xurl and bookmark.
## Songti SC remains the default when installed; Fandol Song is the fallback.
## Set BUILD_DIR to keep all compilation intermediates outside the source tree.

PYTHON ?= python3
XELATEX ?= xelatex
BUILD_DIR ?= .build
SRC := cv.tex
PDF := cv.pdf
DATA := ../data/content.json
GENERATOR := ../scripts/build_content.py

# Make timestamps and trailer identifiers repeatable with the same TeX/fonts.
# The data's review date is UTC; an explicit environment value takes precedence.
SOURCE_DATE_EPOCH ?= $(shell $(PYTHON) -c 'import datetime,json; d=json.load(open("$(DATA)"))["updated"]; print(int(datetime.datetime.fromisoformat(d).replace(tzinfo=datetime.timezone.utc).timestamp()))')
export SOURCE_DATE_EPOCH
export FORCE_SOURCE_DATE := 1

.PHONY: all check clean realclean

all: $(PDF)

$(SRC): $(DATA) template.tex $(GENERATOR)
	$(PYTHON) $(GENERATOR) --cv-only

$(PDF): $(BUILD_DIR)/cv.pdf
	cp "$<" "$@"

$(BUILD_DIR)/cv.pdf: $(SRC) Makefile
	mkdir -p "$(BUILD_DIR)"
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory="$(BUILD_DIR)" $(SRC)
	$(XELATEX) -interaction=nonstopmode -halt-on-error -file-line-error -output-directory="$(BUILD_DIR)" $(SRC)

check:
	$(PYTHON) $(GENERATOR) --cv-only --check

clean:
	rm -f "$(BUILD_DIR)"/cv.aux "$(BUILD_DIR)"/cv.log "$(BUILD_DIR)"/cv.out "$(BUILD_DIR)"/cv.toc "$(BUILD_DIR)"/cv.pdf

realclean: clean
	rm -f $(PDF)
