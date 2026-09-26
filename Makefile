# Build locally; never uploads or publishes.
PYTHON ?= python3

.PHONY: all check
all:
	$(MAKE) -C cv
	$(PYTHON) scripts/build_site.py

check:
	$(PYTHON) scripts/build_site.py --check
