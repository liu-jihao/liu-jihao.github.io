# Jihao Liu — Personal Academic Website

Source for [jihaoliu.org](https://jihaoliu.org/), hosted on GitHub Pages from the `main` branch.

## Site and local preview

The site uses static HTML and CSS, with a local JavaScript enhancement for publication filtering and site search. Crimson Pro is hosted in `assets/fonts/`, with its SIL Open Font License. All publications remain readable without JavaScript; the search page provides a page directory when JavaScript is disabled.

From the repository root:

```sh
python3 -m http.server 8000 --bind 127.0.0.1
```

Open [localhost:8000](http://localhost:8000/).

## Edit and rebuild

Shared content is maintained in:

- `data/content.json`: profile, publications, research topics, teaching, talks, organized events, and CV content.
- `data/navigation.json`: shared primary and secondary navigation.
- `data/seminar-current.json`: current seminar dates, times, venues, and source URLs.
- `cv/template.tex`: CV layout.
- `templates/search.html`, `assets/discovery.js`: site search and publication filters.
- `assets/style.css`: shared website layout and typography.

The generated `assets/search-index.js` contains public content only. Searches run locally in the visitor's browser.

```sh
make
make check
```

Python generation uses only the standard library. Building the CV requires XeLaTeX with EB Garamond, xurl and bookmark. Songti SC is used when installed; Fandol Song is the fallback. `make` builds the CV before the website and search index. `make check` verifies that generated content matches its sources and checks local links and anchors.

To place CV compilation intermediates in a separate directory:

```sh
make BUILD_DIR=/absolute/path/to/cv-build
```

`cv/cv.tex` is generated but remains a self-contained editable LaTeX file. Direct edits to generated pages or CV source will be overwritten by the next build. Update the shared data or templates instead. The explicit `updated` review date and seminar records determine upcoming/past seminar placement.

Notes, AI-result papers, collaborators, conference details and seminar archives retain their established content and files. Specialized content remains in the corresponding HTML unless represented in shared data. Preserve the author's accounts of priority, sharing dates, external verification and publication decisions, including their attribution and qualifications. Neither a private source nor a failed independent search is by itself grounds to omit or downgrade these accounts. Preserve original talk hyperlinks and the separate Northwestern Winter 2024 course/section records.

## Other site assets

`scripts/make_favicon.py` generates favicon and social sharing images; it requires Pillow. `scripts/inject_head_tags.py` maintains the shared metadata blocks. Notes and AI-result entries retain their expandable details and visible status lines.

## Publishing

Review the generated files and `make check` results before committing. Pushing the reviewed commit to `main` triggers the existing GitHub Pages deployment. `CNAME` retains `jihaoliu.org`.

`_config.yml` excludes source data, scripts, templates, LaTeX inputs, build files and this documentation from the deployed website. These files remain versioned in the repository for reproducibility. Only the compiled `cv/cv.pdf` is served as the CV.

## License

Site content © Jihao Liu. Source code (HTML/CSS templates) is free to reuse. Crimson Pro is distributed under its included SIL Open Font License.
