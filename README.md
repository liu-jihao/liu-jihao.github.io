# Jihao Liu — Personal Academic Website

Source for [https://liu-jihao.github.io](https://liu-jihao.github.io).

## Stack

- Pure static HTML + CSS (no build step, no JavaScript dependencies)
- Hosted on GitHub Pages
- Custom academic design system in `assets/style.css`

## Site Structure

```
.
├── index.html              # Home / About
├── publications.html       # Publications and preprints (48 papers)
├── teaching.html           # Teaching experience
├── talks.html              # Invited talks (conferences + seminars)
├── conferences.html        # Conferences and seminars co-organized
├── collaborators.html      # Mathematical collaborators
├── pku-ag-seminar.html     # PKU Algebraic Geometry Seminar (current schedule)
├── ffm-conference.html     # Conference on Foliation, Fibration, and Moduli (2026)
├── cv/                     # CV source + compiled PDF (cv.tex, cv.pdf, Makefile)
└── assets/
    └── style.css           # All site styles
```

## Local Preview

To preview the site locally, you can use any static file server. The simplest options:

```bash
# Using Python (built-in on macOS)
python3 -m http.server 8000

# Or using Node
npx serve .
```

Then open [http://localhost:8000](http://localhost:8000) in your browser.

## How to Edit

All content lives in plain HTML files. To update:

1. Open the relevant `.html` file in any text editor (or in Cursor)
2. Edit the content between the `<main>` tags
3. Save the file
4. Commit and push to GitHub — the site updates automatically within ~1 minute

To change the site-wide style, edit `assets/style.css`.

## Common Updates

### Adding a new publication

Open `publications.html`, find the year heading, and add a new entry block:

```html
<div class="entry">
  <div class="entry-title"><span class="entry-number">49.</span> Title of the paper</div>
  <div class="entry-authors">Author 1, Author 2</div>
  <div class="entry-venue">Journal name, year, pages.</div>
  <div class="entry-links">
    <a href="https://arxiv.org/abs/XXXX.XXXXX">arXiv:XXXX.XXXXX</a>
    <a href="https://doi.org/...">DOI</a>
  </div>
</div>
```

### Adding a seminar talk

Open `talks.html`, find the right section (Conferences or Seminars), and add an `<li>`:

```html
<li>
  <div class="course-name"><span class="entry-number">N.</span> Talk title</div>
  <div class="course-meta">Date · Venue · <a href="...">link</a></div>
</li>
```

## Regenerating favicons / OG image

The favicons (`assets/favicon.{svg,ico}`, `assets/favicon-{32,180,512}.png`)
and the social sharing card (`assets/og-image.png`) are generated from a
single Python script. To rebuild them after changing the design:

```bash
python3 -m venv .venv
.venv/bin/pip install Pillow
.venv/bin/python3 scripts/make_favicon.py
```

## Updating site-wide head tags

The `<!-- META:START --> ... <!-- META:END -->` block in each HTML file
(favicon links, canonical URL, Open Graph + Twitter card meta tags) is
managed by `scripts/inject_head_tags.py`. To add a new page or change
the tags everywhere at once:

1. Edit the `PAGES` list (or the block template) in `scripts/inject_head_tags.py`
2. Run `.venv/bin/python3 scripts/inject_head_tags.py`

The script is idempotent — it replaces an existing block in place if found,
and inserts after the `<meta name="description">` tag otherwise.

## Rebuilding the CV

The CV (`cv/cv.pdf`) is compiled from `cv/cv.tex` with XeLaTeX. After editing
the `.tex`, regenerate the PDF with:

```bash
cd cv && make
```

Two passes are run automatically (the second pass settles the
"page X / Y" footer). The Makefile cleans up `.aux/.log/.out` afterwards
via `make clean`. EB Garamond and Songti SC are required (both ship with
TeX Live / macOS by default).

## Deployment

The site auto-deploys to GitHub Pages from the `main` branch. To deploy:

```bash
git add .
git commit -m "Update content"
git push
```

GitHub Pages serves the site at `https://liu-jihao.github.io` within ~1 minute of pushing.

## License

Site content © Jihao Liu. Source code (HTML/CSS templates) free to reuse.
