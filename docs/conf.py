# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import importlib.metadata

project = "ucla_plha"
author = "Meera L. Kota, Scott J. Brandenberg"
copyright = "2026, Meera L. Kota and Scott J. Brandenberg"
try:
    release = importlib.metadata.version("ucla_plha")
except importlib.metadata.PackageNotFoundError:
    release = "2.0.3"
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",        # pull docstrings from the code
    "sphinx.ext.autosummary",    # generate per-function/class summary pages
    "sphinx.ext.napoleon",       # parse Google-style docstrings (Args:/Returns:)
    "sphinx.ext.viewcode",       # add [source] links
    "sphinx.ext.intersphinx",    # cross-link to numpy/scipy/etc. docs
    "sphinx.ext.mathjax",        # render LaTeX math
    "myst_parser",               # let Sphinx read Markdown (.md) files
    "sphinx_copybutton",         # copy button on code blocks
    "sphinx_design",             # cards / grids for the landing page
]

autosummary_generate = True
autodoc_typehints = "description"
autodoc_member_order = "bysource"
napoleon_google_docstring = True
napoleon_numpy_docstring = False

# Accept both Markdown and reStructuredText sources
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# MyST niceties (dollar-math, AMS math, definition lists, ::: fences)
myst_enable_extensions = ["dollarmath", "amsmath", "deflist", "colon_fence"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Render \( \) and \[ \] math delimiters in addition to $ and $$
mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"], ["$", "$"]],
        "displayMath": [["\\[", "\\]"], ["$$", "$$"]],
    }
}

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_book_theme"
html_title = f"ucla_plha {version}"
html_favicon = "_static/logo.svg"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "repository_url": "https://github.com/meerakota/ucla_plha",
    "use_repository_button": True,
    "use_download_button": True,
    "show_nav_level": 2,
    "logo": {
        "image_light": "_static/logo.svg",
        "image_dark": "_static/logo-dark.svg",
        "text": "ucla_plha",
    },
}
