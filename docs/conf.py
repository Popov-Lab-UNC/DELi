"""configurations for sphinx documentation"""

from pathlib import Path

# Try to read the project version from pyproject.toml (PEP 621). Use the
# stdlib `tomllib` when available (Python 3.11+). If not available, fall
# back to the external `toml` package if present. If both fail, fall back
# to `importlib.metadata` which requires the package to be installed.
try:
    import tomllib as _toml
except Exception:  # pragma: no cover - optional dependency fallback
    try:
        import toml as _toml  # type: ignore
    except Exception:
        _toml = None  # type: ignore

_pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
version = None
if _toml is not None and _pyproject_path.exists():
    try:
        with _pyproject_path.open("rb") as _f:
            _data = _toml.load(_f)
        version = _data.get("project", {}).get("version")
    except Exception:
        version = None

if not version:
    # Final fallback to importlib.metadata (package must be installed).
    from importlib.metadata import version as _get_version

    try:
        # If the installed distribution name differs from the import name
        # (e.g. distribution `deli-chem` but package imports as `deli`), try
        # to resolve the distribution name that provides the `deli` top-level
        # package using `packages_distributions()` when available. Also try
        # common distribution name candidates.
        dist_candidates = []
        try:
            from importlib.metadata import packages_distributions

            dist_candidates = packages_distributions().get("deli", [])
        except Exception:
            dist_candidates = []

        dist_candidates.extend(["deli", "deli-chem"])

        version = None
        for _dist in dist_candidates:
            try:
                version = _get_version(_dist)
                break
            except Exception:
                continue

        if not version:
            version = "0.0.0"
    except Exception:
        version = "0.0.0"

release = version

# Project information
project = "DELi"
copyright = "2025, James Wellnitz"
author = "James Wellnitz"

# Add any Sphinx extension module names here, as strings
extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The theme to use for HTML and HTML Help pages.
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets)
html_static_path = ["_static"]

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Configure autodoc
autodoc_member_order = "bysource"
autoclass_content = "both"

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
}
