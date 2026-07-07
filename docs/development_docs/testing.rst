########
Testing and coverage
########

DELi uses `pytest <https://docs.pytest.org/>`_ for automated tests and `coverage.py <https://coverage.readthedocs.io/>`_ (via `pytest-cov <https://pytest-cov.readthedocs.io/>`_) to measure line coverage of the ``deli`` package.

Run the full suite from a development install:

.. code-block:: shell

   uv sync --dev --extra ml
   uv run pytest

GNN tests are marked ``ml`` and require the optional ML dependencies above. To skip them:

.. code-block:: shell

   uv run pytest -m "not ml"

Generate a coverage report
==========================

.. code-block:: shell

   uv run pytest --cov=deli --cov-report=term-missing --cov-report=html

This runs all tests, measures which lines in ``src/deli`` were executed, and writes:

- a terminal summary
- an HTML report in ``htmlcov/index.html``

Coverage counts **statements** (executable lines) in ``src/deli``. It does not include tests, docs, or example scripts.

Current coverage (DELi 0.2.1)
=============================

Figures below were generated on **2026-07-07** with **199 passing tests** (Python 3.13+, including ML/GNN tests with ``deli-chem[ml]`` installed).

.. list-table:: Total coverage
   :header-rows: 1
   :widths: 30 20 20 30

   * - Metric
     - Covered
     - Total
     - Percent
   * - All ``deli`` modules
     - 4,171
     - 5,385
     - **77.5%**

.. list-table:: Coverage by module
   :header-rows: 1
   :widths: 35 15 15 35

   * - Module
     - Covered
     - Statements
     - Coverage
   * - ``deli.decode``
     - 1,137
     - 1,384
     - 82.2%
   * - ``deli.dels``
     - 859
     - 1,008
     - 85.2%
   * - ``deli.enumeration``
     - 330
     - 409
     - 80.7%
   * - ``deli.analysis``
     - 978
     - 1,329
     - 73.6%
   * - ``deli.dna``
     - 96
     - 99
     - 97.0%
   * - ``deli`` (package root: CLI, configure, etc.)
     - 718
     - 1,044
     - 68.8%
   * - ``deli.utils``
     - 53
     - 73
     - 72.6%
   * - ``deli.design``
     - 0
     - 39
     - 0.0%

.. list-table:: ``deli.analysis`` submodules
   :header-rows: 1
   :widths: 45 15 15 25

   * - File
     - Covered
     - Statements
     - Coverage
   * - ``analysis/cube_class.py``
     - 569
     - 797
     - 71.4%
   * - ``analysis/poly_o.py``
     - 95
     - 99
     - 96.0%
   * - ``analysis/gnn.py``
     - 64
     - 73
     - 87.7%
   * - ``analysis/analysis_parse.py``
     - 125
     - 166
     - 75.3%
   * - ``analysis/analysis_report_gen.py``
     - 62
     - 85
     - 72.9%
   * - ``analysis/analyzers/gnn_analyzer.py``
     - 63
     - 109
     - 57.8%

Interpreting these numbers
==========================

- **Higher coverage** in ``decode`` and ``dels`` reflects unit tests around library definitions, barcodes, and decoding.
- **Analysis coverage** exercises cube enrichment, disynthon plots, PolyO, Random Forest baselines, and GNN helpers; the GNN training loop is partially covered via mocked tests.
- **`deli.design``** is not yet covered by automated tests.

In publications, we cite the total coverage percentage, test count, Python version, and DELi version used when the figures were generated.
