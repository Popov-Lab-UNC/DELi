#!/usr/bin/env python3
"""Export the `docs` dependency group from pyproject.toml to a requirements-style file.

Usage:
    python scripts/export_docs_requirements.py \
        --pyproject pyproject.toml \
        --output docs/requirements-docs.txt

This script prefers the PEP 621 `dependency-groups.docs` table, falling back to
`project.optional-dependencies.docs` if present.
"""

import argparse
import pathlib
import sys


# Prefer stdlib tomllib (3.11+), fall back to tomli package for older runners
try:
    import tomllib  # Python 3.11+

    def loads_toml(s: str):
        return tomllib.loads(s)
except Exception:
    try:
        import tomli

        def loads_toml(s: str):
            return tomli.loads(s)
    except Exception:
        print("tomllib (Python 3.11+) or tomli package is required.", file=sys.stderr)
        raise


def main():
    """run the script"""
    parser = argparse.ArgumentParser(description="Export docs dependency group to a requirements file.")
    parser.add_argument("--pyproject", default="pyproject.toml", help="Path to pyproject.toml")
    parser.add_argument("--output", default="docs/requirements-docs.txt", help="Output requirements file")
    args = parser.parse_args()

    pyproject_path = pathlib.Path(args.pyproject)
    if not pyproject_path.exists():
        print(f"pyproject file not found: {pyproject_path}", file=sys.stderr)
        sys.exit(2)

    data = loads_toml(pyproject_path.read_text(encoding="utf-8"))

    deps = None
    # PEP 621 style: dependency-groups
    dep_groups = data.get("dependency-groups") or {}
    if isinstance(dep_groups, dict) and "docs" in dep_groups:
        deps = dep_groups["docs"]

    # fallback: project.optional-dependencies
    if deps is None:
        project = data.get("project") or {}
        opt = project.get("optional-dependencies") or {}
        if isinstance(opt, dict) and "docs" in opt:
            deps = opt["docs"]

    if not deps:
        print("No `docs` dependency group found in pyproject.toml", file=sys.stderr)
        sys.exit(3)

    # Normalize items to strings
    lines = []
    for item in deps:
        if isinstance(item, str):
            lines.append(item)
        else:
            # fallback conversion
            lines.append(str(item))

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    print(f"Wrote {len(lines)} requirement(s) to {out_path}")


if __name__ == "__main__":
    main()
