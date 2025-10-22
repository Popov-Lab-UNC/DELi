Contributing to DELi Analysis
=============================

The best entry point for contributing to DELi analysis is to review the methods available 
in `cube_class.py`. If you find a method that would be ideal to add into automated analysis/
report generation, follow these steps:

**For methods that should be included in the automated analysis pipeline:**

1. Add the method to `cube_class.py` (if not already present)
2. Add the method call to `analysis_parse.py` with appropriate flag handling
3. Add the method call to `cli.py` in the `analyze` command with appropriate flag handling
4. If the method generates plots/images, update `analysis_report_gen.py` to include the new plot directory
5. If the method generates plots/images, modify `analysis_report.html` template to display the new plots

**Configuration:**
- Add any new configuration flags to the `flags` section of your YAML config files
- Update `analysis_config.rst` documentation to include the new flags

**Development and Testing:**
- `analysis_parse.py` is utilized for development and testing purposes separate from the main CLI
- It provides a standalone entry point for running analysis without the full CLI infrastructure
- Useful for debugging and rapid prototyping of new analysis features

**Custom Reports:**
Users interested in generating custom reports can provide their own `.html` template to 
be used with DELi interfacing with Jinja2.