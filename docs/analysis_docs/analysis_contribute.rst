Contributing to DELi Analysis
=============================

The best entry point for contributing to DELi analysis is to review the methods available 
in `cube_class.py`. If you find a method that would be ideal to add into automated analysis/
report generation, follow these steps:

1. Add the method to `analysis_parse.py`.
2. Update `analysis_report_gen.py` to include the new method.
3. Modify the `template.html` to incorporate the new method.

Users interested in generating custom reports can provide their own `.html` template to 
be used with DELi interfacing with Jinja.