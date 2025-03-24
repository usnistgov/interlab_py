<!-- markdownlint-disable MD041 -->

<!-- prettier-ignore-start -->
[![Repo][repo-badge]][repo-link]
[![Docs][docs-badge]][docs-link]
[![PyPI license][license-badge]][license-link]
[![PyPI version][pypi-badge]][pypi-link]
[![Conda (channel only)][conda-badge]][conda-link]
[![Code style: ruff][ruff-badge]][ruff-link]
[![uv][uv-badge]][uv-link]

<!--
  For more badges, see
  https://shields.io/category/other
  https://naereen.github.io/badges/
  [pypi-badge]: https://badge.fury.io/py/interlab_py
-->

[ruff-badge]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json
[ruff-link]: https://github.com/astral-sh/ruff
[uv-badge]: https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json
[uv-link]: https://github.com/astral-sh/uv
[pypi-badge]: https://img.shields.io/pypi/v/interlab_py
[pypi-link]: https://pypi.org/project/interlab_py
[docs-badge]: https://img.shields.io/badge/docs-sphinx-informational
[docs-link]: https://pages.nist.gov/interlab_py/
[repo-badge]: https://img.shields.io/badge/--181717?logo=github&logoColor=ffffff
[repo-link]: https://github.com/wpk-nist-gov/interlab_py
[conda-badge]: https://img.shields.io/conda/v/wpk-nist/interlab_py
[conda-link]: https://anaconda.org/wpk-nist/interlab_py
[license-badge]: https://img.shields.io/pypi/l/interlab_py?color=informational
[license-link]: https://github.com/wpk-nist-gov/interlab_py/blob/main/LICENSE

<!-- other links -->

<!-- prettier-ignore-end -->

# `interlab`: A python module for consensus analysis in interlaboratory studies

Welcome to the home page for analysis software of `interlab`: A python module
for consensus analysis in interlaboratory studies. This software is a Python
package that will perform consensus analysis on spectral data such as NMR, GC-MS
and LC-MS.

## Overview

`interlab` allows researchers to identify laboratories producing data closest to
the consensus values, thereby ensuring that untargeted studies are using the
most precise data available to them. The software was originally developed for
analyzing NMR data [1, 2] but can be applied to any array data, including Raman
or FTIR spectroscopy and GC-MS or LC-MS. Details on the implementation of the
code can be found in Ref. [1]. The input for the code consists of a set of
sample labels identifying the physical objects measured in the interlaboratory
study, facility labels that identify the facility of origin of the measurements,
and the data themselves. It is the responsibility of the user to format the data
and metadata so that the code can read it. In addition, the user must specify
the distance function that will be used to compare the spectra and the
statistical distribution that these distances will be fit to. One example is
given in the
[example notebook](https://pages.nist.gov/interlab_py/analysis_demo.html). Given
the input data, the code will perform the following tasks:

- Calculates the interspectral distances
- Fits the project’s distribution function to the distance data and calculate
  the corresponding scores.
- Identifies outliers within each spectral population
- Conducts a principal components analysis on the scores and compute the
  projected statistical distance
- Uses the projected statistical distance to determine the data set outliers.

The software cannot be used out of the box. Users must create an interface to
their own software, and that interface will be specific to the user’s
application. The example notebook, demonstrates one such possible interface.

<!-- Quick overview... -->

<!-- ## Status -->

<!-- This package is actively used by the author.  -->
<!-- Please feel free to create a pull request  -->
<!-- for wanted features and suggestions! -->

## Links

<!-- prettier-ignore-start -->
[NIST GitHub Organization](https://github.com/usnistgov)
[NIST MetQual program](https://www.nist.gov/programs-projects/metabolomics-quality-assurance-and-quality-control-materials-metqual-program)
[NIST Chemical Informatics Research Group](https://www.nist.gov/mml/csd/chemical-informatics-research-group)
<!-- prettier-ignore-end -->

<!-- ## Example usage -->

<!-- ```python -->
<!-- import interlab -->
<!-- ``` -->

<!-- end-docs -->

## Installation

<!-- start-installation -->

Use one of the following

```bash
pip install interlab_py
```

or

```bash
conda install -c wpk-nist interlab_py
```

<!-- end-installation -->

## Documentation

See the [documentation][docs-link] for further details.

## License

This is free software. See [LICENSE][license-link].

<!-- ## Related work -->

<!-- Any other stuff to mention.... -->

## Contact

The author can be reached at <david.sheen@nist.gov>.

## Credits

This package was created using
[Cookiecutter](https://github.com/audreyr/cookiecutter) with the
[usnistgov/cookiecutter-nist-python](https://github.com/usnistgov/cookiecutter-nist-python)
template.
