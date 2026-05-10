# Contributing to AbSolution

Thank you for your interest in contributing to AbSolution! Whether you
are reporting a bug, suggesting a feature, improving documentation, or
submitting code, your input is valued.

## How to Contribute

### Reporting Bugs

If you find a bug, please [open an
issue](https://github.com/EDS-Bioinformatics-Laboratory/AbSolution/issues/new)
and include:

- A clear and descriptive title.
- Steps to reproduce the problem.
- Expected vs. actual behavior.
- Your R version, operating system, and AbSolution version
  (`packageVersion("AbSolution")`).
- If possible, a minimal reproducible example using the built-in test
  dataset.

### Suggesting Features or Enhancements

Feature requests are welcome. Please open an issue describing:

- The problem your suggestion would solve.
- Your proposed solution or approach.
- Any alternatives you have considered.

Check the
[ROADMAP.md](https://eds-bioinformatics-laboratory.github.io/AbSolution/ROADMAP.md)
to see if your idea aligns with planned developments.

### Submitting Code

1.  **Fork** the repository and create a new branch from `main`:

    ``` bash
    git checkout -b feature/your-feature-name
    ```

2.  **Follow the existing code style.** AbSolution is built with
    [golem](https://thinkr-open.github.io/golem/) and follows its module
    structure. Shiny modules are located in `R/` and prefixed with
    `mod_`.

3.  **Document your changes.** Use
    [roxygen2](https://roxygen2.r-lib.org/) for function documentation
    and update vignettes if applicable.

4.  **Test your changes.** Run `devtools::check()` and ensure no new
    warnings or errors are introduced. Tests are located in `tests/`.

5.  **Commit** with clear, descriptive messages.

6.  **Open a pull request** against `main`, describing what your changes
    do and why.

### Improving Documentation

Improvements to documentation, vignettes, and tutorials are always
appreciated. If you find something unclear or incomplete, feel free to
submit a pull request or open an issue.

## Development Setup

This section is for contributors who want to modify the AbSolution
source code. If you only want to use AbSolution, see the [Installation
instructions](https://eds-bioinformatics-laboratory.github.io/AbSolution/README.html#installation)
instead.

Clone the repository:

``` bash
git clone https://github.com/EDS-Bioinformatics-Laboratory/AbSolution.git
cd AbSolution
```

Then, in R (with your working directory set to the cloned repository):

``` r

# Install development tools
install.packages("devtools")

# Install all AbSolution dependencies from DESCRIPTION
devtools::install_deps(dependencies = TRUE)

# Install Bioconductor dependencies
BiocManager::install(c("Biostrings", "IRanges", "pwalign"))

# Load AbSolution from local source code (without installing)
# This lets you test your changes immediately
devtools::load_all()

# Launch the app
run_app()
```

## Code of Conduct

Please be respectful and constructive in all interactions. We are
committed to fostering a welcoming and inclusive environment for
everyone, regardless of background or experience level.

## Questions?

If you have questions about contributing, feel free to [open an
issue](https://github.com/EDS-Bioinformatics-Laboratory/AbSolution/issues)
or contact us at <r.garciavaliente@amsterdamumc.nl>.
