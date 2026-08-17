# MIRAGE Reproducible Environment

The MIRAGE package itself is a standard R package whose dependencies are
listed in `DESCRIPTION`. The micromamba environment below provides R together
with the packages needed to install MIRAGE, run the tutorials, and rebuild the
pkgdown site.

```bash
micromamba env create -f envs/mirage-r.yml
micromamba run -n mirage-r R CMD INSTALL .
```

Render the vignettes and pkgdown site:

```bash
micromamba run -n mirage-r R CMD build .
micromamba run -n mirage-r Rscript -e 'pkgdown::build_site(lazy = FALSE)'
```

Every tutorial reads its input tables from `inst/extdata/` (installed with the
package), so no external data or additional environment is required.
