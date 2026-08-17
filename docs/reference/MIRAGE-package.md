# MIRAGE: Mutation-Encoded Inference of RNA Activity

MIRAGE estimates site-level RNA activity from matched treatment and
control read-count tables produced by mutation-encoded RNA assays (for
example BACS, SHAPE-MaP, and PAR-CLIP). It jointly models the on-target
conversion rate, the background conversion rate, sequencing error, and
allele status, and returns a latent activity level `beta_est` for every
site.

## Details

The main entry points are
[`estimate_inference_with_empirical`](https://yangli04.github.io/MIRAGE/reference/estimate_inference_with_empirical.md)
(motif-free inference),
[`estimate_inference_with_prior`](https://yangli04.github.io/MIRAGE/reference/estimate_inference_with_prior.md)
together with
[`compute_prior`](https://yangli04.github.io/MIRAGE/reference/compute_prior.md)
(motif-aware inference), and
[`plot_signal_scatter`](https://yangli04.github.io/MIRAGE/reference/plot_signal_scatter.md)
(diagnostics).

## See also

Useful links:

- <https://github.com/yangli04/MIRAGE>

- <https://yangli04.github.io/MIRAGE/>

- Report bugs at <https://github.com/yangli04/MIRAGE/issues>

## Author

**Maintainer**: Yang Li <yliuchicago@uchicago.edu>

Authors:

- Yang Li <yliuchicago@uchicago.edu>

- Mengjie Chen <mengjiechen@uchicago.edu>

- Shun Liu <shunliu@uchicago.edu>
