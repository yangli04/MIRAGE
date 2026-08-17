# Compute motif prior probabilities

`compute_prior` computes prior probabilities for 5-mer motif sequences.

## Usage

``` r
compute_prior(
  motif_freq_exp,
  motif_freq_bg,
  target_motif = "DRACH",
  target_base_index = 3,
  prob_methylation = 0.007
)
```

## Arguments

- motif_freq_exp:

  A `data.frame` with columns "motif" and "freq". "motif" stores 5-mer
  motif sequences and "freq" stores the frequency of each motif among
  independently defined signal sites.

- motif_freq_bg:

  A `data.frame` with columns "motif", "tx_count", and "genome_count".
  "tx_count" and "genome_count" are motif counts in the transcriptome
  and genome, respectively.

- target_motif:

  A `string` representing the target motif pattern. MIRAGE expands
  supported degenerate bases such as N, R, Y, W, and H. The default is
  "DRACH"; set this explicitly for the assay being analyzed.

- target_base_index:

  A positive integer numeric value indicating the position of the target
  base within `target_motif`. Default is 3.

- prob_methylation:

  A positive numeric value indicating the genome- or transcriptome-wide
  prior probability that the target base carries true signal. This
  argument name is retained for backward compatibility. Default is
  0.007.

## Value

A `data.frame` including the information from `motif_freq_exp` and
`motif_freq_bg`, plus the prior probability for each 5-mer motif
sequence.

## Details

`compute_prior` calculates the prior probability that a site is a true
signal site given its local 5-mer motif sequence. Using the property of
conditional likelihood, we compute the prior probability as follows. Let
\\P(S=1\|D_i)\\ denote the prior probability that site \\i\\ carries
true signal given local motif \\D_i\\, we have \$\$P(S=1\|D_i=j) =
\frac{P(D_i=j\| S=1) P (S=1)}{ P(D_i=j)}\$\$ where \\P(D_i=j\|S=1)\\ is
the likelihood for each motif among true signal sites and can be
calculated using motif frequencies from an independent signal set.
\\P(S=1)\\ is the overall prior probability that the target base carries
true signal. \\P(D_i=j)\\ is the likelihood of observing that motif in
the transcriptome.

## Examples

``` r
# Build motif frequency tables from a PAR-CLIP count table.
count_table <- read.delim(system.file("extdata", "pos.read.count.example.parclip.txt",
                                      package = "MIRAGE"))
motif_counts <- table(count_table$motif)
motif_freq_exp <- data.frame(motif = names(motif_counts),
                             freq = as.numeric(motif_counts) / sum(motif_counts))

motif_freq_bg <- data.frame(motif = names(motif_counts),
                            tx_count = as.numeric(motif_counts),
                            genome_count = as.numeric(motif_counts))

# Compute the prior probability for each 5-mer motif sequence.
motif_prior <- compute_prior(motif_freq_exp, motif_freq_bg, target_motif = "NNUNN",
                             target_base_index = 3, prob_methylation = 0.01)
```
