# Phylogenetics and tree rendering

To be added:

- `phylo_pipeline.sh` — maximum-likelihood reconstruction (Figure 3)
- `beast.sh` — Bayesian phylodynamic analysis, three independent chains (Figure 4)
- `beast_config.xml` — BEAST v1.10.4 run configuration
- `itol/` — iTOL v6 annotation files for Figures 3 and 4 and Supplementary Figure 5

## Analysis summary

**Figure 3.** SNPs called against D39 (CP000410.2) with Bowtie2, SAMtools and BCFtools;
phage regions (Phaster), repeats and recombination (Gubbins) masked; tree built with RAxML.

**Figure 4.** BEAST v1.10.4; TVMe substitution model selected by ModelFinder under the BIC;
uncorrelated exponential relaxed clock; Bayesian skyline coalescent tree prior; year of
isolation as tip dates; asymmetric discrete trait model for country of isolation;
500 million steps sampled every 25,000. Run as three independent chains, assessed for
convergence in Tracer v1.7.1; the chain with the highest posterior (joint ESS 3,631) was
annotated with TreeAnnotator v1.10.
