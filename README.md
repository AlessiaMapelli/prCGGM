# prCGGM

**Prior-informed conditional Gaussian graphical models** for reconstructing
personalised molecular networks, such as protein-protein interaction
networks, from omics data.

prCGGM estimates a Gaussian graphical model whose mean and precision matrix depend
linearly on subject-level covariates,
Ω(X) = Δ₀ + Σₕ Xₕ Δₕ,
and integrates prior knowledge on edges (for example 'STRING' confidence
scores) through a weighted sparse-group lasso penalty applied only to the
population-level network Δ₀. Covariate-specific changes Δₕ stay entirely
data-driven. Personalised networks of new individuals are obtained from
their covariates.

The method is described in:

> Mapelli, A., Massi, M. C., Cuccuru, G., Di Angelantonio, E. and Ieva, F.
> (2026). Prior-informed conditional Gaussian graphical models: an
> application to protein interaction network reconstruction.
> [arXiv:2606.31805](https://arxiv.org/abs/2606.31805).

## Installation

```r
install.packages("prCGGM")
```

## Usage

```r
library(prCGGM)
data(prcggm_example)

set.seed(1)
fit <- prcggm(prcggm_example$x,
              covariates = prcggm_example$covariates,
              prior = prcggm_example$prior)
summary(fit)
plot(fit)                        # baseline network
plot(fit, component = "x1_1")    # change associated with x1 = 1

# Personalised networks of new individuals
nets <- predict(fit, data.frame(x1 = factor(c("0", "1"))))
compare_networks(nets$Subject_1, nets$Subject_2, method = "frobenius")
```

Depending on the inputs, `prcggm()` fits a standard GGM (`x` only), a
prior-weighted GGM (`prior`), a conditional GGM (`covariates`) or the full
prior-informed conditional GGM (both). Node-wise regressions run
sequentially, in parallel (`ncores`), or as SLURM array jobs with
`prcggm_slurm_prepare()`, `prcggm_node()` and `prcggm_slurm_collect()`.

See `vignette("prCGGM")` for a complete walk-through.

## Reproducing the paper

The scripts used for the simulation study and the UK Biobank application are
available at
<https://github.com/AlessiaMapelli/Prior-informed-conditional-GGMs>.
