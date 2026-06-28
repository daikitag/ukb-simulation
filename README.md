# snakemake-quantitative-trait-simulation

This snakemake workflow is used to conduct the following simulations:

1. Simulating genetic trajectories by using underdominance model in SLiM
2. Simulating effect sizes by using the simulated selection coefficient in SLiM simulation
3. Computing genetic values of individuals by using the simulated effect sizes
4. Simulating environmental noise

To run this job, modify `config.yaml` in `profiles` directory to set the correct configuration for slurm to execute the code.

The simulation parameters are specified in `config.yaml` in `config` directory.

Please install the following packages to run snakemake on slurm:
- snakemake
- snakemake-executor-plugin-slurm


## Phenotypic Model

Our multi-dimensional stabilizing selection model is adapted from Simons et al. (2018). The individual phenotypes are determined by

```math
\vec{r}=\sum_{l=1}^L(\vec{a}_l+\vec{a}_{l'}),
```

where $`\vec{a}_l`$ and $`\vec{a}_{l'}`$ are the phenotypic effects of the parents' alleles at site $l$. If we let $`a=||\vec{a}||`$ as the size of mutations in the $n$-dimensional trait space, the selection coefficient of a mutation can be expressed as $`s=a^2/w^2`$. Since we are assuming a diploid setting, the population-scaled selection coefficient can be expressed as

```math
S=2Ns=2N\frac{a^2}{w^2},
```

where $N$ is the population size and $w$ is the parameter that parameterizes the strength of stabilizing selection. In our simulation model, we will be simulating the population-scaled mutation size, $2Na^2$, from a Gamma distribution. The population size, $N$, will be an input parameter, because we are simulating genetic processes from a demographic model with a non-constant population size.
