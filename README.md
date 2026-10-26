# Qualifying Paper 1 

This repository contains the document and R code used for the simulation study conducted as part of my first PhD qualifying paper. The paper I explored was *Large-scale causal discovery using interventional data sheds light on gene network structure in K562 cells* by Brown et al. (2025). This study looks at estimating the network of direct causal effects using a sparse regression method, inspre. My proposed extension investigates whether the sparse direct causal-effect network estimated using inspre can be used to reconstruct total causal effects with greater estimation accuracy and better structural recovery than direct instrumental variable estimation. 

## This repository includes 

1. [Qualifying Paper](Qualifying_Paper_1.pdf): The qualifying paper itself containing the problem definition, significance, limitations, proposed extension, simulation results, and discussion.

2. [Simulation code](simulation.R): The R code used to simulate the gene regulatory networks, estimate total causal effects using direct IV estimation and sparse reconstruction, and evaluate their performance.

3. [Corrected Cross-Validation Function](fixed_inspre_cv.R): A modified version of the inspre package's cross-validation prediction function that uses the estimated direct causal-effect matrix (G_hat) instead of the estimated total causal-effect matrix (R_hat) when calculating structural-model prediction errors.

4. [README](README.md): This file, outlining the contents of the repository! 

