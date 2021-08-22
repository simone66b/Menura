# Menura

Fits user-defined stochastic diffusion models to univariate trait data on a phylogeny. 
Also fits three canned models: Ornstein Uhlenbeck (OU; Gaussian), Cox, Ingersoll, Ross (CIR; Non-Gaussian) and an OU-like 
model with a Beta stationary distribution (Beta, Non-Gaussian). Models are fitted and parameters are estimated using 
a Data Augmentation - Metropolis Hastings algorithm. Output can be analysed using the `coda` package , 
`tidybayes` package or other appropriate packages for Bayesian analysis and visualisation. 

Menura is the genus name for the Australian Superb Lyrebird (Menura novaehollandiae), known for being the world's largest passerine, 
its elaborate tail and its talent for mimicry. 

We wrote menura in order to study the behaviour of SDE models on phylogenies. 
Thus, menura is EXPERIMENTAL SOFTWARE. If you are unsure whether menura is appropriate for your own analyses, it probably isn't.


# Installation

Install directly from using the devtools package

```
library(devtools)
install_github("https://github.com/simone66b/Menura", subdir = "menura")
```
