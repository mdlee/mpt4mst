### Trinity is required

The `trinity` package is required: https://github.com/joachimvandekerckhove/trinity. It is used for the `callbayes` function which connects MATLAB to JAGS.

Python or R users can use `PyJAGS` or `RJAGS` or similar to connect to the JAGS scripts.

### Running the models
Each model is run by a MATLAB `.m` script, and has an associated JAGS `_jags.txt` file. At the top of the MATLAB scripts, there are constants `dataDir`, `storageDir`, and `figDir` and directories defining where the data are stored, where saved JAGS runs should be stored, and where figures should be placed. There are also booleans `preLoad` and `printTrue` controling whether stored results should be used, and whether figures should be printed.

#### Old-new task

The models are:

1. `mptStudyTestON_0`: Basic model of the old-new task. The graphical model is shown in Figure 3 of the paper.
2. `mptStudyTestON_7`: Hierarchical latent-mixture model for two groups using MPT process, plus contamimant groups. The graphical model is shown in Figure 6 of the paper.
3. `mptStudyTestON_8`: Hierarchical latent-mixture model for two groups using MPT process, plus contamimant groups, and logistic model relating lure bin to discriminability. Note that the inferences of `mptStudyTestON_7` are used to remove contaminants.
4. `mptStudyTestON_9`: Recommended first-pass measurement model. This model is not explicitly in the paper, but is the right place to start using the MPT model to measure the performance of a set of individuals. The `rho`, `psi`, `gamma` parameters are independent, and the logistic model relating lure bin to discriminability is included to provide `beta`, `tau` and `delta` parameters for each individual.

#### Old-similar-new task

The models are:

1. `mptStudyTestOSN_0`: Basic model of the old-similar-new task.
2. `mptStudyTestOSN_7`: Hierarchical latent-mixture model for two groups using MPT process, plus contamimant groups.
3. `mptStudyTestOSN_7b`: Hierarchical latent-mixture model for just one group using MPT process, plus contamimant groups.
4. `mptStudyTestON_8`: Hierarchical latent-mixture model for two groups using MPT process, plus contamimant groups, and logistic model relating lure bin to discriminability. Note that the inferences of `mptStudyTestOSN_7` are used to remove contaminants.

### Supporting files

For producing the figures, some additional custom files need to be in the path:
- `Raxes.m` for creating separated axes like R graphics often uses
- `pantoneColors.mat` which provides a palette of colors (shown in `pantoneColors.png`)
- `setFigure.m` for setting the figure window, allowing for changing single and dual monitor setup
- `suplabel.m` for using overall axis labels in multi-panel figures
