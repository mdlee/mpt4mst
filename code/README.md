The `trinity` package is required: https://github.com/joachimvandekerckhove/trinity. It is used for the `callbayes` function which connects MATLAB to JAGS.

Each model is run by a MATLAB `.m` script, and has an associated JAGS `.txt` file. At the top of the MATLAB scripts, there are constants `dataDir`, `storageDir`, and `figDir` and directories defining where the data are stored, where saved JAGS runs should be stored, and where figures should be placed. There are also booleans `preLoad` and `printTrue` controling whether stored results should be used, and whether figures should be printed.

For producing the figures, some additional custom files need to be in the path:
- `Raxes.m` for creating separated axes like R graphics often uses
- `pantoneColors.mat` which provides a palette of colors (shown in `pantoneColors.png`)
- `setFigure.m` for setting the figure window, allowing for changing single and dual monitor setup
- `suplabel.m` for using overall axis labels in multi-panel figures
