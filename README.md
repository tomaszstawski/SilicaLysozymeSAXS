# SilicaLysozymeSAXS

WHAT IS THIS FOR?
This repository constitutes a part of the Supplementary Material for our publication "Mechanism of formation of silica-lysozyme composites followed by in situ fast SAXS" (in preparation). It contains routines in GNU Octave/Matlab used to processes SAXS (small-angle X-ray scattering) patterns. We used the scripts to perform fits and analyze data obtained during the experiment at the BioSAXS beamline P12 of the EMBL at PETRA III (DESY, Germany).

HOW TO USE?
The script "fitting_macro_cluster_SHS.m" was written and tested in GNU Octave >4. It should also work in Matlab. The following packages has to be also installed: NAN, OPTIM.

The input and output directories can be edited within the script. By default we use the following data structure:

  1. The scattering pattern of the unaggregated silica nanoparticles and the size distribution from the Monte Carlo fit are located in:

    /distribution_data/histogram_linear.dat
    /distribution_data/silica_initial.dat

  2. The measured input scattering patterns are located in:

    /ResultBkgSub_selected

    This is only a selected data set. The patterns are bakcground-corrected and the scattered intesities are converted to absolute units. The numbering of the files corresponds to the order of the frames. We can provide a complete data set.

  3. The fits and the corresponding statistics are stored in:

    /fitting_results/selected  

    This directory already contains a complete set of fits to the data from (2). To prevent the overwriting of the results, line 45 in the script should be amended with the new name of the output directory. 