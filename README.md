# SilicaLysozymeSAXS

<b> WHAT IS THIS FOR? </b>

This repository constitutes a part of the supplementary material to our publication 

<i><b>"Mechanism of silica–lysozyme composite formation unravelled by in situ fast SAXS"</b></i>
by Tomasz M. Stawski, Daniela B. van den Heuvel, Rogier Besselink, Dominique J. Tobler and Liane G. Benning

<b> Beilstein J. Nanotechnol. 2019, 10, 182–197. doi:10.3762/bjnano.10.17 </b>

 https://www.beilstein-journals.org/bjnano/articles/10/17

<b> This research was made possible by a Marie Curie grant from the European Commission: the NanoSiAl Individual Fellowship, Project No. 703015. </b>


![graphical_abstract](https://user-images.githubusercontent.com/10513547/51824443-92a8cc80-22e2-11e9-9779-77ae3e3d3763.png)

It contains routines in GNU Octave/Matlab used to processes SAXS (small-angle X-ray scattering) patterns. We used the scripts to perform fits and analyze data obtained during the experiment at the BioSAXS beamline P12 of the EMBL at PETRA III (DESY, Germany).

<b> HOW TO USE? </b>
The script "fitting_macro_cluster_SHS.m" was written and tested in GNU Octave >4. It should also work in Matlab. The following packages have to be installed as well: NAN, OPTIM.

The input and output directories can be edited within the script. By default we use the following data structure:

  1. The scattering pattern of the unaggregated silica nanoparticles and the size distribution from the Monte Carlo fit are located in:

    /distribution_data/histogram_linear.dat
    /distribution_data/silica_initial.dat

  2. The measured input scattering patterns are located in:

    /ResultBkgSub_selected

  This is only a selected data set. The patterns are background-corrected and the scattered intesities are converted to absolute units. The numbering of the files corresponds to the order of frames. We can provide a complete data set.

  3. The fits and the corresponding statistics are stored in:

    /fitting_results/selected  

  This directory already contains a complete set of fits to the data from (2). To prevent overwriting of the results, line 45 in the script should be amended with the new name of the output directory.
