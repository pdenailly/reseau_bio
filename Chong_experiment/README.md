## Chong's experiment

#### Main script:

This folder is for simulating Chong's experiment(Summerized in [Figure 5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4105854/figure/F5/)), `Chong_sim.py` is the main script, it can be used as follows:

    # python Chong_sim.py <nbr_of_sims> <params_file> <output_pre> <output_post>
    python Chong_sim.py 5 params.ini pre post

The script take as parameters :
1. *nbr_of_sims* : Numbner of simuations
2. *params_file* : Path to parameter file
3. *output_pre* : The name of the directory in which the files generated from the first part of the simulation are saved (before adding the gyrase).
4.  *output_post* : The name of the directory in which the files generated from the second part, after gyrase's addition (the resume of the simulation)

NOTE : Both will be combined and plotted in one graph afterwards, in the previous example the results of `pre0`  and `post0` are combined and plotted in the graphe `<plotted_indo>_sim0.pdf` (see `sim_plots` directory, and the description of the output below).

#### Output files:

#### Things to improve:
* Complete the Readme file.
* Group the output files (i.e `pre` and `post`) in one directory.
