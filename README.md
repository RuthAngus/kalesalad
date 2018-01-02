# kalesalad
Measuring K2 rotation periods.

data/download_epic.py: code for downloading the EPIC catalogue.

data/download_fits_files.py: code for downloading k2 fits files.

identify_dwarfs.py: Using the TGAS/EPIC cross-matched catalog, this script
identifies stars that are most likely giants according to their position on
the HR diagram and removes them from the catalog.
It also removes stars cooler than 3200 K and hotter than 8000 K.

kalesalad.py: measures rotation periods of K2 targets using the ACF method.
In order to run kalesalad you must have a file called c<campaign>_targets.txt.
This is a list of epic ids in that campaign, e.g. 201120213.
This file can be created by doing an ls data/c<campaign>/ >
c<campaign>_targets.txt.
kalesalad creates a file called c<campaign>_periods.txt containing a list of
epic ids and a list of rotation periods.
This file is appended to by kalesalad, so you must delete it or rename it
between runs.

simple_acf.py: Simple ACF algorithm with highest peak = period.
This needs some work --- peak detection could be improved!

injection_tests.py: performs some (very crude) injection and recovery tests.

plot_results.py: code for generating results plots.

GPkalesalad.py: measuring rotation periods using GPs (in dev).

GProtation.py: contains functions used by GPkalesalad.

TODO:
=====

1. Crossmatch Gaia and K2 data (at least for campaign 1).

2. Use Gaia parallaxes and colours to cut giants from the sample.

3. Measure rotation periods using ACF and LS for now.

4. Plot Teff vs Prot

Next: Run on all light curves and look at the ones where acf and pgram agree.
