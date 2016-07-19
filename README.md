# kalesalad
Measuring K2 rotation periods.

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
download_epic.py: code for downloading the EPIC catalogue.
download_fits_files.py: code for downloading k2 fits files.
plot_results.py: code for generating results plots.
GPkalesalad.py: measuring rotation periods using GPs (in dev).
GProtation.py: contains functions used by GPkalesalad.
