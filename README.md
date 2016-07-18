# kalesalad
Measuring K2 rotation periods.

kalesalad.py: measures rotation periods of K2 targets using the ACF method.
injection_tests.py: performs some (very crude) injection and recovery tests.
simple_acf.py: Simple ACF algorithm with highest peak = period.
download_epic.py: code for downloading the EPIC catalogue
plot_results.py: code for generating results plots.
GPkalesalad.py: measuring rotation periods using GPs (in dev).
GProtation.py: contains functions used by GPkalesalad.

c01_????_periods.txt: results files. Columns = epicid, period, period_err
(currently all zeros!)

????_ids.txt: text files containing list of epic ids (last 5 digits only).

Now run on campaign 1.
Downloading campaign 2 light curves right now!
