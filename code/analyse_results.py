# Look at the statistics of the ACF and periodogram results.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Load results file
date = "2018-01-01"
df = pd.read_csv("c1_periods_{}.txt".format(date), delimiter=" ")

plt.clf()
plt.plot(df.ACF_period, df.pgram_period, "k.")
plt.xlabel("ACF period")
plt.ylabel("Pgram period")
plt.savefig("ACF_vs_pgram")

percent = 10
lower = df.ACF_period.values - df.ACF_period.values/percent
upper = df.ACF_period.values + df.ACF_period.values/percent
success_inds = (lower < df.pgram_period.values) | (df.pgram_period.values < upper)

print(len(df.ACF_period), sum(success_inds))
