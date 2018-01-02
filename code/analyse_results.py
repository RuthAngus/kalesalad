# Look at the statistics of the ACF and periodogram results.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# Load results file
date = "2018-01-01"
df = pd.read_csv("c1_periods_{}.txt".format(date), delimiter=" ")

xs = np.linspace(min(df.ACF_period), max(df.pgram_period), 100)
plt.clf()
plt.plot(df.ACF_period, df.pgram_period, "k.")
plt.plot(xs, xs, "--", color=".5")
plt.plot(xs, xs*1.1, "--", color="r")
plt.plot(xs, xs*.9, "--", color="r")
plt.xlabel("ACF period")
plt.ylabel("Pgram period")

percent = 10
lower = df.ACF_period.values - df.ACF_period.values/percent
upper = df.ACF_period.values + df.ACF_period.values/percent
success_inds = (lower < df.pgram_period.values) * \
    (df.pgram_period.values < upper)

plt.plot(df.ACF_period[success_inds], df.pgram_period[success_inds], "b.")
plt.savefig("ACF_vs_pgram")

print(df.epic_id.values[success_inds])
print(sum(success_inds), "successful out of", len(df.ACF_period))
