import numpy as np
import matplotlib.pyplot as plt
from kepler_data import load_kepler_data
from simple_acf import simple_acf
import glob

# load lightcurve
kid = "002450729"
path = "/Users/ruthangus/.kplr/data/lightcurves"
fnames = glob.glob("{0}/{1}/kplr*_llc.fits".format(path, kid))
x, y, yerr = load_kepler_data(fnames)

period, acf_smooth, lags = simple_acf(x, y)
print(period)
plt.clf()
plt.plot(lags, acf_smooth)
plt.axvline(period, color="r")
plt.savefig("test")
