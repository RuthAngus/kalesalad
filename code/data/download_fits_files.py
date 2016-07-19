import numpy as np
import matplotlib.pyplot as plt
import wget
from download_epic import get_catalog
import sys

def get_fits(c, df):
    ids = df.epic_number
    print(ids)
    prefixes = [str(i)[:4] for i in ids]
    suffixes = [str(i)[4:] for i in ids]

    for i, id in enumerate(ids):
        url0 = "https://archive.stsci.edu/missions/hlsp/everest/"
        url1 = "c{0}/{1}00000/{2}/".format(c, prefixes[i], suffixes[i])
        url2 = "hlsp_everest_k2_llc_{0}-c{1}_kepler_v1.0_lc.fits".format(id, c)
        wget.download("{0}{1}{2}".format(url0, url1, url2))

if __name__ == "__main__":
    c = str(sys.argv[1])
    df = get_catalog("k2targets")
    get_fits(c, df)
