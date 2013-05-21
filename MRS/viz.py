"""

Visualization functions

"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def stacked_spectra(f, spectra, fig=None):

    """
    Plot stacked spectra

    Parameters
    ----------

    """
    if fig is None:
        fig, ax = plt.subplots(1)
    else:
        fig = fig
        ax = fig.get_axes()[0]

    offset = 0
    for s in spectra:
        ax.plot(f, s+offset, color='k')
        offset += np.var(s)
    ax.invert_xaxis()
    return fig

def plot_spectra(f, spectra, n_boots=1000, fig=None, pct=68):
    """
    Plot average spectra with boot-strapped confidence intervals

    Parameters
    ----------
    pct : how many percentiles to include in the CI (default to 95%)

    """
    if fig is None:
        fig, ax = plt.subplots(1)
    else:
        fig = fig
        ax = fig.get_axes()[0]

    boots = np.zeros((n_boots, spectra.shape[-1]))
    for b in xrange(n_boots):
        boots[b] = np.random.random_integers(0,
                                             spectra.shape[-1]-1,
                                             spectra.shape[-1])

    alpha = (100 - pct)/2.
    CI = np.zeros(spectra.shape[-1])
    for x in xrange(CI.shape[0]):
        CI[x] = (stats.scoreatpercentile(boots[:, x], 100 - alpha) -
                 stats.scoreatpercentile(boots[:, x], alpha))

    y = np.mean(spectra, 0)
    ax.plot(f, y, color='k')
    ax.fill_between(f, y-CI, y+CI, alpha=0.2, color='b')
    ax.invert_xaxis()
    return fig
