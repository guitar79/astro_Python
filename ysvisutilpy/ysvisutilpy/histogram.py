import numpy as np

__all__ = ["hist", "hist2"]


def hist(data, bins, range=None, normed=None, weights=None, density=None):
    """Simplify np.histogram and plt.bar/hbar.
    Parameters
    ----------
    See np.histogram

    Example
    -------
    >>> fig, ax = plt.subplots(1, 1)
    >>> pos, n, width, edges = hist(data, np.arange(-3, 3, 0.5))
    >>> ax.bar(pos, n, width, fill=False, edgecolor='k')
    """
    n, edges = np.histogram(
        data, bins=bins, range=range, normed=normed, weights=weights, density=density
    )
    pos = (edges[1:] + edges[:-1]) / 2

    return pos, n, np.ediff1d(edges), edges


def hist2(x, y, xbins, y_bin_fun=np.nansum):
    """Draw histogram by aggregating y-values according to x-bins."""
    x = np.asarray(x).ravel()  # following np.histogram's _ravel_and_check_weights
    y = np.asarray(y).ravel()
    if x.size != y.size:
        raise ValueError("x, y must have the same size")

    _, x_edges = np.histogram(x, bins=xbins)
    # ^ xdata is used just as a dummy
    pos = (x_edges[1:] + x_edges[:-1]) / 2

    y_binned = []
    for idx, start in enumerate(x_edges[:-1]):
        _end = x_edges[idx + 1]
        _yidxs = np.where((x >= start) & (x < _end))[0]
        if len(_yidxs) == 0:
            y_binned.append(0)
        else:
            y_binned.append(y_bin_fun(y[_yidxs]))
    return pos, np.array(y_binned), np.ediff1d(x_edges), x_edges
