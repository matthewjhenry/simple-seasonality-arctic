import numpy as np


def cfd_old(y, x, axis=-1, inc_bdy='both'):
    """
        Vectorized central finite difference for N-D array Y, 1-D array X,
        whose length matches with the /axis/^th dimension of Y.
        Input:
            y: N-D array
            x: Coordinate array
            axis: Dimension where len(y_axis) == len(x)
            inc_bdy: Include 0 {'low'} or N {'high'} boundary in calculations
    """
    rank = len(y.shape)
    if axis < 0:
        axis = rank+axis
    if axis > rank-1 or axis > 3:
        raise ValueError('Axis {} out of range for {} dimensions'.format(axis, rank))

    dydx = np.zeros(y.shape)

    if inc_bdy in ['low', 'Low', 'L', 'l']:
        inc_l = True
        inc_h = False
    elif inc_bdy in ['high', 'High', 'H', 'h']:
        inc_l = False
        inc_h = True
    else:
        inc_l = True
        inc_h = True

    for i in range(rank):
        if i != axis:
            x = np.expand_dims(x, axis=i)

    if axis == 0:
        dydx[1:-1, ...] = (y[2:, ...]-y[:-2, ...]) / (x[2:, ...]-x[:-2, ...])
        if inc_l:
            dydx[0, ...] = (y[1, ...]-y[0, ...]) / (x[1, ...]-x[0, ...])
        if inc_h:
            dydx[-1, ...] = (y[-1, ...]-y[-2, ...]) / (x[-1, ...]-x[-2, ...])

    elif axis == 1:
        dydx[:, 1:-1, ...] = (y[:, 2:, ...]-y[:, :-2, ...]) /\
                             (x[:, 2:, ...]-x[:, :-2, ...])
        if inc_l:
            dydx[:, 0, ...] = (y[:, 1, ...]-y[:, 0, ...]) /\
                              (x[:, 1, ...]-x[:, 0, ...])
        if inc_h:
            dydx[:, -1, ...] = (y[:, -1, ...]-y[:, -2, ...]) /\
                               (x[:, -1, ...]-x[:, -2, ...])

    elif axis == 2:
        dydx[:, :, 1:-1, ...] = (y[:, :, 2:, ...]-y[:, :, :-2, ...]) /\
                                (x[:, :, 2:, ...]-x[:, :, :-2, ...])
        if inc_l:
            dydx[:, :, 0, ...] = (y[:, :, 1, ...]-y[:, :, 0, ...]) /\
                                 (x[:, :, 1, ...]-x[:, :, 0, ...])
        if inc_h:
            dydx[:, :, -1, ...] = (y[:, :, -1, ...]-y[:, :, -2, ...]) /\
                                  (x[:, :, -1, ...]-x[:, :, -2, ...])

    elif axis == 3:
        dydx[:, :, :, 1:-1] = (y[:, :, :, 2:]-y[:, :, :, :-2]) /\
                              (x[:, :, :, 2:]-x[:, :, :, :-2])
        if inc_l:
            dydx[:, :, :, 0] = (y[:, :, :, 1]-y[:, :, :, 0]) /\
                               (x[:, :, :, 1]-x[:, :, :, 0])
        if inc_h:
            dydx[:, :, :, -1] = (y[:, :, :, -1]-y[:, :, :, -2]) /\
                                (x[:, :, :, -1]-x[:, :, :, -2])

    return dydx


def cfd(F_in, x_in, axis=0, inc_bdy='both', cyclic=False):
    """
    Vectorized central finite difference for N-D array F_in by x_in.
    Parameters
    ----------
    F_in : array_like, N-dimensions
            Array to be differentiated
    x_in : array_like, either 1-D or F_in.shape == x_in.shape
            Coordinate array
    axis : integer
            Dimension where F_in.shape[axis] == x_in.shape[axis]
    inc_bdy : String, optional
            Include 0 'low' or N 'high', or both boundary(ies) in calculations
    cyclic : Boolean, optional
            Data is cyclic on <axis> if true

    Returns
    -------
    Centred
    """
    F = np.swapaxes(F_in, axis, 0)
    if len(x_in.shape) == len(F_in.shape):
        x = np.swapaxes(x_in, axis, 0)
    else:
        x = x_in.copy()
        for i in range(len(F_in.shape) - len(x_in.shape)):
            x = x[:, None]

    if inc_bdy in ['low', 'Low', 'L', 'l']:
        inc_l = True
        inc_h = False
    elif inc_bdy in ['high', 'High', 'H', 'h']:
        inc_l = False
        inc_h = True
    else:
        inc_l = True
        inc_h = True

    dFdx = np.zeros(F.shape)
    if cyclic:
        dFdx[1:-1, ...] = (F[2:, ...] - F[:-2, ...])#/(x[2:, ...] - x[:-2, ...])
        dFdx[0, ...] = (F[1, ...] - F[-1, ...])#/(x[-1, ...] - x[1, ...])
        dFdx[-1, ...] = (F[0, ...] - F[-2, ...])#/(x[-2, ...] - x[0, ...])
    else:
        dFdx[1:-1, ...] = (F[2:, ...] - F[:-2, ...])/(x[2:, ...] - x[:-2, ...])
        if inc_l:
            dFdx[0, ...] = (F[0, ...] - F[1, ...])/(x[0, ...] - x[1, ...])
        if inc_h:
            dFdx[-1, ...] = (F[-2, ...] - F[-1, ...])/(x[-2, ...] - x[-1, ...])

    return np.swapaxes(dFdx, axis, 0)