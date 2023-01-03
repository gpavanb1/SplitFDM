from enum import Enum
from .error import SFDM

# TODO: Add more schemes
Schemes = Enum("Schemes", "CENTRAL")
Directions = Enum("Directions", "WEST EAST")


def Dx(F, cell_sub, scheme):
    """
    Calculate the first derivative of a given stencil.

    Parameters
    ----------
    F : function
        The function.
    cell_sub : list of Cell
        The stencil.
    scheme : Schemes
        The scheme to use.

    Returns
    -------
    numpy.ndarray
        The first derivative.
    """

    if scheme == Schemes.CENTRAL:
        if len(cell_sub) != 3:
            raise SFDM("Improper stencil size for central difference scheme")

        ul = cell_sub[0].values()
        ur = cell_sub[2].values()

        dx = cell_sub[2].x() - cell_sub[0].x()

        return (ur - ul) / dx
    else:
        raise SFDM("Unsupported scheme")


def D2x(D, cell_sub, scheme):
    """
    Calculate the second derivative of a given stencil.

    Parameters
    ----------
    D : function
        The function.
    cell_sub : list of Cell
        The stencil.

    Returns
    -------
    numpy.ndarray
        The second derivative.
    """

    # Only central scheme
    if scheme == Schemes.CENTRAL:
        if len(cell_sub) != 3:
            raise SFDM("Improper stencil size for central difference scheme")

        # West derivative
        ul = cell_sub[0].values()
        uc = cell_sub[1].values()
        Dl = D(ul)
        Dr = D(uc)
        dxw = 0.5 * (cell_sub[1].x() - cell_sub[0].x())
        Dw = (Dr - Dl) / dxw

        # East derivative
        uc = cell_sub[1].values()
        ur = cell_sub[2].values()
        Dl = D(uc)
        Dr = D(ur)
        dxe = 0.5 * (cell_sub[2].x() - cell_sub[1].x())
        De = (Dr - Dl) / dxe

        dx = 0.5 * (dxw + dxe)
        return (De - Dw) / dx
    else:
        raise SFDM("Unsupported scheme")
