from .domain import Domain
from .constants import btype, btype_map


def apply_BC(d: Domain, v: str, bc: str = "periodic", xmin=0.0, xmax=1.0):
    """
    Apply boundary conditions to the given domain.

    Parameters
    ----------
    d : Domain
        The domain to apply the boundary conditions to.
    v : str
        The name of the variable to apply the boundary conditions to.
    bc : str, optional
        The type of boundary condition to apply. Acceptable values are "periodic" and "outflow". Default is "periodic".
    xmin : float, optional
        The minimum x-value of the domain. Default is 0.0.
    xmax : float, optional
        The maximum x-value of the domain. Default is 1.0.

    Raises
    ------
    NotImplementedError
        If an unsupported boundary condition is specified.
    """

    # Get cells
    cells = d.cells()

    # Find index of component
    idx = d.component_index(v)

    # Common values used in all types
    ilo = d.ilo()
    ihi = d.ihi()

    lb, rb = d.boundaries()

    if bc == "periodic":
        # left boundary
        # Ghost cells are (right to left) the rightmost elements (same order)
        for i, b in enumerate(lb):
            shift = xmax - cells[ihi - (i + 1)].x()
            b.set_x(xmin - shift)
            b.set_value(idx, cells[ihi - (i + 1)].value(idx))

        # right boundary
        # Ghost cells are (left to right) the leftmost elements (same order)
        for i, b in enumerate(rb):
            shift = cells[(i + 1) + ilo].x() - xmin
            b.set_x(xmax + shift)
            b.set_value(idx, cells[(i + 1) + ilo].value(idx))

    elif bc == "outflow":
        # left boundary
        for i, b in enumerate(lb):
            # The shift mirrors the interior on the same side
            shift = cells[(i + 1) + ilo].x() - xmin
            b.set_x(xmin - shift)
            # Value same as leftmost interior
            b.set_value(idx, cells[ilo].value(idx))

        # right boundary
        for i, b in enumerate(rb):
            # The shift mirrors the interior on the same side
            shift = xmax - cells[ihi - (i + 1)].x()
            b.set_x(xmax + shift)
            # Value same as extrapolated from interior
            dy = cells[ihi].value(idx) - cells[ihi - 1].value(idx)
            dx = cells[ihi].x() - cells[ihi - 1].x()
            delta_x = cells[ihi + (i + 1)].x() - cells[ihi + i].x()

            b.set_value(idx, cells[ihi + i].value(idx) + (dy/dx) * delta_x)

    # Dictionary-based boundary conditions
    # Dirichlet and Neumann BCs require additional values also
    elif isinstance(bc, dict) and len(bc.keys()) == 1:
        bctype = list(bc.keys())[0]

        if bctype == "neumann":
            # left boundary
            neumann_value = bc[bctype][btype_map[btype.LEFT]]
            for i, b in enumerate(lb):
                # The shift mirrors the interior on the same side
                shift = cells[(i + 1) + ilo].x() - xmin
                b.set_x(xmin - shift)

                # Cell width to the left
                dx = cells[ilo - i].x() - cells[ilo - (i + 1)].x()
                b.set_value(
                    idx, cells[ilo - i].value(idx) - neumann_value * dx)

            # right boundary
            neumann_value = bc[bctype][btype_map[btype.RIGHT]]
            for i, b in enumerate(rb):
                # The shift mirrors the interior on the same side
                shift = xmax - cells[ihi - (i + 1)].x()
                b.set_x(xmax + shift)

                # Cell width to the right
                dx = cells[ihi + (i + 1)].x() - cells[ihi + i].x()
                b.set_value(
                    idx, cells[ihi + i].value(idx) + neumann_value * dx)

        elif bctype == "dirichlet":
            # left boundary
            dirichlet_value = bc[bctype][btype_map[btype.LEFT]]
            for i, b in enumerate(lb):
                # The shift mirrors the interior on the same side
                shift = cells[(i + 1) + ilo].x() - xmin
                b.set_x(xmin - shift)

                b.set_value(idx, dirichlet_value)

            # right boundary
            dirichlet_value = bc[bctype][btype_map[btype.RIGHT]]
            for i, b in enumerate(rb):
                # The shift mirrors the interior on the same side
                shift = xmax - cells[ihi - (i + 1)].x()
                b.set_x(xmax + shift)

                b.set_value(idx, dirichlet_value)
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError
