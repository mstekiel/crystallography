'''Mixin class for adding 2D/3D crystal and lattice visualization methods.'''

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from ..crystal.crystal import Crystal


class CrystalPlotterMixin:
    '''Mixin that adds plotting methods to Crystal or Lattice subclasses.

    Intended for 2D/3D visualization of crystal structures and lattice points.
    Methods delegate to the engine-specific SupercellPlotter implementations.

    TODO
    ----
    - plot_structure_2d(): project crystal structure onto a 2D plane
    - plot_structure_3d(): 3D ball-and-stick view
    - plot_lattice_points(): visualize reciprocal or direct lattice points
    '''
    pass
