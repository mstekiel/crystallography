'''Plotting functionalities for crystal structures and lattices.'''

from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from ..crystal.crystal import Crystal


IMPLEMENTED_SC_ENGINES = ['mpl']
# TODO: implement remaining engines
# PLANNED_SC_ENGINES = ['vispy', 'vispy+', 'qtgraph']


def plot_structure(crystal: 'Crystal',
                   engine: str='mpl',
                   plot_options: dict={}) -> Any:
    '''Render the crystal structure.

    Parameters
    ----------
    crystal: Crystal
        Crystal structure to visualize.
    engine: str
        Rendering library. Currently only 'mpl' (matplotlib) is implemented.
    plot_options: dict
        Additional options passed to the plotter.

    Returns
    -------
        Library-specific objects handling the plot widget.
    '''
    if engine not in IMPLEMENTED_SC_ENGINES:
        raise NotImplementedError(
            f"Engine {engine!r} is not implemented. "
            f"Available engines: {IMPLEMENTED_SC_ENGINES}"
        )

    if engine == 'mpl':
        from .supercell_plotter_mpl import MPLSupercellPlotter as SCPlotter

    plotter = SCPlotter(crystal, plot_options=plot_options)

    return plotter.deploy()
