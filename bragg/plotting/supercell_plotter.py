from abc import ABC, abstractmethod
from typing import Any, Tuple, Union

import numpy as np
from copy import deepcopy


from ..crystal.atom import Atom
from ..crystal.crystal import Crystal

import logging
import logging.config
import traceback
sc_logger = logging.getLogger('SupercellPlotter')

#######################################################################################################
# Tools
def _format_bbox(boundaries: Union[int, list[int]]):
    bbox = np.zeros((3,2))
    if np.shape(boundaries) == ():
        bbox[:,1] = boundaries
    elif np.shape(boundaries) == (3,):
        bbox[:,1] = boundaries
    elif np.shape(boundaries) == (3,2):
        bbox = np.array(boundaries)
    else:
        raise IndexError(f'Unexpected dimension of the boundary box: {np.shape(boundaries)} not in [(),(3,),(3,2)]')

    return np.sort(bbox, axis=1)


#######################################################################################################

class SupercellPlotter(ABC):
    '''Plots the crystal structure (atoms, cell edges, magnetic moments) in 2D or 3D.

    All distances are in Angstroms.
    '''
    def __init__(self, crystal: Crystal):
        self.crystal = crystal
        self.logger = sc_logger

    @abstractmethod
    def plot_balls(self,
                   positions: np.ndarray,
                   sizes: np.ndarray,
                   colors: np.ndarray):
        '''Plot balls.
        Used for drawing atoms.

        Parameters
        ----------
        positions: (n,3)
        sizes: (n,)
        colors: (n,3)

        Returns
        -------
        objects: (n,)
        '''
        pass

    @abstractmethod
    def plot_lines(self,
                   lines: np.ndarray,
                   colors: np.ndarray):
        '''Plot lines.
        Used for drawing cell edges and lattice points.

        Parameters
        ----------
        lines: (n,2,3)
        colors: (n,3)

        Returns
        -------
        objects: (n,)
        '''
        pass

    @abstractmethod
    def plot_arrows(self,
                    positions: np.ndarray,
                    directions: np.ndarray,
                    colors: np.ndarray):
        '''Plot arrows.
        Used for drawing magnetic moments.

        Parameters
        ----------
        positions: (n,3)
        directions: (n,3)
        colors: (n,3)

        Returns
        -------
        objects: (n,)
        '''
        pass

    @abstractmethod
    def plot_labels(self,
                    positions: np.ndarray,
                    labels: np.ndarray,
                    colors: np.ndarray):
        '''Plot text labels.

        Parameters
        ----------
        positions: (n,3)
        labels: (n,)
        colors: (n,3)

        Returns
        -------
        objects: (n,)
        '''
        pass

    @abstractmethod
    def plot_ellipsoids(self,
                        positions: np.ndarray,
                        matrices: np.ndarray,
                        colors: np.ndarray):
        '''Plot ellipsoids.

        Parameters
        ----------
        positions: (n,3)
        matrices: (n,3,3)
        colors: (n,3)

        Returns
        -------
        objects: (n,)
        '''
        pass

    @abstractmethod
    def deploy(self):
        '''Library-specific routines to deploy the plot.

        Returns
        -------
        widget:
            Main widget used to control the plot window.
        '''
        pass

    #######################################################################################################

    def plot(self, plot_options: dict={}):
        '''Main plotting function'''

        boundaries = [[0,1],[0,1],[0,1]]
        if 'boundaries' in plot_options:
            boundaries = plot_options['boundaries']

        self.atom_alpha = plot_options.pop('atom_alpha', 0.8)
        self.atom_scale = plot_options.pop('atom_scale', 1)
        self.spin_scale = plot_options.pop('spin_scale', 1)
        self.arrow_width = plot_options.pop('arrow_width', 0.1)
        self.arrow_head_size = plot_options.pop('arrow_head_size', 3)

        bbox = _format_bbox(boundaries)
        self.logger.info(f"Plotting stuff in bbox: {bbox}")

        atoms, edges = self.get_objects_in_supercell(bbox)

        self.logger.info(f"Plotting atoms: {atoms}")

        # Atoms
        pos = self.crystal.uvw2xyz([atom.r for atom in atoms])
        sizes = np.array([atom.radius for atom in atoms])
        colors = np.array([atom.color for atom in atoms])
        try:
            self.plot_balls(positions=pos, sizes=sizes*self.atom_scale, colors=colors)
        except Exception as e:
            self.logger.error(traceback.format_exc())

        self._structure_center = np.average(pos, axis=0)
        self._largest_distance = np.abs(pos-self._structure_center).max()

        # Magnetic moments
        magnetic_atoms = [atom for atom in atoms if atom.is_mag]
        ma_r = self.crystal.uvw2xyz([atom.r for atom in magnetic_atoms])
        ma_m = np.array([atom.m*atom.s*self.spin_scale for atom in magnetic_atoms])
        ma_colors = np.array([atom.color for atom in magnetic_atoms])
        try:
            self.plot_arrows(positions=ma_r, directions=ma_m, colors=ma_colors)
        except Exception as e:
            self.logger.error(traceback.format_exc())

        # Cell edges
        colors = np.array([[1,1,1] for _ in edges]) # black
        edges = np.array([self.crystal.uvw2xyz(np.array(edge))
                          for edge in edges])
        try:
            self.plot_lines(edges, colors, alpha=0.25, width=1)
        except Exception as e:
            self.logger.error(traceback.format_exc())

        # Zeroth cell has colorful axes
        main_edges = np.array([
            [[0,0,0],[1,0,0]],
            [[0,0,0],[0,1,0]],
            [[0,0,0],[0,0,1]]
        ])
        main_edges = np.array([self.crystal.uvw2xyz(np.array(edge))
                               for edge in main_edges])
        colors = np.eye(3)  # list of R G B
        try:
            self.plot_lines(main_edges, colors, width=4)
        except Exception as e:
            self.logger.error(traceback.format_exc())

        return

    def get_objects_in_supercell(self, boundaries) -> tuple[list[Atom], np.ndarray]:
        '''Find all atoms and cell edges within the given boundaries.

        Returns
        -------
        atoms: list[Atom]
        edges: np.ndarray
        '''
        EPS = 1e-8
        ext_atoms = np.floor(boundaries).astype(int)
        ext_edges = np.trunc(boundaries).astype(int)

        ### EDGES
        edges = []
        it_nx = range(ext_edges[0][0], ext_edges[0][1]+1)
        it_ny = range(ext_edges[1][0], ext_edges[1][1]+1)
        it_nz = range(ext_edges[2][0], ext_edges[2][1]+1)

        for ny in it_ny:
            for nz in it_nz:
                edges.append([[ext_edges[0][0],ny,nz], [ext_edges[0][1],ny,nz]])
        for nx in it_nx:
            for nz in it_nz:
                edges.append([[nx,ext_edges[1][0],nz], [nx,ext_edges[1][1],nz]])
        for nx in it_nx:
            for ny in it_ny:
                edges.append([[nx,ny,ext_edges[2][0]], [nx,ny,ext_edges[2][1]]])

        edges = np.array(edges, dtype=float)

        ### ATOMS
        atoms = []
        it_nx = range(ext_atoms[0][0], ext_atoms[0][1]+1)
        it_ny = range(ext_atoms[1][0], ext_atoms[1][1]+1)
        it_nz = range(ext_atoms[2][0], ext_atoms[2][1]+1)
        low_bound = boundaries[:,0] - EPS
        high_bound = boundaries[:,1] + EPS
        for zcen in it_nz:
            for ycen in it_ny:
                for xcen in it_nx:
                    for atom in self.crystal.atoms_all:
                        atom_candidate = deepcopy(atom)
                        n_uvw = np.array([xcen, ycen, zcen])
                        atom_candidate.r += n_uvw

                        if all(atom_candidate.r > low_bound) and all(atom_candidate.r < high_bound):
                            atoms.append(atom_candidate)

        return atoms, edges
