# -*- coding: utf-8 -*-
r"""Tests lattice math

"""
import numpy as np
from fractions import Fraction
import pytest

from crystallography import Crystal, MSG, Atom
from crystallography.utils.linalg import DMI


class _CouplingStub:
    """Minimal coupling stub for testing crystal.is_respectful_DMI.
    Coupling itself lives in a separate library; this covers only the
    interface that Crystal needs: id1, id2, n_uvw, DMI_vector.
    """
    def __init__(self, label: str, id1: int, id2: int, n_uvw, J: np.ndarray):
        self.label = label
        self.id1 = id1
        self.id2 = id2
        self.n_uvw = np.array(n_uvw, dtype=float)
        self.J = J

    @property
    def DMI_vector(self) -> np.ndarray:
        # Extract DM vector from antisymmetric J using bragg.utils.linalg.DMI convention:
        # DMI([dx,dy,dz]) = [[0, dz,-dy],[-dz,0,dx],[dy,-dx,0]]
        # so dx=J[1,2], dy=J[2,0], dz=J[0,1]
        return np.array([self.J[1, 2], self.J[2, 0], self.J[0, 1]])



def test_crystal_constructor():
    """Test construction of the lattice"""

    atoms = [ Atom(r=[0,0,1/4]), Atom(r=[1/3, 2/3, 3/4])]
    P_194 = MSG.from_xyz_strings(generators=[
        '-y,x-y,z, +1',
        '-x,-y,z+1/2, +1',
        'y,x,-z, +1',
        '-x,-y,-z, +1',
    ])

    print(P_194)
    a, c = 3.6, 12
    graphite = Crystal(lattice_parameters=[a,a,c, 90,90,120],
                      atoms=atoms,
                      MSG=P_194)

    print(graphite)


def test_get_atom_index():
    atoms = [ 
        Atom(label='A1', r=[0,0,1/4], m=[0,0,1], s=1/2),
        Atom(label='A2', r=[1/3, 2/3, 3/4], m=[0,0,1], s=1/2)
        ]
    P_194 = MSG.from_xyz_strings(generators=[
        '-y,x-y,z, +1',
        '-x,-y,z+1/2, -1',
        'y,x,-z, +1',
        '-x,-y,-z, -1',
    ])

    mag_graphite = Crystal(lattice_parameters=[3.6,3.6,12, 90,90,120],
                      atoms=atoms,
                      MSG=P_194)   

    assert mag_graphite.get_magatom_id([0, 0, 0.25]) == 0
    assert mag_graphite.get_magatom_id([0, 0, 0.75]) == 1
    assert mag_graphite.get_magatom_id([1/3, 2/3, 0.75]) == 2
    assert mag_graphite.get_magatom_id([2/3, 1/3, 0.25]) == 3

    with pytest.raises(LookupError) as e_info:
        mag_graphite.get_magatom_id([1/3, 1/3, 0.25]) == 3


def test_constructor_validators():
    # 1. Magnetic moment should obey MSG
    pass

def test_crystal_couplings():
    # Atoms in base plane, with funny moments
    atoms = [
        Atom(label='Dy', r=(0.5,   0.5, 0),   m=(0,0,5), s=2.5),
        Atom(label='Fe', r=(0.1,   0.1, 0),   m=(1,1,0), s=2.5)]
    # atoms = [Atom(label='Fe', r=(0,   0.5, 0),   m=(-1,0,Fz), s=2.5)]
    P4mm = MSG.from_xyz_strings(generators=[
        '-y, x, z, +1', # 4_001
        'y, x, z, +1'   # m_110
    ])

    crystal = Crystal(lattice_parameters=(4,4,10, 90,90,90),
                      MSG=P4mm, atoms=atoms)
    
    print(crystal.MSG.get_point_symmetry([0.1,0.1,0]))
    
    D1 = DMI([2,2,0])
    DMI_appropriate = _CouplingStub(label='D1', id1=0, id2=2, n_uvw=[0,0,0], J=D1)
    assert crystal.is_respectful_DMI(DMI_appropriate)

    D2 = DMI([0.5,0.05,0.005])
    DMI_inappropriate = _CouplingStub(label='D2', id1=0, id2=2, n_uvw=[0,0,0], J=D2)
    inappropriate, D2_symmetrized = crystal.is_respectful_DMI(DMI_inappropriate, return_symmetrized=True)
    
    assert not inappropriate
    assert np.allclose(D2_symmetrized, [0.275, 0.275, 0.005])
    
if __name__ == "__main__":
    # pytest.main()
    test_get_atom_index()


    ### Quick tests
    # import time
    # t_start = time.time()
    # P = MSG.from_xyz_strings(generators=[
    #     # '-y, x, z, +1', # 4_001
    #     # 'z, x, y, +1',    # 3_111
    #     # 'x+1/2, y+1/2, z+1/2, +1', # I centering
    #     '-x, -y, z, +1',    # 2_110
    #     # '-x, -y, -z, +1',   # -1
    #     # 'x, y, z, -1',   # 1'
    #     'y, x, z, +1'   # m_110
    #     ])
    # t_end = time.time()
    # for g in P:
    #     print(g.print())

    # print(P.order)
    # print(P.make_cayley_table())
    # print('Constructor and Cayley table', t_end-t_start)
    # print()
    # t_start = time.time()
    # P_subs = P.get_subgroups()
    # for P_sub in P_subs:
    #     print(P_sub)
    # t_end = time.time()
    # print('subgroups', t_end-t_start)
    # print('orders', [G.order for G in P_subs])

    # pos = [0, 0.5, 0.3]
    # print('point symmetry')
    # print(P.get_point_symmetry(pos))
    # print('orbit')
    # print(P.get_orbit(pos))