import numpy as np
from fractions import Fraction
from crystallography.symmetry import mSymOp, Group, MSG

import logging
logging.getLogger().setLevel(logging.DEBUG)

class MPG(Group[mSymOp]):
    '''Magnetic Point Group'''
    _name: str
    _generators: tuple[mSymOp]
    _operations: tuple[mSymOp]


    ##############################################################################
    # Constructors
    @classmethod
    def from_xyz_strings(cls, base: list[str]):
        '''Construct Magnetic Point Group from a list of xyz_strings,
        that represent MPG operations.
        Any translation in symmetry operation will be ignored.
        '''
        
        generators_sanitized = list()
        for gs in base:
            g = mSymOp.from_string(gs)
            if not np.allclose(g.translation, (0,0,0)):
                # #TODO will be logger in future
                print(f"WRN: Setting translation to 0 for: {gs}")

                g = mSymOp(matrix = g._matrix,
                           translation = np.array([0,0,0], dtype=Fraction),
                           time_reversal=g._time_reversal)
                
            generators_sanitized.append(g)

        return cls(base = generators_sanitized)

def todo():
    Imm3 = MSG.from_xyz_strings(base=[
        'z, x, y, +1',    # 3_111
        'x+1/2, y+1/2, z+1/2, +1', # I centering
        '-x, -y, z, +1',    # 2_110
        '-x, -y, -z, +1',   # -1
        ])
    # This groups has order 48, check only orders of subgroups

    G_i = Imm3.get_subgroups()
    for G_ii in G_i:
        print(G_ii)

    # subgroups_order = [2,4,4,4,6,8,8,8,12,16,24]
    # subgroups = Imm3.get_subgroups()
    # assert sorted(subgroups_order) == sorted([Gsub.order for Gsub in subgroups])


def main():
    m_e = "x, y, z, +1"
    m_3_111 = "z, x, y, +1"
    m_mx = "-x, y, z, +1"
    m_4  = "-y, x, z, +1"
    MPG_mm2 = MPG.from_xyz_strings([m_mx, m_4, m_3_111])

    print(MPG_mm2)
    print(MPG_mm2.operations)
    print("generators? = ", MPG_mm2._adjacency[mSymOp.from_string(m_e)])

    # 3. Time one full closure of the whole group
    # import time
    # t = time.perf_counter()
    # K = MPG.build_group_close([mSymOp.from_string(gs) for gs in [m_mx, m_4, m_3_111]])
    # print("closure of 3 gens:", len(K), "in", time.perf_counter() - t, "s")


    G_i = MPG_mm2.get_subgroups(contains=[mSymOp.from_string(m) for m in [m_e]])
    for n, G_ii in enumerate(G_i):
        print(n, G_ii, "\n")

# mult table 
# Min = 0.0003448000061325729, max=0.11000779998721555, mean = 0.022805130357820808

# build_close
# Min = 0.000450699997600168, max=0.06615349999628961, mean = 0.009198039854752978

# fast_close
# Min = 0.0008243999909609556, max=0.03157410002313554, mean = 0.004614717755390777

# close subset

if __name__ == '__main__':
    main()
    # todo()