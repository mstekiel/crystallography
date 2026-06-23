import pytest

from crystallography.symmetry.spatial import cSymOp, SG

import numpy as np

from . import make_perm_symop, s3_generators, c4_generator



# --------------------------
# Tests
# --------------------------

@pytest.mark.skip(reason="manual")
def hand_test():
    e = cSymOp.identity()
    print(e)
    g1 = cSymOp.from_string('x,y+1/2,z+1/2')
    g2 = cSymOp.from_string('x+1/2,y,z+1/2')
    g3 = cSymOp.from_string('-y,x,z')
    g4 = cSymOp.from_string('-x,-y,-z')
    g5 = cSymOp.from_string('y,x,z')
    g6 = cSymOp.from_string('y,z,x')
    SG = Group([
            cSymOp.from_string('x,y+1/2,z+1/2'),
            # cSymOp.from_string('x+1/2,y,z+1/2'),
            cSymOp.from_string('-y,x,z'),
            cSymOp.from_string('-x,-y,-z'),
            cSymOp.from_string('y,x,z'),
            cSymOp.from_string('y,z,x'),
    ])
    print(SG)
    # print(SG._adjacency_matrix)
    # print(SG._adjacency_tensor)
    print(SG.operations)

    print(g2, SG.index_of(g2))



if __name__ == '__main__':
    pytest.main([__file__, '-v'])
