'''Implementation of basic entities from group theory.
Limited to application of groups to symmetry elements
which allows for certain simplifications and forces conventions.

For example, I call the group elements operations, more natural in crystallography.
Also, each element has its inverse.
'''
import numpy as np
from collections import deque, defaultdict

import logging

from abc import ABC, abstractmethod
from typing import Any, Self, Iterable, Callable, Dict, Generic, List, Tuple, TypeVar

T = TypeVar('T', bound='SymOp')     # Holder for descendants of SymOp
A = TypeVar('A')                    # Holder for arbitrary objects to be symmetrized

logger = logging.getLogger('Symmetry')

MAX_GROUP_ORDER = 1024

class SymOp(ABC):
    '''Abstract class of symmetry operations as members of the crystallographic (SG) and magnetic space groups (MSG).
    Enforces to contain the implementation of functions required to form a SG or MSG:
    - multiplication: `__mul__`
    - string casting: `to_string()`, `from_string()`
    - inverse element: `inv`
    - identity element: `identity`
    - hashing and equality: `__hash__`, `__eq__`
    '''

    @classmethod
    @abstractmethod
    def from_string(cls, xyz_str: str) -> 'Self':
        '''Construct a symmetry operation from an xyz string.

        Contract
        --------
        - Must be the inverse of ``to_string()``: ``cls.from_string(op.to_string()) == op``.
        - Must use ``cls(...)`` so subclass identity is preserved.
        '''

    @abstractmethod
    def to_string(self) -> str:
        '''Serialize the operation to an xyz string.

        Contract
        --------
        - Must round-trip with ``from_string()``: ``cls.from_string(op.to_string()) == op``.
        - Must be a unique representation suitable for hashing and equality checks.
        '''

    @abstractmethod
    def __mul__(self, other: 'Self') -> 'Self':
        '''Compose two symmetry operations.

        Contract
        --------
        - Apply ``other`` first, then ``self`` (standard left-to-right composition).
        - Return ``self.__class__(...)``, not a hardcoded class name.
        '''

    @abstractmethod
    def inv(self) -> 'Self':
        '''Compute the inverse of this symmetry operation.

        Contract
        --------
        - Return ``self.__class__(...)``, not a hardcoded class name.
        - Must satisfy ``op * op.inv() == op.identity()``.
        '''

    @classmethod
    @abstractmethod
    def identity(cls) -> 'Self':
        '''Return the identity element of this symmetry operation class.

        Contract
        --------
        - Must use ``cls(...)`` so subclass identity is preserved.
        - Must satisfy ``op * op.identity() == op`` for all ``op``.
        '''

    @abstractmethod
    def __hash__(self) -> int:
        '''Hash of the symmetry operation.

        Contract
        --------
        - Must be consistent with ``__eq__``: ``a == b`` implies ``hash(a) == hash(b)``.
        - Base on raw internal numerical representation components (matrix ,trnaslation, ...).
        '''

    @abstractmethod
    def __eq__(self, other: object) -> bool:
        '''Equality of two symmetry operations.

        Contract
        --------
        - Must compare the same raw components as ``__hash__``.
        - Return ``NotImplemented`` if ``other`` is not the same class.
        '''


class Group(Generic[T]):
    '''Represent group, holder of group elemets `g` with functionalities.
    
    Attributes
    ----------
    #TODO
    operations: list[g]
        List of elements that make the group
    index_of: Dict[g, int])
        Gives the internal index (integer) which corresponds to the operation.
    mult_table: Dict[Tuple[g, g], g]
        Dictionary that holds the multiplication table of the group.
    adjacency: Dict[g, List[Tuple[g, g]]]
        Adjacency graph of the group. `adjacency[g1]=(g2, g3)`, where g1*g2=g3, 
        so `g2` is the label and `g3` is neighbouring vertex.
    adjacency_matrix: 
    
    TODO
    ----
    [ ] The constructor assumes the generators are unique, i.e. are not repeated
    '''
    _name: str
    _generators: tuple[T]
    _operations: tuple[T]
    _mult_table: Dict[Tuple[T, T], T]
    _adjacency: Dict[T, List[Tuple[str, T]]]
    _index_of:  Dict[T, int]

    def __init__(self, base: list[T], name='Group'):
        """Build the group from the list of base elements.

        1. The `base` is sanitized by finding only unique elements and removing identiy.
        2. The multiplication table is built, as well as adjacency graph and some other entities.
        3. The group generators are taken from the adjacency of identity element.
        """

        # Name
        self._name = name

        # Generate few attributes of the group that can be done in single loop
        # and thus enhances the speed and clarity.
        (operations, mult_table, adjacency, index_of) = self.build_group_cayley(base=base)
        self._operations = operations
        self._mult_table = mult_table
        self._adjacency  = adjacency
        self._index_of   = index_of

        # Take generators as elements adjacent to identity in the adj. graph
        identity = base[0].identity()
        gens = [g for edge, g in adjacency[identity]]
        self._generators = tuple(gens)

        # dict[label] -> n*n numpy arrays
        mats = self.adjacency_tensor(self._generators, self._operations, self._adjacency)
        self._adjacency_tensor = mats

        # Uncolored directed adjacency (counts parallel edges)
        n = len(operations)
        mat_list = [mat for mat in mats.values()]
        self._adjacency_matrix = np.add.reduce(mat_list).astype(int) if mat_list else np.zeros((n, n), dtype=int)

    ##########################################################################################################
    # Properties

    @property
    def generators(self) -> list[T]:
        '''List of all symmetry operations of the group.'''
        return self._generators

    @property
    def operations(self) -> list[T]:
        '''List of all symmetry operations of the group.'''
        return self._operations

    @property
    def order(self) -> int:
        '''Order of the group, that is number of operations in the group.'''
        return len(self._operations)

    ##########################################################################################################
    # Basic methods
    def index_of(self, g: T) -> int:
        '''Index of the group element, as stored in `self.operations`.'''
        return self._index_of[g]
    
    def __iter__(self):
        '''Iterate over the group elements.'''
        return iter(self._operations)

    def __repr__(self) -> str:
        gen_str = ", ".join([g.__repr__() for g in self._generators])
        return f"<{self.__class__.__name__}={self._name}, order={self.order}, gens=[{gen_str}]>"
    
    ##########################################################################################################
    # Advanced methods
    def symmetrize(self, obj: A, transform_func: Callable[[T,A], A], check_attrs: list[str]=[]) -> list[T]:
        '''Symmetrize the `object` of arbitrary type according to the `transform_func`
        within the symmetry of the `Group`.
        In other words, will apply each symmetry operation of the `Group` to the `object`,
        where the symmetry transformation rules are defined within the `transform_func`,
        and return only unique elements from the symmetrized elements.

        >>> class A: r=[0,0,0]
        >>> gen_fun = lambda g, a: A(g.matrix @ a.r)
        >>> MSG.symmetrize(A([0, 0, 0.5]), gen_fun) 
        
        Parameters
        ----------
        obj: T
            Object that will be symmetrized
        transform_obj: Callable[[T], T]
            Recipe how the object transforms under symmetry operations.
        check_attrs: list[str], optional
            If provided, it will check if the attributes are respecting the symmetry conditions of the MSG.

        Returns
        -------
        list[T]
            List of objects created by applying symmetry operations of the `MSG`.
        '''

        # This implementation assumes good hashing of the symmetrized objects
        objs_unique:      list[A]           = []
        objs_symmetrized: list[A]           = []

        # helper dict to assign unique ids to the symmetrized objects
        obj_to_uid:       dict[A, int]      = {}
        # helper dict to group the symmetry operations that produce 
        # the same symmetrized object, used for checking the symmetry conditions
        uid_to_eq_ids:    defaultdict[int, list[int]] = defaultdict(list)   

        for idx, g in enumerate(self.operations):
            obj_new = transform_func(g, obj)
            objs_symmetrized.append(obj_new)

            if obj_new not in obj_to_uid:
                obj_to_uid[obj_new] = len(objs_unique)
                objs_unique.append(obj_new)

            uid = obj_to_uid[obj_new]
            uid_to_eq_ids[uid].append(idx)

        if check_attrs:
            for id_unique, obj_unique in enumerate(objs_unique):
                id_equivalent = uid_to_eq_ids[id_unique]

                for check_field in check_attrs:
                    field_eq      = [getattr(objs_symmetrized[i], check_field) for i in id_equivalent]
                    field_averaged = np.average(field_eq, axis=0)
                    field_unique   = getattr(obj_unique, check_field)

                    if not np.allclose(field_unique, field_averaged):
                        msg  = f'Object property {check_field!r} does not respect symmetry:\n\t{obj_unique}\n'
                        msg += f'\tStored value     : {field_unique}\n'
                        msg += f'\tSymmetrized value: {field_averaged}\n'
                        for i in id_equivalent:
                            msg += f'\t{self.operations[i]}\t-> {check_field}={getattr(objs_symmetrized[i], check_field)}\n'
                        msg += 'You better know what you are doing.'
                        logger.warning(msg)

        return objs_unique

    def close_subset(self, base: list[T]) -> list[T]:
        """Generate a subgroup of `self` by closing the symmetry
         elements from `base` list.
         
         Uses the internal multiplication table for speed,
         thus can close only subgroups of self.

         Parameters
         ----------
         base: list[T]
            List of symmetry operations.
         """
        # DEV: this algorithm is basicaly snippet from the 
        # self.build_group_cayley

        assert set(base) <= set(self._operations)   # `base` operations must be subset of the group

        identity = base[0].identity()
        gens_inv = [g.inv() for g in base]
        
        # use Spm as the generating set for BFS
        steps_bfs = set(base + gens_inv) - {identity}

        # Degenerate case: only the identity exists
        if not steps_bfs:
            return [identity]
    
        # BFS discovery over the Cayley graph
        discovered: set[T]  = {identity}
        operations: list[T] = [identity]
        q = deque([identity])

        while q:
            if len(operations) > MAX_GROUP_ORDER:
                raise ValueError("Group appears larger than max_size (possible infinite subgroup).")

            x = q.popleft()
            for s in steps_bfs:
                # should it be LH or RH multiplicatoin?
                # y = x*s
                y = self._mult_table[(x, s)]

                # if new vertex, register and extend BFS frontier
                if y not in discovered:
                    discovered.add(y)
                    operations.append(y)
                    q.append(y)


        return operations
              
    @staticmethod
    def build_group_cayley(base: list[T]
                           ) -> tuple[tuple[T], 
                                      Dict[tuple[T,T], T], 
                                      Dict[T, list[tuple[T,T]]], 
                                      Dict[T, int]
                                      ]:
        """
        Generate the group from base elements using a Cayley-graph Breadth-First Search algorithm BFS,
        build the full multiplication table, and construct the Cayley graph.

        Parameters
        ----------
        base: tuple[T]
            List of base elements of the group, multiplication of which
            will form the full group.

        Returns
        -------
        tuple(operations, mult_table, adjacency, index_of)
        - operations : list of all operations (discovery order)
        - mult_table : dict mapping (a, b) -> a*b for all pairs a,b in operations
        - adjacency  : dict u -> list of (edge_label, v) for each generator/inverse
        - index_of   : map operation -> index to address rows/cols of the table
        """
        identity = base[0].identity()

        gens_inv = [g.inv() for g in base]
        
        # use Spm as the generating set for BFS
        steps_bfs = set(base + gens_inv) - {identity}
        steps_adj = set(base) - {identity}

        # Degenerate case: only the identity exists
        if not steps_bfs:
            elements = [identity]
            table = {(identity, identity): identity*identity}
            adj = defaultdict(list)
            return (elements, table, adj, {identity: 0})
    
        # BFS discovery over the Cayley graph
        discovered: set[T]  = {identity}
        operations: list[T] = [identity]
        index_of: Dict[T, int] = {identity: 0}
        q = deque([identity])

        # Cayley graph (directed, edge-labeled)
        adjacency: Dict[T, list[tuple[T,T]]] = defaultdict(list)

        # Full multiplication table, filled incrementally
        mult_table: Dict[tuple[T,T], T] = {}
        mult_table[(identity, identity)] = identity*identity

        # Helper: when a new element y is discovered, fill its row/column vs all known elements
        def fill_table_for_new(y: Any) -> None:
            # row y,* and column *,y vs all existing elements (which exclude y)
            for b in operations:
                mult_table[(y, b)] = y*b
                mult_table[(b, y)] = b*y
            # finally y*y
            mult_table[(y, y)] = y*y

        while q:
            if len(operations) > MAX_GROUP_ORDER:
                raise ValueError("Group appears larger than max_size (possible infinite subgroup).")

            x = q.popleft()
            for s in steps_bfs:
                # should it be LH or RH multiplicatoin?
                y = x*s

                # record labeled edge x --(s)--> y in the Cayley graph
                if s in steps_adj:
                    adjacency[x].append((s, y))

                # if new vertex, register and extend BFS frontier
                if y not in discovered:
                    discovered.add(y)
                    fill_table_for_new(y)
                    index_of[y] = len(operations)
                    operations.append(y)
                    q.append(y)


        operations = tuple(operations)
        return (operations, mult_table, adjacency, index_of)
    
    @staticmethod
    def adjacency_tensor(
        generators: List[T],
        operations: List[T],
        adjacency: Dict[T, List[Tuple[str, T]]],
    ) -> Dict[T, np.ndarray]:
        """
        Build a per-generator adjacency tensor for the Cayley diagram of the group.

        Returns a dict keyed by the *generator elements themselves* (not strings):
            { g -> A_g }  where A_g is an n*n numpy array (int)
        and (A_g)[i, j] counts directed edges operations[i] --(g)--> operations[j].

        Parameters
        ----------
        generators : Iterable[T]
            The base generating set S (identity excluded). The tensor is indexed by these objects.
        operations : list[T]
            All group elements, their order will define the row/column order of the matrices.
        adjacency : Dict[T, List[Tuple[str, T]]]
            Directed, labeled edges of the Cayley graph: u -> [(label, v), ...]

        Returns
        -------
        Dict[T, np.ndarray]
            One n*n matrix per generator g in S. Edges labeled by g or g^{-1} are folded into A_g.
        """
        n = len(operations)
        idx = {g: i for i, g in enumerate(operations)}

        mats: Dict[T, np.ndarray] = {g: np.zeros((n, n), dtype=int) 
                                     for g in generators}

        for g1, connections in adjacency.items():
            for gc,g2 in connections:
                if gc in mats:
                    mats[gc][idx[g1],idx[g2]] += 1

        return mats
    
    def get_subgroups(self, contains: Iterable[T] = None):
        """Find all subgroups of space group, by creating cyclic subgroups
        of all group operations and then extending them with leftover operations."""  

        # set needs to be freezed to allow hashing it
        cyclic_subgroups = set(frozenset(self.close_subset([g])) 
                               for g in self.operations)
        
        # TODO here filter only those cyclic cubgroups 
        # that contain the required symmetry element
        if contains is None:
            subgroups = cyclic_subgroups
        else:
            subgroups = set(G for G in cyclic_subgroups if set(contains) <= G)


        logger.debug(subgroups)

        n = self.order
        q = list(subgroups) # creates a safe copy
        while q:
            H = q.pop()
            if len(H) == n:
                continue

            # Record which elements were already covered in extension,
            # to avoid rediscovery
            covered = set(H)
            for g in self.operations:
                if g in covered:
                    continue

                new_gens = list(H | {g})
                K = frozenset(self.close_subset(new_gens))

                # Add extended group to the covere elements set
                covered |= K
                if K not in subgroups:
                    subgroups.add(K)
                    q.append(K)

        ### TODO Return
        # Right now, the group constructor creates multiplication
        # table from scratch, which is relatively slow.
        # For now, keep returning just set of symops
        # but i nthe future, the constructor should allow setting all
        # internal fields.

        # 1.
        return sorted(subgroups, key=len)

        # 2. TODO
        # s = []
        # for H in subgroups:
        #     print("New group from: ", H)
        #     new_H = deepcopy(self)
        #     new_H.__init__(base=list(new_H), name=self._name+"_sub")
        #     s.append(new_H)


        # return s
        # return [self.__init__(base=list(Gi), name=f"self._name}_sub{n}")
        #         for n,Gi in enumerate(subgroups)]
    




def plot_network(group, layout: str='spring', seed:int=0, node_size: float=100, font_size: int=8):
    """Export Cayley digraph as a NetworkX MultiDiGraph with string node ids."""
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError("Either `matplotlib` or `networkx` are not installed.") from e


    G = nx.DiGraph(name=group._name)

    # nodes: use string ids; keep useful attrs
    for op in group._operations:
        nid = op.to_string()
        G.add_node(nid, index=group._index_of[op], op_str=nid, op=op)

    # edges: label may be SymOp or str; store as 'generator'
    for u, edges in group._adjacency.items():
        u_id = u.to_string()
        for label, v in edges:
            if label in group._generators:
                v_id = v.to_string()
                G.add_edge(u_id, v_id, generator=label.to_string())



    # positions
    if layout == "spring":
        pos = nx.spring_layout(G, seed=seed)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = nx.spring_layout(G, seed=seed)

    # nodes
    nx.draw_networkx_nodes(G, pos, node_size=node_size)
    nx.draw_networkx_labels(G, pos, font_size=font_size)

    # assign distinct colors per generator label
    gens = sorted({d["generator"] for _, _, d in G.edges(data=True)})
    import itertools
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    color_of = {g: c for g, c in zip(gens, itertools.cycle(palette))}

    # draw edges per generator with its color
    for g in gens:
        edgelist = [(u, v) for u, v, d in G.edges(data=True) if d["generator"] == g]
        nx.draw_networkx_edges(G, pos, edgelist=edgelist, edge_color=color_of[g], arrows=True, width=3)

    # edge labels (works for DiGraph)
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels={(u, v): d["generator"] for u, v, d in G.edges(data=True)},
        font_size=font_size, alpha=0.8
    )
    plt.axis("off")
    plt.show()

