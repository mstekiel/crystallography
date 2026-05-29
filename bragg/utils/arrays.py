import inspect
from functools import wraps
from typing import Iterable
import numpy as np

##################################################################################################
# Shape matching
def match_shape(arr_shape: tuple[int, ...], expected_shape: tuple[int, ...]) -> bool:
    """
    Recursively match `arr_shape` against `expected_shape` with ... wildcards.
    Supports multiple ellipses anywhere.
    """
    # print(f'DEBUG: matching: {arr_shape=} with {expected_shape=}')
    if not expected_shape:
        return not arr_shape

    head, *tail = expected_shape

    if head is ...:
        # try all possible splits: let ... absorb 0,1,2,... dimensions
        for k in range(len(arr_shape) + 1):
            if match_shape(arr_shape[k:], tuple(tail)):
                return True
        return False

    if not arr_shape or arr_shape[0] != head:
        return False

    return match_shape(arr_shape[1:], tuple(tail))

def ensure_shape(**shapes):
    """
    Decorator to validate numpy array arguments by name and shape.

    Usage
    -----
    @ensure_shape(
        a=(..., 3),          # last dimension must be 3
        b=(2, ...),          # first dimension must be 2
        c=(2, 3, ..., 3, 3)  # fixed prefix and suffix, flexible middle
    )
    def func(a, b, c): ...
    """
    def decorator(func):
        sig = inspect.signature(func)

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            for name, expected_shape in shapes.items():
                if name not in bound.arguments:
                    raise ValueError(f"Parameter '{name}' not found in function arguments")
                
                val = bound.arguments.get(name)
                if val is None:
                    continue

                arr_shape = np.shape(val)
                if not match_shape(arr_shape, expected_shape):
                    raise ValueError(
                        f"[{func.__name__}] Parameter '{name}': "
                        f"expected shape {expected_shape}, got {arr_shape}"
                    )

            return func(*args, **kwargs)

        return wrapper
    
    return decorator

##################################################################################################
# Array creation

def create_mesh(*vs: Iterable) -> np.ndarray:
    '''Create a coordinate mesh from N arrays or scalars.

    Each input can be a 1D array or a scalar. Returns an array of shape
    `(*sizes, N)` where `result[i, j, ..., :]` holds the coordinates
    of that grid point. Scalar inputs are broadcast and squeezed out of
    the spatial dimensions.
    
    >>> m = create_mesh([1,2], [3,4], 0)
    array([[[1, 3, 0],
            [1, 4, 0]],

           [[2, 3, 0],
            [2, 4, 0]]])
    >>> m.shape
    (2, 2, 3)
    '''

    arrays = [np.atleast_1d(v) for v in vs]
    grids = np.meshgrid(*arrays, indexing='ij')
    return np.stack(grids, axis=-1).squeeze()