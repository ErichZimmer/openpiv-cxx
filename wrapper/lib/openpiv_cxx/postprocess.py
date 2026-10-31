from __future__ import annotations

from numbers import Integral, Real
from typing import TypeAlias

import numpy as np

from . import pyopenpivcore as cpiv
from .parsing import convert_vector_to_grid_data_type


__all__ = [
    'substitute_nans', 
    'replace_nans'
]


CppOutput: TypeAlias = tuple[cpiv.grid_data, int | bool]
NumpyOutput: TypeAlias = tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    int | bool,
]
PostprocessOutput: TypeAlias = CppOutput | NumpyOutput
SubstituteOutput: TypeAlias = (
    tuple[cpiv.grid_data, int] |
    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]
)


def _prepare_field(
    u: np.ndarray | cpiv.grid_data | None,
    v: np.ndarray | None,
    field_data: cpiv.grid_data | None,
    arrays: dict[str, np.ndarray | None],
    require_secondary: bool=False,
) -> cpiv.grid_data:
    if isinstance(u, cpiv.grid_data):
        if field_data is not None:
            raise ValueError('field_data was provided more than once')
        
        field_data = u
        u = None

    if field_data is not None:
        supplied = [
            name for name, value in {'u': u, 'v': v, **arrays}.items()
            if value is not None
        ]

        if supplied:
            names = ', '.join(supplied)

            msg = 'provide either field_data or NumPy field arrays; do not ' +\
                 f'mix them (received arrays: {names})'
            raise ValueError(msg)
        
        if not isinstance(field_data, cpiv.grid_data):
            raise TypeError('field_data must be of type grid_data')
        
        return field_data

    required = {'u': u, 'v': v}
    if require_secondary:
        required.update(u2=arrays.get('u2'), v2=arrays.get('v2'))

    missing = [name for name, value in required.items() if value is None]
    if missing:
        msg = f'NumPy field input requires {", ".join(required)} arrays '+\
              f'(missing: {", ".join(missing)})'
        raise ValueError(msg)

    data = convert_vector_to_grid_data_type(u, v, **arrays)
    for name in ('s2n', 'p2p', 'peak'):
        if arrays.get(name) is None:
            np.array(getattr(data, name), copy=False)[:] = np.nan

    if arrays.get('flag') is None:
        flag = np.array(data.flag, copy=False)
        flag[:] = int(cpiv.FLAG.VALID)
        invalid = (
            ~np.isfinite(np.array(data.u, copy=False)) |
            ~np.isfinite(np.array(data.v, copy=False))
        )
        flag[invalid] = int(cpiv.FLAG.INVALID)

    return data


def substitute_nans(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    s2n: np.ndarray | None=None,
    peak: np.ndarray | None=None,
    u2: np.ndarray | None=None,
    v2: np.ndarray | None=None,
    s2n2: np.ndarray | None=None,
    peak2: np.ndarray | None=None,
    flag: np.ndarray | None=None,
    field_data: cpiv.grid_data | None=None,
    parse_output: bool=True,
) -> SubstituteOutput:
    '''Replace invalid vectors with the next stored candidate.

    Substitute invalid vectors with secondary peak information. This allows
    for possibly correct vectors to be recovered.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array. A ``cpiv.grid_data`` object may be passed here
        as the first positional argument instead.
    v : np.ndarray or None
        V-displacement array. Required when ``u`` is a NumPy array.
    s2n : np.ndarray or None
        Optional signal-to-noise ratio array. Omitted values become NaN.
    peak : np.ndarray or None
        Optional correlation-peak value array. Omitted values become NaN.
    u2 : np.ndarray or None
        Secondary u-displacement array. Required for NumPy inputs.
    v2 : np.ndarray or None
        Secondary v-displacement array. Required for NumPy inputs.
    s2n2 : np.ndarray or None
        Optional secondary signal-to-noise ratio array. Omitted values
        remain NaN.
    peak2 : np.ndarray or None
        Optional secondary correlation-peak value array. Omitted values
        remain NaN.
    flag : np.ndarray or None
        Optional vector-status flag array, including MASKED locations. If
        omitted, finite u/v are VALID and nonfinite u/v are INVALID.
    field_data : cpiv.grid_data or None
        Reusable C++ field, updated in place. Cannot be mixed with arrays.
    parse_output : bool
        If True, return independent NumPy copies. If False, return the same
        grid_data object with the unsuccessful-substitution count.

    Returns
    -------
    u : np.ndarray
        U-displacement copy, returned when parse_output is True.
    v : np.ndarray
        V-displacement copy, returned when parse_output is True.
    s2n : np.ndarray
        Signal-to-noise ratio copy, returned when parse_output is True.
    peak : np.ndarray
        Correlation-peak value copy, returned when parse_output is True.
    flag : np.ndarray
        Independent uint8 copy of the vector-status flags, returned when
        parse_output is True: 0 valid, 1 invalid, 2 substituted,
        3 interpolated, and 4 masked.
    field_data : cpiv.grid_data
        Updated field when parse_output is False, instead of the five arrays.
    failed_substitutions : int
        Number of INVALID locations not replaced on this call. This is the
        last returned item in either output mode.

    Raises
    ------
    TypeError
        If an input has an unsupported type or parse_output is not a bool.
    ValueError
        If input modes are mixed, required arrays are missing, or array shapes
        differ or are not two-dimensional.
    RuntimeError
        If the underlying substitution routine fails.

    '''
    if not isinstance(parse_output, bool):
        raise TypeError('parse_output must be a bool')
    
    data = _prepare_field(
        u, v, field_data,
        {
            's2n': s2n, 
            'peak': peak, 
            'flag': flag,
            'u2': u2, 
            'v2': v2, 
            's2n2': s2n2, 
            'peak2': peak2,
        },
        require_secondary=True,
    )
    failed_substitutions = cpiv.substitute_nans(data)

    if not parse_output:
        return data, failed_substitutions

    return (
        np.array(data.u, copy=True),
        np.array(data.v, copy=True),
        np.array(data.s2n, copy=True),
        np.array(data.peak, copy=True),
        np.array(data.flag, copy=True),
        failed_substitutions,
    )


def replace_nans(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    s2n: np.ndarray | None=None,
    p2p: np.ndarray | None=None,
    peak: np.ndarray | None=None,
    flag: np.ndarray | None=None,
    peak_ind: np.ndarray | None=None,
    field_data: cpiv.grid_data | None=None,
    include_corners: bool=False,
    tolerance: float=0.0001,
    max_iterations: int=10,
    parse_output: bool=True,
) -> PostprocessOutput:
    '''Impute invalid vectors using the spring-leaf analogy.

    Interpolate missing vectors using the leaf-spring analogy through an
    iterative solver. Isolated vectors are effectively replaced via the mean
    of their neighbors. Large clusters of invalid vectors are replaced by
    iteratively distributing information from valid vectors to invalid ones.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array, or grid_data supplied positionally.
    v : np.ndarray or None
        V-displacement array. Required for NumPy inputs.
    s2n, peak : np.ndarray or None
        Optional measurement arrays. Omitted arrays become NaN.
        Imputed locations receive NaN in both arrays.
    p2p : np.ndarray or None
        Optional measurement peak ratio, preserved during imputation. If
        omitted in NumPy mode, it is initialized to NaN.
    flag : np.ndarray or None
        Optional status codes, including MASKED locations. In NumPy mode, if
        omitted, finite u/v are VALID and nonfinite u/v are INVALID.
    peak_ind : np.ndarray or None
        Optional zero-based candidate indices, preserved during imputation.
    field_data : cpiv.grid_data or None
        Reusable C++ field, updated in place. Cannot be mixed with arrays.
    include_corners : bool
        Include diagonal neighbors when solving for outlier values.
    tolerance : float
        Finite, nonnegative convergence tolerance.
    max_iterations : int
        Positive iteration limit. Higher is needed for large outlier clusters.
    parse_output : bool
        If True, return independent NumPy copies. If False, return the same
        grid_data object with the convergence result.

    Returns
    -------
    u, v, s2n, p2p, peak, flag, peak_ind : np.ndarray
        Independent array copies, in this order, when parse_output is True.
    field_data : cpiv.grid_data
        Updated field when parse_output is False, instead of the seven arrays.
    converged : bool
        Whether the convergence criterion was met. This is the last
        returned item in either output mode.

    Raises
    ------
    TypeError
        If an array, field, or option has an unsupported type.
    ValueError
        If input modes are mixed, required arrays are missing, or array shapes
        differ or are not two-dimensional, or imputation settings are
        invalid.
    RuntimeError
        If the imputation routine has no finite valid input vectors.

    Notes
    -----

    '''
    if not isinstance(parse_output, bool):
        raise TypeError('parse_output must be a bool')
    
    if not isinstance(include_corners, bool):
        raise TypeError('include_corners must be a bool')
    
    if isinstance(tolerance, bool) or not isinstance(tolerance, Real):
        raise TypeError('tolerance must be a real number')
    
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral):
        raise TypeError('max_iterations must be an integer')
    
    data = _prepare_field(
        u, v, field_data,
        {
            's2n': s2n, 
            'p2p': p2p, 
            'peak': peak, 
            'flag': flag,
            'peak_ind': peak_ind,
        },
    )

    converged = cpiv.impute_nans(
        data, 
        include_corners, 
        float(tolerance), 
        int(max_iterations)
    )

    if not parse_output:
        return field_data, converged
        
    return (
        np.array(field_data.u, copy=True),
        np.array(field_data.v, copy=True),
        np.array(field_data.s2n, copy=True),
        np.array(field_data.p2p, copy=True),
        np.array(field_data.peak, copy=True),
        np.array(field_data.flag, copy=True),
        np.array(field_data.peak_ind, copy=True),
        converged
    )
