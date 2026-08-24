from __future__ import annotations

from numbers import Real
from typing import Any
from typing import TypeAlias

import numpy as np

from . import pyopenpivcore as cpiv
from .parsing import (
    convert_stats_to_grid_data_type,
    convert_vector_to_grid_data_type,
)

ValidationResult: TypeAlias = np.ndarray | cpiv.image_g_u8


__all__ = [
    'range_val',
    'global_std_val',
    'sig2noise_val',
    'peak2peak_val',
    'local_difference_val',
    'local_median_val',
    'local_normalized_median_val',
]


def _require_grid_data(field_data: Any) -> cpiv.grid_data:
    '''Validate an object as OpenPIV C++ grid data.

    Parameters
    ----------
    field_data : object
        Object to validate.

    Returns
    -------
    field_data : cpiv.grid_data
        The original C++ grid-data object.

    Raises
    ------
    TypeError
        If ``field_data`` is not ``cpiv.grid_data``.

    '''
    if not isinstance(field_data, cpiv.grid_data):
        raise TypeError(
            'field_data must be of type grid_data '
            f'(got {type(field_data).__name__})'
        )
    return field_data


def _resolve_vector_field(
    u: np.ndarray | cpiv.grid_data | None,
    v: np.ndarray | None,
    field_data: cpiv.grid_data | None,
) -> cpiv.grid_data:
    '''Resolve vector input into one C++ grid-data object.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array or a C++ grid-data object supplied positionally.
    v : np.ndarray or None
        V-displacement array.
    field_data : cpiv.grid_data or None
        C++ grid data supplied through the explicit parameter.

    Returns
    -------
    data : cpiv.grid_data
        Reused C++ data or a newly allocated object containing copies of u and
        v.

    Raises
    ------
    TypeError
        If a supplied object has an unsupported type.
    ValueError
        If NumPy and C++ representations are mixed, required inputs are
        missing, or array dimensions differ.

    '''
    # Allow the C++ object as the first positional argument as well as through
    # the explicit field_data keyword.
    if isinstance(u, cpiv.grid_data):
        if field_data is not None:
            raise ValueError('field_data was provided more than once')
        if v is not None:
            raise ValueError('v must be None when u is a grid_data object')
        field_data = u
        u = None

    if field_data is not None:
        if u is not None or v is not None:
            raise ValueError(
                'u and v must be None when field_data is provided'
            )
        return _require_grid_data(field_data)

    if u is None or v is None:
        raise ValueError('provide both u and v arrays, or provide field_data')
    return convert_vector_to_grid_data_type(u, v)


def _resolve_stat_field(
    values: np.ndarray | cpiv.grid_data | None,
    field_data: cpiv.grid_data | None,
    statistic: str,
) -> cpiv.grid_data:
    '''Resolve one statistic array into C++ grid data.

    Parameters
    ----------
    values : np.ndarray or cpiv.grid_data or None
        Statistic array or a C++ grid-data object supplied positionally.
    field_data : cpiv.grid_data or None
        C++ grid data supplied through the explicit parameter.
    statistic : str
        Grid-data member to populate, either ``'s2n'`` or ``'p2p'``.

    Returns
    -------
    data : cpiv.grid_data
        Reused C++ data or a newly allocated object containing a copy of the
        statistic array.

    Raises
    ------
    TypeError
        If a supplied object has an unsupported type.
    ValueError
        If NumPy and C++ representations are mixed or no field is supplied.

    '''
    if isinstance(values, cpiv.grid_data):
        if field_data is not None:
            raise ValueError('field_data was provided more than once')
        field_data = values
        values = None

    if field_data is not None:
        if values is not None:
            raise ValueError(
                f'{statistic} must be None when field_data is provided'
            )
        return _require_grid_data(field_data)

    if values is None:
        raise ValueError(
            f'provide a {statistic} array, or provide field_data'
        )
    if statistic == 's2n':
        return convert_stats_to_grid_data_type(s2n=values)
    return convert_stats_to_grid_data_type(p2p=values)


def _threshold(value: float | None, name: str) -> float:
    '''Convert a required threshold into a native float.

    Parameters
    ----------
    value : float or None
        Threshold value to validate.
    name : str
        Parameter name to include in validation errors.

    Returns
    -------
    threshold : float
        Validated native floating-point threshold.

    Raises
    ------
    TypeError
        If ``value`` is not a real number.
    ValueError
        If ``value`` is None.

    '''
    if value is None:
        raise ValueError(f'{name} must be specified')
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f'{name} must be a real number')
    return float(value)


def _nonnegative_threshold(value: float | None, name: str) -> float:
    '''Convert a required nonnegative threshold into a float.

    Parameters
    ----------
    value : float or None
        Threshold value to validate.
    name : str
        Parameter name to include in validation errors.

    Returns
    -------
    threshold : float
        Validated nonnegative threshold.

    Raises
    ------
    TypeError
        If ``value`` is not a real number.
    ValueError
        If ``value`` is None or negative.

    '''
    result = _threshold(value, name)
    if result < 0:
        raise ValueError(f'{name} must not be negative')
    return result


def _threshold_pair(
    value: tuple[float, float] | list[float] | None,
    name: str,
) -> tuple[float, float]:
    '''Validate an ordered pair of lower and upper thresholds.

    Parameters
    ----------
    value : tuple[float, float] or list[float] or None
        Minimum and maximum accepted values.
    name : str
        Parameter name to include in validation errors.

    Returns
    -------
    threshold : tuple[float, float]
        Validated ``(minimum, maximum)`` threshold pair.

    Raises
    ------
    TypeError
        If ``value`` is not a suitable sequence or contains non-real values.
    ValueError
        If the pair is missing, has the wrong length, or has reversed bounds.

    '''
    if value is None:
        raise ValueError(f'{name} must be specified')
    if isinstance(value, (str, bytes)):
        raise TypeError(f'{name} must be a two-element sequence')

    try:
        components = tuple(value)
    except TypeError as exc:
        raise TypeError(f'{name} must be a two-element sequence') from exc
    if len(components) != 2:
        raise ValueError(f'{name} must contain exactly two elements')

    lower = _threshold(components[0], f'{name}[0]')
    upper = _threshold(components[1], f'{name}[1]')
    if lower > upper:
        raise ValueError(f'{name} lower bound must not exceed its upper bound')
    return lower, upper


def _format_flag(flag: Any, return_raw: bool) -> Any:
    '''Format a C++ validation flag image for Python callers.

    Parameters
    ----------
    flag : cpiv.image_g_8
        C++ validation flag image.
    return_raw : bool
        If True, return ``flag`` unchanged. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Boolean NumPy mask or the original C++ flag image.

    Raises
    ------
    TypeError
        If ``return_raw`` is not a bool.

    '''
    if not isinstance(return_raw, bool):
        raise TypeError('return_raw must be a bool')
    if return_raw:
        return flag
    return np.array(flag, dtype=np.bool_, copy=True)


def range_val(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    u_threshold: tuple[float, float] | list[float] | None=None,
    v_threshold: tuple[float, float] | list[float] | None=None,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors against global component ranges.

    Flag a vector when its u or v component lies outside the corresponding
    inclusive threshold interval. Supply either separate NumPy u and v arrays
    or one reusable ``cpiv.grid_data`` object.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array. A ``cpiv.grid_data`` object may be passed here as
        the first positional argument instead.
    v : np.ndarray or None
        V-displacement array. Required when ``u`` is a NumPy array.
    u_threshold : tuple[float, float] or list[float] or None
        Inclusive ``(minimum, maximum)`` range for u displacements.
    v_threshold : tuple[float, float] or list[float] or None
        Inclusive ``(minimum, maximum)`` range for v displacements.
    field_data : cpiv.grid_data or None
        Reusable C++ vector data. It cannot be combined with u or v arrays.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values at invalid vectors.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If the field representation is incomplete or mixed, vector shapes
        differ, or a threshold interval is invalid.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_vector_field(u, v, field_data)
    threshold_u = _threshold_pair(u_threshold, 'u_threshold')
    threshold_v = _threshold_pair(v_threshold, 'v_threshold')
    flag = cpiv.validate_range(data, threshold_u, threshold_v)
    return _format_flag(flag, return_raw)


def global_std_val(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    std_threshold: float=3.0,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors against global standard-deviation limits.

    Compute the global mean and standard deviation of each displacement
    component and flag vectors farther than ``std_threshold`` standard
    deviations from either component mean. Nonfinite vectors are also flagged.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array. A ``cpiv.grid_data`` object may be passed here as
        the first positional argument instead.
    v : np.ndarray or None
        V-displacement array. Required when ``u`` is a NumPy array.
    std_threshold : float
        Nonnegative number of standard deviations allowed from each global
        component mean.
    field_data : cpiv.grid_data or None
        Reusable C++ vector data. It cannot be combined with u or v arrays.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values at invalid vectors.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If the field representation is incomplete or mixed, vector shapes
        differ, or ``std_threshold`` is negative.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_vector_field(u, v, field_data)
    threshold = _nonnegative_threshold(std_threshold, 'std_threshold')
    flag = cpiv.validate_z_score(data, threshold)
    return _format_flag(flag, return_raw)


def sig2noise_val(
    sig: np.ndarray | cpiv.grid_data | None=None,
    threshold: float | None=None,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors using their signal-to-noise ratios.

    Flag vectors whose signal-to-noise ratio is less than ``threshold``.
    Supply either a NumPy signal-to-noise array or a reusable C++ grid-data
    object containing its ``s2n`` member.

    Parameters
    ----------
    sig : np.ndarray or cpiv.grid_data or None
        Signal-to-noise ratio array. A ``cpiv.grid_data`` object may be passed
        here as the first positional argument instead.
    threshold : float or None
        Minimum accepted signal-to-noise ratio.
    field_data : cpiv.grid_data or None
        Reusable C++ grid data. It cannot be combined with ``sig``.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values below the threshold.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If no field is supplied or NumPy and C++ field inputs are mixed.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_stat_field(sig, field_data, 's2n')
    checked_threshold = _threshold(threshold, 'threshold')
    flag = cpiv.validate_s2n(data, checked_threshold)
    return _format_flag(flag, return_raw)


def peak2peak_val(
    sig: np.ndarray | cpiv.grid_data | None=None,
    threshold: float | None=None,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors using their primary-to-secondary peak ratios.

    Flag vectors whose primary-to-secondary correlation-peak ratio is less
    than ``threshold``. Supply either a NumPy ratio array or a reusable C++
    grid-data object containing its ``p2p`` member.

    Parameters
    ----------
    sig : np.ndarray or cpiv.grid_data or None
        Primary-to-secondary peak ratio array. A ``cpiv.grid_data`` object may
        be passed here as the first positional argument instead.
    threshold : float or None
        Minimum accepted primary-to-secondary peak ratio.
    field_data : cpiv.grid_data or None
        Reusable C++ grid data. It cannot be combined with ``sig``.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values below the threshold.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If no field is supplied or NumPy and C++ field inputs are mixed.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_stat_field(sig, field_data, 'p2p')
    checked_threshold = _threshold(threshold, 'threshold')
    flag = cpiv.validate_p2p(data, checked_threshold)
    return _format_flag(flag, return_raw)


def local_difference_val(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    u_threshold: float | None=None,
    v_threshold: float | None=None,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors against differences from neighboring vectors.

    Compare each vector with its local 3-by-3 neighborhood and flag vectors
    that differ from too many neighbors by more than the selected u or v
    component threshold. Nonfinite vectors are also flagged.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array. A ``cpiv.grid_data`` object may be passed here as
        the first positional argument instead.
    v : np.ndarray or None
        V-displacement array. Required when ``u`` is a NumPy array.
    u_threshold : float or None
        Maximum allowed absolute u-component difference.
    v_threshold : float or None
        Maximum allowed absolute v-component difference.
    field_data : cpiv.grid_data or None
        Reusable C++ vector data. It cannot be combined with u or v arrays.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values at invalid vectors.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If the field representation is incomplete or mixed, vector shapes
        differ, or a threshold is negative.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_vector_field(u, v, field_data)
    threshold_u = _nonnegative_threshold(u_threshold, 'u_threshold')
    threshold_v = _nonnegative_threshold(v_threshold, 'v_threshold')
    flag = cpiv.validate_difference(data, threshold_u, threshold_v)
    return _format_flag(flag, return_raw)


def local_median_val(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    u_threshold: float | None=None,
    v_threshold: float | None=None,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors against local component medians.

    Compute u and v medians over each local 3-by-3 neighborhood and flag a
    vector when either absolute component residual exceeds its threshold.
    Nonfinite vectors and neighborhoods with insufficient data are flagged.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array. A ``cpiv.grid_data`` object may be passed here as
        the first positional argument instead.
    v : np.ndarray or None
        V-displacement array. Required when ``u`` is a NumPy array.
    u_threshold : float or None
        Maximum allowed absolute residual from the local u median.
    v_threshold : float or None
        Maximum allowed absolute residual from the local v median.
    field_data : cpiv.grid_data or None
        Reusable C++ vector data. It cannot be combined with u or v arrays.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values at invalid vectors.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If the field representation is incomplete or mixed, vector shapes
        differ, or a threshold is negative.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_vector_field(u, v, field_data)
    threshold_u = _nonnegative_threshold(u_threshold, 'u_threshold')
    threshold_v = _nonnegative_threshold(v_threshold, 'v_threshold')
    flag = cpiv.validate_median(data, threshold_u, threshold_v)
    return _format_flag(flag, return_raw)


def local_normalized_median_val(
    u: np.ndarray | cpiv.grid_data | None=None,
    v: np.ndarray | None=None,
    threshold: float | None=None,
    field_data: cpiv.grid_data | None=None,
    return_raw: bool=False,
) -> ValidationResult:
    '''Validate vectors using the normalized local-median test.

    Compute the normalized residual of each vector from the local u and v
    medians and flag vectors whose combined residual exceeds ``threshold``.
    Nonfinite vectors and neighborhoods with insufficient data are flagged.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data or None
        U-displacement array. A ``cpiv.grid_data`` object may be passed here as
        the first positional argument instead.
    v : np.ndarray or None
        V-displacement array. Required when ``u`` is a NumPy array.
    threshold : float or None
        Nonnegative maximum normalized-median residual.
    field_data : cpiv.grid_data or None
        Reusable C++ vector data. It cannot be combined with u or v arrays.
    return_raw : bool
        If True, return the C++ flag image. Otherwise, return a boolean NumPy
        copy.

    Returns
    -------
    invalid : np.ndarray or cpiv.image_g_8
        Flag image with nonzero or True values at invalid vectors.

    Raises
    ------
    TypeError
        If a field, threshold, or option has an unsupported type.
    ValueError
        If the field representation is incomplete or mixed, vector shapes
        differ, or ``threshold`` is negative.
    RuntimeError
        If the underlying C++ validation routine fails.

    '''
    data = _resolve_vector_field(u, v, field_data)
    checked_threshold = _nonnegative_threshold(threshold, 'threshold')
    flag = cpiv.validate_normalized_median(data, checked_threshold)
    return _format_flag(flag, return_raw)