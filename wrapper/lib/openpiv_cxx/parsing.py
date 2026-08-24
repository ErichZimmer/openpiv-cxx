from __future__ import annotations

from numbers import Integral
from typing import TypeAlias

import numpy as np

from . import pyopenpivcore as cpiv


__all__ = [
    'check_nd',
    'convert_to_image_type',
    'convert_to_list_type',
    'convert_cxx_data_to_numpy',
    'convert_coords_to_grid_coords_type',
    'convert_vector_to_grid_data_type',
    'convert_stats_to_grid_data_type',
]


ImageInput: TypeAlias = np.ndarray | cpiv.image_g_f32 | cpiv.image_g_f64
WindowInput: TypeAlias = int | tuple[int, int] | list[int]

SUPPORTED_IMAGE_TYPES = (
    cpiv.image_g_f32,
    cpiv.image_g_f64,
    np.ndarray,
)

# Retain the misspelled historical constant for callers that imported it.
SUPPPORTED_IMAGE_TYPES = SUPPORTED_IMAGE_TYPES


def check_nd(
    ndim: int=2, 
    **arrays: np.ndarray
) -> None:
    '''Validate the dimensionality of one or more arrays.

    Check that every supplied array has the requested number of dimensions.
    The keyword used for each array is included in any resulting error.

    Parameters
    ----------
    ndim : int
        The number of dimensions each supplied array must have.
    **arrays : np.ndarray
        Arrays to validate, supplied as named keyword arguments.

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If ``ndim`` is not an integer.
    ValueError
        If ``ndim`` is negative or an array has a different dimensionality.

    '''
    if isinstance(ndim, bool) or not isinstance(ndim, Integral):
        raise TypeError('ndim must be an integer')
    if ndim < 0:
        raise ValueError('ndim must not be negative')

    for name, array in arrays.items():
        actual = np.ndim(array)
        if actual != ndim:
            raise ValueError(
                f'{name} must be a {ndim}D array (got {actual} dimensions)'
            )


def convert_to_image_type(
    image: ImageInput,
    double: bool=False,
) -> cpiv.image_g_f32 | cpiv.image_g_f64:
    '''Convert an array into an OpenPIV C++ grayscale image.

    Convert either a NumPy array or a C++ ``core::image`` into the requested
    data type. C++ images of the requested type are returned without further
    processing.

    Parameters
    ----------
    image : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        The image to convert into a supported C++ image type.
    double : bool
        If True, return ``cpiv.image_g_f64``. Otherwise, return
        ``cpiv.image_g_f32``.

    Returns
    -------
    image : cpiv.image_g_f32 or cpiv.image_g_f64
        The converted image. A matching C++ image is returned unchanged;
        other supported inputs are copied.

    Raises
    ------
    TypeError
        If ``double`` is not a bool or ``image`` has an unsupported type.
    ValueError
        If a NumPy input is not two-dimensional.

    '''
    if not isinstance(double, bool):
        raise TypeError('double must be a bool')
    image_type = cpiv.image_g_f64 if double else cpiv.image_g_f32

    if not isinstance(image, SUPPORTED_IMAGE_TYPES):
        raise TypeError(
            f'image has unsupported type {type(image).__name__}; expected a '
            '2D NumPy array, image_g_f32, or image_g_f64'
        )

    if isinstance(image, np.ndarray):
        check_nd(ndim=2, image=image)
        converted = image_type(list(image.shape[::-1]))
        np.array(converted, copy=False)[:] = image
        return converted

    if isinstance(image, image_type):
        return image

    converted = image_type(image.size())
    np.array(converted, copy=False)[:] = np.array(image, copy=False)
    return converted


def convert_to_list_type(window: WindowInput) -> list[int]:
    '''Convert a window or overlap size into a two-element list.

    Expand an integer into equal x and y components, or copy the two
    components from a tuple or list. The returned values are native Python
    integers suitable for conversion to the C++ size arguments.

    Parameters
    ----------
    window : int or tuple[int, int] or list[int]
        Scalar or two-component window or overlap size.

    Returns
    -------
    window : list[int]
        The normalized ``[x, y]`` components.

    Raises
    ------
    TypeError
        If the input is not a supported container or a component is not an
        integer.
    ValueError
        If a container does not have two elements or a component is negative.

    '''
    if isinstance(window, bool):
        raise TypeError('window/overlap size must contain integers, not bools')

    if isinstance(window, Integral):
        components = [int(window), int(window)]
    elif isinstance(window, (tuple, list)):
        if len(window) != 2:
            raise ValueError(
                'window/overlap size must contain exactly two elements '
                f'(got {len(window)})'
            )
        components = list(window)
    else:
        raise TypeError(
            'window/overlap size must be an integer, tuple, or list '
            f'(got {type(window).__name__})'
        )

    if any(
        isinstance(value, bool) or not isinstance(value, Integral)
        for value in components
    ):
        raise TypeError('window/overlap components must be integers')

    normalized = [int(value) for value in components]
    if any(value < 0 for value in normalized):
        raise ValueError('window/overlap components must not be negative')
    return normalized


def convert_cxx_data_to_numpy(
    field_coords: cpiv.grid_coords,
    field_data: cpiv.grid_data,
    copy: bool=True,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    '''Convert OpenPIV C++ vector-field types into NumPy arrays.

    Extract the coordinate, displacement, and correlation-statistic buffers
    from ``field_coords`` and ``field_data``. By default the returned arrays
    are independent copies. With ``copy=False``, they reference memory owned
    by the C++ wrapper objects, which must remain alive while the views are in
    use.

    Parameters
    ----------
    field_coords : cpiv.grid_coords
        C++ vector-grid coordinates containing x and y values.
    field_data : cpiv.grid_data
        C++ vector data containing u, v, s2n, p2p, and related members.
    copy : bool
        If True, copy the C++ buffers. If False, return NumPy views where the
        binding permits them.

    Returns
    -------
    x : np.ndarray
        Two-dimensional x-coordinate array.
    y : np.ndarray
        Two-dimensional y-coordinate array.
    u : np.ndarray
        Two-dimensional u-displacement array.
    v : np.ndarray
        Two-dimensional v-displacement array.
    s2n : np.ndarray
        Two-dimensional signal-to-noise ratio array.
    p2p : np.ndarray
        Two-dimensional primary-to-secondary peak ratio array.

    Raises
    ------
    TypeError
        If the wrapper objects have the wrong types or ``copy`` is not a bool.
    ValueError
        If the coordinate buffer does not have shape ``(height, width, 2)``.

    '''
    if not isinstance(field_coords, cpiv.grid_coords):
        raise TypeError(
            'field_coords must be of type grid_coords '
            f'(got {type(field_coords).__name__})'
        )
    if not isinstance(field_data, cpiv.grid_data):
        raise TypeError(
            'field_data must be of type grid_data '
            f'(got {type(field_data).__name__})'
        )
    if not isinstance(copy, bool):
        raise TypeError('copy must be a bool')

    coord_array = np.array(field_coords, copy=copy)
    if coord_array.ndim != 3 or coord_array.shape[-1] != 2:
        raise ValueError(
            'grid_coords buffer must have shape (height, width, 2); '
            f'got {coord_array.shape}'
        )

    x = coord_array[:, :, 0]
    y = coord_array[:, :, 1]
    u = np.array(field_data.u, copy=copy)
    v = np.array(field_data.v, copy=copy)
    s2n = np.array(field_data.s2n, copy=copy)
    p2p = np.array(field_data.p2p, copy=copy)
    return x, y, u, v, s2n, p2p


def _require_2d_array(
    array: np.ndarray, 
    name: str
) -> np.ndarray:
    '''Validate one input as a two-dimensional NumPy array.

    Parameters
    ----------
    array : np.ndarray
        Array to validate.
    name : str
        Parameter name to include in validation errors.

    Returns
    -------
    array : np.ndarray
        The original array, returned without copying.

    Raises
    ------
    TypeError
        If ``array`` is not a NumPy array.
    ValueError
        If ``array`` is not two-dimensional.

    '''
    if not isinstance(array, np.ndarray):
        raise TypeError(f'{name} must be a NumPy ndarray')
    check_nd(ndim=2, **{name: array})
    return array


def _require_same_shape(
    left: np.ndarray,
    right: np.ndarray,
    left_name: str,
    right_name: str,
) -> None:
    '''Validate that two arrays have identical shapes.

    Parameters
    ----------
    left : np.ndarray
        First array to compare.
    right : np.ndarray
        Second array to compare.
    left_name : str
        Name of the first array for error reporting.
    right_name : str
        Name of the second array for error reporting.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If the array shapes differ.

    '''
    if left.shape != right.shape:
        raise ValueError(
            f'{left_name} and {right_name} must have the same shape '
            f'(got {left.shape} and {right.shape})'
        )


def convert_coords_to_grid_coords_type(
    x: np.ndarray,
    y: np.ndarray,
) -> cpiv.grid_coords:
    '''Convert coordinate arrays into OpenPIV C++ grid coordinates.

    Allocate one ``cpiv.grid_coords`` object and copy the x and y arrays into
    its two coordinate components.

    Parameters
    ----------
    x : np.ndarray
        Two-dimensional x-coordinate array.
    y : np.ndarray
        Two-dimensional y-coordinate array with the same shape as ``x``.

    Returns
    -------
    field_coords : cpiv.grid_coords
        C++ grid coordinates containing copies of ``x`` and ``y``.

    Raises
    ------
    TypeError
        If either input is not a NumPy array.
    ValueError
        If an input is not two-dimensional or their shapes differ.

    '''
    x = _require_2d_array(x, 'x')
    y = _require_2d_array(y, 'y')
    _require_same_shape(x, y, 'x', 'y')

    field_coords = cpiv.grid_coords(list(x.shape[::-1]))
    coord_buffer = np.array(field_coords, copy=False)
    coord_buffer[:, :, 0] = x
    coord_buffer[:, :, 1] = y
    return field_coords


def convert_vector_to_grid_data_type(
    u: np.ndarray,
    v: np.ndarray,
    s2n: np.ndarray | None=None,
    p2p: np.ndarray | None=None,
    peak: np.ndarray | None=None,
    flag: np.ndarray | None=None,
) -> cpiv.grid_data:
    '''Convert vector-field arrays into OpenPIV C++ grid data.

    Allocate one ``cpiv.grid_data`` object and copy each supplied array into
    its corresponding member. The u and v arrays are required. Correlation
    statistics and vector flags may be omitted when they are not needed.

    Parameters
    ----------
    u : np.ndarray
        Two-dimensional u-displacement array.
    v : np.ndarray
        Two-dimensional v-displacement array.
    s2n : np.ndarray or None
        Optional signal-to-noise ratio array.
    p2p : np.ndarray or None
        Optional primary-to-secondary peak ratio array.
    peak : np.ndarray or None
        Optional primary correlation-peak value array.
    flag : np.ndarray or None
        Optional vector-status flag array.

    Returns
    -------
    field_data : cpiv.grid_data
        C++ grid data containing copies of all supplied arrays.

    Raises
    ------
    TypeError
        If a supplied value is not a NumPy array.
    ValueError
        If an array is not two-dimensional or the supplied shapes differ.

    '''
    u = _require_2d_array(u, 'u')
    v = _require_2d_array(v, 'v')
    _require_same_shape(u, v, 'u', 'v')

    optional_arrays = {
        's2n': s2n,
        'p2p': p2p,
        'peak': peak,
        'flag': flag,
    }
    checked_arrays: dict[str, np.ndarray] = {}
    for name, array in optional_arrays.items():
        if array is None:
            continue
        checked = _require_2d_array(array, name)
        _require_same_shape(u, checked, 'u', name)
        checked_arrays[name] = checked

    field_data = cpiv.grid_data(list(u.shape[::-1]))
    np.array(field_data.u, copy=False)[:] = u
    np.array(field_data.v, copy=False)[:] = v
    for name, array in checked_arrays.items():
        np.array(getattr(field_data, name), copy=False)[:] = array
    return field_data


def convert_stats_to_grid_data_type(
    s2n: np.ndarray | None=None,
    p2p: np.ndarray | None=None,
) -> cpiv.grid_data:
    '''Convert correlation-statistic arrays into C++ grid data.

    Allocate one ``cpiv.grid_data`` object and copy the supplied signal-to-noise
    or peak-ratio arrays into their corresponding members. This is intended for
    validation routines that only consume one correlation statistic.

    Parameters
    ----------
    s2n : np.ndarray or None
        Optional two-dimensional signal-to-noise ratio array.
    p2p : np.ndarray or None
        Optional two-dimensional primary-to-secondary peak ratio array.

    Returns
    -------
    field_data : cpiv.grid_data
        C++ grid data containing the supplied statistic arrays.

    Raises
    ------
    TypeError
        If a supplied value is not a NumPy array.
    ValueError
        If neither array is supplied, an array is not two-dimensional, or the
        two supplied arrays have different shapes.

    '''
    if s2n is None and p2p is None:
        raise ValueError('either s2n or p2p must be specified')

    if s2n is not None:
        s2n = _require_2d_array(s2n, 's2n')
    if p2p is not None:
        p2p = _require_2d_array(p2p, 'p2p')
    if s2n is not None and p2p is not None:
        _require_same_shape(s2n, p2p, 's2n', 'p2p')

    source = s2n if s2n is not None else p2p
    assert source is not None
    field_data = cpiv.grid_data(list(source.shape[::-1]))

    if s2n is not None:
        np.array(field_data.s2n, copy=False)[:] = s2n
    if p2p is not None:
        np.array(field_data.p2p, copy=False)[:] = p2p
    return field_data