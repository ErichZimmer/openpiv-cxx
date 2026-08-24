from __future__ import annotations

from numbers import Integral
from typing import TypeAlias

import numpy as np

from . import pyopenpivcore as cpiv
from .parsing import (
    convert_coords_to_grid_coords_type,
    convert_cxx_data_to_numpy,
    convert_to_image_type,
    convert_to_list_type,
    convert_vector_to_grid_data_type,
)


__all__ = [
    'first_pass',
    'multi_pass'
]


ImageInput: TypeAlias = np.ndarray | cpiv.image_g_f32 | cpiv.image_g_f64
WindowInput: TypeAlias = int | tuple[int, int] | list[int]
CppOutput: TypeAlias = tuple[cpiv.grid_coords, cpiv.grid_data]
NumpyOutput: TypeAlias = tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]
ProcessOutput: TypeAlias = CppOutput | NumpyOutput


def _prepare_images(
    image_a: ImageInput,
    image_b: ImageInput,
) -> tuple[cpiv.image_g_f32, cpiv.image_g_f32]:
    '''Convert and validate a pair of particle images.

    Parameters
    ----------
    image_a : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        First particle image.
    image_b : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        Second particle image.

    Returns
    -------
    image_a_cxx : cpiv.image_g_f32
        First image in the C++ processing type.
    image_b_cxx : cpiv.image_g_f32
        Second image in the C++ processing type.

    Raises
    ------
    TypeError
        If either image has an unsupported type.
    ValueError
        If an image is not two-dimensional or the image shapes differ.

    '''
    image_a_cxx = convert_to_image_type(image_a)
    image_b_cxx = convert_to_image_type(image_b)

    shape_a = (image_a_cxx.height(), image_a_cxx.width())
    shape_b = (image_b_cxx.height(), image_b_cxx.width())
    if shape_a != shape_b:
        raise ValueError(
            'image_a and image_b must have the same shape '
            f'(got {shape_a} and {shape_b})'
        )

    return image_a_cxx, image_b_cxx


def _prepare_windows(
    window_size: WindowInput,
    overlap: WindowInput,
) -> tuple[list[int], list[int]]:
    '''Normalize and validate window and overlap components.

    Parameters
    ----------
    window_size : int or tuple[int, int] or list[int]
        Interrogation-window dimensions.
    overlap : int or tuple[int, int] or list[int]
        Interrogation-window overlap dimensions.

    Returns
    -------
    window : list[int]
        Normalized x and y window dimensions.
    overlap_size : list[int]
        Normalized x and y overlap dimensions.

    Raises
    ------
    TypeError
        If either input has unsupported components.
    ValueError
        If window components are zero or overlap is not smaller than the
        corresponding window component.

    '''
    window = convert_to_list_type(window_size)
    overlap_size = convert_to_list_type(overlap)

    if any(component == 0 for component in window):
        raise ValueError('window_size components must be greater than zero')
    if any(o >= w for o, w in zip(overlap_size, window)):
        raise ValueError(
            'overlap components must be smaller than the corresponding '
            'window_size components'
        )

    return window, overlap_size


def _validate_execution_options(
    threads: int,
    simd: bool,
    robust: bool,
    zero_pad: bool,
) -> None:
    '''Validate boolean and thread processing options.

    Parameters
    ----------
    threads : int
        Requested processing thread count.
    simd : bool
        Whether to enable the SIMD processing path.
    robust : bool
        Whether to enable robust phase correlation.
    zero_pad : bool
        Whether to zero-pad interrogation windows before correlation.

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If ``threads`` is not an integer or a boolean option is not a bool.

    '''
    if isinstance(threads, bool) or not isinstance(threads, Integral):
        raise TypeError('threads must be an integer')
    if not isinstance(simd, bool):
        raise TypeError('simd must be a bool')
    if not isinstance(robust, bool):
        raise TypeError('robust must be a bool')
    if not isinstance(zero_pad, bool):
        raise TypeError('zero_pad must be a bool')


def _format_output(
    field_coords: cpiv.grid_coords,
    field_data: cpiv.grid_data,
    parse_output: bool,
) -> ProcessOutput:
    '''Format a C++ PIV result for the requested Python return mode.

    Parameters
    ----------
    field_coords : cpiv.grid_coords
        C++ vector-grid coordinates.
    field_data : cpiv.grid_data
        C++ vector data.
    parse_output : bool
        If True, return NumPy copies. If False, return the C++ objects.

    Returns
    -------
    output : tuple[cpiv.grid_coords, cpiv.grid_data] or tuple[np.ndarray, ...]
        C++ objects or the ``(x, y, u, v, s2n, p2p)`` NumPy arrays.

    Raises
    ------
    TypeError
        If ``parse_output`` is not a bool or the C++ result types are invalid.
    ValueError
        If the coordinate buffer has an invalid shape.

    '''
    if not isinstance(parse_output, bool):
        raise TypeError('parse_output must be a bool')
    if not parse_output:
        return field_coords, field_data

    return convert_cxx_data_to_numpy(field_coords, field_data, copy=True)


def _deformation_method(value: str | cpiv.deform_method) -> cpiv.deform_method:
    '''Convert a deformation-method option into its C++ enum value.

    Parameters
    ----------
    value : str or cpiv.deform_method
        Deformation interpolation method.

    Returns
    -------
    method : cpiv.deform_method
        Corresponding C++ deformation-method enum.

    Raises
    ------
    TypeError
        If ``value`` is neither a string nor ``cpiv.deform_method``.
    ValueError
        If a string does not name a supported method.

    '''
    if isinstance(value, cpiv.deform_method):
        return value
    if not isinstance(value, str):
        raise TypeError('method must be a deform_method or string')

    methods = {
        'sinc': cpiv.deform_method.SINC,
        'lagrange': cpiv.deform_method.LAGRANGE,
        'lanczos': cpiv.deform_method.LANCZOS,
    }
    try:
        return methods[value.casefold()]
    except KeyError as exc:
        raise ValueError(
            f'deformation method {value!r} is not supported; '
            'expected sinc, lagrange, or lanczos'
        ) from exc


def _deformation_order(value: str | cpiv.deform_order) -> cpiv.deform_order:
    '''Convert a deformation-order option into its C++ enum value.

    Parameters
    ----------
    value : str or cpiv.deform_order
        Deformation order.

    Returns
    -------
    order : cpiv.deform_order
        Corresponding C++ deformation-order enum.

    Raises
    ------
    TypeError
        If ``value`` is neither a string nor ``cpiv.deform_order``.
    ValueError
        If a string does not name a supported order.

    '''
    if isinstance(value, cpiv.deform_order):
        return value
    if not isinstance(value, str):
        raise TypeError('order must be a deform_order or string')

    orders = {
        'forward': cpiv.deform_order.FORWARD,
        'symmetric': cpiv.deform_order.SYMMETRIC,
    }
    try:
        return orders[value.casefold()]
    except KeyError as exc:
        raise ValueError(
            f'deformation order {value!r} is not supported; '
            'expected forward or symmetric'
        ) from exc


def _prepare_previous_field(
    x: np.ndarray | None,
    y: np.ndarray | None,
    u: np.ndarray | None,
    v: np.ndarray | None,
    s2n: np.ndarray | None,
    p2p: np.ndarray | None,
    peak: np.ndarray | None,
    flag: np.ndarray | None,
    field_coords: cpiv.grid_coords | None,
    field_data: cpiv.grid_data | None,
) -> tuple[cpiv.grid_coords, cpiv.grid_data]:
    '''Build or reuse one complete previous-pass vector field.

    Select exactly one input mode. C++ mode requires both ``field_coords`` and
    ``field_data`` and forwards them unchanged. NumPy mode requires x, y, u,
    and v arrays and copies them, together with any optional field members,
    into new C++ grid objects.

    Parameters
    ----------
    x : np.ndarray or None
        Previous x-coordinate array.
    y : np.ndarray or None
        Previous y-coordinate array.
    u : np.ndarray or None
        Previous u-displacement array.
    v : np.ndarray or None
        Previous v-displacement array.
    s2n : np.ndarray or None
        Optional signal-to-noise ratio array.
    p2p : np.ndarray or None
        Optional primary-to-secondary peak ratio array.
    peak : np.ndarray or None
        Optional primary correlation-peak array.
    flag : np.ndarray or None
        Optional vector-status flag array.
    field_coords : cpiv.grid_coords or None
        Reusable C++ grid coordinates.
    field_data : cpiv.grid_data or None
        Reusable C++ vector data.

    Returns
    -------
    coords : cpiv.grid_coords
        Reused or newly constructed C++ grid coordinates.
    data : cpiv.grid_data
        Reused or newly constructed C++ vector data.

    Raises
    ------
    TypeError
        If a supplied array or C++ wrapper has an unsupported type.
    ValueError
        If the two representations are mixed, a representation is incomplete,
        or component shapes differ.

    '''
    array_inputs = {
        'x': x,
        'y': y,
        'u': u,
        'v': v,
        's2n': s2n,
        'p2p': p2p,
        'peak': peak,
        'flag': flag,
    }
    supplied_arrays = [
        name for name, value in array_inputs.items() if value is not None
    ]
    wrappers_supplied = field_coords is not None or field_data is not None

    if wrappers_supplied and supplied_arrays:
        names = ', '.join(supplied_arrays)
        raise ValueError(
            'provide either field_coords and field_data, or NumPy field '
            f'arrays; do not mix them (received arrays: {names})'
        )

    if wrappers_supplied:
        if field_coords is None or field_data is None:
            raise ValueError(
                'field_coords and field_data must be provided together'
            )
        if not isinstance(field_coords, cpiv.grid_coords):
            raise TypeError('field_coords must be of type grid_coords')
        coords = field_coords
        if not isinstance(field_data, cpiv.grid_data):
            raise TypeError('field_data must be of type grid_data')
        data = field_data
    else:
        required_arrays = {'x': x, 'y': y, 'u': u, 'v': v}
        missing = [
            name for name, value in required_arrays.items() if value is None
        ]
        if missing:
            names = ', '.join(missing)
            raise ValueError(
                'NumPy field input requires x, y, u, and v arrays '
                f'(missing: {names})'
            )

        assert x is not None and y is not None
        assert u is not None and v is not None
        coords = convert_coords_to_grid_coords_type(x, y)
        data = convert_vector_to_grid_data_type(
            u,
            v,
            s2n=s2n,
            p2p=p2p,
            peak=peak,
            flag=flag,
        )

    coord_shape = (coords.height(), coords.width())
    data_shape = (data.u.height(), data.u.width())
    v_shape = (data.v.height(), data.v.width())
    if data_shape != v_shape:
        raise ValueError(
            'field_data.u and field_data.v must have the same shape '
            f'(got {data_shape} and {v_shape})'
        )
    if coord_shape != data_shape:
        raise ValueError(
            'field_coords and field_data must have the same shape '
            f'(got {coord_shape} and {data_shape})'
        )

    return coords, data


def first_pass(
    image_a: ImageInput,
    image_b: ImageInput,
    window_size: WindowInput=32,
    overlap: WindowInput=16,
    robust: bool=False,
    zero_pad: bool=False,
    centered: bool=False,
    limit_search: bool=False,
    simd: bool=False,
    threads: int=1,
    parse_output: bool=True
) -> ProcessOutput:
    '''Perform a standard first-pass PIV evaluation.

    Convert the input images into OpenPIV C++ grayscale images, evaluate the
    vector field with the requested interrogation-window settings, and return
    either NumPy arrays or the reusable C++ grid objects.

    Parameters
    ----------
    image_a : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        First particle image. NumPy and double-precision C++ images are copied
        into the single-precision image type used by the PIV implementation.
    image_b : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        Second particle image. It must have the same shape as ``image_a``.
    window_size : int or tuple[int, int] or list[int]
        Interrogation-window size in pixels. A scalar applies to both axes.
    overlap : int or tuple[int, int] or list[int]
        Number of shared pixels between adjacent interrogation windows. A
        scalar applies to both axes, and each component must be smaller than
        the corresponding ``window_size`` component.
    robust : bool
        If True, enable robust minimum quadratic differences template matching
    zero_pad : bool
        If True, zero-pad interrogation windows before correlation.
    centered : bool
        If True, center the interrogation grid within the images.
    limit_search : bool
        If True, restrict the peak search to the central search region.
    simd : bool
            If True, enable the SIMD processing path where supported.
    threads : int
        Requested processing thread count. Nonpositive values defer thread
        selection to the C++ implementation.
    parse_output : bool
            If True, return NumPy copies. If False, return ``cpiv.grid_coords``
            and ``cpiv.grid_data`` objects that can be reused by ``multi_pass``.

    Returns
    -------
    field_coords : cpiv.grid_coords
        C++ vector-grid coordinates, returned when ``parse_output`` is False.
    field_data : cpiv.grid_data
        C++ vector data, returned when ``parse_output`` is False.
    x : np.ndarray
        X-coordinate array, returned when ``parse_output`` is True.
    y : np.ndarray
        Y-coordinate array, returned when ``parse_output`` is True.
    u : np.ndarray
        U-displacement array, returned when ``parse_output`` is True.
    v : np.ndarray
        V-displacement array, returned when ``parse_output`` is True.
    s2n : np.ndarray
        Signal-to-noise ratio array, returned when ``parse_output`` is True.
    p2p : np.ndarray
        Primary-to-secondary peak ratio array, returned when ``parse_output``
        is True.

    Raises
    ------
    TypeError
        If an image or processing option has an unsupported type.
    ValueError
        If image shapes differ or a window or overlap option is invalid.
    RuntimeError
        If the underlying C++ processing routine fails.

    '''
    image_a_cxx, image_b_cxx = _prepare_images(image_a, image_b)
    window, overlap_size = _prepare_windows(window_size, overlap)
    _validate_execution_options(threads, simd, robust, zero_pad)

    if robust:
        field_coords, field_data = cpiv.process_images_standard(
            image_a_cxx,
            image_b_cxx,
            window,
            overlap_size,
            step=False,
            zero_pad=zero_pad,
            centered=centered,
            limit_search=limit_search,
            simd=simd,
            threads=int(threads),
        )
    else:
        field_coords, field_data = cpiv.process_images_robust(
            image_a_cxx,
            image_b_cxx,
            window,
            overlap_size,
            step=False,
            zero_pad=zero_pad,
            centered=centered,
            limit_search=limit_search,
            simd=simd,
            threads=int(threads),
        )
    
    return _format_output(field_coords, field_data, parse_output)


def multi_pass(
    image_a: ImageInput,
    image_b: ImageInput,
    field_coords: cpiv.grid_coords | None=None,
    field_data: cpiv.grid_data | None=None,
    x: np.ndarray | None=None,
    y: np.ndarray | None=None,
    u: np.ndarray | None=None,
    v: np.ndarray | None=None,
    window_size: WindowInput=32,
    overlap: WindowInput=16,
    robust: bool=False,
    method: str | cpiv.deform_method='lagrange',
    order: str | cpiv.deform_order='forward',
    k: int=3,
    zero_pad: bool=False,
    centered: bool=False,
    limit_search: bool=False,
    simd: bool=False,
    threads: int=1,
    parse_output: bool=True,
) -> ProcessOutput:
    '''Perform a window-deformation PIV refinement pass.

    Refine a previous vector field by deforming the particle images and
    repeating the PIV evaluation with the requested interrogation-window
    settings. Supply the previous field either as a complete pair of C++ grid
    objects or as NumPy arrays. The two representations cannot be mixed.

    Parameters
    ----------
    image_a : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        First particle image.
    image_b : np.ndarray or cpiv.image_g_f32 or cpiv.image_g_f64
        Second particle image. It must have the same shape as ``image_a``.
    field_coords : cpiv.grid_coords or None
        Previous C++ grid coordinates. It must be supplied together with
        ``field_data`` and without any NumPy field arrays.
    field_data : cpiv.grid_data or None
        Previous C++ vector data. It must be supplied together with
        ``field_coords`` and without any NumPy field arrays.
    x : np.ndarray or None
            Previous x-coordinate array. Required for NumPy input mode.
    y : np.ndarray or None
        Previous y-coordinate array. Required for NumPy input mode.
    u : np.ndarray or None
        Previous u-displacement array. Required for NumPy input mode.
    v : np.ndarray or None
        Previous v-displacement array. Required for NumPy input mode.
    window_size : int or tuple[int, int] or list[int]
        Interrogation-window size in pixels. A scalar applies to both axes.
    overlap : int or tuple[int, int] or list[int]
        Number of shared pixels between adjacent interrogation windows. Each
        component must be smaller than the corresponding window component.
    robust : bool
        If True, enable robust phase correlation.
    method : str or cpiv.deform_method
        Image-deformation interpolation method. Supported strings are
        ``'sinc'``, ``'lagrange'``, and ``'lanczos'``.
    order : str or cpiv.deform_order
        Deformation order. Supported strings are ``'forward'`` and
        ``'symmetric'``.
    k : int
        Positive interpolation-kernel parameter passed to the selected
        deformation method.
    zero_pad : bool
        If True, zero-pad interrogation windows before correlation.
    centered : bool
        If True, center the interrogation grid within the images.
    limit_search : bool
        If True, restrict the peak search to the central search region.
    simd : bool
        If True, enable the SIMD processing path where supported.
    threads : int
        Requested processing thread count. Nonpositive values defer thread
        selection to the C++ implementation.
    parse_output : bool
        If True, return NumPy copies. If False, return reusable C++ grid
        objects.

    Returns
    -------
    field_coords : cpiv.grid_coords
        Refined C++ vector-grid coordinates, returned when ``parse_output`` is
        False.
    field_data : cpiv.grid_data
        Refined C++ vector data, returned when ``parse_output`` is False.
    x : np.ndarray
        Refined x-coordinate array, returned when ``parse_output`` is True.
    y : np.ndarray
        Refined y-coordinate array, returned when ``parse_output`` is True.
    u : np.ndarray
        Refined u-displacement array, returned when ``parse_output`` is True.
    v : np.ndarray
        Refined v-displacement array, returned when ``parse_output`` is True.
    s2n : np.ndarray
        Refined signal-to-noise ratio array, returned when ``parse_output`` is
        True.
    p2p : np.ndarray
        Refined primary-to-secondary peak ratio array, returned when
        ``parse_output`` is True.

    Raises
    ------
    TypeError
        If an image, previous-field object, or processing option has an
        unsupported type.
    ValueError
        If the previous-field modes are mixed or incomplete, array or image
        shapes differ, or a processing option is invalid.
    RuntimeError
        If the underlying C++ processing routine fails.

    '''
    image_a_cxx, image_b_cxx = _prepare_images(image_a, image_b)
    window, overlap_size = _prepare_windows(window_size, overlap)
    field_coords, field_data = _prepare_previous_field(
        field_coords=field_coords,
        field_data=field_data,
        x=x,
        y=y,
        u=u,
        v=v,
    )
    deform_method = _deformation_method(method)
    deform_order = _deformation_order(order)
    _validate_execution_options(threads, simd, robust, zero_pad)

    if isinstance(k, bool) or not isinstance(k, Integral):
        raise TypeError('k must be an integer')
    if k <= 0:
        raise ValueError('k must be greater than zero')

    field_coords, field_data = cpiv.process_images_multipass(
        image_a_cxx,
        image_b_cxx,
        field_coords,
        field_data,
        window,
        overlap_size,
        robust=robust,
        method=deform_method,
        order=deform_order,
        k=int(k),
        step=False,
        zero_pad=zero_pad,
        centered=centered,
        limit_search=limit_search,
        simd=simd,
        threads=int(threads),
    )

    return _format_output(field_coords, field_data, parse_output)