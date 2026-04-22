from typing import Tuple
from .parsing import check_nd, convert_stats_to_grid_data_type, convert_vector_to_grid_data_type

import numpy as np
from . import pyopenpivcore as cpiv


__all__ = [
    'range_val',
    'global_std_val',
    'sig2noise_val',
    'peak2peak_val',
    'local_difference_val',
    'local_median_val',
    'local_normalized_median_val'
]


def range_val(
    u: np.ndarray,
    v: np.ndarray,
    u_threshold: Tuple[float, float],
    v_threshold: Tuple[float, float],
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
)-> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_vector_to_grid_data_type(
            u,
            v
        )

    # Call the c++ function
    c_flag = cpiv.validate_range(
        field_data,
        u_threshold,
        v_threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag


def global_std_val(
    u: np.ndarray,
    v: np.ndarray,
    std_threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
)-> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_vector_to_grid_data_type(
            u,
            v
        )

    # Call the c++ function
    c_flag = cpiv.validate_z_score(
        field_data,
        std_threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag


def sig2noise_val(
    sig: np.ndarray,
    threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
) -> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if sig is not None:
            msg = 'sig must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_stats_to_grid_data_type(
            s2n=sig
        )

    # Call the c++ function
    c_flag = cpiv.validate_s2n(
        field_data,
        threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag


def peak2peak_val(
    sig: np.ndarray,
    threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
) -> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if sig is not None:
            msg = 'sig must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_stats_to_grid_data_type(
            p2p=sig
        )

    # Call the c++ function
    c_flag = cpiv.validate_p2p(
        field_data,
        threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag


def local_difference_val(
    u: np.ndarray,
    v: np.ndarray,
    u_threshold: float,
    v_threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
)-> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_vector_to_grid_data_type(
            u,
            v
        )

    # Call the c++ function
    c_flag = cpiv.validate_difference(
        field_data,
        u_threshold,
        v_threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag


def local_median_val(
    u: np.ndarray,
    v: np.ndarray,
    u_threshold: float,
    v_threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
)-> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_vector_to_grid_data_type(
            u,
            v
        )

    # Call the c++ function
    c_flag = cpiv.validate_median(
        field_data,
        u_threshold,
        v_threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag


def local_normalized_median_val(
    u: np.ndarray,
    v: np.ndarray,
    threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
)-> np.ndarray:
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        if not isinstance(cpiv.grid_data):
            msg = 'field_data must be of type grid_data'
            return TypeError(msg)
    else:
        field_data = convert_vector_to_grid_data_type(
            u,
            v
        )

    # Call the c++ function
    c_flag = cpiv.validate_normalized_median(
        field_data,
        threshold
    )

    if return_raw:
        return c_flag
    else:
        flag = np.array(c_flag).astype(bool)
        return flag