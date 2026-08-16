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


def _range_val_cpiv(
    field_data: cpiv.grid_data,
    u_threshold: Tuple[float, float],
    v_threshold: Tuple[float, float]
)-> np.ndarray:
    if not isinstance(field_data, cpiv.grid_data):
        msg = 'field_data must be of type grid_data'
        return TypeError(msg)
        
    # Call the c++ function
    c_flag = cpiv.validate_range(
        field_data,
        u_threshold,
        v_threshold
    )

    flag = np.array(c_flag).astype(bool)
    
    return flag 


def _range_val_array(
    u: np.ndarray,
    v: np.ndarray,
    u_threshold: Tuple[float, float],
    v_threshold: Tuple[float, float]
)-> np.ndarray:
    field_data = convert_vector_to_grid_data_type(
        u,
        v
    )

    return _range_val_cpiv(
        field_data,
        u_threshold,
        v_threshold
    )


def range_val(
    u: np.ndarray,
    v: np.ndarray,
    u_threshold: Tuple[float, float],
    v_threshold: Tuple[float, float],
    field_data: cpiv.grid_data=None
) -> np.ndarray:
    """Eliminate spurious vectors with a global threshold.

    This validation method tests for the spatial consistency of the data
    and outliers vector are flagged as 1 (True) if at least one of the two
    velocity components is out of a specified global range.

    Parameters
    ----------
    u : np.ndarray or cpiv.grid_data
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    u_thresholds : tuple of floats
        u_thresholds = (u_min, u_max). If ``u<u_min`` or ``u>u_max``
        the vector is treated as an outlier.
    v_thresholds : tuple of floats
        ``v_thresholds = (v_min, v_max)``. If ``v<v_min`` or ``v>v_max``
        the vector is treated as an outlier.
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, u and v must be set to None.
    
    """
    # If field_data is passed, we can skip making copies of the original data and us e it directly
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        return _range_val_cpiv(
            field_data,
            u_threshold,
            v_threshold
        )
            
    else:
        return _range_val_array(
            u,
            v,
            u_threshold,
            v_threshold
        )
        

def _global_std_val_cpiv(
    field_data: cpiv.grid_data,
    std_threshold: float
)-> np.ndarray:
    if not isinstance(field_data, cpiv.grid_data):
        msg = 'field_data must be of type grid_data'
        return TypeError(msg)
        
    # Call the c++ function
    c_flag = cpiv.validate_z_score(
        field_data,
        std_threshold
    )

    flag = np.array(c_flag).astype(bool)
    
    return flag 


def _global_std_val_array(
    u: np.ndarray,
    v: np.ndarray,
    std_threshold: float
)-> np.ndarray:
    field_data = convert_vector_to_grid_data_type(
        u,
        v
    )

    return _global_std_val_cpiv(
        field_data,
        std_threshold
    )
    
def global_std_val(
    u: np.ndarray,
    v: np.ndarray,
    std_threshold: float=3.0,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
)-> np.ndarray:
    """Eliminate spurious vectors with a global threshold defined by the
    standard deviation.

    This validation method tests for the spatial consistency of the data
    and outliers vector are flagges as 1 (True) if at least one of the two
    velocity components is out of a specified global range.

    Parameters
    ----------
    u : np.ndarray
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    std_threshold: float
        If the length of the vector (actually the sum of squared components) is
        larger than std_threshold times standard deviation of the flow field,
        then the vector is treated as an outlier. [default = 3]
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, u and v must be set to None.
    return_raw : bool
        If true, return the raw c++ data. Otherwise, copy the data as a numpy
        ndarray of dtype bool.
    
    """
    if field_data is not None:
        if u is not None or v is not None:
            msg = 'u and v must be set to None when using field_data'
            raise ValueError(msg)

        return _global_std_val_cpiv(
            field_data,
            u_threshold,
            v_threshold
        )
            
    else:
        return _global_std_val_array(
            u,
            v,
            u_threshold,
            v_threshold
        )


def sig2noise_val(
    sig: np.ndarray,
    threshold: float,
    field_data: cpiv.grid_data=None,
    return_raw: bool=False
) -> np.ndarray:
    """Marks spurious vectors if signal to noise ratio is below a specified threshold.

    This function validates velocity vectors based on the signal-to-noise ratio
    from the cross-correlation function. Vectors with a signal-to-noise ratio
    below the specified threshold are marked as outliers.

    Parameters
    ----------
    sig : ndarray
        A 2D numpy array containing signal to noise ratios of a vector field.
    threshold : float, optional
        The signal to noise ratio threshold value. Vectors with s2n < threshold
        will be marked as outliers.
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, sig must be set to None.
    return_raw : bool
        If true, return the raw c++ data. Otherwise, copy the data as a numpy
        ndarray of dtype bool.
    
    """
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
    """Marks spurious vectors if peak to peak ratio is below a specified threshold.

    This function validates velocity vectors based on the signal-to-noise ratio
    from the cross-correlation function. Vectors with a signal-to-noise ratio
    below the specified threshold are marked as outliers.

    Parameters
    ----------
    sig : ndarray
        A 2D numpy array containing peak to peak ratios of a vector field.
    threshold : float, optional
        The peak to peak ratio threshold value. Vectors with p2p < threshold
        will be marked as outliers.
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, sig must be set to None.
    return_raw : bool
        If true, return the raw c++ data. Otherwise, copy the data as a numpy
        ndarray of dtype bool.
    
    """
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
    """Eliminate spurious vectors with a local difference threshold.

    This validation method tests for the spatial consistency of the data
    through calculating the difference between a vector and its neighbors.
    Outliers are flagged as 1 (True) when at least half of the neighborly
    differences are above the specified threshold.

    Parameters
    ----------
    u : np.ndarray
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    u_thresholds : float
        The threshold value for the u-component of a vector field.
    v_thresholds : float
        The threshold value for the v-component of a vector field.
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, u and v must be set to None.
    return_raw : bool
        If true, return the raw c++ data. Otherwise, copy the data as a numpy
        ndarray of dtype bool.
    
    """
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
    """Eliminate spurious vectors with a local median threshold.

    This validation method tests for the spatial consistency of the data
    through calculating the median of the neighboring vectors and subtracting
    the vector of interest from that median. Outliers are flagged as 1 (True) 
    when either component thresholds are exceeded.

    Parameters
    ----------
    u : np.ndarray
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    u_thresholds : float
        The threshold value for the u-component of a vector field.
    v_thresholds : float
        The threshold value for the v-component of a vector field.
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, u and v must be set to None.
    return_raw : bool
        If true, return the raw c++ data. Otherwise, copy the data as a numpy
        ndarray of dtype bool.
    
    """
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
    """Eliminate spurious vectors with a local normalized median threshold.

    This validation method tests for the spatial consistency of the data
    through calculating the normalized median test which basically normalized
    the residuals to make the filter more robust. Outliers are flagged as 1
    (True) values are above the specified threshold.

    Parameters
    ----------
    u : np.ndarray
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    threshold : float
        The threshold value for the displacement vector of a vector field.
    field_data : grid_data, optional
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form. If specified, u and v must be set to None.
    return_raw : bool
        If true, return the raw c++ data. Otherwise, copy the data as a numpy
        ndarray of dtype bool.
    
    """
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