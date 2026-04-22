from typing import Tuple, Union, List
from copy import deepcopy

import numpy as np

from . import pyopenpivcore as cpiv


__all__ = [
    'check_nd',
    'convert_to_image_type',
    'convert_to_list_type',
    'convert_cxx_data_to_numpy',
    'convert_vector_to_grid_data_type',
    'convert_stats_to_grid_data_type'
]


SUPPPORTED_IMAGE_TYPES = (
    cpiv.image_g_f32,
    cpiv.image_g_f64,
    np.ndarray
)

SUPPORTED_WINDOW_TYPES = (
    tuple,
    list,
    int
)


def check_nd(
    ndim: int=2, 
    **kwargs
) -> None:
    for arg in kwargs:
        arr = kwargs[arg]

        if np.ndim(arr) != ndim:
            raise ValueError(f"{arg} is not a {ndim}D array")


def convert_to_image_type(
    image: Union[np.ndarray, cpiv.image_g_f32, cpiv.image_g_f64],
    double: bool=False
) -> Union[cpiv.image_g_f32, cpiv.image_g_f64]:
    if double:
        image_type = cpiv.image_g_f64
    else:
        image_type = cpiv.image_g_f32
        

    if not isinstance(image, SUPPPORTED_IMAGE_TYPES):
        msg = f'Image is not in a supported type (type = {type(image)})'
        raise TypeError(msg)
        
    if isinstance(image, np.ndarray):
        # We can only convert 2D numpy arrays to an image
        check_nd(
            ndim=2, 
            image=image
        )
        
        image_reformatted = image_type(list(image.shape[::-1]))
        image_buf = np.array(image_reformatted, copy=False)
        image_buf[:] = image

        return image_reformatted

    # Cast the image to new image type if it is not a Numpy array
    if not isinstance(image, image_type):
        image_reformatted = image_type(image.size())
        image_buf_old = np.array(image, copy=False)
        image_buf_new = np.array(image_reformatted, copy=False)
        image_buf_new[:] = image_buf_old

        return image_reformatted

    # If no conversions are needed, simply return the original image
    return image


def convert_to_list_type(
    window: Union[Tuple[int, int], int]
) -> List[[int]*2]:
    MAX_LEN = 2
    
    # window must be able to be converted to two-element list for c++ to accept
    if not isinstance(window, SUPPORTED_WINDOW_TYPES):
        msg = f'window/overlap size is not in a supported type (type = {type(window)})'
        raise TypeError(msg)

    # Deep copy so we don't accidentaly change the original variable
    window_fixed = deepcopy(window)

    # If window is already list, make sure it is valid
    if isinstance(window, list):
        if len(window) != MAX_LEN:
            msg = f'list must contain two elements. Got {len(window)}'
            return ValueError(msg)

    #  For tuples, check to see if we can convert to two-element list
    if isinstance(window, tuple):
        if len(window) != MAX_LEN:
            msg = f'tuple must contain two elements. Got {len(window)}'
            return ValueError(msg)

            window_fixed = list(window)

    # If scalar, go ahead and make a list
    if isinstance(window, int):
        window_fixed = [window, window]
        
    return window_fixed


def convert_cxx_data_to_numpy(
    field_coords: cpiv.grid_coords,
    field_data: cpiv.grid_data,
    copy: bool=True
) -> List[[np.ndarray]*6]:
    if not isinstance(field_coords, cpiv.grid_coords):
        msg = f'field_coords must be of type grid_coords. Got {type(field_coords)}'
        raise TypeError(msg)

    if not isinstance(field_data, cpiv.grid_data):
        msg = f'field_data must be of type grid_data. Got {type(field_data)}'
        raise TypeError(msg)

    # In c++, grid coords is a 2D matrix containing [x,y] pairs. We can unpack them
    # using a simple non-coping conversion by exploiting the underlying py::buffer_protocol.
    coord_arr = np.array(field_coords, copy=copy)
    x = coord_arr[:,:,0]
    y = coord_arr[:,:,1]

    # For u, v, s2n, ... data, we have to select the member access specifier to the struct
    # of 2D matrices and then convert to numpy array following similar py::buffer_protocol.
    # Note: The c++ software can store many peaks and subsequent u,v, ... data. We are only
    # interested in the first peak and its associated data.
    u = np.array(field_data.u, copy=copy)
    v = np.array(field_data.v, copy=copy)
    s2n = np.array(field_data.s2n, copy=copy)
    p2p = np.array(field_data.p2p, copy=copy)
    
    return x, y, u, v, s2n, p2p


def convert_vector_to_grid_data_type(
    u: np.ndarray,
    v: np.ndarray
) -> cpiv.grid_data:
    # Make sure we are using numpy ndarrays so things don't break
    if not isinstance(u, np.ndarray) or \
        not isinstance(v, np.ndarray):
        msg = 'u and v must be numpy ndarrays'
        raise TypeError(msg)
    
    # Make sure u and v are both 2D so we an get vector field shape
    check_nd(
        ndim=2,
        u=u,
        v=v,
    )

    common_shape = u.shape

    if v.shape != common_shape:
        msg = f'u and v must be of same shape (got {u.shape} and {v.shape})'
        raise ValueError(msg)

    field_data = cpiv.grid_data(list(u.shape)[::-1])
    u_buf = np.array(field_data.u, copy=False)
    v_buf = np.array(field_data.v, copy=False)

    u_buf[:] = u
    v_buf[:] = v

    return field_data


def convert_stats_to_grid_data_type(
    s2n: np.ndarray=None,
    p2p: np.ndarray=None
) -> cpiv.grid_data:
    data_shape = None
    
    # Make sure we are using numpy ndarrays so things don't break
    if s2n is not None:
        if not isinstance(s2n, np.ndarray):
            msg = 's2n must be numpy ndarray'
            raise TypeError(msg)
        
        # Make sure s2n is a 2D array since core::grid_data's shape would be inferred from it
        check_nd(
            ndim=2,
            s2n=s2n
        )

        data_shape = list(s2n.shape)[::-1]
            
    if p2p is not None:
        if not isinstance(p2p, np.ndarray):
            msg = 'p2p must be numpy ndarray'
            raise TypeError(msg)
        
        # Make sure s2n is a 2D array since core::grid_data's shape would be inferred from it
        check_nd(
            ndim=2,
            p2p=p2p
        )

        data_shape = list(p2p.shape)[::-1]

    # This happens when neither s2n or p2p is sent (e.g., empty function call)
    if data_shape is None:
        msg = 'Either s2n or p2p must be specified'
        raise ValueError(msg)

    # Now build vector correlation statistics
    field_data = cpiv.grid_data(data_shape)

    if s2n is not None:
        s2n_buf = np.array(field_data.s2n, copy=False)
        s2n_buf[:] = s2n
    
    if p2p is not None:
        p2p_buf = np.array(field_data.p2p, copy=False)
        p2p_buf[:] = p2p

    return field_data