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
    """Validate the dimensions of a numpy array.

    Parameters
    ----------
    ndim : int
        The number of dimensions the input arrays are expected to have.
    **kwargs : np.ndarray
        The arrays to be validated
    
    Returns
    -------
    None

    Raises
    ------
    ValueError
        if there is a dimension mismatch between expected and input dimensions.
    
    Examples
    --------
    import numpy as np
    arr1 = np.zeros([256, 256], dtype=float)
    arr2 = np.zeros([256, 256], dtype=float)

    check_nd(
        ndims=2,
        arr1=arr1,
        arr2=arr2
    )

    """
    for arg in kwargs:
        arr = kwargs[arg]

        if np.ndim(arr) != ndim:
            raise ValueError(f"{arg} is not a {ndim}D array")


def convert_to_image_type(
    image: Union[np.ndarray, cpiv.image_g_f32, cpiv.image_g_f64],
    double: bool=False
) -> Union[cpiv.image_g_f32, cpiv.image_g_f64]:
    """Convert array to c++ image type.

    Make sure the input image is the expected c++ data type before further
    use. For numpy arrays, the array is copied into the buffer of a c++
    image type. A c++ image type of a different floating point type may also
    be converted into the requested image type.

    Parameters
    ----------
    image : np.ndarray, image_g_f32, image_g_f64
        An array of pixel intensities representing a 2D image.
    double : bool
        If true, use 64-bit floating precision instead of the standard 32 bit
        precision.
    
    Returns
    -------
    image : image_g_f32, image_g_f64
        An array of pixel intensities representing a 2D image.
    
    Raises
    ------
    ValueError
        if there is a dimension mismatch between expected and input dimensions.
    TypeError
        if the image type is not of a supported conversion type.

    """
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
    """Convert window/overlap size to list type.

    Make sure the input window or overlap size is the expected type
    and form before use in the c++ l=library.

    Parameters
    ----------
    window : int, (int, int)

    Returns
    -------
    window : list[int, int]
        A list of window or overlap sizes
    
    Raises
    ------
    ValueError
        if there is a length mismatch between expected and input lengths.
    TypeError
        if the window type is not of a supported conversion type.

    """
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
    """Convert c++ data type to numpy arrays.

    Convert c++ grid_coords and grid_data types into numpy arrays for native
    Python use. Data can be chosen to copy the original data or simply
    reference the underlying data.

    Parameters
    ----------
    field_coords : grid_coords
        c++ data type which holds x and y coordinates in floating point form.
    field_data : grid_data
        c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form.
    copy : bool
        If true, copy the underlying data without providing a reference to
        the original data.
    
    Returns
    -------
    x : np.ndarray
        A 2D numpy array containing x-coordinates of a vector field.
    y : np.ndarray
        A 2D numpy array containing y-coordinates of a vector field.
    u : np.ndarray
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    s2n : np.ndarray
        A 2D numpy array containing signal to noise ratios of a vector field.
    p2p : np.ndarray
        A 2D numpy array containing peak to peak ratios of a vector field.
    
    Raises
    ------
    TypeError
        if the field_coords or field_data types are not of a supported conversion type.

    """
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
    """Convert u and v arrays to c++ grid_data type.

    Convert u and v 2D numpy arrays to c++ grid_data type for use in c++
    functions. 

    Parameters
    ----------
    u : np.ndarray
        A 2D numpy array containing u-components of a vector field.
    v : np.ndarray
        A 2D numpy array containing v-components of a vector field.
    
    Returns
    -------
    field_data : grid_data
        The c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form.
    
    Raises
    ------
    ValueError
        if there is a dimension mismatch between expected and input dimensions.
    TypeError
        if the u or v types are not of a supported conversion type.

    """
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
    """Convert s2n or p2p arrays to c++ grid_data type.

    Convert u and v 2D numpy arrays to c++ grid_data type for use in c++
    functions. 

    Parameters
    ----------
    s2n : np.ndarray
        A 2D numpy array containing signal to noise ratios of a vector field.
    p2p : np.ndarray
        A 2D numpy array containing peak to peak ratios of a vector field.
    
    Returns
    -------
    field_data : grid_data
        c++ data type which holds u, v, s2n, p2p, and other data in floating
        point form.
    
    Raises
    ------
    ValueError
        if there is a dimension mismatch between expected and input dimensions.
    TypeError
        if the s2n or p2p types are not of a supported conversion type.

    """

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