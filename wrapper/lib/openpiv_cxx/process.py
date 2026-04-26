from typing import Union, Tuple, List
from . import pyopenpivcore as cpiv
from .parsing import convert_to_image_type, convert_to_list_type, convert_cxx_data_to_numpy

import numpy as np


__all__ = [
    'first_pass',
    'multi_pass'
]


def first_pass(
    image_a: Union[np.ndarray, cpiv.image_g_f32],
    image_b: Union[np.ndarray, cpiv.image_g_f32],
    window_size: Union[int, Tuple[int,int]]=32,
    overlap: Union[int, Tuple[int,int]]=16,
    correlation_method: str='circular',
    centered: bool=False,
    limit_search: bool=False,
    parse_output: bool=True,
    threads: int=1
)-> Union[Tuple[cpiv.grid_coords, cpiv.grid_data], List[[np.ndarray]*6]]: 
    """Standard PIV cross correlation algorithm.

    Standard PIV cross correlation algorithm where the mean is substracted
    from each interrogation window and the correlation plane is divided by
    the autocorrelation of a weighting function (currently uniform) to 
    normalize the correlation peaks into a -1 to 1 range.

    Parameters
    ----------
    frame_a : ndarray
        An array of pixel intensities representing a 2D image.
    frame_b : ndarray
        An array of pixel intensities representing a 2D image.
    window_size : int or list of int
        The size of the interrogation window. [default: 32 pix]
    overlap : int or list of int
        The number of pixels by which two adjacent windows overlap.
        [default: 16 pix]
    correlation_method : string
        The type of FFT cross correlation to perform
        'circular'
            The FFT signals are assumed periodic. This is not typically
            very robust against image background noise and intensity changes
            since frequencies may wrap around the correlation matrix.
        
        'linear'
            The interrogation windows is padded by two times the window size
            to effectively remove periodic signals. For noisy image data,
            this makes the cross correlation peak more well defined, but
            this mode typically takes twice as long to perform calculations.
    centered : bool, optional
        If true, center the interrogation windows in the image such that there
        is an equal offset on all four corners of the vector field from the
        image borders.
    limit_search : bool, optional
        If true, only search for the correlation peak in half of the
        correlation matrix. This aligns with the 1/4 rule where 
        displacements should never really exceed 1/4 the interrogation window
        due to signal loss from periodic tendencies of FFTs. For window
        deformation algorithms, it makes sense to always enable this to
        increase performance during processing.
    parse_output : bool, optional
        If true, parse the c++ data types into numpy arrays to return the
        typical x, y, u, v data returned by other PIV software. If false,
        return the raw grid_coords and grid_data c++ types containing the
        vector field and grid coordinates. [defualt: True]
    threads : int, optional
        The amound of threads to use during the cross-correlation process.
        Thread counts less than one automatically set the thread count to the
        maximum supported threads. Thread counts greater than one allow for
        manual selection of the amount of threads to use during processing.
        Thread counts of one sets the processing state to serial.
        [defualt: 1]

    Returns
    -------
    field_coords : grid_coords
        Returned if ``parse_output`` is ``False``. The c++ data type which holds
        x and y coordinates in floating point form.
    field_data : grid_data
        Returned if ``parse_output`` is ``False``. The c++ data type which holds
        u, v, s2n, p2p, and other vector data in floating point form.
    x : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        x-coordinates of a vector field.
    y : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        y-coordinates of a vector field.
    u : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        u-components of a vector field.
    v : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        v-components of a vector field.
    s2n : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        signal to noise ratios of a vector field.
    p2p : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        peak to peak ratios of a vector field.
    
    Raises
    ------
    ValueError
        if there is a dimension mismatch between expected and input dimensions.
    TypeError
        if some data is not a supported conversion type.
    RuntimeError
        if any error is encountered in the c++ library.

    """ 
    # Images must be of core::image c++ types
    image_a = convert_to_image_type(image_a)
    image_b = convert_to_image_type(image_b)

    if correlation_method == 'circular':
        zero_pad = False
    elif correlation_method == 'linear':
        zero_pad = True
    else:
        msg = f'Correlation method {correlation_method} is not supported'
        raise ValueError(msg)

    # In the c++ library, step toggles between OpenPIV-style overlap and
    # the actual spacing of the interrogation windows. For instance,
    # a window size of 32 and overlap of 8 would result in a grid spacing
    # of 24 if step = False, and 8 if step = True.
    step = False

    # Make sure window_size and overlap are lists of ints.
    # List of ints can be converted to core::size c++ data types
    window_size = convert_to_list_type(window_size)
    overlap = convert_to_list_type(overlap)

    # Call the c++ function
    field_coords, field_data = cpiv.process_images_standard(
        image_a,
        image_b,
        window_size,
        overlap,
        step=step,
        zero_pad=zero_pad,
        centered=centered,
        limit_search=limit_search,
        threads=threads
    )

    if not parse_output:
        return field_coords, field_data
    
    # Parse the c++ data types and make numpy array copies of the data
    return convert_cxx_data_to_numpy(
        field_coords,
        field_data,
        copy=True
    )


def multi_pass(
    image_a: Union[np.ndarray, cpiv.image_g_f32],
    image_b: Union[np.ndarray, cpiv.image_g_f32],
    window_size: Union[int, Tuple[int,int]]=32,
    overlap: Union[int, Tuple[int,int]]=16,
    correlation_method: str='circular',
    centered: bool=False,
    limit_search: bool=False,
    parse_output: bool=True,
    threads: int=1
)-> List[[np.ndarray]*6]: 
    """Standard PIV cross correlation algorithm.

    Standard PIV cross correlation algorithm where the mean is substracted
    from each interrogation window and the correlation plane is divided by
    the autocorrelation of a weighting function (currently uniform) to 
    normalize the correlation peaks into a -1 to 1 range.

    Parameters
    ----------
    frame_a : ndarray
        An array of pixel intensities representing a 2D image.
    frame_b : ndarray
        An array of pixel intensities representing a 2D image.
    window_size : int or list of int
        The size of the interrogation window. [default: 32 pix]
    overlap : int or list of int
        The number of pixels by which two adjacent windows overlap.
        [default: 16 pix]
    correlation_method : string
        The type of FFT cross correlation to perform
        'circular'
            The FFT signals are assumed periodic. This is not typically
            very robust against image background noise and intensity changes
            since frequencies may wrap around the correlation matrix.
        
        'linear'
            The interrogation windows is padded by two times the window size
            to effectively remove periodic signals. For noisy image data,
            this makes the cross correlation peak more well defined, but
            this mode typically takes twice as long to perform calculations.
    centered : bool, optional
        If true, center the interrogation windows in the image such that there
        is an equal offset on all four corners of the vector field from the
        image borders.
    limit_search : bool, optional
        If true, only search for the correlation peak in half of the
        correlation matrix. This aligns with the 1/4 rule where 
        displacements should never really exceed 1/4 the interrogation window
        due to signal loss from periodic tendencies of FFTs. For window
        deformation algorithms, it makes sense to always enable this to
        increase performance during processing.
    parse_output : bool, optional
        If true, parse the c++ data types into numpy arrays to return the
        typical x, y, u, v data returned by other PIV software. If false,
        return the raw grid_coords and grid_data c++ types containing the
        vector field and grid coordinates. [defualt: True]
    threads : int, optional
        The amound of threads to use during the cross-correlation process.
        Thread counts less than one automatically set the thread count to the
        maximum supported threads. Thread counts greater than one allow for
        manual selection of the amount of threads to use during processing.
        Thread counts of one sets the processing state to serial.
        [defualt: 1]

    Returns
    -------
    field_coords : grid_coords
        Returned if ``parse_output`` is ``False``. The c++ data type which holds
        x and y coordinates in floating point form.
    field_data : grid_data
        Returned if ``parse_output`` is ``False``. The c++ data type which holds
        u, v, s2n, p2p, and other vector data in floating point form.
    x : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        x-coordinates of a vector field.
    y : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        y-coordinates of a vector field.
    u : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        u-components of a vector field.
    v : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        v-components of a vector field.
    s2n : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        signal to noise ratios of a vector field.
    p2p : ndarray
        Returned if ``parse_output`` is ``True``. A 2D numpy array containing
        peak to peak ratios of a vector field.
    
    Raises
    ------
    ValueError
        if there is a dimension mismatch between expected and input dimensions.
    TypeError
        if some data is not a supported conversion type.
    RuntimeError
        if any error is encountered in the c++ library.

    """ 
    # Images must be of core::image c++ types
    image_a = convert_to_image_type(image_a)
    image_b = convert_to_image_type(image_b)

    if correlation_method == 'circular':
        zero_pad = False
    elif correlation_method == 'linear':
        zero_pad = True
    else:
        msg = f'Correlation method {correlation_method} is not supported'
        raise ValueError(msg)

    # In the c++ library, step toggles between OpenPIV-style overlap and
    # the actual spacing of the interrogation windows. For instance,
    # a window size of 32 and overlap of 8 would result in a grid spacing
    # of 24 if step = False, and 8 if step = True.
    step = False

    # Make sure window_size and overlap are lists of ints.
    # List of ints can be converted to core::size c++ data types
    window_size = convert_to_list_type(window_size)
    overlap = convert_to_list_type(overlap)

    # Call the c++ function
    field_coords, field_data = cpiv.process_images_standard(
        image_a,
        image_b,
        window_size,
        overlap,
        step=step,
        zero_pad=zero_pad,
        centered=centered,
        limit_search=limit_search,
        threads=threads
    )

    if not parse_output:
        return field_coords, field_data
    
    # Parse the c++ data types and make numpy array copies of the data
    return convert_cxx_data_to_numpy(
        field_coords,
        field_data,
        copy=True
    )