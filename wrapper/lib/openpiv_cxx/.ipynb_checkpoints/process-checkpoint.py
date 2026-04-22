from . import pyopenpivcore as cpiv
from .parsing import convert_to_image_type, convert_to_list_type, convert_cxx_data_to_numpy


__all__ = [
    'first_pass'
]


def first_pass(
    image_a,
    image_b,
    window_size,
    overlap,
    correlation_method='circular',
    centered=False,
    parse_output: bool=True,
    threads=1
)-> List[[np.ndarray]*6]:  
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

    # limit_search is a performance optimization where the correlation peaks are only
    # searched for 1/2 of the correlation plane. This aligns with the 1/4 rule where 
    # displacements should never really exceet 1/4 the interrogation window due to
    # signal loss from periodic tendencies of the Fast Fourier Transform.
    limit_search=False

    # Make sure window_size and overlap are lists of ints so we don't get a overload error
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