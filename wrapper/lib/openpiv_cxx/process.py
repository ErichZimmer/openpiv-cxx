from .pyopenpivcore import process_images_standard, image_g_f32

def first_pass(
    image_a,
    image_b,
    window_size,
    overlap,
    correlation_method='circular',
    centered=False,
    threads=1
):  
    if not isinstance(image_a, image_g_f32):
        image_a_temp = image_a.copy()
        image_a = image_g_f32(list(image_a_temp.shape[::-1]))
        image_a_buf = np.array(image_a, copy=False)
        image_a_buf[:] = image_a_temp

    if not isinstance(image_b, image_g_f32):
        image_b_temp = image_b.copy()
        image_b = image_g_f32(list(image_b_temp.shape[::-1]))
        image_b_buf = np.array(image_b, copy=False)
        image_b_buf[:] = image_b_temp

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
    # searched for 1/2 of the correlation plane.
    limit_search=False
    
    field_coords, field_data = process_images_standard(
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

    # Convert to NumPy and unpack coords
    coords = np.array(field_coords, copy=False)

    x = coords[:,:, 0]
    y = coords[:,:, 1]

    # Unpack u, v, and s2n data for first peak
    u   = np.array(field_data.u, copy=False)
    v   = np.array(field_data.v, copy=False)
    s2n = np.array(field_data.s2n, copy=False)

    return x, y, u, v, s2n