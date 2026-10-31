import numpy as np
import pyopenpivcore as cpiv
import pytest


CANDIDATES = (
    'u1', 'v1', 's2n1', 'peak1',
    'u2', 'v2', 's2n2', 'peak2',
    'u3', 'v3', 's2n3', 'peak3',
)


@pytest.mark.parametrize('use_size', [False, True])
def test_initialization(use_size):
    field = cpiv.grid_data(cpiv.size(4, 3)) if use_size else cpiv.grid_data(4, 3)

    for name in CANDIDATES:
        values = np.asarray(getattr(field, name))
        assert values.shape == (3, 4)
        assert values.dtype == np.float64
        assert np.isnan(values).all()

    for name in ('u', 'v', 's2n', 'peak', 'p2p', 'flag', 'peak_ind'):
        values = np.asarray(getattr(field, name))
        assert values.shape == (3, 4)
        assert (values == 0).all()
    assert np.asarray(field.peak_ind).dtype == np.uint8


def test_independent_storage_and_copy():
    field = cpiv.grid_data(4, 3)
    for value, name in enumerate(CANDIDATES, 1):
        np.asarray(getattr(field, name))[0, 0] = value
    np.asarray(field.peak_ind)[0, 0] = 2
    np.asarray(field.flag)[0, 0] = 2

    copied = cpiv.grid_data(field)
    for name in ('u', 'v', 's2n', 'peak', 'p2p'):
        np.asarray(getattr(field, name))[:] = 100.0
        assert (np.asarray(getattr(copied, name)) == 0).all()

    for value, name in enumerate(CANDIDATES, 1):
        original = np.asarray(getattr(field, name))
        duplicate = np.asarray(getattr(copied, name))
        assert original[0, 0] == value
        assert duplicate[0, 0] == value
        assert np.isnan(duplicate[0, 1])
        assert not np.shares_memory(original, duplicate)
        duplicate[0, 0] = -value
        assert original[0, 0] == value

    assert np.asarray(copied.peak_ind)[0, 0] == 2
    assert np.asarray(copied.flag)[0, 0] == 2
    np.asarray(copied.peak_ind)[0, 0] = 3
    np.asarray(copied.flag)[0, 0] = 1
    assert np.asarray(field.peak_ind)[0, 0] == 2
    assert np.asarray(field.flag)[0, 0] == 2
