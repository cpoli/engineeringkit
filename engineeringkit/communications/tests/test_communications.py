import numpy as np

import engineeringkit.communications as m


def test_capacity_and_db():
    assert np.isclose(m.from_db(m.db(10)), 10)
    assert np.isclose(m.shannon_capacity(1, 1), 1)


def test_ber():
    assert m.ber_bpsk_awgn(10) < m.ber_bpsk_awgn(1)
