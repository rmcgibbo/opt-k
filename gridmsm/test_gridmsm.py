import numpy as np
from gridmsm import _log_transmat_evidence

def test_log_transmat_evidence():
    "Example from Diaconis and Rolles (2006)"
    countsmat = np.array([
        [91, 160, 261, 108],
        [213, 351, 161, 249],
        [251, 224, 388, 201],
        [66, 239, 254, 152]])
    result = _log_transmat_evidence(countsmat, starting_state=3, prior=1)
    expected = np.log(2.166939224648291) - 1961 * np.log(10)
    np.testing.assert_almost_equal(expected, result)
