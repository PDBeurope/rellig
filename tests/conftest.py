import logging
from unittest.mock import MagicMock

import pytest


@pytest.fixture
def logger():
    """A logger double that accepts any of the log calls used across the pipelines."""
    return MagicMock(spec=logging.Logger)
