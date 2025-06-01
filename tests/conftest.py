"""configure testing environment for deli"""

import os

import pytest

from deli.configure import set_deli_data_dir


@pytest.fixture(autouse=True, scope="function")
def set_up_test_env():
    """Set up the test environment by setting the necessary environment variables."""
    os.environ["DELI_CONFIG"] = os.path.join(os.path.dirname(__file__), "test_data", ".deli")
    os.environ["DELI_DATA_DIR"] = os.path.join(
        os.path.dirname(__file__), "test_data", "test_deli_data_dir"
    )
    set_deli_data_dir(os.path.join(os.path.dirname(__file__), "test_data", "test_deli_data_dir"))
