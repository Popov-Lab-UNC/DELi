"""configure testing environment for deli"""

import os

import pytest


@pytest.fixture(autouse=True, scope="session")
def set_up_test_env():
    """Set up the test environment by setting the necessary environment variables."""
    os.environ["DELI_CONFIG"] = os.path.join(os.path.dirname(__file__), "test_data", ".deli")
    os.environ["DELI_DATA_DIR"] = os.path.join(
        os.path.dirname(__file__), "test_data", "test_deli_data_dir"
    )
