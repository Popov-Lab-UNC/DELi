"""configure testing environment for deli"""

import os
from pathlib import Path

import pytest

from deli.configure import init_deli_config, set_deli_data_dir


@pytest.fixture(autouse=True, scope="function")
def set_up_test_env():
    """Set up the test environment by setting the necessary environment variables."""
    os.environ["DELI_CONFIG"] = os.path.join(os.path.dirname(__file__), "test_data", ".deli")
    os.environ["DELI_DATA_DIR"] = os.path.join(os.path.dirname(__file__), "test_data", "test_deli_data_dir")

    init_deli_config(Path(os.path.join(os.path.dirname(__file__), "test_data", ".deli")), fail_on_exist=False)
    set_deli_data_dir(os.path.join(os.path.dirname(__file__), "test_data", "test_deli_data_dir"))
