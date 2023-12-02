import os
from pathlib import Path


def pytest_sessionstart():
    file_dir = Path(__file__).parent
    os.environ["QUACC_CONFIG_FILE"] = str(file_dir / ".quacc.yaml")
