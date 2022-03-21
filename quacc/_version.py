"""Version handling"""
from pkg_resources import get_distribution

__version__ = get_distribution("quacc").version
