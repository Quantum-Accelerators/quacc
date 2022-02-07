from pkg_resources import DistributionNotFound, get_distribution

try:
    __version__ = get_distribution("quacc").version
except DistributionNotFound:
    __version__ = ""
