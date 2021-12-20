from htase.schemas.vasp.calculation import get_calculation_summary
import os


def get_results(dir_name=None, **kwargs):
    if dir_name is None:
        dir_name = os.getcwd()

    results = get_calculation_summary(dir_name, **kwargs)
    # We might want to do more later, but for now this is simple

    return results
