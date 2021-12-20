from atomate.vasp.schemas.calculation import Calculation


# TODO
def get_calculation_summary(dir_name, task_name, **kwargs):
    input_output_metadata = Calculation.from_vasp_files(task_name=task_name, **kwargs)
