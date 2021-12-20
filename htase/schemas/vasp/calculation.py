from atomate2.vasp.schemas.task import TaskDocument
import os


def get_calculation_summary(dir_name=None, **kwargs):
    if dir_name is None:
        dir_name = os.getcwd()
    taskdoc = TaskDocument.from_directory(dir_name, **kwargs)

    return taskdoc
