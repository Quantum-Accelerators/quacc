"""Recipes for VASP."""

from quacc.recipes.vasp.core import relax_job, static_job
from quacc.recipes.vasp.md import md_job

__all__ = ["relax_job", "static_job", "md_job"]

