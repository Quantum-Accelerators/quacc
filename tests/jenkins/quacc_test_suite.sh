#!/bin/bash
pytest tests/local/recipes/gaussian_recipes tests/local/recipes/gulp_recipes tests/local/recipes/orca_recipes --noconftest --cov=quacc --cov-report=xml
