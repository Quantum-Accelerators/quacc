# Overview

## Approving Jenkins Builds

When a new pull request is opened in the project and the author of the pull request isn't white-listed, the builder will ask "Can one of the admins verify this patch?".

- "ok to test" to accept this pull request for testing
- "add to whitelist" to add the author to the whitelist
- "run jenkins" to run the Jenkins build once approved

## Manual Upload to PyPi

To upload to PyPI manually in the case of a GitHub Actions failure:

```python
pip install -e .
python -m build
twine check dist/*
twine upload dist/*
```
