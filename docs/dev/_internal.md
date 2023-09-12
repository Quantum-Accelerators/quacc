To upload to PyPI:

```python
pip install -e .
python -m build
twine check dist/*
twine upload dist/*
```
