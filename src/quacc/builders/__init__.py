"""
Modules related to builders, which take in an output from at least one calculation and
collate/analyze the data to produce valuable derived properties.

Builders do not call a compute-intensive (e.g. DFT) executable. Rather, they are
manipulations on pre-computed data. Most commonly, this is a many-to-one transformation,
and the builders can be run independently of the recipes if desired.
"""
