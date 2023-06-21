# Parsl Setup

## Introduction

If you would prefer to use [Parsl](https://github.com/Parsl/parsl) to write and dispatch your workflows instead of Covalent, follow the instructions below. For additional details, refer to the full [Parsl documentation](https://parsl.readthedocs.io/en/stable/).

## Installation

In your activated Python environment, install Parsl (with support for a web dashboard) as follows:

```bash
pip install parsl[monitoring]
```

## Executor Configuration

```{note}
If you are just starting out, try running some test calculations locally first. Then come back and set up the relevant configuration files for your desired machines.
```

Out-of-the-box, Parsl will run on your local machine. To configure Parsl for the high-performance computing environment of your choice, refer to the [Configuration](https://parsl.readthedocs.io/en/stable/userguide/configuring.html) page in the documentation.
