# Prefect Setup

## Introduction

If you would prefer to use [Prefect](https://www.prefect.io/) to write and dispatch your workflows instead of Covalent, follow the instructions below. For additional details, refer to the full [Prefect documentation](https://docs.prefect.io/).

## Installation

1. In your activated Python environment, install Prefect via `pip install prefect prefect-dask dask-jobqueue`.
2. Make an account on [Prefect Cloud](https://app.prefect.cloud/).
3. Make an [API Keý](https://docs.prefect.io/cloud/users/api-keys/) and (optionally) store it in a `PREFECT_API_KEY` environment variable (e.g. in your `~/.bashrc`).
4. Run `prefect cloud login` from the command-line and enter your API key (or use the browser, if possible).
