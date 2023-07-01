# Comparing Syntax

## Introduction

All of the recommended workflow solutions have a similar decorator-based syntax for compute tasks and workflows. Here, we highlight these differences.

## Simple Workflow

### No Workflow Engine

```python
def add(a, b):
    return a + b

def mult(a, b):
    return a * b

def workflow(a, b, c):
    return mult(add(a, b), c)

result = workflow(1, 2, 3) # 9
```

### Covalent

```python
import covalent as ct

@ct.electron
def add(a, b):
    return a + b

@ct.electron
def mult(a, b):
    return a * b

@ct.lattice
def workflow(a, b, c):
    return mult(add(a, b), c)

# Locally
result = workflow(1, 2, 3) # 9

# Dispatched
dispatch_id = ct.dispatch(workflow)(1, 2, 3)
result = ct.get_result(dispatch_id, wait=True) # 9
```

### Parsl

```python
import parsl
from parsl import python_app

@python_app
def add(a, b):
    return a + b

@python_app
def mult(a, b):
    return a * b

def workflow(a, b, c):
    future1 = add(a, b)
    return mult(future1.result(), c)

result = workflow(1, 2, 3).result() # 9
```

### Prefect

```python
from prefect import task, flow

@task
def add(a, b):
    return a + b

@task
def mult(a, b):
    return a * b

@flow
def workflow(a, b, c):
    future1 = add.submit(a, b)
    return mult.submit(future1.result(), c)

result = workflow(1, 2, 3).result() # 9
```

### Jobflow

```python
from jobflow import job, Flow, run_locally

@job
def add(a, b):
    return a + b

@job
def mult(a, b):
    return a * b

job1 = add(1, 2)
job2 = mult(job1.output, 3)
flow = Flow([job1, job2], output=job2.output)

responses = run_locally(flow)
result = responses[job2.uuid][1].output # 9
```

## Dynamic Workflow

### No Workflow Engine

```python
def add(a, b):
    return a + b

def make_more(val):
    return [val] * 3

def workflow(a, b, c):
    result1 = add(a, b)
    result2 = make_more(result1)
    return [add(val, c) for val in result2]

result = workflow(1, 2, 3) # [6, 6, 6]
```

### Covalent

```python
import covalent as ct

@ct.electron
def add(a, b):
    return a + b

@ct.electron
def make_more(val):
    return [val] * 3

@ct.lattice
def workflow(a, b, c):

    @ct.electron
    @ct.lattice
    def _add_distributed(vals):
        return [add(val, c) for val in vals]

    result1 = add(a, b)
    result2 = make_more(result1)
    return _add_distributed(result2)

# Locally
result = workflow(1, 2, 3) # [6, 6, 6]

# Dispatched
dispatch_id = ct.dispatch(workflow)(1, 2, 3)
result = ct.get_result(dispatch_id, wait=True) # [6, 6, 6]
```

### Parsl

```python
import parsl
from parsl import join_app, python_app

@python_app
def add(a, b):
    return a + b

@python_app
def make_more(val):
    return [val] * 3

@join_app
def workflow(a, b, c):

    def _add_distributed(vals):
        return [add(val, c) for val in vals]

    result1 = add(a, b)
    result2 = make_more(result1.result())
    return _add_distributed(result2.result())

result = workflow(1, 2, 3).result() # [6, 6, 6]
```

### Prefect

### Jobflow

```python
from jobflow import job, Flow, Response, run_locally

@job
def add(a, b):
    return a + b

@job
def make_more(val):
    return [val] * 3

@job
def add_distributed(vals, c):
    jobs = []
    for val in vals:
        jobs.append(add(val, c))

    flow = Flow(jobs)
    return Response(detour=flow)

job1 = add(1, 2)
job2 = make_more(job1.output)
job3 = add_distributed(job2.output, 3)
flow = Flow([job1, job2, job3])

responses = run_locally(flow) # [6, 6, 6] in final 3 jobs
```
