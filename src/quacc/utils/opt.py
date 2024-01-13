class ASEOptimizerSuperCharged:
    def __init__(self):
        self.tasks = {}
        self.every = {}
        self.results = {}

    def step(self):
        super().step()
        for task in self.tasks:
            if self.nstep % self.every[task] != 0:
                continue
            self.results[task] = self.tasks[task]()

    def attach_task(self, task, every = 1):
        self.task[task.__name__] = task
        self.every[task.__name__] = every