"""Write slurm job submission scripts."""

from dataclasses import dataclass

from .io_base import PathLike


@dataclass
class SlurmWriter:
    """
    Simple dataclass wrapper to write slurm sbatch submission scripts for HPC computations.
    
    Attributes
    ----------
    shebang: str, optional
        Path to the shell executable. Defaults to "/bin/bash". The "#!" shebang is
        added automatically.

    jobname: str, optional
        Name of the job in the job queue. Defaults to "jobname".

    output: str, optional
        Name of the slurm output file capturing stdout output.
        Defaults to "slurm_%j.out", where %j is replaced with Job ID at execution time.
        See Slurm documentation for more details about '%' shortcuts behaviors.

    error: str, optional
        Name of the slurm error file capturing stderr output.
        Defaults to "slurm_%j.err", where %j is replaced with Job ID at execution time.
        See Slurm documentation for more details about '%' shortcuts behaviors.

    nodes: int, optional
        Number of reserved computation nodes. Must be at least 1.
        Defaults to 1.

    ntasks: int, optional
        Max number of tasks to do in parallel in total. If not given,
        `ntasks_per_node` must be given.

    ntasks_per_node: int, optional
        Max number of tasks to do in parallel in each node. If not given,
        `ntasks` must be given.

    n_gpus: int, optional
        Number of reserved GPUs in each node. Defaults to 0 for CPU only execution.

    cpus_per_task: int, optional
        Number of reserved CPU cores for each task. Defaults to 1 for sequential execution.
        
    time: str, optional
        Max timeout for the job before stopping. Format is "d-hh:mm:ss",
        where d stands for days, h for hours, m for minutes and s for seconds.
        If the job is less than a day, the "d-" part is optional.

    constraint: str, optional
        Apply a constraint to the job parameters, such as reserving only specific nodes.
        Must be a valid constraint defined for the slurm installation. See documentation
        for more details.

    partition: str, optional
        Name of a specific partition to target for the run.

    hint: str, optional
        Apply a specific property to the run, such as deactivationg hyperthreading. See
        documentation for more details.

    qos: str, optional
        Apply a defined Quality of Service. See documentation for more details.

    account: str, optional
        Define an account to debit computation hours from, if applicable.

    content_lines: list[str] | tuple[str], optional
        List or tuple of the command lines to write in the sbatch script for execution.
        Newline characters are added between each element automatically, and the
        resulting content is written as-is after all sbatch parameters.
    """
    shebang: str = "/bin/bash"
    jobname: str = "jobname"
    output: str = "slurm_%j.out"
    error: str = "slurm_%j.err"
    nodes: int = 1
    ntasks: int = 0
    ntasks_per_node: int = 0
    n_gpus: int = 0
    cpus_per_task: int = 1
    time: str | None = None
    constraint: str | None = None
    partition: str | None = None
    hint: str | None = None
    qos: str | None = None
    account: str | None = None
    content_lines: list[str] | tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        assert isinstance(self.shebang, PathLike)
        assert isinstance(self.jobname, str)
        assert isinstance(self.output, str)
        assert isinstance(self.error, str)
        assert self.nodes > 0
        assert bool(self.ntasks) ^ bool(self.ntasks_per_node)
        assert self.ntasks >= 0
        assert self.ntasks_per_node >= 0
        assert self.n_gpus >= 0
        assert self.cpus_per_task > 0
        assert isinstance(self.time, (type(None), str))
        assert isinstance(self.constraint, (type(None), str))
        assert isinstance(self.partition, (type(None), str))
        assert isinstance(self.hint, (type(None), str))
        assert isinstance(self.qos, (type(None), str))
        assert isinstance(self.account, (type(None), str))
        assert isinstance(self.content_lines, (type(None), list, tuple))

    def write_script(self, filename: PathLike) -> None:
        """Generate the script corresponding to attributes values."""
        sbatch_header = "#SBATCH "
        lines = f"#!{self.shebang}\n"
        lines += f"{sbatch_header} --job-name={self.jobname}\n"
        lines += f"{sbatch_header} --output={self.output}\n"
        lines += f"{sbatch_header} --error={self.error}\n"

        if self.time:
            lines += f"{sbatch_header} --time={self.time}\n"

        lines += f"{sbatch_header} --nodes={self.nodes}\n"

        if self.ntasks:
            lines += f"{sbatch_header} --ntasks={self.ntasks_per_node}\n"
        if self.ntasks_per_node:
            lines += f"{sbatch_header} --ntasks_per_node={self.ntasks_per_node}\n"
        if self.n_gpus:
            lines += f"{sbatch_header} --gres=gpu:{self.n_gpus}\n"

        lines += f"{sbatch_header} --cpus_per_task={self.cpus_per_task}\n"

        if self.constraint:
            lines += f"{sbatch_header} --constraint={self.constraint}\n"
        if self.partition:
            lines += f"{sbatch_header} --partition={self.partition}\n"
        if self.hint:
            lines += f"{sbatch_header} --hint={self.hint}\n"
        if self.qos:
            lines += f"{sbatch_header} --qos={self.qos}\n"
        if self.account:
            lines += f"{sbatch_header} --account={self.account}\n"

        lines += "\n"

        if self.content_lines:
            lines += "\n".join(self.content_lines)

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write(lines)
