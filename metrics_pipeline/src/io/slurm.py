"""Write slurm job submission scripts."""

import typing as tp
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
        Defaults to 'slurm_%j.out', or 'slurm_%A_%a.out' if 'array' is given.
        See Slurm documentation for more details about '%' shortcuts behaviors.

    error: str, optional
        Name of the slurm error file capturing stderr output.
        Defaults to 'slurm_%j.err', or 'slurm_%A_%a.err' if 'array' is given.
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
        Must be a valid constraint defined for the slurm installation. See cluster documentation
        for more details.

    partition: str, optional
        Name of a specific partition to target for the run.

    hint: str, optional
        Apply a specific property to the run, such as deactivationg hyperthreading.
        See cluster documentation for more details.

    qos: str, optional
        Apply a defined Quality of Service. See cluster documentation for more details.

    account: str, optional
        Define an account to debit computation hours from, if applicable.

    array: Iterable[int] | str, optional
        Define a queue of jobs using the same script with indices.
        Use the slurm environment variable `$SLURM_ARRAY_TASK_ID` in `content_lines` to use
        the array index of each queued job and distinguish between each job similar actions.
        - If a string is given, it will be written as-is in the script, and `array_max_exec`
        is ignored.
        - If an iterable of indices is given, consecutive indices will be concatenated
        automatically before writing for conciseness.

    array_max_exec: int, optional
        If `array` is given as an iterable of indices, define max number of jobs in the array
        to execute at the same time. Ignored if `array` is not given or given as a string.

    content_lines: list[str] | tuple[str], optional
        List or tuple of the command lines to write in the sbatch script for execution.
        Newline characters are added between each element automatically, and the
        resulting content is written as-is after all sbatch parameters.
    """
    _sbatch_header: str = "#SBATCH"
    shebang: str = "/bin/bash"
    jobname: str = "jobname"
    output: str | None = None
    error: str | None = None
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
    array: tp.Iterable[int] | str | None = None
    array_max_exec: int | None = None
    content_lines: list[str] | tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        assert isinstance(self.shebang, PathLike)
        assert isinstance(self.jobname, str)
        assert isinstance(self.output, (str, type(None)))
        if self.output is None:
            self.output = f"slurm_{'%j' if self.array is None else '%A_%a'}.out"
        assert isinstance(self.error, (str, type(None)))
        if self.error is None:
            self.error = f"slurm_{'%j' if self.array is None else '%A_%a'}.err"
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
        assert isinstance(self.array, (type(None), tp.Iterable, str))
        assert isinstance(self.array_max_exec, (type(None), int))
        assert isinstance(self.content_lines, (type(None), list, tuple))

    def _process_array_indices(self) -> str | None:
        """Process job array indices into the most compact valid string for the script."""
        if isinstance(self.array, (type(None), str)):
            return self.array

        # Avoid any duplication and sort
        indices = sorted(set(self.array))
        consecutives: list[list[int]] = [[]]

        # Group consecutive indices
        for idx, array_idx in enumerate(indices[1:], start=1):
            if array_idx == indices[idx -1] + 1:
                consecutives[-1].append(array_idx)
            else:
                consecutives.append([array_idx])

        array_str_list = []

        # Process each group appropriately
        for group in consecutives:
            if len(group) > 1:
                array_str_list.append(f"{group[0]}-{group[-1]}")
            else:
                array_str_list.append(f"{group[0]}")
        array_str = ",".join(array_str_list)

        # Add max parallel execution if defined
        if self.array_max_exec:
            array_str += f"%{self.array_max_exec}"

        return array_str

    def write_script(self, filename: PathLike) -> None:
        """Generate the script corresponding to attributes values."""
        HEAD = self._sbatch_header
        lines = [f"#!{self.shebang}"]
        lines.append(f"{HEAD} --job-name={self.jobname}")
        lines.append(f"{HEAD} --output={self.output}")
        lines.append(f"{HEAD} --error={self.error}")

        if self.time:
            lines.append(f"{HEAD} --time={self.time}")

        lines.append(f"{HEAD} --nodes={self.nodes}")

        if self.ntasks:
            lines.append(f"{HEAD} --ntasks={self.ntasks_per_node}")
        if self.ntasks_per_node:
            lines.append(f"{HEAD} --ntasks_per_node={self.ntasks_per_node}")
        if self.n_gpus:
            lines.append(f"{HEAD} --gres=gpu:{self.n_gpus}")

        lines.append(f"{HEAD} --cpus_per_task={self.cpus_per_task}")

        if self.constraint:
            lines.append(f"{HEAD} --constraint={self.constraint}")
        if self.partition:
            lines.append(f"{HEAD} --partition={self.partition}")
        if self.hint:
            lines.append(f"{HEAD} --hint={self.hint}")
        if self.qos:
            lines.append(f"{HEAD} --qos={self.qos}")
        if self.account:
            lines.append(f"{HEAD} --account={self.account}")
        if self.array:
            array_str = self._process_array_indices()
            lines.append(f"{HEAD} --array={array_str}")

        lines.append("")

        if self.content_lines:
            lines.extend(self.content_lines)

        with open(filename, "wt", encoding="utf-8") as fp:
            fp.write("\n".join(lines))
