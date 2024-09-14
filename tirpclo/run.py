from typing.io import TextIO
import math
import time
from tirpclo.data_types import SequenceDB
from tirpclo.tiep_index import TiepIndex
from tirpclo import main_algorithm
from tirpclo import stis2seq
from tirpclo import utils


def run_tirpclo(run_config: utils.RunConfig) -> None:
    """
    runs TIRPClo to discover the entire set of frequent TIRPs
    :param run_config: (utils.RunConfig) run configuration
    :return: (None)
    """

    out_file: TextIO = utils.out_file_set_up(run_config.out_file_path)
    min_support: int = math.ceil(run_config.num_entities * run_config.min_support_percentage)
    print(
        f'***START***\nrunning TIRPClo on dataset: {run_config.in_file_path}, '
        f'support - {run_config.min_support_percentage}, gap - {run_config.maximal_gap}.'
    )

    start_time: float = time.time()

    index: TiepIndex = TiepIndex()
    initial_seq_db: SequenceDB = stis2seq.transform_input_file_to_seq_db(run_config.in_file_path, index)
    main_algorithm.discover_tirps(
        index, initial_seq_db, min_support, run_config.maximal_gap, out_file, run_config.is_closed_tirp_mining
    )

    end_time: float = time.time()

    runtime_sec: float = end_time - start_time
    print(f'***FINISHED***\nruntime - {runtime_sec} (sec).')
    out_file.close()
    utils.generate_sorted_output_file(run_config.out_file_path)
    utils.generate_stats_output_file(run_config.out_file_path, runtime_sec)


if __name__ == '__main__':
    parsed_args = utils.parse_arguments()
    input_run_config = utils.RunConfig(**parsed_args)
    run_tirpclo(input_run_config)
    # run_tirpclo(
    #     num_entities=65, min_support_percentage=0.5, maximal_gap=30, in_file_path='datasets/asl/asl.csv'
    # )
