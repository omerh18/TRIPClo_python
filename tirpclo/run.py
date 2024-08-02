from typing.io import TextIO
import math
import time
from tirpclo.data_types import SequenceDB
from tirpclo.tiep_index import TiepIndex
from tirpclo import main_algorithm
from tirpclo import stis2seq
from tirpclo import utils


def run_tirpclo(
        num_entities: int,
        min_support_percentage: float,
        maximal_gap: int,
        in_file_path: str,
        out_file_path: str = None
) -> None:
    """
    runs TIRPClo to discover the entire set of frequent TIRPs
    :param num_entities: (int) number of dataset's entities
    :param min_support_percentage: (float) minimum vertical support threshold
    :param maximal_gap: (int) maximal gap
    :param in_file_path: (str) path to input file
    :param out_file_path: (str) path to TIRPs output file
    :return: (None)
    """

    if out_file_path is None:
        out_file_path: str = f'{in_file_path[: -4]}-support-{min_support_percentage}-gap-{maximal_gap}.txt'

    out_file: TextIO = utils.out_file_set_up(out_file_path)
    min_support: int = math.ceil(num_entities * min_support_percentage)
    print(f'***START***\nrunning TIRPClo on dataset: {in_file_path}, support - {min_support_percentage}, gap - {maximal_gap}.')

    start_time: float = time.time()

    index: TiepIndex = TiepIndex()
    initial_seq_db: SequenceDB = stis2seq.transform_input_file_to_seq_db(in_file_path, index)
    main_algorithm.discover_tirps(index, initial_seq_db, min_support, maximal_gap, out_file)

    end_time: float = time.time()

    runtime_sec: float = end_time - start_time
    print(f'***FINISHED***\nruntime - {runtime_sec} (sec).')
    out_file.close()
    utils.generate_sorted_output_file(out_file_path)
    utils.generate_stats_output_file(out_file_path, runtime_sec)


if __name__ == '__main__':
    parsed_args = utils.parse_arguments()
    run_tirpclo(**parsed_args)
    # run_tirpclo(
    #     num_entities=65, min_support_percentage=0.5, maximal_gap=30, in_file_path='datasets/asl/asl.csv'
    # )
