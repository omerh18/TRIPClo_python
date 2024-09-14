from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from typing.io import TextIO
from pathlib import Path
import distutils.util
import argparse
from tirpclo.data_types import Tiep, STI
from tirpclo import constants


@dataclass
class RunConfig:
    """this class represents a run configuration

    Attributes:  # noqa
        is_closed: (bool) whether mining only closed TIRPs or not
        num_entities: (int) number of entities
        min_support_percentage: (float) minimum vertical support percentage
        maximal_gap: (int) maximal gap
        in_file_path: (str) path to input file
        out_file_path: (str) path to output file
    """
    is_closed_tirp_mining: bool
    num_entities: int
    min_support_percentage: float
    maximal_gap: int
    in_file_path: str
    out_file_path: Optional[str] = None

    def __post_init__(self):
        if self.out_file_path is None:
            self.out_file_path = f'{self.in_file_path[: -4]}-support-{self.min_support_percentage}-gap-{self.maximal_gap}.txt'


def parse_arguments(
) -> Dict[str, Any]:
    """
    parses input arguments
    :return: (Dict[str, Any]) parsed arguments
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--is_closed_tirp_mining', type=lambda x: bool(distutils.util.strtobool(x)), required=True
    )
    parser.add_argument(
        '-n', '--num_entities', type=int, required=True
    )
    parser.add_argument(
        '-s', '--min_support_percentage', type=float, required=True
    )
    parser.add_argument(
        '-g', '--maximal_gap', type=int, required=True
    )
    parser.add_argument(
        '-f', '--in_file_path', type=str, required=True
    )

    parsed_args = parser.parse_args()

    return parsed_args.__dict__


def out_file_set_up(
        out_file_path: str
) -> TextIO:
    """
    sets-up an output file, in append mode
    :param out_file_path: (str) output file path
    :return: (TextIO) output file
    """

    if Path(out_file_path).is_file():
        raise Exception(constants.OUT_FILE_EXISTS_ERR)

    out_file: TextIO = open(out_file_path, 'a+')
    return out_file


def max_gap_holds(
        pattern_minimal_finish_time: float,
        candidate_tiep: Tiep,
        maximal_gap: int
) -> bool:
    """
    checks whether the maximal gap constraint holds for a given pattern w.r.t an additional tiep
    :param pattern_minimal_finish_time: (float) a given pattern's minimal finish time
    :param candidate_tiep: (Tiep) additional tiep
    :param maximal_gap: (int) maximal gap
    :return: (bool) whether the maximal gap constraint holds for a given pattern w.r.t an additional tiep
    """

    candidate_sti: STI = candidate_tiep.sti
    return maximal_gap > candidate_sti.start_time - pattern_minimal_finish_time


def get_sorted_output_file_name(
        out_file_path: str
) -> str:
    """
    returns sorted output file path
    :param out_file_path: (str) output file path
    :return: (str) sorted output file path
    """

    return f'{out_file_path[:-4]}_sorted.txt'


def generate_sorted_output_file(
        out_file_path: str
) -> None:
    """
    generates sorted output file
    :param out_file_path: (str) output file path
    :return: (None)
    """

    with open(out_file_path, 'r') as out_file:
        output_lines: List[str] = out_file.readlines()

    output_lines.sort()
    sorted_out_file_path = get_sorted_output_file_name(out_file_path)

    open(sorted_out_file_path, 'w').writelines(output_lines)


def get_stats_output_file_name(
        out_file_path: str
) -> str:
    """
    returns stats output file path
    :param out_file_path: (str) output file path
    :return: (str) stats output file path
    """

    return f'{out_file_path[:-4]}_stats.txt'


def generate_stats_output_file(
        out_file_path: str,
        runtime_sec: float
) -> None:
    """
    generates stats output file
    :param out_file_path: (str) output file path
    :param runtime_sec: (float) runtime in seconds
    :return: (None)
    """

    open(get_stats_output_file_name(out_file_path), 'w').writelines([f'{runtime_sec}'])
