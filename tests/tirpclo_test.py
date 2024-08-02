import os
from main import run_tirpclo
import utils


def test_tirpclo_on_asl(
):
    """runs TIRPClo on the ASL dataset and verifies a correct output is produced"""

    num_entities: int = 65
    min_support_percentage: float = 0.5
    maximal_gap: int = 30
    in_file_path: str = '../datasets/asl/asl.csv'
    out_file_path: str = f'asl-support-{min_support_percentage}-gap-{maximal_gap}.txt'

    run_tirpclo(
        num_entities=num_entities,
        min_support_percentage=min_support_percentage,
        maximal_gap=maximal_gap,
        in_file_path=in_file_path,
        out_file_path=out_file_path
    )

    verified_output_file = f'verified-asl-support-{min_support_percentage}-gap-{maximal_gap}.txt'
    output_lines = open(utils.get_sorted_output_file_name(out_file_path), 'r').readlines()
    verified_lines = open(utils.get_sorted_output_file_name(verified_output_file), 'r').readlines()

    os.remove(out_file_path)
    os.remove(utils.get_sorted_output_file_name(out_file_path))
    os.remove(utils.get_stats_output_file_name(out_file_path))

    assert '\n'.join(output_lines) == '\n'.join(verified_lines)
