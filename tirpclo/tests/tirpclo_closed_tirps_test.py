import os
from tirpclo.run import run_tirpclo
from tirpclo import utils


def test_tirpclo_closed_tirps_on_asl(
):
    """runs TIRPClo, mining closed frequent TIRPs, on the ASL dataset and verifies a correct output is produced"""

    is_closed_tirp_mining = True
    num_entities = 65
    min_support_percentage = 0.1
    maximal_gap = 30
    in_file_path = '../../datasets/asl/asl.csv'
    out_file_path = f'asl-closed-support-{min_support_percentage}-gap-{maximal_gap}.txt'

    run_config: utils.RunConfig = utils.RunConfig(
        is_closed_tirp_mining=is_closed_tirp_mining,
        num_entities=num_entities,
        min_support_percentage=min_support_percentage,
        maximal_gap=maximal_gap,
        in_file_path=in_file_path,
        out_file_path=out_file_path
    )

    run_tirpclo(run_config)

    verified_output_file = f'verified-asl-closed-support-{min_support_percentage}-gap-{maximal_gap}.txt'
    output_lines = open(utils.get_sorted_output_file_name(out_file_path), 'r').readlines()
    verified_lines = open(utils.get_sorted_output_file_name(verified_output_file), 'r').readlines()

    os.remove(out_file_path)
    os.remove(utils.get_sorted_output_file_name(out_file_path))
    os.remove(utils.get_stats_output_file_name(out_file_path))

    assert '\n'.join(output_lines) == '\n'.join(verified_lines)
