from typing.io import TextIO
from typing import List
from tirpclo.data_types import STI, SequenceDB
from tirpclo import constants


def write_tirp(
        seq_db: SequenceDB,
        out_file: TextIO
) -> None:
    """
    writes a TIRP represented by the sequence database to the output file
    :param seq_db: (SequenceDB) sequence database representing a TIRP
    :param out_file: (TextIOWrapper) output file
    :return: (None)
    """

    support: int = seq_db.support
    length: int = len(seq_db.db[0][1].tieps) // 2
    stis: List[STI] = [tiep.sti for tiep in seq_db.db[0][1].tieps if tiep.type == constants.START_REP]
    stis.sort()

    tirp_text_representation = f"{length} "
    tirp_text_representation += "-".join([f"{sti.symbol}" for sti in stis]) + " "

    if length == 1:
        tirp_text_representation += "-."
    else:
        for i in range(length):
            for j in range(i + 1, length):
                tirp_text_representation += __get_relation(stis[i], stis[j]) + "."

    tirp_text_representation += f" {support} "
    tirp_text_representation += f"{support if length == 1 else round(len(seq_db.db) / support, 2)} "
    tirp_text_representation += f"{seq_db.db[0][0].entity} {__get_stis_as_str(stis)} "

    for i in range(1, len(seq_db.db)):
        stis = [tiep.sti for tiep in seq_db.db[i][1].tieps if tiep.type == constants.START_REP]
        stis.sort()
        tirp_text_representation += f"{seq_db.db[i][0].entity} {__get_stis_as_str(stis)} "

    out_file.write(f'{tirp_text_representation[: -1]}\n')


def __get_stis_as_str(
        stis: List[STI]
) -> str:
    """
    returns string representation of stis
    :param stis: (List[STI]) list of stis
    :return: (str) string representation of stis
    """

    return "".join([repr(sti) for sti in stis])


def __get_relation(
        sti1: STI,
        sti2: STI
) -> str:
    """
    returns temporal relation between two STIs
    :param sti1: (STI) first STI
    :param sti2: (STI) second STI
    :return: (str) temporal relation
    """

    if sti1.finish_time < sti2.start_time:
        return constants.ALLEN_BEFORE
    if sti1.finish_time == sti2.start_time:
        return constants.ALLEN_MEET
    if sti1.start_time == sti2.start_time and sti1.finish_time == sti2.finish_time:
        return constants.ALLEN_EQUAL
    if sti1.start_time < sti2.start_time and sti1.finish_time > sti2.finish_time:
        return constants.ALLEN_CONTAIN
    if sti1.start_time == sti2.start_time and sti1.finish_time < sti2.finish_time:
        return constants.ALLEN_STARTS
    if sti1.start_time < sti2.start_time and sti1.finish_time == sti2.finish_time:
        return constants.ALLEN_FINISHBY
    return constants.ALLEN_OVERLAP
