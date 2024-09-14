from typing import List, Tuple, Optional
from pathlib import Path
from tirpclo.data_types import STI, Tiep, Coincidence, CoincidenceSequence, PatternInstance, SequenceDB
from tirpclo.tiep_index import TiepIndex
from tirpclo import constants


class EndTime:
    """this class represents an End-Time, i.e., collection of co-occurring STI end-points having the same type

    Attributes:  # noqa
        stis: (List[STI]) list of stis for which one of the end-points has the given type and occurs in the given time
        entry_time: (int) time
        entry_type: (bool) type
    """
    def __init__(self, sti: STI, entry_time: int, entry_type: bool):
        self.stis: List[STI] = [sti]
        self.entry_time: int = entry_time
        self.entry_type: bool = entry_type

    def add_sti(self, new_sti: STI) -> None:
        """
        adds an STI to the end-time
        :param new_sti: (STI) STI to add
        :return: (None)
        """

        min_idx: int = 0
        max_idx: int = len(self.stis) - 1

        while min_idx <= max_idx:
            mid: int = (min_idx + max_idx) // 2
            diff: int = new_sti.symbol - self.stis[mid].symbol
            if diff < 0:
                max_idx = mid - 1
            else:
                min_idx = mid + 1

        self.stis.insert(min_idx, new_sti)

    def compare(self, other: 'EndTime') -> int:
        """
        compares current end-time to other end-time
        :param other: (EndTime) end-time to compare
        :return: (int) comparison result of current end-time and other end-time
        """

        time_diff: int = self.entry_time - other.entry_time
        if time_diff != 0:
            return time_diff

        type_diff: bool = self.entry_type == other.entry_type
        if not type_diff:
            if self.entry_type == constants.FINISH:
                return -1
            else:
                return 1

        return 0


def transform_input_file_to_seq_db(
        in_file_path: str,
        tiep_index: TiepIndex
) -> SequenceDB:
    """
    transforms an input STIs series file into a sequence database, while populating the tiep index
    :param in_file_path: (str) input file path
    :param tiep_index: (TiepIndex) main tiep index to populate
    :return: (SequenceDB) initial sequence database
    """

    if not Path(in_file_path).is_file():
        raise Exception(constants.IN_FILE_NOT_EXISTS_ERR)

    with open(in_file_path, 'r') as in_file:

        # move on until the significant start
        while line := in_file.readline():
            if line.startswith(constants.FILE_START):
                break
        if not (line.startswith(constants.FILE_START) and in_file.readline().startswith(constants.FILE_NUM)):
            raise Exception(constants.IN_FILE_FORMAT_ERR)

        # start reading the entities
        seq_db: List[Tuple[CoincidenceSequence, PatternInstance]] = []
        while line := in_file.readline().strip():

            entity_id: str = line.split(';')[0].split(',')[0]

            line = in_file.readline().strip()
            stis_line_components = line.split(';')
            end_time_list: List[EndTime] = []

            for i in range(len(stis_line_components) - 1):
                sti_components: List[str] = stis_line_components[i].split(',')
                sti: STI = STI(
                    start_time=int(sti_components[constants.STI_START_INDEX]),
                    finish_time=int(sti_components[constants.STI_FINISH_INDEX]),
                    symbol=int(sti_components[constants.STI_SYMBOL_INDEX])
                )
                __add_sti_to_end_times(sti, end_time_list)

            coincidence_seq: CoincidenceSequence = __convert_event_seq_to_coincidence_seq(
                entity_id, end_time_list, tiep_index
            )
            pattern_instance: PatternInstance = PatternInstance()
            seq_db.append((coincidence_seq, pattern_instance))

    return SequenceDB(seq_db, None, 0, None)


def __add_sti_to_end_times(
        sti: STI,
        end_time_list: List[EndTime]
) -> None:
    """
    adds the end-points of an STI to the end-time list
    :param sti: (STI) input STI to add
    :param end_time_list: (List[EndTime]) cumulative end-time list
    :return: (None)
    """

    start_end_time: EndTime = EndTime(sti, sti.start_time, constants.START)
    finish_end_time: EndTime = EndTime(sti, sti.finish_time, constants.FINISH)
    __add_point_to_end_times(start_end_time, end_time_list)
    __add_point_to_end_times(finish_end_time, end_time_list)


def __add_point_to_end_times(
        new_et: EndTime,
        end_time_list: List[EndTime]
) -> None:
    """
    adds an end-point to the end-time list in an ordered manner
    :param new_et: (EndTime) end-time of additional end-point
    :param end_time_list: (List[EndTime]) cumulative end-time list
    :return: (None)
    """

    place, end_time = __get_end_time_place(new_et, end_time_list)
    if place >= 0:
        end_time_list.insert(place, new_et)
    else:
        end_time.add_sti(new_et.stis[0])


def __get_end_time_place(
        end_time: EndTime,
        end_time_list: List[EndTime]
) -> Tuple[int, Optional[EndTime]]:
    """
    finds the correct place within the ordered end-time list for the additional end-time
    :param end_time: (EndTime) end-time of additional end-point
    :param end_time_list: (List[EndTime]) cumulative end-time list
    :return: (Tuple[int, Optional[EndTime]]) a pair of <position to insert (if not found), matching end-time (if found)>
    """

    min_index: int = 0
    max_index: int = len(end_time_list) - 1

    while min_index <= max_index:
        mid_index = (min_index + max_index) // 2
        diff: int = end_time.compare(end_time_list[mid_index])
        if diff == 0:
            return -1, end_time_list[mid_index]
        elif diff < 0:
            max_index = mid_index - 1
        else:
            min_index = mid_index + 1

    return min_index, None


def __convert_event_seq_to_coincidence_seq(
        entity: str,
        end_time_list: List[EndTime],
        tiep_index: TiepIndex
) -> CoincidenceSequence:
    """
    converts an event sequence represented as an end-time list into a coincidence sequence,
        while populating the tiep index
    :param entity: (str) entity ID
    :param end_time_list: (List[EndTime]) cumulative end-time list
    :param tiep_index: tiep_index: (TiepIndex) main tiep index to populate
    :return: (CoincidenceSequence) converted coincidence sequence
    """

    curr_coincidence: Optional[Coincidence] = None
    first_coincidence: Optional[Coincidence] = None
    last_end_time: Optional[EndTime] = None

    for index, end_time in enumerate(end_time_list):
        is_meet: bool = False
        entry_rep: str = constants.FINISH_REP

        if end_time.entry_type == constants.START:
            entry_rep = constants.START_REP
            if last_end_time is not None and last_end_time.entry_time == end_time.entry_time:
                is_meet = True

        coincidence: Coincidence = Coincidence(index, is_meet)
        __generate_and_index_coincidence_tieps(
            tiep_index, coincidence, end_time.stis, end_time.entry_time, entry_rep, entity
        )

        if index == 0:
            first_coincidence = coincidence
        else:
            curr_coincidence.next = coincidence

        index += 1
        curr_coincidence = coincidence
        last_end_time = end_time

    return CoincidenceSequence(entity, first_coincidence)


def __generate_and_index_coincidence_tieps(
        tiep_index: TiepIndex,
        coincidence: Coincidence,
        stis: List[STI],
        time: int,
        tieps_type: str,
        entity: str
) -> None:
    """
    generates a coincidence's tieps and populates the main tiep index with them
    :param tiep_index: (TiepIndex) main tiep index to populate
    :param coincidence: (Coincidence) coincidence to populate
    :param stis: (List[STI]) list of stis whose end-points belong to the coincidence
    :param time: (int) time-point of coincidence
    :param tieps_type: (str) whether start (+) or finish (-)
    :param entity: (str) entity ID
    :return: (None)
    """

    for sti in stis:
        tiep: Tiep = Tiep(time, sti, coincidence, tieps_type)
        coincidence.tieps.append(tiep)
        sti.entity_sti_index = tiep_index.add_tiep_occurrence(tiep.primitive_rep, entity, tiep)
