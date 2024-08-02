from typing import List, Optional, Dict, Tuple
from dataclasses import dataclass, field
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from tiep_index import TiepIndex
import constants


@dataclass(order=True)
class STI:
    """this class represents an STI, i.e., Symbolic Time Interval

    Attributes:  # noqa
        start_time: (int) start-time
        finish_time: (int) finish-time
        symbol: (int) symbol
    """
    start_time: int
    finish_time: int
    symbol: int

    def __repr__(self):
        return f"[{self.start_time}-{self.finish_time}]"


@dataclass
class Tiep:
    """this class represents a Tiep, i.e., Time Interval End-Point

    Attributes:  # noqa
        symbol: (int) symbol of tiep
        time: (int) time of occurrence
        sti: (STI) STI of the tiep
        coincidence: (Coincidence) coincidence to which the tiep belongs
        type: (str) type of tiep
        primitive_rep: (str) tiep primitive representation, e.g., A+ or B-
        orig_tiep: (Optional[Tiep]) original tiep object from which current tiep is derived
            (relevant for meet / co-occurrence tieps only)
        entity_tiep_index: (int) index of tiep within ordered list of tieps having the same primitive_rep
            in the respective entity
    """
    symbol: int = field(init=False)
    time: int
    sti: STI
    coincidence: 'Coincidence'
    type: str
    primitive_rep: str = field(init=False)
    orig_tiep: Optional['Tiep'] = None
    entity_tiep_index: int = -1

    def __post_init__(self):
        self.symbol = self.sti.symbol
        self.primitive_rep = f'{self.symbol}{self.type}'


@dataclass
class Coincidence:
    """this class represents a Coincidence, i.e., list of coinciding tieps within a sequence

    Attributes:  # noqa
        index: (int) index of coincidence within sequence
        is_meet: (bool) whether it is a 'meet' coincidence
            (i.e., includes start tieps co-occurring with other finish tieps)
        is_co: (bool) whether it is a 'co-occurrence' coincidence (relevant for partially projected coincidences only)
        tieps: (List[tiep]) list of coinciding tieps
        next: (Optional[Coincidence]) next coincidence in sequence
    """
    index: int
    is_meet: bool = False
    is_co: bool = False
    tieps: List[Tiep] = field(default_factory=lambda: [])
    next: Optional['Coincidence'] = None


@dataclass
class CoincidenceSequence:
    """this class represents a Coincidence Sequence, i.e., ordered list of coincidence objects representing an entity

    Attributes:  # noqa
        entity: (str) entity represented by the sequence
        first_co: (Coincidence) first coincidence of sequence.
        partial_co: (Optional[Coincidence]) partially projected coincidence with which sequence currently begins,
            if exists
    """
    entity: str
    first_co: Coincidence
    partial_co: Optional[Coincidence] = None


@dataclass
class PatternInstance:
    """this class represents a Pattern Instance, i.e., a specific instance of a specific pattern (TIRP)

    Attributes:  # noqa
        tieps: (List[Tiep]) ordered list of tieps
        symbol_db_indices: (Dict[int, int]) index of latest entity occurrence of each symbol of the pattern instance
        minimal_finish_time: (float) minimal finish time of an STI within the pattern instance
        pre_matched: (List[STI]) list of STI for which only the start tieps are included within the pattern instances
        first_expected_finish_time: (float) earliest finish time expected to match the current STIs of the
            pattern instance
    """
    tieps: List[Tiep] = field(default_factory=lambda: [])
    symbol_db_indices: Dict[int, int] = field(default_factory=lambda: {})
    minimal_finish_time: float = float('inf')
    pre_matched: List[STI] = field(default_factory=lambda: [])
    first_expected_finish_time: float = float('inf')

    def pre_extend_copy(self, current_pattern_instance: 'PatternInstance') -> None:
        """
        copies (deep) contents of input pattern instance prior to extension
        :param current_pattern_instance: (PatternInstance) input pattern instance to duplicate
        :return: (None)
        """

        for t in current_pattern_instance.tieps:
            self.tieps.append(t)

        self.symbol_db_indices = dict(current_pattern_instance.symbol_db_indices)
        self.minimal_finish_time = current_pattern_instance.minimal_finish_time

        for sti in current_pattern_instance.pre_matched:
            self.pre_matched.append(sti)

        self.first_expected_finish_time = current_pattern_instance.first_expected_finish_time

    def extend_pattern_instance(self, new_tiep: Tiep) -> None:
        """
        extends the current pattern instance with an additional tiep
        :param new_tiep: (Tiep) additional tiep
        :return: (None)
        """

        self.tieps.append(new_tiep)

        if new_tiep.sti in self.pre_matched:
            self.pre_matched.remove(new_tiep.sti)
            if len(self.pre_matched) == 0:
                self.first_expected_finish_time = float('inf')
            else:
                self.first_expected_finish_time = min([sti.finish_time for sti in self.pre_matched])

        else:
            self.symbol_db_indices[new_tiep.symbol] = new_tiep.entity_tiep_index
            self.pre_matched.append(new_tiep.sti)
            self.first_expected_finish_time = min(self.first_expected_finish_time, new_tiep.sti.finish_time)

        if new_tiep.type == constants.START_REP:
            self.minimal_finish_time = min(self.minimal_finish_time, new_tiep.sti.finish_time)


@dataclass
class SequenceDB:
    """this class represents a Sequence Database

    Attributes:  # noqa
        db: (List[Tuple[CoincidenceSequence, PatternInstance]]) list of coincidence sequences & respective discovered
            pattern instances
        entries_prev_indices: (Optional[List[int]]) indices of db entries within previous sequence database,
            from which the current sequences database has been projected
        support: (int) vertical support of pattern represented by this sequence database
    """
    db: List[Tuple[CoincidenceSequence, PatternInstance]]
    entries_prev_indices: Optional[List[int]]
    support: int

    def filter_infrequent_tieps_from_initial_seq_db(self, index: 'TiepIndex') -> None:
        """
        filters infrequent tieps from this sequence database
        :param index: (TiepIndex) main tiep index
        :return: (None)
        """
        for coincidence_seq, _ in self.db:
            current_coincidence: Coincidence = coincidence_seq.first_co
            previous_coincidence: Optional[Coincidence] = None
            removed_coincidences: int = 0
            removed_recent: bool = False

            while current_coincidence is not None:
                tieps: List[Tiep] = current_coincidence.tieps
                removed_co_tieps: int = 0
                running_index: int = 0

                while running_index < len(tieps) + removed_co_tieps:
                    if tieps[running_index - removed_co_tieps].primitive_rep not in index.master_tieps:
                        tieps.remove(tieps[running_index - removed_co_tieps])
                        removed_co_tieps += 1
                    running_index += 1

                if len(tieps) == 0:
                    removed_coincidences += 1
                    removed_recent = True
                    if previous_coincidence is not None:
                        previous_coincidence.next = current_coincidence.next

                else:
                    current_coincidence.index -= removed_coincidences
                    if removed_recent:
                        current_coincidence.is_meet = False
                    if previous_coincidence is None:
                        coincidence_seq.first_co = current_coincidence
                    removed_recent = False
                    previous_coincidence = current_coincidence

                current_coincidence = current_coincidence.next


@dataclass
class TiepProjector:
    """this class represents a Tiep Projector, i.e., a data structure used to efficiently project a tiep from
        a given sequence database

    Attributes:  # noqa
        supporting_entities: (List[str]) list of supporting entities of the tiep-projector's tiep
        first_indices: (Dict[int, int]) first index of the tiep-projector's tiep within each record
            of a sequence database
    """
    supporting_entities: List[str] = field(default_factory=lambda: [])
    first_indices: Dict[int, int] = field(default_factory=lambda: {})
