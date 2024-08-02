from typing import List, Optional, Tuple
import copy
from data_types import SequenceDB, CoincidenceSequence, PatternInstance, Tiep, TiepProjector, Coincidence
from tiep_index import TiepIndex, MasterTiep
import constants
import utils


def project_initial_seq_db(
		initial_seq_db: SequenceDB,
		tiep: str,
		supporting_entities: List[str],
		index: TiepIndex
) -> SequenceDB:
	"""
	projects initial sequence database by a tiep
	:param initial_seq_db: (SequenceDB) initial sequence database
	:param tiep: (str) tiep for projection
	:param supporting_entities: (List[str]) list of supporting entities of the tiep
	:param index: (TiepIndex) main tiep index
	:return: (SequenceDB) projected sequence database
	"""

	projected_db: List[Tuple[CoincidenceSequence, PatternInstance]] = []
	master_tiep: MasterTiep = index.master_tieps[tiep]

	for coincidence_seq, pattern_instance in initial_seq_db.db:
		entity_id: str = coincidence_seq.entity
		if entity_id not in supporting_entities:
			continue

		entity_tiep_instances: List[Tiep] = master_tiep.tiep_occurrences[entity_id]
		for i in range(len(entity_tiep_instances)):
			tiep_instance: Tiep = entity_tiep_instances[i]
			projected_record: Optional[CoincidenceSequence] = __project_seq_by_tiep_instance(
				tiep_instance, tiep, coincidence_seq, pattern_instance
			)

			if projected_record is not None:
				extended_pattern_instance: PatternInstance = PatternInstance()
				extended_pattern_instance.extend_pattern_instance(tiep_instance)
				projected_db.append((projected_record, extended_pattern_instance))

	return SequenceDB(projected_db, None, len(master_tiep.supporting_entities))


def project_projected_seq_db(
		seq_db: SequenceDB,
		tiep: str,
		tiep_projector: TiepProjector,
		index: TiepIndex,
		maximal_gap: int
) -> SequenceDB:
	"""
	projects a projected sequence database by a tiep
	:param seq_db: (SequenceDB) initial sequence database
	:param tiep: (str) tiep for projection
	:param tiep_projector: (TiepProjector) tiep-projector of the tiep
	:param index: (TiepIndex) main tiep index
	:param maximal_gap: (int) maximal gap
	:return: (SequenceDB) projected sequence database
	"""

	projected_db: List[Tuple[CoincidenceSequence, PatternInstance]] = []
	projected_indices: List[int] = []

	base_tiep_form: str = tiep
	is_meet: bool = False
	is_co: bool = False

	# get base primitive form of tiep
	if base_tiep_form[0] == constants.MEET_REP:
		is_meet = True
		base_tiep_form = base_tiep_form[1:]
	elif base_tiep_form[0] == constants.CO_REP:
		is_co = True
		base_tiep_form = base_tiep_form[1:]

	master_tiep: MasterTiep = index.master_tieps[base_tiep_form]
	is_start_tiep: bool = base_tiep_form[-1] == constants.START_REP
	supporting_entities: List[str] = []

	for db_entry_index, first_index in tiep_projector.first_indices.items():
		coincidence_seq, pattern_instance = seq_db.db[db_entry_index]
		entity_id: str = coincidence_seq.entity
		entity_tiep_instances: List[Tiep] = master_tiep.tiep_occurrences[entity_id]

		for i in range(first_index, len(entity_tiep_instances)):
			if entity_tiep_instances[i].time > pattern_instance.first_expected_finish_time:
				continue
			if is_start_tiep and not utils.max_gap_holds(
					pattern_instance.minimal_finish_time, entity_tiep_instances[i], maximal_gap
			):
				break

			tiep_instance: Tiep = entity_tiep_instances[i]
			projected_record: Optional[CoincidenceSequence] = __project_seq_by_tiep_instance(
				tiep_instance, tiep, coincidence_seq, pattern_instance
			)

			if projected_record is not None:
				if entity_id not in supporting_entities:
					supporting_entities.append(entity_id)
				extended_pattern_instance: PatternInstance = PatternInstance()
				extended_pattern_instance.pre_extend_copy(pattern_instance)
				extended_pattern_instance.extend_pattern_instance(tiep_instance)
				projected_db.append((projected_record, extended_pattern_instance))
				projected_indices.append(db_entry_index)

			if is_co or is_meet or not is_start_tiep:
				break

	return SequenceDB(projected_db, projected_indices, len(supporting_entities))


def __project_seq_by_tiep_instance(
		tiep_instance: Tiep,
		tiep: str,
		coincidence_seq: CoincidenceSequence,
		pattern_instance: PatternInstance
) -> Optional[CoincidenceSequence]:
	"""
	projects a coincidence sequence by a tiep and returns the projected sequence
	:param tiep_instance: (Tiep) specific tiep instance
	:param tiep: (str) tiep representation
	:param coincidence_seq: (CoincidenceSequence) coincidence sequence to project
	:param pattern_instance: (PatternInstance) respective pattern instance
	:return: (Optional[CoincidenceSequence]) projected coincidence sequence, if succeeded
	"""

	first_projected_co, succeeded = __get_projected_coincidence_seq(
		tiep_instance, tiep, coincidence_seq, pattern_instance
	)
	if not succeeded:
		return None

	cs: CoincidenceSequence = CoincidenceSequence(
		coincidence_seq.entity,
		first_projected_co,
		first_projected_co if first_projected_co is not None and first_projected_co.index == tiep_instance.coincidence.index else None
	)
	return cs


def __get_projected_coincidence_seq(
		tiep_instance: Tiep,
		tiep: str,
		coincidence_seq: CoincidenceSequence,
		pattern_instance: PatternInstance
) -> Tuple[Optional[Coincidence], bool]:
	"""
	returns first coincidence from which the projected coincidence sequence begins after projection by a tiep,
		if projection is possible
	:param tiep_instance: (Tiep) specific tiep instance
	:param tiep: (str) tiep representation
	:param coincidence_seq: (CoincidenceSequence) coincidence sequence to project
	:param pattern_instance: (PatternInstance) respective pattern instance
	:return: (Tuple[Optional[Coincidence], bool]) tuple of
		<first coincidence of projected sequence, projection success status>
	"""

	if coincidence_seq.partial_co is not None and coincidence_seq.partial_co.index == tiep_instance.coincidence.index:
		current_coincidence: Coincidence = coincidence_seq.partial_co
	else:
		current_coincidence: Coincidence = tiep_instance.coincidence

	coincidence_tieps: List[Tiep] = current_coincidence.tieps
	for i in range(len(coincidence_tieps)):
		found: bool = tiep_instance == coincidence_tieps[i].orig_tiep if tiep[0] == constants.CO_REP else coincidence_tieps[i] == tiep_instance
		if found:
			current_tiep: Tiep = coincidence_tieps[i]
			if current_tiep.type == constants.FINISH_REP and not __is_tiep_valid_for_extension(
					current_tiep, pattern_instance
			):
				return None, False

			projected_seq_first_co: Coincidence = Coincidence(current_coincidence.index, is_co=True)
			if i < len(coincidence_tieps) - 1:
				for k in range(i + 1, len(coincidence_tieps)):
					if current_coincidence.is_co:
						current_tiep = coincidence_tieps[k]
					else:
						current_tiep = copy.copy(coincidence_tieps[k])
						current_tiep.orig_tiep = coincidence_tieps[k]
					projected_seq_first_co.tieps.append(current_tiep)

			projected_seq_first_co.next = current_coincidence.next

			if len(projected_seq_first_co.tieps) == 0:
				projected_seq_first_co = projected_seq_first_co.next

			return projected_seq_first_co, True

	return None, False


def __is_tiep_valid_for_extension(
		tiep: Tiep,
		pattern_instance: PatternInstance
) -> bool:
	"""
	return whether it is valid to extend the current pattern instance by the given tiep or not
	:param tiep: (Tiep) tiep for extension
	:param pattern_instance: (PatternInstance) current pattern instance
	:return: (bool) whether it is valid to extend the current pattern instance by the given tiep or not
	"""

	return tiep.sti in pattern_instance.pre_matched
