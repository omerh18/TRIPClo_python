from typing import List, Optional, Tuple, Dict
import copy
from tirpclo.data_types import SequenceDB, CoincidenceSequence, PatternInstance, Tiep, \
	TiepProjector, Coincidence, BackwardExtensionTiep
from tirpclo.tiep_index import TiepIndex, MasterTiep
from tirpclo import closure_checking
from tirpclo import constants
from tirpclo import utils


def project_initial_seq_db(
		initial_seq_db: SequenceDB,
		tiep_primitive_rep: str,
		supporting_entities: List[str],
		index: TiepIndex,
		maximal_gap: int,
		is_closed_tirp_mining: bool
) -> Tuple[SequenceDB, Optional[bool], Optional[Dict[str, List[BackwardExtensionTiep]]]]:
	"""
	projects initial sequence database by a tiep
	:param initial_seq_db: (SequenceDB) initial sequence database
	:param tiep_primitive_rep: (str) tiep for projection
	:param supporting_entities: (List[str]) list of supporting entities of the tiep
	:param index: (TiepIndex) main tiep index
	:param maximal_gap: (int) maximal gap
	:param is_closed_tirp_mining: (bool) whether mining only closed TIRPs or not
	:return: (Tuple[SequenceDB, bool, Dict[str, List[BackwardExtensionTiep]]]) projected sequence database, as well as
		whether the projected pattern has the potential of being closed and its backward-extension tieps
	"""

	projected_db: List[Tuple[CoincidenceSequence, PatternInstance]] = []
	master_tiep: MasterTiep = index.master_tieps[tiep_primitive_rep]
	cumulative_be_tieps: Optional[Dict[str, BackwardExtensionTiep]] = None
	entry_index: int = 0

	for coincidence_seq, pattern_instance in initial_seq_db.db:
		entity_id: str = coincidence_seq.entity
		if entity_id not in supporting_entities:
			continue

		entity_tiep_instances: List[Tiep] = master_tiep.tiep_occurrences[entity_id]
		entity_be_tieps: Dict[str, BackwardExtensionTiep] = {}
		for i, tiep_instance in enumerate(entity_tiep_instances):
			projected_record: Optional[CoincidenceSequence] = __project_seq_by_tiep_instance(
				tiep_instance, tiep_primitive_rep, coincidence_seq, pattern_instance
			)

			if projected_record is not None:
				extended_pattern_instance: PatternInstance = PatternInstance()
				if is_closed_tirp_mining:
					extended_pattern_instance.next_coincidences.append(coincidence_seq.first_co)
				extended_pattern_instance.extend_pattern_instance(
					tiep_instance, projected_record.first_co, is_closed_tirp_mining
				)
				projected_db.append((projected_record, extended_pattern_instance))

				if is_closed_tirp_mining:
					closure_checking.collect_be_tieps_wrt_tiep_instance(
						tiep_instance, coincidence_seq.first_co if i == 0 else entity_tiep_instances[i - 1].coincidence,
						entry_index, entity_be_tieps, cumulative_be_tieps, maximal_gap
					)

				entry_index += 1

		cumulative_be_tieps = entity_be_tieps

	be_tieps_lists: Optional[Dict[str, List[BackwardExtensionTiep]]] = None
	may_be_closed: Optional[bool] = None
	pre_matched: Optional[List[str]] = None
	if is_closed_tirp_mining:
		may_be_closed, be_tieps_lists = closure_checking.finalize_initial_be_tieps(cumulative_be_tieps)
		pre_matched = [tiep_primitive_rep.replace(constants.START_REP, constants.FINISH_REP)]

	return SequenceDB(
		projected_db, None, len(master_tiep.supporting_entities), pre_matched
	), may_be_closed, be_tieps_lists


def project_projected_seq_db(
		seq_db: SequenceDB,
		tiep: str,
		tiep_projector: TiepProjector,
		index: TiepIndex,
		maximal_gap: int,
		is_closed_tirp_mining: bool
) -> SequenceDB:
	"""
	projects a projected sequence database by a tiep
	:param seq_db: (SequenceDB) initial sequence database
	:param tiep: (str) tiep for projection
	:param tiep_projector: (TiepProjector) tiep-projector of the tiep
	:param index: (TiepIndex) main tiep index
	:param maximal_gap: (int) maximal gap
	:param is_closed_tirp_mining: (bool) whether mining only closed TIRPs or not
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
				extended_pattern_instance.pre_extend_copy(pattern_instance, is_closed_tirp_mining)
				extended_pattern_instance.extend_pattern_instance(
					tiep_instance, projected_record.first_co, is_closed_tirp_mining
				)
				projected_db.append((projected_record, extended_pattern_instance))
				projected_indices.append(db_entry_index)

			if is_co or is_meet or not is_start_tiep:
				break

	pre_matched: Optional[List[str]] = None
	if is_closed_tirp_mining:
		pre_matched = seq_db.pre_matched.copy()
		if base_tiep_form[-1] == constants.START_REP:
			pre_matched.append(base_tiep_form.replace(constants.START_REP, constants.FINISH_REP))
		else:
			pre_matched.remove(base_tiep_form)

	return SequenceDB(projected_db, projected_indices, len(supporting_entities), pre_matched)


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
