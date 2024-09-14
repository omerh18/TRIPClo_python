from typing import List, Dict, Optional, Tuple
from tirpclo.data_types import SequenceDB, TiepProjector, BackwardExtensionTiep, Coincidence, Tiep
from tirpclo import constants
from tirpclo import utils


def may_tirp_be_closed(
		pattern_seq_db: SequenceDB,
		tiep_projectors: Dict[str, TiepProjector],
		be_tieps_lists: Dict[str, List[BackwardExtensionTiep]]
) -> bool:
	"""
	returns whether the TIRP represented by the sequence database may be a closed TIRP or not,
		based on backward-extension and forward-extension tieps
	:param pattern_seq_db: (SequenceDB) projected sequence database
	:param tiep_projectors: (Dict[str, TiepProjector]) tiep-projectors serving as forward-extension tieps
	:param be_tieps_lists: (Dict[str, List[BackwardExtensionTiep]]) backward-extension tieps
	:return: (bool) whether the TIRP represented by the sequence database may be a closed TIRP or not
	"""

	for tiep, tiep_projector in tiep_projectors.items():

		if pattern_seq_db.support == len(tiep_projector.supporting_entities):

			if tiep[-1] == constants.START_REP:
				return False

			tiep_primitive_rep: str = tiep[1:] if tiep[0] == constants.CO_REP else tiep
			complement_start_tiep_rep: str = tiep_primitive_rep.replace(constants.FINISH_REP, constants.START_REP)
			if complement_start_tiep_rep in be_tieps_lists:
				if __do_be_fe_match_in_all_entities(
						be_tieps_lists[complement_start_tiep_rep], tiep_projector, pattern_seq_db
				):
					return False

	return True


def __do_be_fe_match_in_all_entities(
		start_tiep_be_tieps: List[BackwardExtensionTiep],
		finish_fe_tiep_projector: TiepProjector,
		pattern_seq_db: SequenceDB
) -> bool:
	"""
	returns whether all entities include an STI which matches the backward-extension start tiep
		and the forward-extension finish tiep having the same symbol
	:param start_tiep_be_tieps: (List[BackwardExtensionTiep]) backward-extension tieps of a start tiep
	:param finish_fe_tiep_projector: (TiepProjector) tiep-projector representing a forward-extension tiep
		of a finish tiep having the same symbol of the start tiep
	:param pattern_seq_db: (SequenceDB) projected sequence database
	:return: (bool) whether all entities include an STI which matches the backward and forward-extension tieps
	"""

	for be_tiep in start_tiep_be_tieps:
		matching_entities: List[str] = []

		for db_entry_index, finish_tiep_first_index in finish_fe_tiep_projector.first_indices.items():
			if db_entry_index in be_tiep.stis_per_entry:
				entity_id: str = pattern_seq_db.db[db_entry_index][0].entity
				if entity_id not in matching_entities:
					for be_tiep_sti in be_tiep.stis_per_entry[db_entry_index]:
						if be_tiep_sti.entity_sti_index >= finish_tiep_first_index:
							matching_entities.append(entity_id)
							break

		if pattern_seq_db.support == len(matching_entities):
			return True

	return False


def __do_be_be_match_in_all_entities(
		start_tiep_be_tieps: List[BackwardExtensionTiep],
		finish_be_tiep: BackwardExtensionTiep,
		pattern_seq_db: SequenceDB
) -> bool:
	"""
	returns whether all entities include an STI which matches the backward-extension start tiep
		and the backward-extension finish tiep having the same symbol
	:param start_tiep_be_tieps: (List[BackwardExtensionTiep]) backward-extension tieps of a start tiep
	:param finish_be_tiep: (BackwardExtensionTiep) backward-extension tiep of a finish tiep
		having the same symbol of the start tiep
	:param pattern_seq_db: (SequenceDB) projected sequence database
	:return: (bool) whether all entities include an STI which matches the backward-extension tieps
	"""

	for be_tiep in start_tiep_be_tieps:
		matching_entities: List[str] = []

		for db_entry_index, finish_tiep_stis in finish_be_tiep.stis_per_entry.items():
			if db_entry_index in be_tiep.stis_per_entry:
				entity_id: str = pattern_seq_db.db[db_entry_index][0].entity
				if entity_id not in matching_entities:
					for finish_tiep_sti in finish_tiep_stis:
						if finish_tiep_sti in be_tiep.stis_per_entry[db_entry_index]:
							matching_entities.append(entity_id)
							break

		if pattern_seq_db.support == len(matching_entities):
			return True

	return False


def back_scan(
		pattern_seq_db: SequenceDB,
		maximal_gap: int
) -> Tuple[bool, Dict[str, List[BackwardExtensionTiep]]]:
	"""
	checks whether a projected pattern has the potential of being closed or not, based on its backward extension tieps
	:param pattern_seq_db: (SequenceDB) projected sequence database
	:param maximal_gap: (int) maximal gap
	:return: (bool) whether the projected pattern has the potential of being closed or not,
		based on its backward extension tieps
	"""

	be_tieps_lists: Dict[str, List[BackwardExtensionTiep]] = {}
	for i in range(len(pattern_seq_db.db[0][1].tieps)):
		cumulative_ith_before_be_tieps: Optional[Dict[str, BackwardExtensionTiep]] = None
		entity_ith_before_be_tieps: Optional[Dict[str, BackwardExtensionTiep]] = None
		entry_index: int = 0
		for coincidence_seq, pattern_instance in pattern_seq_db.db:

			if entry_index == 0 or coincidence_seq.entity != pattern_seq_db.db[entry_index - 1][0].entity:
				cumulative_ith_before_be_tieps = entity_ith_before_be_tieps
				if cumulative_ith_before_be_tieps is not None and len(cumulative_ith_before_be_tieps) == 0:
					break
				entity_ith_before_be_tieps = {}

			tiep_instance: Tiep = pattern_instance.tieps[i]
			current_coincidence: Coincidence = pattern_instance.next_coincidences[i]

			coincidence_prefix: str = constants.CO_REP if current_coincidence.is_co else \
				(constants.MEET_REP if current_coincidence.is_meet else '*')
			while current_coincidence.index != tiep_instance.coincidence.index:

				if current_coincidence.index == tiep_instance.coincidence.index - 1 and tiep_instance.coincidence.is_meet:
					for current_tiep in current_coincidence.tieps:
						tiep_full_rep: str = coincidence_prefix + constants.MEET_REP + current_tiep.primitive_rep
						__add_current_tiep_to_entity_be_tieps(
							tiep_instance, current_tiep, tiep_full_rep, entry_index,
							entity_ith_before_be_tieps, cumulative_ith_before_be_tieps, maximal_gap, check_gap=False
						)

				else:
					for current_tiep in current_coincidence.tieps:
						tiep_full_rep: str = coincidence_prefix + '*' + current_tiep.primitive_rep
						__add_current_tiep_to_entity_be_tieps(
							tiep_instance, current_tiep, tiep_full_rep, entry_index,
							entity_ith_before_be_tieps, cumulative_ith_before_be_tieps, maximal_gap, check_gap=True
						)

				coincidence_prefix = constants.MEET_REP if current_coincidence.is_co and current_coincidence.next.is_meet else '*'
				current_coincidence = current_coincidence.next

			for current_tiep in current_coincidence.tieps:
				if current_tiep == tiep_instance or tiep_instance == current_tiep.orig_tiep:
					break
				tiep_full_rep: str = coincidence_prefix + constants.CO_REP + current_tiep.primitive_rep
				__add_current_tiep_to_entity_be_tieps(
					tiep_instance, current_tiep, tiep_full_rep, entry_index,
					entity_ith_before_be_tieps, cumulative_ith_before_be_tieps, maximal_gap, check_gap=False
				)
			entry_index += 1

		cumulative_ith_before_be_tieps = entity_ith_before_be_tieps
		if len(cumulative_ith_before_be_tieps) == 0:
			continue

		if not finalize_ith_before_be_tieps(cumulative_ith_before_be_tieps, be_tieps_lists, pattern_seq_db):
			return False, be_tieps_lists

	return True, be_tieps_lists


def finalize_ith_before_be_tieps(
		cumulative_ith_before_be_tieps: Dict[str, BackwardExtensionTiep],
		be_tieps_lists: Dict[str, List[BackwardExtensionTiep]],
		pattern_seq_db: SequenceDB
) -> bool:
	"""
	checks whether a given projected pattern has the potential of being closed or not, and also returns
		its final ith-before backward-extension tieps, for some i
	:param cumulative_ith_before_be_tieps: (Dict[str, BackwardExtensionTiep]) cumulative ith-before backward-extension
		tieps over all entities
	:param be_tieps_lists: (Dict[str, List[BackwardExtensionTiep]]) backward-extension tieps
	:param pattern_seq_db: (SequenceDB) projected sequence database
	:return: (bool) whether a given projected pattern
		has the potential of being closed or not based on the ith-before backward-extension tieps
	"""

	for be_tiep_full_rep, be_tiep in cumulative_ith_before_be_tieps.items():
		be_tiep_primitive_rep: str = be_tiep_full_rep[2:]

		if be_tiep_primitive_rep[-1] == constants.START_REP:
			if be_tiep_primitive_rep not in be_tieps_lists:
				be_tieps_lists[be_tiep_primitive_rep] = []
			be_tieps_lists[be_tiep_primitive_rep].append(be_tiep)

	for be_tiep_full_rep, be_tiep in cumulative_ith_before_be_tieps.items():
		be_tiep_primitive_rep: str = be_tiep_full_rep[2:]

		if be_tiep_primitive_rep[-1] == constants.FINISH_REP:
			be_tiep_primitive_rep = be_tiep_primitive_rep.replace(constants.FINISH_REP, constants.START_REP)
			if be_tiep_primitive_rep in be_tieps_lists and \
				__do_be_be_match_in_all_entities(be_tieps_lists[be_tiep_primitive_rep], be_tiep, pattern_seq_db):
				return False

	return True


def collect_be_tieps_wrt_tiep_instance(
		tiep_instance: Tiep,
		current_coincidence: Coincidence,
		entry_index: int,
		entity_be_tieps: Dict[str, BackwardExtensionTiep],
		cumulative_be_tieps: Optional[Dict[str, BackwardExtensionTiep]],
		maximal_gap: int
) -> None:
	"""
	collects backward-extension tieps w.r.t the projected tiep instance
	:param tiep_instance: (Tiep) tiep instance w.r.t which sequence database has been projected
	:param current_coincidence: (Coincidence) coincidence from which to look for backward-extension tieps
	:param entry_index: (int) index of entry in projected sequence database
	:param entity_be_tieps: (Dict[str, BackwardExtensionTiep]) entity backward-extension tieps
	:param cumulative_be_tieps: (Optional[Dict[str, BackwardExtensionTiep]]) cumulative backward-extension
		tieps over all entities
	:param maximal_gap: (int) maximal gap
	:return: (None)
	"""

	for tiep_rep, be_tiep in entity_be_tieps.items():
		if tiep_rep[0] != constants.CO_REP and tiep_rep[0] != constants.MEET_REP:
			if entry_index - 1 not in be_tiep.stis_per_entry:
				continue
			for sti in be_tiep.stis_per_entry[entry_index - 1]:
				if utils.max_gap_holds(sti.finish_time, tiep_instance, maximal_gap):
					entity_be_tieps[tiep_rep].add_sti_in_entry(entry_index, sti)

	while current_coincidence.index != tiep_instance.coincidence.index:
		if current_coincidence.index == tiep_instance.coincidence.index - 1 and tiep_instance.coincidence.is_meet:
			for current_tiep in current_coincidence.tieps:
				tiep_full_rep: str = constants.MEET_REP + current_tiep.primitive_rep
				__add_current_tiep_to_entity_be_tieps(
					tiep_instance, current_tiep, tiep_full_rep, entry_index,
					entity_be_tieps, cumulative_be_tieps, maximal_gap, check_gap=False
				)
		else:
			for current_tiep in current_coincidence.tieps:
				tiep_full_rep: str = '*' + current_tiep.primitive_rep
				__add_current_tiep_to_entity_be_tieps(
					tiep_instance, current_tiep, tiep_full_rep, entry_index,
					entity_be_tieps, cumulative_be_tieps, maximal_gap, check_gap=True
				)
		current_coincidence = current_coincidence.next

	for current_tiep in current_coincidence.tieps:
		if current_tiep == tiep_instance:
			break
		tiep_full_rep: str = constants.CO_REP + current_tiep.primitive_rep
		__add_current_tiep_to_entity_be_tieps(
			tiep_instance, current_tiep, tiep_full_rep, entry_index,
			entity_be_tieps, cumulative_be_tieps, maximal_gap, check_gap=False
		)


def __add_current_tiep_to_entity_be_tieps(
		projected_tiep_instance: Tiep,
		current_tiep: Tiep,
		tiep_full_rep: str,
		entry_index: int,
		entity_be_tieps: Dict[str, BackwardExtensionTiep],
		cumulative_be_tieps: Optional[Dict[str, BackwardExtensionTiep]],
		maximal_gap: int,
		check_gap: bool
) -> None:
	"""
	adds a current tiep as a backward-extension tiep w.r.t the projected tiep instance, if necessary
	:param projected_tiep_instance: (Tiep) tiep instance w.r.t which sequence database has been projected
	:param current_tiep: (Tiep) tiep to add to a backward-extension tiep
	:param tiep_full_rep: (str) tiep string representation
	:param entry_index: (int) index of entry in projected sequence database
	:param entity_be_tieps: (Dict[str, BackwardExtensionTiep]) entity backward-extension tieps
	:param cumulative_be_tieps: (Optional[Dict[str, BackwardExtensionTiep]]) cumulative backward-extension
		tieps over all entities
	:param maximal_gap: (int) maximal gap
	:param check_gap: (bool) whether maximal gap has to be tested or not
	:return: (None)
	"""

	if (cumulative_be_tieps is None or tiep_full_rep in cumulative_be_tieps) and \
		(not check_gap or utils.max_gap_holds(current_tiep.sti.finish_time, projected_tiep_instance, maximal_gap)):
		if tiep_full_rep not in entity_be_tieps:
			if cumulative_be_tieps is None:
				entity_be_tieps[tiep_full_rep] = BackwardExtensionTiep()
			else:
				entity_be_tieps[tiep_full_rep] = cumulative_be_tieps[tiep_full_rep]
		entity_be_tieps[tiep_full_rep].add_sti_in_entry(entry_index, current_tiep.sti)


def finalize_initial_be_tieps(
		cumulative_be_tieps: Dict[str, BackwardExtensionTiep]
) -> Tuple[bool, Dict[str, List[BackwardExtensionTiep]]]:
	"""
	checks whether a given initial (one-tiep) pattern has the potential of being closed or not, and also returns
		its final backward-extension tieps
	:param cumulative_be_tieps: (Dict[str, BackwardExtensionTiep]) cumulative backward-extension
		tieps over all entities
	:return: (Tuple[bool, Dict[str, List[BackwardExtensionTiep]]]) whether a given initial (one-tiep) pattern
		has the potential of being closed or not, and its final backward-extension tieps
	"""

	be_tieps_lists: Dict[str, List[BackwardExtensionTiep]] = {}

	for be_tiep_full_rep, be_tiep in cumulative_be_tieps.items():
		be_tiep_primitive_rep: str = be_tiep_full_rep[1:]

		if be_tiep_primitive_rep[-1] == constants.START_REP:
			if be_tiep_primitive_rep not in be_tieps_lists:
				be_tieps_lists[be_tiep_primitive_rep] = []
			be_tieps_lists[be_tiep_primitive_rep].append(be_tiep)

		else:
			return False, be_tieps_lists

	return True, be_tieps_lists
