from typing import List, Dict, Optional
from tirpclo.data_types import TiepProjector, SequenceDB, Coincidence, Tiep, PatternInstance
from tirpclo.tiep_index import TiepIndex, MasterTiep
from tirpclo import constants
from tirpclo import utils


def get_tiep_projectors(
		seq_db: SequenceDB,
		pattern_last_tiep: str,
		previous_tiep_projectors: Dict[str, TiepProjector],
		index: TiepIndex,
		min_support: int,
		maximal_gap: int
) -> Dict[str, TiepProjector]:
	"""
	generates and returns new tiep-projectors based on the previous ones, if exist
	:param seq_db: (SequenceDB) projected sequence database
	:param pattern_last_tiep: (str) last tiep of current pattern represented by the sequence database
	:param previous_tiep_projectors: (Dict[str, TiepProjector]) mapping of all previous tiep-projectors,
		based on which new tiep-projectors are created
	:param index: (TiepIndex) main tiep-index
	:param min_support: (int) minimum vertical support threshold
	:param maximal_gap: (int) maximal gap
	:return: (Dict[str, TiepProjector]) new tiep-projectors
	"""

	if previous_tiep_projectors is None:
		# initial tiep-projectors
		return __get_initial_tiep_projectors(seq_db, pattern_last_tiep, maximal_gap)

	# get base primitive form of last tiep
	if pattern_last_tiep[0] == constants.CO_REP or pattern_last_tiep[0] == constants.MEET_REP:
		pattern_last_tiep = pattern_last_tiep[1:]

	tiep_projectors: Dict[str, TiepProjector] = {}
	allowed_non_supporting_records: int = len(seq_db.db) - min_support

	# populate tiep-projectors based on recent ones
	__populate_tiep_projectors_based_on_recent(
		seq_db, pattern_last_tiep, previous_tiep_projectors, index,
		min_support, maximal_gap, tiep_projectors, allowed_non_supporting_records
	)

	# if the last tiep was a start tiep, add complement to tiep-projectors
	if pattern_last_tiep[-1] == constants.START_REP:
		__add_complement_finish_tiep_to_tiep_projectors(
			seq_db, pattern_last_tiep, index, tiep_projectors, allowed_non_supporting_records
		)

	# if the last tiep was a finish tiep, potentially add complement to tiep-projectors
	if pattern_last_tiep[-1] == constants.FINISH_REP:
		__add_complement_start_tiep_to_tiep_projectors(
			seq_db, pattern_last_tiep, index, maximal_gap, tiep_projectors, allowed_non_supporting_records
		)

	# add meet & co-occurrence tieps to tiep-projectors
	entry_index: int = 0
	for coincidence_seq, pattern_instance in seq_db.db:
		entity_id: str = coincidence_seq.entity
		current_coincidence: Optional[Coincidence] = coincidence_seq.first_co
		if current_coincidence is None:
			entry_index += 1
			continue
		__add_relevant_meet_co_tieps_to_tiep_projectors(
			current_coincidence, entity_id, entry_index, tiep_projectors, pattern_instance
		)
		entry_index += 1

	return tiep_projectors


def __populate_tiep_projectors_based_on_recent(
		seq_db: SequenceDB,
		pattern_last_tiep: str,
		previous_tiep_projectors: Dict[str, TiepProjector],
		index: TiepIndex,
		min_support: int,
		maximal_gap: int,
		tiep_projectors: Dict[str, TiepProjector],
		allowed_non_supporting_records: int
) -> None:
	"""
	populates newly created tiep-projectors based on recent ones
	:param seq_db: (SequenceDB) projected sequence database
	:param pattern_last_tiep: (str) last tiep of current pattern represented by the sequence database
	:param previous_tiep_projectors: (Dict[str, TiepProjector]) mapping of all previous tiep-projectors,
		based on which new tiep-projectors are created
	:param index: (TiepIndex) main tiep-index
	:param min_support: (int) minimum vertical support threshold
	:param maximal_gap: (int) maximal gap
	:param tiep_projectors: (Dict[str, TiepProjector]) incrementally populated new tiep-projectors
	:param (int) allowed_non_supporting_records: maximal allowed number of non-supporting records for
		a tiep to be surely concluded as infrequent
	:return: (None)
	"""

	for tiep, previous_tiep_projector in previous_tiep_projectors.items():

		# make sure the tiep has been recently frequent, and does not represent a special case handled separately
		# (i.e., not a co-occurrence / meet tiep as well as not equal to the pattern last tiep)
		if len(previous_tiep_projector.supporting_entities) < min_support:
			continue
		if tiep[0] == constants.CO_REP or tiep[0] == constants.MEET_REP:
			continue
		if pattern_last_tiep == tiep:
			continue

		master_tiep: MasterTiep = index.master_tieps[tiep]
		is_finish_tiep: bool = tiep[-1] == constants.FINISH_REP
		entry_index: int = 0
		non_supporting_records: int = 0

		for coincidence_seq, pattern_instance in seq_db.db:

			if non_supporting_records > allowed_non_supporting_records:
				break

			# reach first non special coincidence
			entity_id: str = coincidence_seq.entity
			current_coincidence: Optional[Coincidence] = coincidence_seq.first_co
			if current_coincidence is None:
				non_supporting_records += 1
				entry_index += 1
				continue
			if current_coincidence.is_co:
				current_coincidence = current_coincidence.next
				if current_coincidence is None:
					non_supporting_records += 1
					entry_index += 1
					continue
			if current_coincidence.is_meet:
				current_coincidence = current_coincidence.next
				if current_coincidence is None:
					non_supporting_records += 1
					entry_index += 1
					continue

			start_co_index: int = current_coincidence.index
			previous_entry_index: int = seq_db.entries_prev_indices[entry_index]
			# if entry is not in previous tiep projector
			if previous_entry_index not in previous_tiep_projector.first_indices:
				non_supporting_records += 1
				entry_index += 1
				continue

			tiep_instances: List[Tiep] = master_tiep.tiep_occurrences[entity_id]
			# if it is a finish tiep, check if the projected coincidence sequence includes the specific instance
			# of the start tiep which complements the current tiep within the pattern instance
			if is_finish_tiep:
				tiep_index: int = pattern_instance.symbol_db_indices[tiep_instances[0].symbol]
				if tiep_instances[tiep_index].coincidence.index >= start_co_index:
					__add_tiep_instance_to_tiep_projectors(
						tiep, entity_id, entry_index, tiep_projectors, tiep_index, validate_first=False
					)
				else:
					non_supporting_records += 1
				entry_index += 1
				continue

			# for a start tiep, reach the first instance within the projected coincidence sequence, if exist
			prev_start_index: int = previous_tiep_projector.first_indices[previous_entry_index]
			found: bool = False
			for i in range(prev_start_index, len(tiep_instances)):
				if not utils.max_gap_holds(pattern_instance.minimal_finish_time, tiep_instances[i], maximal_gap):
					break
				if tiep_instances[i].coincidence.index >= start_co_index:
					__add_tiep_instance_to_tiep_projectors(
						tiep, entity_id, entry_index, tiep_projectors, i, validate_first=False
					)
					found = True
					break
			if not found:
				non_supporting_records += 1

			entry_index += 1


def __add_complement_finish_tiep_to_tiep_projectors(
		seq_db: SequenceDB,
		pattern_last_tiep: str,
		index: TiepIndex,
		tiep_projectors: Dict[str, TiepProjector],
		allowed_non_supporting_records: int
) -> None:
	"""
	populates newly created tiep-projectors with complementing finish tiep of the pattern last tiep
	:param seq_db: (SequenceDB) projected sequence database
	:param pattern_last_tiep: (str) last tiep of current pattern represented by the sequence database
	:param index: (TiepIndex) main tiep-index
	:param tiep_projectors: (Dict[str, TiepProjector]) incrementally populated new tiep-projectors
	:param (int) allowed_non_supporting_records: maximal allowed number of non-supporting records for
		a tiep to be surely concluded as infrequent
	:return: (None)
	"""

	finish_tiep_rep: str = pattern_last_tiep.replace(constants.START_REP, constants.FINISH_REP)
	master_tiep: MasterTiep = index.master_tieps[finish_tiep_rep]
	entry_index: int = 0
	non_supporting_records: int = 0

	for coincidence_seq, pattern_instance in seq_db.db:
		if non_supporting_records > allowed_non_supporting_records:
			break

		# reach first non special coincidence
		entity_id: str = coincidence_seq.entity
		current_coincidence: Optional[Coincidence] = coincidence_seq.first_co
		if current_coincidence is None:
			non_supporting_records += 1
			entry_index += 1
			continue
		if current_coincidence.is_co:
			current_coincidence = current_coincidence.next
			if current_coincidence is None:
				non_supporting_records += 1
				entry_index += 1
				continue

		# add the specific finish tiep complementing the recently added start tiep for each record
		start_co_index: int = current_coincidence.index
		tiep_instances: List[Tiep] = master_tiep.tiep_occurrences[entity_id]
		tiep_index: int = pattern_instance.symbol_db_indices[tiep_instances[0].symbol]
		if tiep_instances[tiep_index].coincidence.index >= start_co_index:
			__add_tiep_instance_to_tiep_projectors(
				finish_tiep_rep, entity_id, entry_index, tiep_projectors, tiep_index, validate_first=False
			)
		else:
			non_supporting_records += 1

		entry_index += 1


def __add_complement_start_tiep_to_tiep_projectors(
		seq_db: SequenceDB,
		pattern_last_tiep: str,
		index: TiepIndex,
		maximal_gap: int,
		tiep_projectors: Dict[str, TiepProjector],
		allowed_non_supporting_records: int
) -> None:
	"""
	populates newly created tiep-projectors with complementing start tiep of the pattern last tiep
	:param seq_db: (SequenceDB) projected sequence database
	:param pattern_last_tiep: (str) last tiep of current pattern represented by the sequence database
	:param index: (TiepIndex) main tiep-index
	:param maximal_gap: (int) maximal gap
	:param tiep_projectors: (Dict[str, TiepProjector]) incrementally populated new tiep-projectors
	:param (int) allowed_non_supporting_records: maximal allowed number of non-supporting records for
		a tiep to be surely concluded as infrequent
	:return: (None)
	"""

	start_tiep: str = pattern_last_tiep.replace(constants.FINISH_REP, constants.START_REP)
	master_tiep: MasterTiep = index.master_tieps[start_tiep]
	non_supporting_records: int = 0
	entry_index: int = 0

	for coincidence_seq, pattern_instance in seq_db.db:
		if non_supporting_records > allowed_non_supporting_records:
			break

		# reach first non special coincidence
		entity_id: str = coincidence_seq.entity
		current_coincidence: Optional[Coincidence] = coincidence_seq.first_co
		if current_coincidence is None:
			non_supporting_records += 1
			entry_index += 1
			continue
		if current_coincidence.is_co:
			current_coincidence = current_coincidence.next
			if current_coincidence is None:
				non_supporting_records += 1
				entry_index += 1
				continue
		if current_coincidence.is_meet:
			current_coincidence = current_coincidence.next
			if current_coincidence is None:
				non_supporting_records += 1
				entry_index += 1
				continue

		start_co_index: int = current_coincidence.index
		tiep_instances: List[Tiep] = master_tiep.tiep_occurrences[entity_id]
		tiep_index: int = pattern_instance.tieps[-1].entity_tiep_index + 1
		found: bool = False

		# reach the first instance within the projected coincidence sequence, if exist
		for i in range(tiep_index, len(tiep_instances)):
			if not utils.max_gap_holds(pattern_instance.minimal_finish_time, tiep_instances[i], maximal_gap):
				break
			if tiep_instances[i].coincidence.index >= start_co_index:
				__add_tiep_instance_to_tiep_projectors(
					start_tiep, entity_id, entry_index, tiep_projectors, i, validate_first=False
				)
				found = True
				break
		if not found:
			non_supporting_records += 1

		entry_index += 1


def __get_initial_tiep_projectors(
		seq_db: SequenceDB,
		pattern_last_tiep: str,
		maximal_gap: int
) -> Dict[str, TiepProjector]:
	"""
	generates and returns new tiep-projectors when no previous ones exist
	:param seq_db: (SequenceDB) projected sequence database
	:param pattern_last_tiep: (str) last tiep of current pattern represented by the sequence database
	:param maximal_gap: (int) maximal gap
	:return: (Dict[str, TiepProjector]) new tiep-projectors
	"""

	tiep_projectors: Dict[str, TiepProjector] = {}
	entry_index: int = 0

	for coincidence_seq, pattern_instance in seq_db.db:
		current_coincidence: Optional[Coincidence] = coincidence_seq.first_co
		entity_id: str = coincidence_seq.entity
		found_complement: bool = False
		beyond_gap: bool = False

		while current_coincidence is not None:

			if beyond_gap and found_complement:
				break

			tieps: List[Tiep] = current_coincidence.tieps
			is_finish_tieps_coincidence: bool = tieps[0].type == constants.FINISH_REP
			# skip finish tieps after the complement of last tiep has been found & start tieps
			# after violating maximal gap
			if (found_complement and is_finish_tieps_coincidence) or (beyond_gap and not is_finish_tieps_coincidence):
				current_coincidence = current_coincidence.next
				continue

			for current_tiep in tieps:
				# for a finish tiep, add only if complements last tiep
				if is_finish_tieps_coincidence:
					if pattern_last_tiep == current_tiep.primitive_rep.replace(constants.FINISH_REP, constants.START_REP):
						__add_tiep_instance_to_tiep_projectors(
							current_tiep.primitive_rep, entity_id, entry_index,
							tiep_projectors, current_tiep.entity_tiep_index, validate_first=False
						)
						found_complement = True
						break
					continue

				# for a start-tiep which differs from the last tiep, add only if does not violate maximal gap
				if pattern_last_tiep == current_tiep.primitive_rep:
					continue
				if not utils.max_gap_holds(pattern_instance.minimal_finish_time, current_tiep, maximal_gap):
					beyond_gap = True
					break

				current_tiep_rep: str = (constants.CO_REP if current_coincidence.is_co else '') + current_tiep.primitive_rep
				current_tiep_orig_tiep: Tiep = current_tiep if current_tiep.orig_tiep is None else current_tiep.orig_tiep
				__add_tiep_instance_to_tiep_projectors(
					current_tiep_rep, entity_id, entry_index, tiep_projectors,
					current_tiep_orig_tiep.entity_tiep_index, validate_first=True
				)

			current_coincidence = current_coincidence.next

		entry_index += 1

	return tiep_projectors


def __add_tiep_instance_to_tiep_projectors(
		tiep_rep: str,
		entity_id: str,
		entry_index: int,
		tiep_projectors: Dict[str, TiepProjector],
		first_index: int,
		validate_first: bool = False
) -> None:
	"""
	adds a new tiep instance to the tiep-projector of a tiep
	:param tiep_rep: (str) tiep representation
	:param entity_id: (str) entity ID
	:param entry_index: (int) index of entry in sequences database
	:param tiep_projectors: (Dict[str, TiepProjector]) tiep-projectors
	:param first_index: (int) tiep first index within entry's coincidence sequence
	:param validate_first: (bool) whether to check or not for an already recorded first index
	:return:
	"""

	if tiep_rep not in tiep_projectors:
		tiep_projectors[tiep_rep] = TiepProjector()

	if entity_id not in tiep_projectors[tiep_rep].supporting_entities:
		tiep_projectors[tiep_rep].supporting_entities.append(entity_id)

	if not validate_first or entry_index not in tiep_projectors[tiep_rep].first_indices:
		tiep_projectors[tiep_rep].first_indices[entry_index] = first_index


def __add_relevant_meet_co_tieps_to_tiep_projectors(
		current_coincidence: Coincidence,
		entity_id: str,
		entry_index: int,
		tieps_projectors: Dict[str, TiepProjector],
		pattern_instance: PatternInstance
) -> None:
	"""
	add meet & co-occurrence tieps to tiep-projectors
	:param current_coincidence: (Coincidence) current coincidence sequence
	:param entity_id: (str) entity ID
	:param entry_index: (int) index of entry in sequences database
	:param tieps_projectors: (Dict[str, TiepProjector]) tiep-projectors
	:param pattern_instance: (PatternInstance) current pattern instance
	:return:
	"""

	tieps: List[Tiep] = current_coincidence.tieps

	if current_coincidence.is_co:
		is_finish_tieps_coincidence: bool = current_coincidence.tieps[0].type == constants.FINISH_REP
		for i in range(len(tieps)):
			tiep_rep: str = constants.CO_REP + '' + tieps[i].primitive_rep
			if is_finish_tieps_coincidence and tieps[i].sti not in pattern_instance.pre_matched:
				continue
			__add_tiep_instance_to_tiep_projectors(
				tiep_rep, entity_id, entry_index, tieps_projectors, tieps[i].orig_tiep.entity_tiep_index
			)

		if current_coincidence.next is not None and current_coincidence.next.is_meet:
			current_coincidence = current_coincidence.next
			tieps = current_coincidence.tieps
			for i in range(len(tieps)):
				tiep_rep: str = constants.MEET_REP + '' + tieps[i].primitive_rep
				__add_tiep_instance_to_tiep_projectors(
					tiep_rep, entity_id, entry_index, tieps_projectors, tieps[i].entity_tiep_index
				)

	elif current_coincidence.is_meet:
		for i in range(len(tieps)):
			tiep_rep: str = constants.MEET_REP + '' + tieps[i].primitive_rep
			__add_tiep_instance_to_tiep_projectors(
				tiep_rep, entity_id, entry_index, tieps_projectors, tieps[i].entity_tiep_index
			)
