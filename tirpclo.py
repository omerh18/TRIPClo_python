from typing import List, Dict, Optional, Tuple
from data_types import SequenceDB, TiepProjector, CoincidenceSequence, PatternInstance
from typing.io import TextIO
from tiep_index import TiepIndex
import candidate_generation
import tirp_writing
import projection
import constants


def discover_tirps(
		index: TiepIndex,
		initial_seq_db: SequenceDB,
		min_support: int,
		maximal_gap: int,
		out_file: TextIO
) -> None:
	"""
	discovers all frequent TIRPs
	:param index: (TiepIndex) main tiep index
	:param initial_seq_db: (SequenceDB) initial sequence database
	:param min_support: (int) minimum vertical support threshold
	:param maximal_gap: (int) maximal gap
	:param out_file: (TextIOWrapper) output file
	:return: (None)
	"""

	infrequent_tieps: List[str] = []

	for tiep, master_tiep in index.master_tieps.items():
		if len(master_tiep.supporting_entities) < min_support:
			infrequent_tieps.append(tiep)

	for tiep in infrequent_tieps:
		index.master_tieps.pop(tiep)

	initial_seq_db.filter_infrequent_tieps_from_initial_seq_db(index)
	for tiep, master_tiep in index.master_tieps.items():
		if tiep[-1] == constants.START_REP:
			projDB: SequenceDB = projection.project_initial_seq_db(
				initial_seq_db, tiep, master_tiep.supporting_entities, index
			)
			__extend_tirp(index, projDB, tiep, None, min_support, maximal_gap, out_file)


def __extend_tirp(
		index: TiepIndex,
		pattern_seq_db: SequenceDB,
		pattern_last_tiep: str,
		previous_tiep_projectors: Optional[Dict[str, TiepProjector]],
		min_support: int,
		maximal_gap: int,
		out_file
) -> None:
	"""
	recursively extends a current pattern
	:param index: (TiepIndex) main tiep index
	:param pattern_seq_db: (SequenceDB) projected sequence database
	:param pattern_last_tiep: (str) last tiep of current pattern represented by the sequence database
	:param previous_tiep_projectors: (Optional[Dict[str, TiepProjector]]) mapping of all previous tiep-projectors,
		based on which new tiep-projectors are created
	:param min_support: (int) minimum vertical support threshold
	:param maximal_gap: (int) maximal gap
	:param out_file: (TextIOWrapper) output file
	:return: (None)
	"""

	tiep_projectors: Dict[str, TiepProjector] = candidate_generation.get_tiep_projectors(
		pattern_seq_db, pattern_last_tiep, previous_tiep_projectors, index, min_support, maximal_gap
	)

	if __all_in_pairs(pattern_seq_db.db):
		tirp_writing.write_tirp(pattern_seq_db, out_file)

	for tiep, tiep_projector in tiep_projectors.items():
		if len(tiep_projector.supporting_entities) < min_support:
			continue

		projected_seq_db: SequenceDB = projection.project_projected_seq_db(
			pattern_seq_db, tiep, tiep_projector, index, maximal_gap
		)
		if projected_seq_db.support >= min_support:
			__extend_tirp(index, projected_seq_db, tiep, tiep_projectors, min_support, maximal_gap, out_file)


def __all_in_pairs(
		pattern_db: List[Tuple[CoincidenceSequence, PatternInstance]]
) -> bool:
	"""
	returns whether current pattern does not include start-tieps without their complementing finish-tieps,
		and thus represents a TIRP
	:param pattern_db: (List[Tuple[CoincidenceSequence, PatternInstance]]) records of a current pattern's
		sequence database
	:return: (bool) whether current pattern does not include start-tieps without their complementing finish-tieps,
		and thus represents a TIRP
	"""

	return len(pattern_db[0][1].pre_matched) == 0
