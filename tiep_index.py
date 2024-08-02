from typing import List, Dict
from data_types import Tiep


class MasterTiep:
    """this class represents a Master Tiep, i.e., a data structure

    Attributes:  # noqa
        tiep_occurrences: (Dict[str, List[Tiep]]) all tiep instances by entity ID
        supporting_entities: (List[str]) list of supporting entities
    """
    def __init__(self):
        self.tiep_occurrences: Dict[str, List[Tiep]] = {}
        self.supporting_entities: List[str] = []

    def add_occurrence(self, entity: str, tiep: Tiep) -> None:
        """
        indexes a tiep instance within a specific entity
        :param entity: (str) entity ID
        :param tiep: (Tiep) additional tiep instance
        :return: (None)
        """

        if entity not in self.tiep_occurrences:
            self.tiep_occurrences[entity] = []
            self.supporting_entities.append(entity)

        tiep.entity_tiep_index = len(self.tiep_occurrences[entity])
        self.tiep_occurrences[entity].append(tiep)


class TiepIndex:
    def __init__(self):
        self.master_tieps: Dict[str, MasterTiep] = {}

    def add_tiep_occurrence(self, tiep_rep: str, entity: str, tiep: Tiep) -> None:
        """
        adds a tiep instance to the index
        :param tiep_rep: (ste) tiep representation
        :param entity: (str) entity ID
        :param tiep: (Tiep) tiep instance
        :return:
        """

        if tiep_rep not in self.master_tieps:
            self.master_tieps[tiep_rep] = MasterTiep()

        self.master_tieps[tiep_rep].add_occurrence(entity, tiep)
