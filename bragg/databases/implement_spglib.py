'''
Load space group and magnetic space group databases.

TODO: This module wraps spglib to expose SG/MSG symmetry operations.
      Once data_tables/SG_generators.txt and data_tables/MSG_generators.txt
      are generated from spglib, this wrapper can be completed.
'''
from dataclasses import dataclass
from pathlib import Path

from .database import db_entry, Database

# TODO: uncomment once SG_generators.txt / MSG_generators.txt are generated from spglib
# import spglib

SW_PATH = Path(__file__).parent

####################################################################################################

@dataclass
class SG_entry(db_entry):
    number: int
    international_short: str
    international_full: str
    international: str
    schoenflies: str
    hall_number: int
    hall_symbol: str
    choice: str
    pointgroup_international: str
    pointgroup_schoenflies: str
    arithmetic_crystal_class_number: int
    arithmetic_crystal_class_symbol: str

    units = dict(number=1, international_short=None, international_full=None, international=None,
                 schoenflies=None, hall_number=None, hall_symbol=None, choice=None,
                 pointgroup_international=None, pointgroup_schoenflies=None,
                 arithmetic_crystal_class_number=None, arithmetic_crystal_class_symbol=None)

    @property
    def unique_label(self):
        return self.hall_number

    def get_symmetry_ops(self):
        # TODO: implement once spglib integration is complete
        # return spglib.get_symmetry_from_database(self.hall_number)
        raise NotImplementedError("spglib integration not yet complete; generate SG_generators.txt first")


class Database_SG(Database):
    '''Database of crystallographic space groups in all settings.

    TODO: Wrapper over the `spglib` library. Requires SG_generators.txt to be generated.
    '''

    source_filename = SW_PATH / 'data_tables/SG_generators.txt'
    entry_type = SG_entry
    header = ''
    entries = {}

    def __init__(self):
        with open(self.source_filename, 'r') as ff:
            lines = ff.readlines()

            self.header = ''.join(lines[:1])

            for line in lines[2:]:
                fields = [s.strip() for s in line.split('  ') if s.strip()]
                entry = self.entry_type(*fields)
                self.entries[entry.unique_label] = entry


# TODO: instantiate once SG_generators.txt is generated from spglib
# SG_data = Database_SG()
# SG_data.__doc__ = Database_SG.__doc__

####################################################################################################

@dataclass
class MSG_entry(db_entry):
    # uni_number litvin_number bns_number og_number number type

    uni_number: int
    litvin_number: int
    bns_number: str
    og_number: str
    number: int
    type: int

    units = dict(uni_number=None, litvin_number=None, bns_number=None,
                 og_number=None, number=None, type=None)

    @property
    def unique_label(self):
        return self.uni_number-1

    def get_symmetry_ops(self):
        '''Get symmetry operations of the Magnetic Space Group'''
        # TODO: implement once spglib integration is complete
        # return spglib.get_magnetic_symmetry_from_database(self.uni_number)
        raise NotImplementedError("spglib integration not yet complete; generate MSG_generators.txt first")


class Database_MSG(Database):
    '''Database of magnetic space groups.

    TODO: Wrapper over the `spglib` library. Requires MSG_generators.txt to be generated.
    '''

    source_filename = SW_PATH / 'data_tables/MSG_generators.txt'
    entry_type = MSG_entry
    header = ''
    entries = {}

    def __init__(self):
        with open(self.source_filename, 'r') as ff:
            lines = ff.readlines()

            self.header = ''.join(lines[:1])

            for line in lines[2:]:
                fields = [s.strip() for s in line.split('  ') if s.strip()]
                entry = self.entry_type(*fields)
                self.entries[entry.unique_label] = entry


# TODO: instantiate once MSG_generators.txt is generated from spglib
# MSG_data = Database_MSG()
# MSG_data.__doc__ = Database_MSG.__doc__
