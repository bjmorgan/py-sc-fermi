import unittest
from unittest.mock import Mock

from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.defect_charge_state import DefectChargeState

class TestDefectSpeciesInit(unittest.TestCase):

    def test_defect_species_is_initialised(self):
        name = 'foo'
        nsites = 2
        mock_charge_states = [ Mock(spec=DefectChargeState), Mock(spec=DefectChargeState) ]
        mock_charge_states[0].charge = 0
        mock_charge_states[1].charge = 1
        defect_species = DefectSpecies( name=name,
                                        nsites=nsites,
                                        charge_states=mock_charge_states )
        self.assertEqual( defect_species._name, name )
        self.assertEqual( defect_species._nsites, nsites )
        self.assertEqual( defect_species._charge_states[0], mock_charge_states[0] )
        self.assertEqual( defect_species._charge_states[1], mock_charge_states[1] )
        self.assertEqual( defect_species._fixed_concentration, None )

if __name__ == '__main__':
    unittest.main()
