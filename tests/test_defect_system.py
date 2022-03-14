import unittest
from unittest.mock import Mock, patch

from py_sc_fermi.defect_system import DefectSystem
from py_sc_fermi.defect_species import DefectSpecies
from py_sc_fermi.dos import DOS


class TestDefectSystemInit(unittest.TestCase):
    def test_defect_system_is_initialised(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        dos = Mock(spec=DOS)
        spin_pol = 1
        temperature = 298
        defect_system = DefectSystem(
            defect_species=mock_defect_species,
            volume=volume,
            dos=dos,
            spin_pol=spin_pol,
            temperature=temperature,
        )
        self.assertEqual(defect_system.volume, volume)
        self.assertEqual(defect_system.spin_pol, spin_pol)
        self.assertEqual(defect_system.dos, dos)
        self.assertEqual(defect_system.temperature, temperature)
        self.assertEqual(defect_system.defect_species[0], mock_defect_species[0])
        self.assertEqual(defect_system.defect_species[1], mock_defect_species[1])

class TestDefectSystem(unittest.TestCase):
    def setUp(self):
        volume = 100
        mock_defect_species = [Mock(spec=DefectSpecies), Mock(spec=DefectSpecies)]
        mock_defect_species[0].name = 'v_O'
        dos = Mock(spec=DOS)
        spin_pol = 1
        temperature = 298
        self.defect_system = DefectSystem(
            defect_species=mock_defect_species,
            volume=volume,
            dos=dos,
            spin_pol=spin_pol,
            temperature=temperature,
        )

    def test_defect_species_by_name(self):
        self.assertEqual(self.defect_system.defect_species_by_name('v_O'),self.defect_system.defect_species[0])
        
        
# if __name__ == "__main__":
#     unittest.main()
