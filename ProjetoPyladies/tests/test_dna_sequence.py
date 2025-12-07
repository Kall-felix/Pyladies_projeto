import unittest
from dna_sequence import DNASequence, InvalidSequenceError


class TestDNASequence(unittest.TestCase):
    
    def setUp(self):
        # Sequência base para testes
        self.seq = DNASequence("ATGCATGC", seq_id="test1")
    
    def test_create_valid(self):
        self.assertEqual(self.seq.sequence, "ATGCATGC")
        self.assertEqual(len(self.seq), 8)
    
    def test_invalid_sequence(self):
        # Deve levantar exceção para bases inválidas
        with self.assertRaises(InvalidSequenceError):
            DNASequence("ATGCUX")
    
    def test_gc_content(self):
        self.assertAlmostEqual(self.seq.gc_content(), 50.0)
        self.assertEqual(DNASequence("AATT").gc_content(), 0.0)
    
    def test_complement(self):
        comp = self.seq.complement()
        self.assertEqual(comp.sequence, "TACGTACG")
    
    def test_reverse_complement(self):
        rc = self.seq.reverse_complement()
        self.assertEqual(rc.sequence, "GCATGCAT")
    
    def test_transcribe(self):
        rna = self.seq.transcribe()
        self.assertEqual(rna, "AUGCAUGC")
    
    def test_find_all_occurrences(self):
        seq = DNASequence("ATGATGATG")
        positions = seq.find_all_occurrences("ATG")
        self.assertEqual(positions, [0, 3, 6])
    
    def test_base_composition(self):
        comp = self.seq.get_base_composition()
        self.assertEqual(comp['A'], 2)
        self.assertEqual(comp['T'], 2)
        self.assertEqual(comp['G'], 2)
        self.assertEqual(comp['C'], 2)
        self.assertEqual(comp['N'], 0)
    
    def test_slicing(self):
        sub = self.seq[0:4]
        self.assertTrue(isinstance(sub, DNASequence))
        self.assertEqual(sub.sequence, "ATGC")
        self.assertEqual(sub.id, "test1_slice_4")


if __name__ == '__main__':
    unittest.main()