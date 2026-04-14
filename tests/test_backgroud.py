"""
test_background.py
Unit tests to verify the mathematical correctness of the Markov Background Model.
"""
import sys
import os
import unittest

# Ensure the src directory is in the path for importing
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from background import train_background_model

class TestBackgroundModel(unittest.TestCase):

    def setUp(self):
        # A simple synthetic sequence for testing
        self.test_seqs = ["ACGTACGTACGT", "AAAAACCCCCGGGGGTTTTT"]

    def test_zero_order_probabilities_sum_to_one(self):
        """Ensures that 0-order background probabilities sum to exactly 1.0"""
        model = train_background_model(self.test_seqs, order=0)
        total_prob = sum(model[""].values())
        self.assertAlmostEqual(total_prob, 1.0, places=5, 
                               msg="0-order probabilities do not sum to 1.0")

    def test_first_order_probabilities_sum_to_one(self):
        """Ensures that for any given base (A, C, G, T), the transition probs sum to 1.0"""
        model = train_background_model(self.test_seqs, order=1)
        bases = ['A', 'C', 'G', 'T']
        for base in bases:
            total_prob = sum(model[base].values())
            self.assertAlmostEqual(total_prob, 1.0, places=5, 
                                   msg=f"1st-order transition probabilities for context '{base}' do not sum to 1.0")

if __name__ == '__main__':
    unittest.main()
