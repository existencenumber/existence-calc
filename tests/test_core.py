import unittest
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from app import compute

class TestCore(unittest.TestCase):
    def test_polynomial(self):
        result = compute("Sum(n**2,(n,1,oo))")
        self.assertEqual(result['status'], 'success')
        self.assertAlmostEqual(result['value'], 0.0, places=4)

    def test_zeta(self):
        result = compute("Sum(n,(n,1,oo))")
        self.assertAlmostEqual(result['value'], -1/12, places=4)

    def test_geometric(self):
        result = compute("Sum(2**n,(n,0,oo))")
        self.assertAlmostEqual(result['value'], -1.0, places=4)

    def test_alternating_harmonic(self):
        result = compute("Sum((-1)**(n+1)/n,(n,1,oo))")
        self.assertAlmostEqual(result['value'], 0.693147, places=4)

if __name__ == '__main__':
    unittest.main()
