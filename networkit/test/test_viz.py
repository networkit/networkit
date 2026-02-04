#!/usr/bin/env python3

import unittest
import networkit as nk


class TestViz(unittest.TestCase):

    def getSmallGraph(self, weighted=False, directed=False):
        """Create a small test graph"""
        G = nk.Graph(5, weighted, directed)
        G.addEdge(0, 1, 1.0)
        G.addEdge(0, 2, 2.0)
        G.addEdge(1, 2, 3.0)
        G.addEdge(1, 3, 4.0)
        G.addEdge(2, 3, 5.0)
        G.addEdge(3, 4, 6.0)
        return G

    def testMaxentStress2DCoordinates(self):
        """Test that 2D start coordinates work (backward compatibility, warm start)"""
        G = self.getSmallGraph(weighted=True)

        # Create 2D start coordinates for all nodes (initial positions)
        coords_2d = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.5, 0.5)]

        # Create MaxentStress with 2D start coordinates
        maxent = nk.viz.MaxentStress(
            G,
            dim=2,
            k=2,
            coordinates=coords_2d,
            tolerance=1e-3
        )

        # Run the algorithm (will optimize from the start coordinates)
        maxent.run()

        # Get the resulting coordinates (after optimization)
        result_coords = maxent.getCoordinates()

        # Verify we got coordinates for all nodes
        self.assertEqual(len(result_coords), G.numberOfNodes())

        # Verify each coordinate is 2D
        for coord in result_coords:
            self.assertEqual(len(coord), 2)

    def testMaxentStress3DCoordinates(self):
        """Test that 3D start coordinates work correctly (warm start)"""
        G = self.getSmallGraph(weighted=True)

        # Create 3D start coordinates for all nodes (these are initial positions)
        coords_3d = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.5, 0.5, 1.0)
        ]

        # Create MaxentStress with 3D start coordinates
        maxent = nk.viz.MaxentStress(
            G,
            dim=3,
            k=2,
            coordinates=coords_3d,
            tolerance=1e-3
        )

        # Run the algorithm (will optimize from the start coordinates)
        maxent.run()

        # Get the resulting coordinates (after optimization)
        result_coords = maxent.getCoordinates()

        # Verify we got coordinates for all nodes
        self.assertEqual(len(result_coords), G.numberOfNodes())

        # Verify each coordinate is 3D
        for coord in result_coords:
            self.assertEqual(len(coord), 3)


    def testMaxentStressInvalidCoordinateDimensions(self):
        """Test that invalid coordinate dimensions raise ValueError"""
        G = self.getSmallGraph(weighted=True)

        # Test 1D coordinates (invalid)
        coords_1d = [(0.0,), (1.0,), (2.0,), (3.0,), (4.0,)]
        with self.assertRaises(ValueError) as context:
            nk.viz.MaxentStress(G, dim=2, k=2, coordinates=coords_1d, tolerance=1e-3)
        self.assertIn("Coordinates must be 2D or 3D", str(context.exception))
        self.assertIn("got 1 dimensions", str(context.exception))

        # Test 4D coordinates (invalid)
        coords_4d = [
            (0.0, 0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
            (1.0, 1.0, 0.0, 0.0),
            (0.5, 0.5, 1.0, 1.0)
        ]
        with self.assertRaises(ValueError) as context:
            nk.viz.MaxentStress(G, dim=3, k=2, coordinates=coords_4d, tolerance=1e-3)
        self.assertIn("Coordinates must be 2D or 3D", str(context.exception))
        self.assertIn("got 4 dimensions", str(context.exception))

    def testMaxentStressWithoutCoordinates(self):
        """Test that MaxentStress works without providing start coordinates"""
        G = self.getSmallGraph(weighted=True)

        # Test 2D without coordinates
        maxent_2d = nk.viz.MaxentStress(G, dim=2, k=2, tolerance=1e-3)
        maxent_2d.run()
        result_2d = maxent_2d.getCoordinates()
        self.assertEqual(len(result_2d), G.numberOfNodes())
        for coord in result_2d:
            self.assertEqual(len(coord), 2)

        # Test 3D without coordinates
        maxent_3d = nk.viz.MaxentStress(G, dim=3, k=2, tolerance=1e-3)
        maxent_3d.run()
        result_3d = maxent_3d.getCoordinates()
        self.assertEqual(len(result_3d), G.numberOfNodes())
        for coord in result_3d:
            self.assertEqual(len(coord), 3)

    def testMaxentStressMethods(self):
        """Test that various MaxentStress methods work correctly"""
        G = self.getSmallGraph(weighted=True)

        # Test with 3D coordinates
        coords_3d = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (1.0, 1.0, 0.0),
            (0.5, 0.5, 1.0)
        ]

        maxent = nk.viz.MaxentStress(G, dim=3, k=2, coordinates=coords_3d, tolerance=1e-3)
        maxent.run()

        # Test various methods don't crash
        maxent.scaleLayout()
        scaling_factor = maxent.computeScalingFactor()
        self.assertIsInstance(scaling_factor, float)

        stress = maxent.fullStressMeasure()
        self.assertIsInstance(stress, float)

        maxent_measure = maxent.maxentMeasure()
        self.assertIsInstance(maxent_measure, float)

        mde = maxent.meanDistanceError()
        self.assertIsInstance(mde, float)

        ldme = maxent.ldme()
        self.assertIsInstance(ldme, float)

if __name__ == "__main__":
    unittest.main()
