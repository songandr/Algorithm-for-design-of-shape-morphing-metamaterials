# Algorithm-for-design-of-shape-morphing-metamaterials
Elastic energy minimization algorithm discussed in the manuscript "Algorithmic design framework for shape-morphing metamaterials" section 3.

The files in this repository are all that are necessary for running the algorithm discussed in section 3 of "Algorithmic design framework for shape-morphing metamaterials." To generate designs of desired topology:

1. Set up initial conditions for x in 2D, y in 3D, R = Id in R^{3x3}, the rigidity constraint matrices L_0 and L, and panel assignments Pj.
2. Set the desired tolerance (e-5 recommended for max strain ~0.1%). 
3. Run the minimizationAlgorithm script with the specified inputs from the above steps.
4. Use plot4vectors3D for visualization. If a new topology (not Miura, Miura4x4, rotating squares, helical waterbomb), the necessary edge graphs must be written into the plot1vector.m, plot1vector3D.m files beforehand.

See Tests directory for examples. Post-processing results is recommended for all design generation, done effectively by running minimizationAlgorithm again on outputs perturbed by RNG vectors sized 2N and 3N (where N is the number of nodes) both of which sampled by a normal random distribution ~ (0, epsilon). This ensures results indeed sit at local minimum in the energy landscape (opposed to saddle points, etc).
