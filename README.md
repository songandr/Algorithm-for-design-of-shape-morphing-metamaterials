# Algorithm for design of shape morphing metamaterials
Elastic energy minimization algorithm discussed in the manuscript "Algorithmic design framework for shape-morphing metamaterials" section 3.

The files in the /MATLABscript directory are all that are necessary for running the algorithm discussed in section 3 of "Algorithmic design framework for shape-morphing metamaterials." To generate designs of desired topology:

1. Set up initial conditions for $\mathbf{x}\in\mathbb{R}^{2I}$, $\mathbf{y}\in\mathbb{R}^{3I}$, $\mathbf{R}_j=I\in\mathbb{R}^{3x3}$, the rigidity constraint matrices $\mathbf{L}_0$ and $mathbf{L}$, and panel assignments $P_j$. $\mathbf{L}$ must be ordered such that the first $3|B_1|$ rows correspond to $\mathbf{d}_1$ and the next $3|B_2|$ rows correspond to $\mathbf{d}_2$. The last row is reserved for $\boldsymbol{\chi}_1$ to fix the first node to the origin. Any redundant constraints must be removed to successfully run post-processing geometric stiffness calculations e.g. for an abelian group of isometries $\mathbf{g}_1$ and $\mathbf{g}_2$, if $\mathbf{y}_{i'} = \mathbf{g}_1(\mathbf{g}_2(\mathbf{y}_i) = \mathbf{g}_1(\mathbf{y}_{i''}) = \mathbf{g}_2(\mathbf{y}_{i'''}), then only one of the constraints among \mathbf{y}_{i'} = \mathbf{g}_1(\mathbf{y}_{i''}) and \mathbf{y}_{i'} = \mathbf{g}_2(\mathbf{y}_{i'''}) should be embedded into $\mathbf{L}$.
2. Set the desired tolerance (e-5 recommended for max strain ~0.1%). 
3. Run the minimizationAlgorithm script with the specified inputs from the above steps.
4. Use plot4vectors3D for visualization within Mathematica. If a new topology (not Miura, Miura4x4, rotating squares, helical waterbomb), the necessary edge graphs must be written into the plot1vector.m, plot1vector3D.m files beforehand.
   <br> Alternatively, for better visualization in Mathematica, use makeMathematicaCoords.m to output Mathematica compatible coordinate data. See /Mathematica_Files directory for Mathematica visualization examples. Geo_to_STL.m can also be used to convert coordinate data with a connectivity matrix to an STL file for visualization in other software.

See /MATLABscript/Tests directory for examples. Post-processing results is recommended for all design generation, done effectively by running minimizationAlgorithm again on outputs perturbed by RNG vectors sized 2N and 3N (where N is the number of nodes) both of which sampled by a normal random distribution ~ (0, epsilon). This ensures results indeed sit at local minimum in the energy landscape (opposed to saddle points, etc).

Geometric stiffness can be calculated for any solution through geometricStiffness.m given additional information: $|B_1|$, $|B_2|$, and whether the pattern is planar (dim = 2) or not (dim = 3).

Figures generated for use in the manuscript are contained in the /Figures directory.
