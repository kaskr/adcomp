/* $Id: cppad_colpack.cpp 3240 2014-04-29 13:12:19Z bradbell $ */
/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2003-14 Bradley M. Bell

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    GNU General Public License Version 3.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

# include <vector>
# include <cppad/vector.hpp>
# include <ColPack/ColPackHeaders.h>

namespace CppAD {
	void cppad_colpack_general(
		      CppAD::vector<size_t>&         color         ,
		size_t                               m             ,
		size_t                               n             ,
		const CppAD::vector<unsigned int*>&  adolc_pattern )
	{	size_t i, k;
		CPPAD_ASSERT_UNKNOWN( adolc_pattern.size() == m );
	
		// Use adolc sparsity pattern to create corresponding bipartite graph
		ColPack::BipartiteGraphPartialColoringInterface graph(
				SRC_MEM_ADOLC,
				adolc_pattern.data(),
				m,
				n
		);
	
		// row ordered Partial-Distance-Two-Coloring of the bipartite graph 
		graph.PartialDistanceTwoColoring(
			"SMALLEST_LAST", "ROW_PARTIAL_DISTANCE_TWO"
		);
	
		// Use coloring information to create seed matrix
		int n_seed_row;
		int n_seed_col;
		double** seed_matrix = graph.GetSeedMatrix(&n_seed_row, &n_seed_col);
		CPPAD_ASSERT_UNKNOWN( size_t(n_seed_col) == m );
	
		// now return coloring in format required by CppAD
		for(i = 0; i < m; i++)
			color[i] = m;
		for(k = 0; k < size_t(n_seed_row); k++)
		{	for(i = 0; i < m; i++)
			{	if( seed_matrix[k][i] != 0.0 ) 
				{	// check that no row appears twice in the coloring
					CPPAD_ASSERT_UNKNOWN( color[i] == m );
					color[i] = k;
				}
			}
		}
# ifndef NDEBUG
		// check that all non-zero rows appear in the coloring
		for(i = 0; i < m; i++)
			CPPAD_ASSERT_UNKNOWN(color[i] < m || adolc_pattern[i][0] == 0);

		// check that no rows with the same color have overlapping entries
		CppAD::vector<bool> found(n);
		for(k = 0; k < size_t(n_seed_row); k++)
		{	size_t j, ell;
			for(j = 0; j < n; j++)
				found[j] = false;
			for(i = 0; i < m; i++) if( color[i] == k )
			{	for(ell = 0; ell < adolc_pattern[i][0]; ell++)
				{	j = adolc_pattern[i][1 + ell];
					CPPAD_ASSERT_UNKNOWN( ! found[j] );
					found[j] = true;
				}
			}
		}
# endif
		return;
	}
	// The following routine is not yet used or tested.
	void cppad_colpack_symmetric(
		      CppAD::vector<size_t>&         color         ,
		size_t                               n             ,
		const CppAD::vector<unsigned int*>&  adolc_pattern )
	{	size_t i, k;
		CPPAD_ASSERT_UNKNOWN( adolc_pattern.size() == n );
	
		// Use adolc sparsity pattern to create corresponding bipartite graph
		ColPack::GraphColoringInterface graph(
				SRC_MEM_ADOLC,
				adolc_pattern.data(),
				n
		);
	
		// Color the graph with the speciied ordering
		// graph.Coloring("SMALLEST_LAST", "STAR") is slower in adolc testing
		graph.Coloring("SMALLEST_LAST", "ACYCLIC_FOR_INDIRECT_RECOVERY");
	
		// Use coloring information to create seed matrix
		int n_seed_row;
		int n_seed_col;
		double** seed_matrix = graph.GetSeedMatrix(&n_seed_row, &n_seed_col);
		CPPAD_ASSERT_UNKNOWN( size_t(n_seed_col) == n );
	
		// now return coloring in format required by CppAD
		for(i = 0; i < n; i++)
			color[i] = n;
		for(k = 0; k < size_t(n_seed_row); k++)
		{	for(i = 0; i < n; i++)
			{	if( seed_matrix[k][i] != 0.0 ) 
				{	CPPAD_ASSERT_UNKNOWN( color[i] == n );
					color[i] = k;
				}
			}
		}
# ifndef NDEBUG
		for(i = 0; i < n; i++)
			CPPAD_ASSERT_UNKNOWN(color[i] < n || adolc_pattern[i][0] == 0);

		// The coloring above will probably fail this  test.
		// Check that no rows with the same color have overlapping entries:
		CppAD::vector<bool> found(n);
		for(k = 0; k < size_t(n_seed_row); k++)
		{	size_t j, ell;
			for(j = 0; j < n; j++)
				found[j] = false;
			for(i = 0; i < n; i++) if( color[i] == k )
			{	for(ell = 0; ell < adolc_pattern[i][0]; ell++)
				{	j = adolc_pattern[i][1 + ell];
					CPPAD_ASSERT_UNKNOWN( ! found[j] );
					found[j] = true;
				}
			}
		}
# endif
		return;
	}
}
