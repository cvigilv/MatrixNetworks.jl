module MatrixNetworks
	using LinearAlgebra
	using DelimitedFiles

	export random_adjmat, symmetrize, neighbours, degree, load_matrix_network

	"""
		symmetrize(M::Matrix) -> SymM::Symmetric
	
		
	Convert rectangular (`n_L` × `n_T`) matrix into a symmetric (`n_L` + `n_T`) order 
	matrix.

	# Examples
	```julia-repl
	julia> A = [ 1 0 1 ; 0 0 1 ]
	2×3 Array{Int64,2}:
	1  0  1
	0  0  1
	julia> symmetrize(A)
	5×5 Array{Int64,2}:
	0  0  1  0  1
	0  0  0  0  1
	1  0  0  0  0
	0  0  0  0  0
	1  1  0  0  0
	```

	# Arguments
	- `M::Matrix`: rectangular (`n_L` × `n_T`) matrix.

	"""
	function symmetrize(M::Matrix)::Symmetric
		# Get dimensions of non-symmetric matrix
		R, C = size(M)
		
		# Create empty matrices to fill spaces between our matrix
		Rzeros = zeros(Int64, R, R)
		Czeros = zeros(Int64, C, C)
		
		# Create symmetric matrix
		M	 = [ Rzeros M ; M' Czeros ]
		SymM = Symmetric(M)

		return SymM
	end

	"""
		neighbours(M::Matrix, node::Int, orientation::String = "row") -> neighbours_array::Array{Int}
	
		
	Get neighbours of node in adjacency matrix.
	
	# Examples
	```julia-repl
	julia> A = [ 0 0 1 0 1 ; 0 0 0 0 1 ; 1 0 0 0 0 ; 0 0 0 0 0 ; 1 1 0 0 0 ]
	5×5 Array{Int64,2}:
	0  0  1  0  1
	0  0  0  0  1
	1  0  0  0  0
	0  0  0  0  0
	1  1  0  0  0
	julia> neighbours(A, 1)
	2-element Array{Int64,1}:
	3
	5
	```
	
	```julia-repl
	julia> B = [ 0 1 0 0 ; 1 0 1 0 ] 
	2×4 Array{Int64,2}:
	0  1  0  0
	1  0  1  0
	julia> neighbours(B, 2, "row")
	2-element Array{Int64,1}:
	1
	3
	julia> neighbours(B, 2, "column")
	1-element Array{Int64,1}:
	2
	```
	
	# Arguments
	- `M::Matrix`: adjacency matrix.
	- `node::Int`: idx of node.
	- `orientation::String="row"` : orientation of search (column or row-wise).

	"""
	function neighbours(M::Matrix, node::Int, orientation::String = "row")::Array{Int}
		if orientation == "row"
			node_adjacency	 = M[node,:]
		elseif orientation == "column"
			node_adjacency	 = M[:, node]
		else
			error("Orientation value must be either \"row\" of \"column\"")
		end

		neighbours_array = [idx for (idx,neighbour) in enumerate(node_adjacency) 
							if neighbour != 0]

		return neighbours_array
	end

	"""
		degree(M::Matrix, node::Int, orientation::String = "row") -> node_degree::Int
	
		
	Get degree of node in symmetric adjacency matrix.
	
	# Examples
	```julia-repl
	julia> A = [ 0 0 1 0 1 ; 0 0 0 0 1 ; 1 0 0 0 0 ; 0 0 0 0 0 ; 1 1 0 0 0 ]
	5×5 Array{Int64,2}:
	0  0  1  0  1
	0  0  0  0  1
	1  0  0  0  0
	0  0  0  0  0
	1  1  0  0  0
	julia> degree(A, 1)
	2
	```
	
	```julia-repl
	julia> B = [ 0 1 0 0 ; 1 0 1 0 ] 
	2×4 Array{Int64,2}:
	0  1  0  0
	1  0  1  0
	julia> degree(B, 2, "row")
	2
	julia> degree(B, 2, "column")
	1
	```
	
	# Arguments
	- `M::Matrix`: adjacency matrix.
	- `node::Int`: idx of node.
	- `orientation::String="row"` : orientation of search (column or row-wise).

	"""
	function degree(M::Matrix, node::Int, orientation::String = "row")::Int
		if orientation == "row"
			node_neighbours = neighbours(M, node, "row")
		elseif orientation == "column"
			node_neighbours = neighbours(M, node, "column")
		else
			error("Orientation value must be either \"row\" of \"column\"")
		end

		node_degree = length(node_neighbours)

		return node_degree
	end

	"""
		load_matrix_network(filepath::String, separator::Char = ' ') -> A::Matrix{Int64}
	"""
	function load_matrix_network(filepath::String, separator::Char = ' ')::Matrix{Int64}
		A = readdlm(filepath, separator, Int64, '\n')
		return A
	end

end
