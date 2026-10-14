rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "P" ]
		node [ id 5 label "L" ]
		node [ id 6 label "C" ]
		node [ id 7 label "P" ]
		node [ id 8 label "L" ]
		node [ id 9 label "C" ]
		node [ id 10 label "C" ]
		node [ id 11 label "P" ]
		node [ id 12 label "L" ]
		node [ id 13 label "P" ]
		node [ id 14 label "L" ]
		node [ id 15 label "L" ]
		node [ id 16 label "C" ]
		node [ id 17 label "P" ]
		node [ id 18 label "C" ]
		node [ id 19 label "L" ]
		node [ id 20 label "P" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 6 label "next" ]
		edge [ source 0 target 9 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 6 label "next" ]
		edge [ source 3 target 14 label "cl" ]
		edge [ source 3 target 16 label "next" ]
		edge [ source 4 target 14 label "bound" ]
		edge [ source 5 target 18 label "cl" ]
		edge [ source 5 target 20 label "bound" ]
		edge [ source 6 target 7 label "cp" ]
		edge [ source 6 target 19 label "cl" ]
		edge [ source 8 target 16 label "cl" ]
		edge [ source 9 target 13 label "cp" ]
		edge [ source 9 target 15 label "cl" ]
		edge [ source 9 target 18 label "next" ]
		edge [ source 10 target 11 label "cp" ]
		edge [ source 10 target 12 label "cl" ]
		edge [ source 10 target 16 label "next" ]
		edge [ source 11 target 12 label "bound" ]
		edge [ source 16 target 17 label "cp" ]
		edge [ source 18 target 20 label "cp" ]
	]
	constrainNoEdge [
		source 7
		target 19
	]
	constrainNoEdge [
		source 8
		target 17
	]
	constrainNoEdge [
		source 13
		target 15
	]
]