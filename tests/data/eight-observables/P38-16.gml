rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "P" ]
		node [ id 5 label "C" ]
		node [ id 6 label "P" ]
		node [ id 7 label "C" ]
		node [ id 8 label "C" ]
		node [ id 9 label "P" ]
		node [ id 10 label "L" ]
		node [ id 11 label "C" ]
		node [ id 12 label "P" ]
		node [ id 13 label "P" ]
		node [ id 14 label "L" ]
		node [ id 15 label "L" ]
		node [ id 16 label "C" ]
		node [ id 17 label "L" ]
		node [ id 18 label "P" ]
		node [ id 19 label "L" ]
		node [ id 20 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 5 label "next" ]
		edge [ source 0 target 7 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 8 label "next" ]
		edge [ source 3 target 11 label "next" ]
		edge [ source 3 target 14 label "cl" ]
		edge [ source 4 target 14 label "bound" ]
		edge [ source 5 target 6 label "cp" ]
		edge [ source 5 target 11 label "next" ]
		edge [ source 5 target 17 label "cl" ]
		edge [ source 7 target 13 label "cp" ]
		edge [ source 7 target 15 label "cl" ]
		edge [ source 7 target 16 label "next" ]
		edge [ source 8 target 9 label "cp" ]
		edge [ source 8 target 10 label "cl" ]
		edge [ source 9 target 10 label "bound" ]
		edge [ source 11 target 12 label "cp" ]
		edge [ source 11 target 20 label "cl" ]
		edge [ source 16 target 18 label "cp" ]
		edge [ source 16 target 19 label "cl" ]
		edge [ source 18 target 19 label "bound" ]
	]
	constrainNoEdge [
		source 6
		target 17
	]
	constrainNoEdge [
		source 12
		target 20
	]
	constrainNoEdge [
		source 13
		target 15
	]
]