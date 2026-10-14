rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "L" ]
		node [ id 5 label "C" ]
		node [ id 6 label "P" ]
		node [ id 7 label "P" ]
		node [ id 8 label "L" ]
		node [ id 9 label "L" ]
		node [ id 10 label "C" ]
		node [ id 11 label "P" ]
		node [ id 12 label "L" ]
		node [ id 13 label "C" ]
		node [ id 14 label "C" ]
		node [ id 15 label "P" ]
		node [ id 16 label "L" ]
		node [ id 17 label "C" ]
		node [ id 18 label "P" ]
		node [ id 19 label "P" ]
		node [ id 20 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 14 label "next" ]
		edge [ source 0 target 17 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 7 label "cp" ]
		edge [ source 3 target 8 label "cl" ]
		edge [ source 3 target 13 label "next" ]
		edge [ source 4 target 10 label "cl" ]
		edge [ source 5 target 6 label "cp" ]
		edge [ source 5 target 10 label "next" ]
		edge [ source 5 target 12 label "cl" ]
		edge [ source 5 target 13 label "next" ]
		edge [ source 7 target 8 label "bound" ]
		edge [ source 9 target 17 label "cl" ]
		edge [ source 10 target 11 label "cp" ]
		edge [ source 10 target 17 label "next" ]
		edge [ source 13 target 19 label "cp" ]
		edge [ source 13 target 20 label "cl" ]
		edge [ source 14 target 15 label "cp" ]
		edge [ source 14 target 16 label "cl" ]
		edge [ source 15 target 16 label "bound" ]
		edge [ source 17 target 18 label "cp" ]
	]
	constrainNoEdge [
		source 4
		target 11
	]
	constrainNoEdge [
		source 6
		target 12
	]
	constrainNoEdge [
		source 9
		target 18
	]
	constrainNoEdge [
		source 19
		target 20
	]
]