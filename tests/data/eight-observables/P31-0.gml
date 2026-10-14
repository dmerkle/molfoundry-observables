rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "P" ]
		node [ id 5 label "C" ]
		node [ id 6 label "L" ]
		node [ id 7 label "C" ]
		node [ id 8 label "P" ]
		node [ id 9 label "P" ]
		node [ id 10 label "L" ]
		node [ id 11 label "C" ]
		node [ id 12 label "P" ]
		node [ id 13 label "L" ]
		node [ id 14 label "L" ]
		node [ id 15 label "C" ]
		node [ id 16 label "C" ]
		node [ id 17 label "P" ]
		node [ id 18 label "L" ]
		node [ id 19 label "P" ]
		node [ id 20 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 11 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 14 label "cl" ]
		edge [ source 3 target 16 label "next" ]
		edge [ source 5 target 9 label "cp" ]
		edge [ source 5 target 10 label "cl" ]
		edge [ source 5 target 15 label "next" ]
		edge [ source 6 target 11 label "cl" ]
		edge [ source 7 target 8 label "cp" ]
		edge [ source 7 target 11 label "next" ]
		edge [ source 7 target 13 label "cl" ]
		edge [ source 7 target 15 label "next" ]
		edge [ source 9 target 10 label "bound" ]
		edge [ source 11 target 12 label "cp" ]
		edge [ source 15 target 19 label "cp" ]
		edge [ source 15 target 20 label "cl" ]
		edge [ source 16 target 17 label "cp" ]
		edge [ source 16 target 18 label "cl" ]
		edge [ source 17 target 18 label "bound" ]
	]
	constrainNoEdge [
		source 4
		target 14
	]
	constrainNoEdge [
		source 6
		target 12
	]
	constrainNoEdge [
		source 8
		target 13
	]
	constrainNoEdge [
		source 19
		target 20
	]
]