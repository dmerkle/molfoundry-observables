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
		node [ id 8 label "P" ]
		node [ id 9 label "L" ]
		node [ id 10 label "C" ]
		node [ id 11 label "L" ]
		node [ id 12 label "C" ]
		node [ id 13 label "P" ]
		node [ id 14 label "P" ]
		node [ id 15 label "L" ]
		node [ id 16 label "L" ]
		node [ id 17 label "C" ]
		node [ id 18 label "P" ]
		node [ id 19 label "L" ]
		node [ id 20 label "C" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 6 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 5 label "cl" ]
		edge [ source 4 target 5 label "bound" ]
		edge [ source 6 target 7 label "cp" ]
		edge [ source 6 target 16 label "cl" ]
		edge [ source 6 target 17 label "next" ]
		edge [ source 8 target 20 label "cp" ]
		edge [ source 9 target 20 label "cl" ]
		edge [ source 10 target 14 label "cp" ]
		edge [ source 10 target 15 label "cl" ]
		edge [ source 10 target 20 label "next" ]
		edge [ source 11 target 17 label "cl" ]
		edge [ source 12 target 13 label "cp" ]
		edge [ source 12 target 17 label "next" ]
		edge [ source 12 target 19 label "cl" ]
		edge [ source 12 target 20 label "next" ]
		edge [ source 14 target 15 label "bound" ]
		edge [ source 17 target 18 label "cp" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
	constrainNoEdge [
		source 7
		target 16
	]
	constrainNoEdge [
		source 8
		target 9
	]
	constrainNoEdge [
		source 11
		target 18
	]
	constrainNoEdge [
		source 13
		target 19
	]
]