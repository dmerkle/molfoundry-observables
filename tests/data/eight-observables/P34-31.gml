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
		node [ id 8 label "C" ]
		node [ id 9 label "L" ]
		node [ id 10 label "C" ]
		node [ id 11 label "P" ]
		node [ id 12 label "P" ]
		node [ id 13 label "L" ]
		node [ id 14 label "L" ]
		node [ id 15 label "C" ]
		node [ id 16 label "P" ]
		node [ id 17 label "L" ]
		node [ id 18 label "L" ]
		node [ id 19 label "C" ]
		node [ id 20 label "C" ]
		node [ id 21 label "P" ]
		node [ id 22 label "P" ]
		node [ id 23 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 6 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 5 label "cl" ]
		edge [ source 4 target 5 label "bound" ]
		edge [ source 6 target 7 label "cp" ]
		edge [ source 6 target 18 label "cl" ]
		edge [ source 6 target 20 label "next" ]
		edge [ source 8 target 12 label "cp" ]
		edge [ source 8 target 13 label "cl" ]
		edge [ source 8 target 19 label "next" ]
		edge [ source 9 target 15 label "cl" ]
		edge [ source 10 target 11 label "cp" ]
		edge [ source 10 target 15 label "next" ]
		edge [ source 10 target 17 label "cl" ]
		edge [ source 10 target 19 label "next" ]
		edge [ source 12 target 13 label "bound" ]
		edge [ source 14 target 20 label "cl" ]
		edge [ source 15 target 16 label "cp" ]
		edge [ source 15 target 20 label "next" ]
		edge [ source 19 target 22 label "cp" ]
		edge [ source 19 target 23 label "cl" ]
		edge [ source 20 target 21 label "cp" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
	constrainNoEdge [
		source 7
		target 18
	]
	constrainNoEdge [
		source 9
		target 16
	]
	constrainNoEdge [
		source 11
		target 17
	]
	constrainNoEdge [
		source 14
		target 21
	]
	constrainNoEdge [
		source 22
		target 23
	]
]