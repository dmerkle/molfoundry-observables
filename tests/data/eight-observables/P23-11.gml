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
		node [ id 10 label "P" ]
		node [ id 11 label "L" ]
		node [ id 12 label "C" ]
		node [ id 13 label "L" ]
		node [ id 14 label "C" ]
		node [ id 15 label "P" ]
		node [ id 16 label "P" ]
		node [ id 17 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 6 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 5 label "cl" ]
		edge [ source 4 target 5 label "bound" ]
		edge [ source 6 target 7 label "cp" ]
		edge [ source 6 target 13 label "cl" ]
		edge [ source 6 target 14 label "next" ]
		edge [ source 8 target 14 label "cl" ]
		edge [ source 9 target 10 label "cp" ]
		edge [ source 9 target 11 label "cl" ]
		edge [ source 9 target 12 label "next" ]
		edge [ source 9 target 14 label "next" ]
		edge [ source 12 target 16 label "cp" ]
		edge [ source 12 target 17 label "cl" ]
		edge [ source 14 target 15 label "cp" ]
		edge [ source 16 target 17 label "bound" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
	constrainNoEdge [
		source 7
		target 13
	]
	constrainNoEdge [
		source 8
		target 15
	]
	constrainNoEdge [
		source 10
		target 11
	]
]