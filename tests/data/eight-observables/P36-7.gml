rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "P" ]
		node [ id 5 label "P" ]
		node [ id 6 label "L" ]
		node [ id 7 label "C" ]
		node [ id 8 label "P" ]
		node [ id 9 label "L" ]
		node [ id 10 label "C" ]
		node [ id 11 label "P" ]
		node [ id 12 label "L" ]
		node [ id 13 label "L" ]
		node [ id 14 label "C" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 7 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 9 label "cl" ]
		edge [ source 3 target 10 label "next" ]
		edge [ source 5 target 14 label "cp" ]
		edge [ source 6 target 14 label "cl" ]
		edge [ source 7 target 8 label "cp" ]
		edge [ source 7 target 13 label "cl" ]
		edge [ source 7 target 14 label "next" ]
		edge [ source 10 target 11 label "cp" ]
		edge [ source 10 target 12 label "cl" ]
		edge [ source 10 target 14 label "next" ]
		edge [ source 11 target 12 label "bound" ]
	]
	constrainNoEdge [
		source 4
		target 9
	]
	constrainNoEdge [
		source 5
		target 6
	]
	constrainNoEdge [
		source 8
		target 13
	]
]