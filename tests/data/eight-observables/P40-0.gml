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
		node [ id 12 label "L" ]
		node [ id 13 label "C" ]
		node [ id 14 label "P" ]
		node [ id 15 label "P" ]
		node [ id 16 label "L" ]
		node [ id 17 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 7 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 5 label "next" ]
		edge [ source 3 target 12 label "cl" ]
		edge [ source 3 target 13 label "next" ]
		edge [ source 4 target 12 label "bound" ]
		edge [ source 5 target 9 label "cp" ]
		edge [ source 5 target 10 label "cl" ]
		edge [ source 5 target 11 label "next" ]
		edge [ source 6 target 13 label "cl" ]
		edge [ source 6 target 14 label "bound" ]
		edge [ source 7 target 8 label "cp" ]
		edge [ source 7 target 13 label "next" ]
		edge [ source 7 target 17 label "cl" ]
		edge [ source 11 target 15 label "cp" ]
		edge [ source 11 target 16 label "cl" ]
		edge [ source 13 target 14 label "cp" ]
		edge [ source 15 target 16 label "bound" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
	constrainNoEdge [
		source 8
		target 17
	]
	constrainNoEdge [
		source 9
		target 10
	]
]