rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "L" ]
		node [ id 4 label "C" ]
		node [ id 5 label "P" ]
		node [ id 6 label "L" ]
		node [ id 7 label "C" ]
		node [ id 8 label "C" ]
		node [ id 9 label "P" ]
		node [ id 10 label "L" ]
		node [ id 11 label "C" ]
		node [ id 12 label "P" ]
		node [ id 13 label "P" ]
		node [ id 14 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 8 label "next" ]
		edge [ source 0 target 11 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 11 label "cl" ]
		edge [ source 4 target 5 label "cp" ]
		edge [ source 4 target 6 label "cl" ]
		edge [ source 4 target 7 label "next" ]
		edge [ source 4 target 11 label "next" ]
		edge [ source 7 target 13 label "cp" ]
		edge [ source 7 target 14 label "cl" ]
		edge [ source 8 target 9 label "cp" ]
		edge [ source 8 target 10 label "cl" ]
		edge [ source 9 target 10 label "bound" ]
		edge [ source 11 target 12 label "cp" ]
		edge [ source 13 target 14 label "bound" ]
	]
	constrainNoEdge [
		source 3
		target 12
	]
	constrainNoEdge [
		source 5
		target 6
	]
]