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
		node [ id 7 label "L" ]
		node [ id 8 label "C" ]
		node [ id 9 label "P" ]
		node [ id 10 label "L" ]
		node [ id 11 label "C" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 11 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 7 label "cl" ]
		edge [ source 3 target 8 label "next" ]
		edge [ source 5 target 11 label "cp" ]
		edge [ source 6 target 11 label "cl" ]
		edge [ source 8 target 9 label "cp" ]
		edge [ source 8 target 10 label "cl" ]
		edge [ source 8 target 11 label "next" ]
		edge [ source 9 target 10 label "bound" ]
	]
	constrainNoEdge [
		source 4
		target 7
	]
	constrainNoEdge [
		source 5
		target 6
	]
]