rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "P" ]
		node [ id 5 label "L" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 5 label "cl" ]
		edge [ source 4 target 5 label "bound" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
]