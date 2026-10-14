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
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 6 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 5 label "cl" ]
		edge [ source 3 target 9 label "next" ]
		edge [ source 4 target 5 label "bound" ]
		edge [ source 6 target 7 label "cp" ]
		edge [ source 6 target 8 label "cl" ]
		edge [ source 6 target 9 label "next" ]
		edge [ source 9 target 10 label "cp" ]
		edge [ source 9 target 11 label "cl" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
	constrainNoEdge [
		source 7
		target 8
	]
	constrainNoEdge [
		source 10
		target 11
	]
]