rule [
	ruleID "locally-normalized-count"
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "P" ]
		node [ id 2 label "L" ]
		node [ id 3 label "C" ]
		node [ id 4 label "P" ]
		node [ id 5 label "L" ]
		node [ id 6 label "L" ]
		node [ id 7 label "C" ]
		node [ id 8 label "L" ]
		node [ id 9 label "P" ]
		node [ id 10 label "L" ]
		node [ id 11 label "C" ]
		node [ id 12 label "P" ]
		node [ id 13 label "C" ]
		node [ id 14 label "P" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 11 label "next" ]
		edge [ source 1 target 2 label "bound" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 8 label "cl" ]
		edge [ source 3 target 13 label "next" ]
		edge [ source 4 target 8 label "bound" ]
		edge [ source 5 target 13 label "cl" ]
		edge [ source 5 target 14 label "bound" ]
		edge [ source 6 target 11 label "cl" ]
		edge [ source 7 target 9 label "cp" ]
		edge [ source 7 target 10 label "cl" ]
		edge [ source 7 target 11 label "next" ]
		edge [ source 7 target 13 label "next" ]
		edge [ source 11 target 12 label "cp" ]
		edge [ source 13 target 14 label "cp" ]
	]
	constrainNoEdge [
		source 6
		target 12
	]
	constrainNoEdge [
		source 9
		target 10
	]
]