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
		node [ id 12 label "P" ]
		node [ id 13 label "L" ]
		node [ id 14 label "C" ]
		edge [ source 0 target 1 label "cp" ]
		edge [ source 0 target 2 label "cl" ]
		edge [ source 0 target 3 label "next" ]
		edge [ source 0 target 14 label "next" ]
		edge [ source 3 target 4 label "cp" ]
		edge [ source 3 target 10 label "cl" ]
		edge [ source 3 target 11 label "next" ]
		edge [ source 5 target 14 label "cp" ]
		edge [ source 6 target 14 label "cl" ]
		edge [ source 7 target 11 label "cl" ]
		edge [ source 8 target 9 label "cp" ]
		edge [ source 8 target 11 label "next" ]
		edge [ source 8 target 13 label "cl" ]
		edge [ source 8 target 14 label "next" ]
		edge [ source 11 target 12 label "cp" ]
	]
	constrainNoEdge [
		source 1
		target 2
	]
	constrainNoEdge [
		source 4
		target 10
	]
	constrainNoEdge [
		source 5
		target 6
	]
	constrainNoEdge [
		source 7
		target 12
	]
	constrainNoEdge [
		source 9
		target 13
	]
]