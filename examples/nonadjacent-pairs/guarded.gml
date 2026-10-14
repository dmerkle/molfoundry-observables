rule [
	ruleID "isolated-alkene-pair-net-photocycloaddition"
	left [
		edge [ source 0 target 1 label "=" ]
		edge [ source 2 target 3 label "=" ]
	]
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
	]
	right [
		edge [ source 0 target 1 label "-" ]
		edge [ source 0 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
	]
	constrainNoEdge [
		source 0
		target 2
	]
	constrainNoEdge [
		source 1
		target 3
	]
]