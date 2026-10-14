rule [
	ruleID " *common* "
	left [
		edge [ source 0 target 1 label "=" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 4 target 5 label "=" ]
	]
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "H" ]
		node [ id 4 label "O" ]
		node [ id 5 label "C" ]
		node [ id 6 label "C" ]
		edge [ source 6 target 5 label "=" ]
	]
	right [
		edge [ source 0 target 1 label "-" ]
		edge [ source 1 target 2 label "=" ]
		edge [ source 4 target 5 label "-" ]
		edge [ source 3 target 4 label "-" ]
		edge [ source 0 target 5 label "-" ]
	]
]