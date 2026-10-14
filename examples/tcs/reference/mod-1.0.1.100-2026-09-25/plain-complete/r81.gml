rule [
	ruleID "r_{81}"
	labelType "string"
	left [
		edge [ source 0 target 1 label "=" ]
		edge [ source 0 target 6 label "-" ]
		edge [ source 6 target 7 label "-" ]
		edge [ source 8 target 9 label "=" ]
	]
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		node [ id 2 label "O" ]
		node [ id 3 label "H" ]
		node [ id 4 label "O" ]
		node [ id 5 label "C" ]
		node [ id 6 label "O" ]
		node [ id 7 label "H" ]
		node [ id 8 label "O" ]
		node [ id 9 label "C" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
		edge [ source 4 target 5 label "=" ]
	]
	right [
		edge [ source 0 target 1 label "-" ]
		edge [ source 0 target 6 label "=" ]
		edge [ source 1 target 9 label "-" ]
		edge [ source 7 target 8 label "-" ]
		edge [ source 8 target 9 label "-" ]
	]
]