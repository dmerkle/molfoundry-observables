rule [
	ruleID "tethered-context-guard-ablated-photocycloaddition"
	left [
		edge [ source 0 target 1 label "=" ]
		edge [ source 2 target 3 label "=" ]
	]
	context [
		node [ id 0 label "C" ]
		node [ id 1 label "C" ]
		node [ id 2 label "C" ]
		node [ id 3 label "C" ]
		node [ id 4 label "C" ]
		node [ id 5 label "C" ]
		edge [ source 0 target 4 label "-" ]
		edge [ source 2 target 5 label "-" ]
		edge [ source 4 target 5 label "-" ]
	]
	right [
		edge [ source 0 target 1 label "-" ]
		edge [ source 0 target 3 label "-" ]
		edge [ source 1 target 2 label "-" ]
		edge [ source 2 target 3 label "-" ]
	]
]