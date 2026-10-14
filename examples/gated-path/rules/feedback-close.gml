rule [
  ruleID "one-closing-match-per-path-vertex"
  context [
    node [ id 0 label "U" ] node [ id 1 label "W" ] node [ id 2 label "V" ]
    node [ id 3 label "A" ]
    edge [ source 0 target 1 label "-" ] edge [ source 1 target 2 label "-" ]
  ]
  right [ edge [ source 0 target 2 label "-" ] ]
  constrainNoEdge [ source 0 target 2 ]
]
