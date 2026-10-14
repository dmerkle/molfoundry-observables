rule [
  ruleID "close-controller"
  context [
    node [ id 0 label "U" ] node [ id 1 label "W" ] node [ id 2 label "V" ]
    edge [ source 0 target 1 label "-" ] edge [ source 1 target 2 label "-" ]
  ]
  right [ edge [ source 0 target 2 label "-" ] ]
  constrainNoEdge [ source 0 target 2 ]
]
