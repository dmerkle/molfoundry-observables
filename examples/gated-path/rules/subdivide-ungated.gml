rule [
  ruleID "subdivide-any-path-edge-without-gate"
  left [ edge [ source 0 target 1 label "-" ] ]
  context [
    node [ id 0 label "A" ] node [ id 1 label "A" ]
    node [ id 3 label "U" ] node [ id 4 label "W" ] node [ id 5 label "V" ]
    edge [ source 3 target 4 label "-" ] edge [ source 4 target 5 label "-" ]
  ]
  right [
    node [ id 2 label "A" ]
    edge [ source 0 target 2 label "-" ] edge [ source 2 target 1 label "-" ]
  ]
]
