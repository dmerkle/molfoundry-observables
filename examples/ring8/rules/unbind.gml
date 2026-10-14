rule [
            ruleID "unbind"
            left [ edge [ source 1 target 2 label "bound" ] ]
            context [ 
        node [ id 0 label "C" ]
        node [ id 1 label "P" ]
        node [ id 2 label "L" ]
        edge [ source 0 target 1 label "cp" ]
        edge [ source 0 target 2 label "cl" ]
     ]
        ]
