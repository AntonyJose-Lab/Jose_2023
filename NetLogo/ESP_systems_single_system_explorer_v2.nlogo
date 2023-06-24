directed-link-breed [active-links active-link]  ; These are links that are visible.
directed-link-breed [inactive-links inactive-link] ; These are links that are hidden.

turtles-own [
  val ; an entity/sensor's current 'activity' (e.g., concentration, % conformational change, sequence, etc.) represented visually as a scaled *relative* size of the node. The command 'ask turtles [ show val ]' will show the actual absolute numerical values.
  norm-val ; normalized value of each entity/sensor used for plotting that monitors the change across cell division cycles or equivalently across generations.
  property ; The steps of change by which activity changes if a positive or negative interaction crosses the threshold for change. Properties are characteristic of the entities/sensors themselves and do not change as the system evolves through interactions.
           ; Note that this specifies that every sensor sees the same property of a given entity/sensor.
  inactive-fraction ; Fraction of entity/sensor not available for regulatory interactions at each tick because of processes like protein folding, compartmentalization, diffusion, etc. This is a characteristic of each entity/sensor.
]

links-own [ weight ] ; weight indicates the number of sensors needed to change one unit of property for each entity (characteristic of each 'interaction').

globals [
  ; system-id ; this is the global variable that will be used to run the random generator and hence is a 'system-id' for the ESP system. Moved this to Interface. This system-id variable was obtained for one behavior-space-run of one particular molecule-kind, perturb-kind, and perturb-phase. So, they have all been made input variables.
  ; perturb-phase ; phase of perturbation for robustness tests. Moved this to Interface.
  ; reporter-inactive-fraction ; is typically different from any entity/sensor inactive-fraction. Moved to interface
  ; node-to-change ; The node that should be changed for interacting with a robust system. Moved to interface
  ; from-node ; starting entity/sensor of link to change. Moved to Interface.
  ; to-node ; ending entity/sensor of link to change. Moved to interface.

  ; Additional globals on the interface are cycle-time, molecule-kinds, link-chance, positive-interactions, max-molecules, stasis-level, max-ever-molecules, perturb-kind

  max-val ; maximum value among those of each entity/sensor in the system. Useful for scaling all entity/sensor sizes without getting them to be too big.
  sum-val ; total value of all entities/sensors in the system. This is equal to the total of the number of molecules in the system.
  min-val ; minimum value of all entities/sensors in the system. This is useful for setting the limits of the loss-of-function perturbations.
  stable-gen ; The number of generations a regulatory architecture has been stable.
  stability-gens ; Number of generations of stability for a system to be considered stable.
  test-node ; Entity/Sensor being dynamically picked for purturbation tests.
  perturb-time ; This is the duration (in ticks) for which a system is perturbed each time.
  perturb-freq ; This is the frequency with which the system will be perturbed.
  perturb-value ; The value an entity/sensor is to be held to during a perturbation.
  time-now ; Keeps track of ticks.
  positive-links ; Keeps track of positive regulation in final surviving regulatory architectures.
  negative-links ; Keeps track of negative regulation in final surviving regulatory architectures.
  system-change ; A measure of susceptibility to heritable epigenetic change (number of times the regulatory architecture changes).
  pre-turtles ; Variable that helps with measuring system-change
  post-turtles ; Variable that helps with measuring system-change
  turtle-numbers ; Count of total numbers of turtles in the system in the beginning. Useful for figuring out the 'who' of a newly added turtle half-way through the simulation.
  random-donor-number ; Variable for picking random sensor as donor of link to new turtle/entity/sensor added.
  new-neighbor-nodes ; Neighbors to the new entity/sensor added.
  report-val ; The val of the entity/sensor for which a reporter is added.
  report-property ; The property of the entity/sensor for which a reporter is made
  reporter-link-thickness ; The thickness of the link to reporter when the same regulatory relationships are recreated.
  reporter-neighbor-turtle-whos ; Getting the who numbers of the neighbors for each reporter to recreate the regulation of an entity/sensor for the reporter.
  reporter-link-color ; a variable for positve vs. negative regulation.
  added-reporter ; a variable to keep track of reporters  added.
  in-report-neighbors ; Variable to keep track of the sensors with regulatory input into a reporter.
  reporter-from-neighbor-nodes ; Variable to keep track of the sensors with regulatory input into a reporter.
]

to setup
  clear-all
  reset-timer ; Note this command in the set up procedure clears all global variables. If you want to preserve any value through this process, first assign that variable to another local variable using let before the clear-all and then assign the local variable again to the global after the clear-all.

 ; SET RANDOM SEED FROM INPUT
 random-seed system-id

  ; INITIALIZE PARAMETERS
  set-default-shape turtles "circle"
  set-default-shape links "small-arrow-link"
  set stable-gen 0 ; initializes the variable for storing the number of generations of stability that a regulatory architecture has remained stable.
  set max-val 0 ; initializes the variable for storing the total value summed across all entities/sensors.
  set turtle-numbers 0 ; initializes the number of entities/sensors
  set positive-links 0 ; gets positive links in the final system.
  set negative-links 0 ; gets negative links in the final system.
  set system-change -1 ; initializes for getting number of heritable epigenetic changes detected.
  set pre-turtles 0 ; helps with measuring system-change
  set post-turtles 0 ; helps with measuring system-change
  set perturb-time 5 ; sets the duration in ticks of each perturbation.
  set perturb-freq 100 ; sets the frequency (in ticks = 2x generations) with which the system is perturbed.
  set stability-gens 50 ; sets the numebr of ticks (i.e., 2x generations) that a regulatory architecture has to be stable for before it is considered stable.

  ; COLOR WORLD
  ask patches [ set pcolor white ]

  ; HATCH TURTLES
  while [turtle-numbers != molecule-kinds] [
    let patch-x random-pxcor * 0.9
    let patch-y random-pycor * 0.9
    ask patch patch-x patch-y
    [ if not any? turtles-on patch patch-x patch-y
     [ sprout 1
       set turtle-numbers turtle-numbers + 1
      ]  ; creates random entities throughout the world for representing each molecule kind and the relative number of molecules of that kind.
  ]
  ]
  layout-circle turtles (world-width / 2 - 2) ; spreads the turtles into a circle.

  ; BUILD ENTITY-SENSOR-PROPERTY SYSTEM
  ask turtles [
    let per-turtle-molecules (round (max-molecules / molecule-kinds)) ; scaling the number of molecules to the total molecules considered for the simulation.
    set val 1 + random (per-turtle-molecules - 1) ; start with a random value between 1 and per-turtle-molecules of molecules for each entity/sensor
    set inactive-fraction random-float 1 ; set a random value between 0 and 1 for the fraction of each entity/sensor that remains inactive.
    let neighbor-nodes turtle-set other turtles
      ; setting neighbor-nodes as the set of all other turtles. (this can be restricted to von Neumann neighborhood using neighbors4)
    create-active-links-to neighbor-nodes
      ; create a directed network such that each node has a link-chance % chance of forming a link with one of its neighbors
    [
      set weight random-float 10 ; number of sensors between 1 and 10 needed to change one unit of property for each entity to be depicted as the thickness of the link (characteristic of each 'interaction').
      ifelse random-float 100 > link-chance ; Using link chance to make a set of links visible and thus 'active' and keep the others hidden as 'inactive'.
      [
        set breed inactive-links
        hide-link
      ]
      [ ; Among the active links, using the % of positive and negative interactions to color the positive and negative interactions differently in the model.
        ifelse (random-float 100 < positive-interactions) [
        set color gray
        set thickness (0.5 - weight / 20)
    ]
        [
          set color black
          set thickness (0.5 - weight / 20)
        ]
  ]
  ]
  ]

  set max-val max [ val ] of turtles ; After the system has been setup, make sure that the max-val is not set to zero through the run to avoid division by zero errors later.

  ; set property values and specify if measured entity/sensor or unmeasured entity.
   ask turtles [
    ifelse any? out-active-link-neighbors [
      ifelse any? in-active-link-neighbors [
        set color red
        set property 1 + random 9
      ]
      [
        set color brown
        set property 1 + random 9
      ]; Colors the measured entities/sensors red. Colors measuring receptors (a special subset of sensors) brown. Such measuring receptors cannot last without production from outside the system, but they can be  part of the initial system.
    ]
    [
      ifelse any? in-active-link-neighbors [
      set color blue ; Colors the unmeasured entities blue.
      set property 1 + random 9
      ]
      [
      set color grey ; Colors the unregulated entities grey.
      set val 0 ; Kills all unregulated entities.
      ]
    ]
  ]
  layout-circle turtles (world-width / 2 - 2)

  reset-ticks
end

to go
  set time-now ticks ; Keep track of time.
  ifelse any? turtles [ ; Check in the beginning if there are turtles or not and stop if there are none becasue they have been lost through dilution or regulation.
    if (sum-val > max-ever-molecules) [
      stop ] ; The maximum number of molecules there can ever be - puts an upper bound on the speed with which number of molecules can grow per tick. This simulates the cell death through explosive production of some molecule.
    if (sum-val < stasis-level) [ ; Wait for molecules that may have accrued above a large number (i.e., carrying capacity of the cell) in previous cycle to drop down by dilution/regulation. This is equivalent to the system running out of 'food' molecules.
    ask turtles [
        if (added-reporter != 0) and turtle node-to-report != nobody and turtle (turtle-numbers - 1) != nobody [
          if turtle node-to-report != nobody [ ask turtle node-to-report [ set val val / 2 ] ]
          if turtle (turtle-numbers - 1) != nobody [ ask turtle (turtle-numbers - 1) [set val val / 2 ] ]
        ]
    let donors turtle-set [turtles-here] of in-active-link-neighbors ; Calculation of the input regulators for all entities/sensors, which automatically takes care of all outputs.
    let current-recipient who ; Capture name of current recipient to make all the changes that need to be made before moving to the next recipient.
    let change 0 ; Initializing the variable that will be used to change the value of each entity/sensor at each tick.
    if name? [ ; can use this switch to see labels if needed.
          set label who
          set label-color black
        ]

    ;; VALUE CHANGE CALCULATION
    ;; Overall logic implemented below:
    ;; number of entities increased/decreased at one time = property (characteristic of an entity/sensor).
    ;; number of sensors needed to change one unit of property with each entity = weight of link (characteristic of each 'interaction').
    ;; only a fraction of the sensors are considered inactive to account for processes like protein folding, diffusion, etc.
    ;; that can limit the number of interactions between the sensor and the entity.

    if any? donors [ ; need to calculate property values to add or subtract from all donors to the recipient entity/sensor.
       ask donors [
            ask turtle who [
            let donor-link-weight [ weight ] of link-with turtle current-recipient
            if (donor-link-weight != 0) [  ; Check if a turtle that was already accepted as donor dies as a result of a previous operation, this avoids death of program.
            let donor-link-kind [ color ] of link-with turtle current-recipient
            let donor-val (val - (val * inactive-fraction))
                ifelse (donor-link-kind = gray) [
            set change change + round (donor-val / donor-link-weight)
            ]
            [
                set change change - round (donor-val / donor-link-weight)
            ]
              ]
        ]
      ]
  set val (val + (change * property))
          ; This adds the amount accrued. Note that only values are changing.


          show change
          show property
          show val
  if (val < 0 ) [ set val 0 ]
        ]
  ]
    ]

     ;; CELL DIVISION
     ;; After every two ticks, only about half the number of each entity/sensor, selected using a random number generator, is kept to simulate dilution upon cell division.
     if ((ticks mod cycle-time) = 0) [
    ask turtles [
        let old-val-here round (val)
        let new-val 0
        let counter 0
        while [ counter != old-val-here ] [
        if ((random-float 1) > 0.5) [
      set new-val (new-val + 1)
      ]
          set counter (counter + 1)
        ]
        set val new-val
    if (val < 1) [ set val 0 ] ; This prevents infinitissmally small values from never disappearing.
  ]
            ]

  ;; UPDATE DISPLAYS AND OTHER VARIABLES BEFORE END OF EACH TICK
  update-globals-and-visuals
  ]
  [
    stop ; To end prematurely if all turtles are lost.
  ]
  tick
end

to count-links
  if any? active-links [
    ask active-links [
      ifelse (color = gray) [ set positive-links positive-links + 1 ] [ set negative-links negative-links + 1 ] ]
  ]
end

to update-globals-and-visuals
  set max-val max [ val ] of turtles
  set sum-val sum [ val ] of turtles

  ; Keeping track of the stability of regulatory architectures.
  ifelse any? (turtles with [ val = 0 ])
  [
     set stable-gen 0 ; Reset value to detect next instability followed by stability.
  ]
  [
    set stable-gen stable-gen + 1 ; Update the stable-gen variable every time there has not been a death of a turtle.
    if (stable-gen = stability-gens) [
    spread-nodes ; Spread the nodes to see the stable regulatory architecture better.
    set pre-turtles count turtles ; get occurances of heritable epigenetic change.
    ]

  if max-val != 0 [
      ask turtles with-min [ val ] [ set min-val val ]; Set the lowest non-zero value as the minimum value.
      ask turtles [ set size 0.1 + 2 * sqrt (val / max-val) ]; scale the size to be between 0.1 and 2.0
  ]

    ; ROBUSTNESS TEST: If this test fails, a heritable epigenetic change has been induced.
    if ((perturb-kind != "none") and (ticks > perturb-time)) [ test-robustness-of-nodes ]
    ; changes in entity/sensor values tested only if robustness switch is on. The specific ticks number is needed to ensure that the purturbation does not occur in the first few ticks because the mod function used to set the duration of puturbation (purturb-time) will be satisfied during the initial ticks.
  ]

  ; Kill turtles with zero value.
  ask turtles with [ val = 0 ] [ die ] ; the associated links are automatically killed!

  ; Change in regulatory architecture that include heritable epigenetic change through loss of a node and links.
  if (stable-gen = stability-gens - 1)[
    set post-turtles count turtles ; to get occurances of heritable epigenetic change.
    if ( pre-turtles != post-turtles) [
      set system-change system-change + 1 ]
  ]

  ; Reassign measured and unmeasured entities in the new regulatory architectures.
  ask turtles [
    ifelse any? out-active-link-neighbors [
      ifelse any? in-active-link-neighbors [ set color red ] [set color brown]
      ; colors the measured entities/sensors red. Colors measuring receptors (a special subset of sensors) brown.
    ]
    [
      ifelse any? in-active-link-neighbors [
      set color blue
        ; colors the unmeasured entities blue.
      ]
      [
        set color grey ; colors unlinked entities grey.
      ]
    ]
  ]
end

to spread-nodes
  layout-circle turtles (world-width / 2 - 2)
end

to test-robustness-of-nodes
  ; test every perturb-freq ticks for perturb-time number of ticks.
  if ((stable-gen + perturb-phase) mod perturb-freq = 0) [
    ask one-of turtles [ set test-node who ]; Pick one entity/sensor to purturb for this cycle.
    if (perturb-kind = "gof") [ set perturb-value (max-val * 2.0) ] ; Pick value to increase to for this cycle.
    if (perturb-kind = "lof") [ set perturb-value (min-val * 0.5) ] ; Pick value to decrease to for this cycle.
  ]
  if (((stable-gen + perturb-phase) mod perturb-freq < perturb-time) and (perturb-value != 0)) [ change-node-x ]
  ; Need the second condition 'perturb-value != 0' because its value is initially zero before the first stable ESP systems emerge.
end

to change-node-x
  ifelse any? turtles [
  ask turtles [
     if ((who = test-node) and (val != 0)) [ ; Change values only if the value is not zero.
    ask turtle test-node [
      let old-val val
      if (perturb-kind = "gof") [
            ifelse (val < perturb-value) [ ; if gain of function, hold value high.
            set val perturb-value
            if (max-val != min-val) [ set size 0.1 + 2 * sqrt (val / (max-val - old-val + val)) ]
            set color green
          ]
          [
            set val val
              if (max-val != min-val) [ set size 0.1 + 2 * sqrt (val / (max-val - old-val + val)) ]
            set color green
          ]
            ]
          if (perturb-kind = "lof") [
              ifelse (val > perturb-value) [ ; if loss of function, hold value low.
            set val perturb-value
            if (max-val != min-val) [ set size 0.1 + 2 * sqrt (val / (max-val - old-val + val)) ]
            set color green
          ]
          [
            set val val
              if (max-val != min-val) [ set size 0.1 + 2 * sqrt (val / (max-val - old-val + val)) ]
            set color green
          ]
            ]
    ]
  ]
    ]
  ]
  [
    stop
  ]
end

to remove-a-link
  if any? active-links [
    ask one-of active-links [
      set breed inactive-links
      hide-link
    ]
  ]
end

to remove-link-x-y

  if active-link node-x node-y != nobody [
    ask active-link node-x node-y [
      set breed inactive-links
      hide-link
    ]
  ]
end

to remove-a-node
  ask one-of turtles [ die ]
end

to change-a-node
   ; Kill turtles with zero value.
  ask turtles with [ val = 0 ] [ die ] ; the associated links are automatically killed!
  ifelse turtle node-to-change != nobody [
    ask turtle node-to-change [
    let old-val val
    set val val * -x-  ; -x- indicates fold by which the value of a particular entity/sensor is lowered.
    set color green
      ifelse val != 0  [ set size 0.1 + 2 * sqrt (val / (max-val - old-val + val)) ] [ die ]
  ]
  ]
  [
    stop
  ]
end

to threshold-hold
  ; changes a chosen link to a particular fraction of the current value. Takes values from 0 to 1. This could immitate the introduction of modified enzymes that are more active or less active.
  let target (active-link from-node to-node)
  ifelse (target != nobody) [ ; Check to make sure that the turtles or link have not been killed by past activity.
      ask target [
      set weight weight * -n- ; -n- indicates fold by which the value of the weight of a link is lowered.
      set thickness (0.5 - weight / 20)
    ]
  ]
  [
    stop
  ]
end

to add-a-node
    ; HATCH TURTLES
  set-default-shape turtles "circle"
  set-default-shape links "small-arrow-link"
  create-turtles 1 ; Using the built-in feature of continuous addition of turtle numbers. 'who' goes from 0 to turtle-numbers - 1.
  set turtle-numbers turtle-numbers + 1
  ask turtle (turtle-numbers - 1) [
    let per-turtle-molecules (round (max-molecules / (molecule-kinds + 1)))
    set val 1 + random (per-turtle-molecules - 1) ; start with a random value between 1 and per-turtle-molecules of molecules for the new entity/sensor.
    set inactive-fraction random-float 1 ; set a random value between 0 and 1 for the fraction of each entity/sensor that is not available for interacting with the entity it regulates per tick.
        if name? [ ; can use this to see labels if needed. This is useful for holding chosen entities/sensors at chosen values.
          set label who
          set label-color black
        ]
    set new-neighbor-nodes turtle-set other turtles
    create-active-links-to new-neighbor-nodes
      ; create a directed network such that each node has a link-chance % chance of forming a link with one of its neighbors
    [
      set weight random-float 10 ; setting number of sensors needed to change one unit of property for each entity (characteristic of each 'interaction').
      ifelse random-float 100 > link-chance ; Using link chance to make a set of links visible and thus 'active' and keeping the others 'inactive'.
      [
        set breed inactive-links
        hide-link
      ]
      [ ; Among the active links, using the % of positive and negative interactions to differently color the positive and negative interactions in the model.
        ifelse (random-float 100 < positive-interactions) [
        set color gray
        set thickness (0.5 - weight / 20)
    ]
        [
          set color black
          set thickness (0.5 - weight / 20)
        ]
  ]
  ]
  ]

  ; Kill turtles with zero value.

  ask turtles with [ val = 0 ] [ die ] ; the associated links are automatically killed.

  ; pick a random entity/sensor to interact with the new entity.
  if (count turtles != 0) [
   ask one-of new-neighbor-nodes [ create-active-link-to turtle (turtle-numbers - 1)
    set random-donor-number who
    ]
  let recipient-number turtle-numbers - 1
    ask active-link random-donor-number recipient-number [
       set weight random-float 10
      ifelse (random-float 100 < positive-interactions) [
        set color gray
         set thickness (0.5 - weight / 20)
        set thickness (0.5 - weight / 20)
    ]
        [
          set color black
        ]
  ]
  ]

  set max-val max [ val ] of turtles ; After the system has been setup, make sure that the max-val is not set to zero throughout the run.

  ; set property values and specify if the new turtle is a measured entity/sensor or unmeasured entity.
   ask turtles [
    ifelse any? out-active-link-neighbors [
      ifelse any? in-active-link-neighbors [
        set color red
        set property 1 + random 9
      ]
      [
        set color brown
        set property 1 + random 9
      ]; Colors the measured entities/sensors red. Colors measuring receptors (a special subset of sensors) brown.
    ]
    [
      ifelse any? in-active-link-neighbors [
      set color blue ; Colors the unmeasured entities blue.
      set property 1 + random 9
      ]
      [
      set color grey ; Colors the unmeasured entities blue.
        set val 0 ; killing off all unconnected molecules.
      ]
    ]
  ]
    layout-circle turtles (world-width / 2 - 2) ; spreads the turtles into a circle.
end

to add-a-reporter ; Note all ideal reporters are entities and not sensors because they interact with the system ideally without disrupting anything within the system.
  ; Kill turtles with zero value.
  ask turtles with [ val = 0 ] [ die ] ; the associated links are automatically killed.
  set added-reporter added-reporter + 1
if turtle node-to-report != nobody [
   set-default-shape turtles "circle"
  set-default-shape links "small-arrow-link"
  create-turtles 1 ; Using the built-in feature of continuous addition of turtle numbers to discover the who number of reporter (goes from 0 to turtle-numbers - 1).
 set turtle-numbers turtle-numbers + 1
    ask turtle node-to-report [
    set report-val val
    set report-property property
    set in-report-neighbors in-active-link-neighbors
  ]
  ask turtle (turtle-numbers - 1) [
    set val report-val
    set property report-property
    set inactive-fraction reporter-inactive-fraction
    set color grey
    set reporter-from-neighbor-nodes turtle-set in-report-neighbors
    create-active-links-from reporter-from-neighbor-nodes
      [
      ifelse perfect? [
        ask turtle-set reporter-from-neighbor-nodes [
          set reporter-neighbor-turtle-whos who
          ask active-link reporter-neighbor-turtle-whos node-to-report [
              set reporter-link-thickness thickness
              set reporter-link-color color ]
          ask active-link reporter-neighbor-turtle-whos (turtle-numbers - 1) [
              set thickness reporter-link-thickness
              set weight (0.5 - reporter-link-thickness) * 20
              set color reporter-link-color ]
        ]
      ]
      [
      set weight random-float 10 ; setting number of sensors needed to change one unit of property for the reporter (characteristic of each 'interaction').
       ; Among the active links, using the % of positive and negative interactions to differently color the positive and negative interactions in the model.
        ifelse (random-float 100 < positive-interactions) [
        set color gray
        set thickness (0.5 - weight / 20)
    ]
        [
          set color black
          set thickness (0.5 - weight / 20)
        ]
  ]
  ]
  ]
  ]
  layout-circle turtles (world-width / 2 - 2) ; spreads the turtles into a circle.
end


;; Copyright 2020 Antony Jose for this model and 2008 Uri Wilensky for NetLogo.
;; See Info tab for full copyright and license.
@#$#@#$#@
GRAPHICS-WINDOW
460
10
888
439
-1
-1
20.0
1
20
1
1
1
0
0
0
1
-10
10
-10
10
1
1
1
ticks
30.0

SLIDER
140
45
295
78
link-chance
link-chance
0
100
50.0
1
1
%
HORIZONTAL

PLOT
890
10
1390
440
phenotype
time
relative [Entity] or [Sensor]
0.0
500.0
0.0
1.0
true
true
"ask turtles [\n  create-temporary-plot-pen (word who) ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color one-of base-colors ; this assigns each pen a random one of the base colors in NetLogo [5 15 25 35 45 55 65 75 85 95 105 115 125 135]\n  ]" " ask turtles [\n  set norm-val val / max-val\n if who < (molecule-kinds - 1) [\n set-current-plot-pen (word who)\n ]\n plot norm-val\n]\n"
PENS

SLIDER
295
45
460
78
positive-interactions
positive-interactions
0
100
50.0
1
1
%
HORIZONTAL

SLIDER
140
80
295
113
max-molecules
max-molecules
0
1000
500.0
1
1
NIL
HORIZONTAL

MONITOR
625
440
720
485
total molecules
sum-val
0
1
11

SWITCH
365
150
455
183
name?
name?
0
1
-1000

SLIDER
295
10
460
43
cycle-time
cycle-time
1
100
2.0
1
1
NIL
HORIZONTAL

SLIDER
295
80
460
113
stasis-level
stasis-level
max-molecules
100 * max-molecules
5000.0
1
1
NIL
HORIZONTAL

SLIDER
140
115
460
148
max-ever-molecules
max-ever-molecules
stasis-level
100 * stasis-level
500000.0
1
1
NIL
HORIZONTAL

MONITOR
720
440
830
485
stable generations
stable-gen / 2
17
1
11

MONITOR
1175
440
1270
485
perturbed node
test-node
17
1
11

MONITOR
550
440
627
485
generation
ticks / 2
17
1
11

MONITOR
915
440
1012
485
NIL
perturb-value
17
1
11

BUTTON
195
150
250
183
NIL
stop
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
140
150
195
183
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
305
150
360
183
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1010
440
1095
485
perturb-time
perturb-time
17
1
11

MONITOR
1095
440
1175
485
perturb-freq
perturb-freq
17
1
11

MONITOR
830
440
917
485
stability gen
stability-gens / 2
17
1
11

CHOOSER
365
190
460
235
perturb-kind
perturb-kind
"none" "lof" "gof"
0

INPUTBOX
140
185
210
245
system-id
11.0
1
0
Number

INPUTBOX
210
185
290
245
molecule-kinds
2.0
1
0
Number

INPUTBOX
290
185
365
245
perturb-phase
4.0
1
0
Number

BUTTON
240
450
340
483
NIL
remove-a-node
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
340
450
435
483
NIL
remove-a-link
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
160
250
287
283
change-a-node
change-a-node
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
410
285
460
345
-n-
0.5
1
0
Number

BUTTON
315
250
437
283
link-hold
threshold-hold
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
290
285
355
345
from-node
0.0
1
0
Number

INPUTBOX
355
285
410
345
to-node
3.0
1
0
Number

BUTTON
160
450
240
483
NIL
add-a-node
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
160
10
280
43
NIL
clear-all
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
250
150
305
183
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
240
285
290
345
-x-
0.0
1
0
Number

INPUTBOX
140
285
240
345
node-to-change
2.0
1
0
Number

BUTTON
140
350
240
383
NIL
remove-link-x-y
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
140
385
190
445
node-x
1.0
1
0
Number

INPUTBOX
190
385
245
445
node-y
2.0
1
0
Number

BUTTON
260
350
355
383
NIL
add-a-reporter
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
260
385
355
445
node-to-report
2.0
1
0
Number

INPUTBOX
355
385
455
445
reporter-inactive-fraction
0.5
1
0
Number

SWITCH
355
350
455
383
perfect?
perfect?
0
1
-1000

@#$#@#$#@
## ABOUT THE MODEL

This Netlogo model simulates entity-sensor-property systems (Jose, AM, J. R. Soc. Interface, 2020) and their evolution across generations to explore regulatory architectures.

Particular regulatory architectures identified as being stable for 250 generations using the related ESP_systems_explorer_v1 (78285 out of 225000 tested) can be explored in detail using this agent-based model.

The initial features of a model are recreated using the system-id, which is the behaviorspace-run-number used to discover stable regulatory arhitectures. The other parameters needed for faithful recreation of the previous system behavior are:
- cycle-time = 2; This is the timing in ticks for each generation.
- link-chance = 50%; This gives the probability that any two entities/sensors will interact as part of the system.
- positive-interactions = 50%; This gives the probability that an interaction is positive.
- max-molecules = 500; This is the maximum number of total molecules at the start of the simulation.
- stasis-level = 5000; This is the number of molecules that arrests growth until molecules get diluted upon cell division.
- max-ever-molecules = 500000; This is the most number of molecules of all kinds put together that can be within any system at any time. Simulates living systems existing in finite space.
- molecule-kinds; This is the number of entities/sensors that are part of the regulatory architecture.
- perturb-kind (none, lof, or gof); This is a chooser for perturbing a random entity/sensor every ~50 generations for 2.5 generations by increasing (gof) or decreasing (lof) its value (i.e., concentration/number) by two fold of the maximal or minimal values, respectively, of all the entities/sensors.
- perturb-phase; This is the precise timing for starting the periodic perturbations (0 = starting @ tick 100; 1 = starting @ tick 101; 2 = starting @ tick 102; 3 = starting @ tick 103; 4 = starting @ tick 104)

Once the ESP system has been recreated, the following aspects can be changed:
1. The value (number of molecules) of an entity using the 'change-a-node' button. Any entity can be chosen using its 'who' number (easily discoverable by toggling 'name' switch on or off). Its value can be changed using a multiplier, -x- (e.g., -x- = 0.5 will halve the value and -x- = 0 will remove the entity).
2. The strength of regulatory input using the 'link-hold' button. The weight of any link, specified using 'from-node' and 'to-node', can be changed using a multiplier -n- ( e.g., -n- = 0.5 will halve the weight. Note, setting -n- to 0 does not remove the link. To remove a particular link, the remove-link-x-y button can be used.)
3. A new entity with random value, property, inactive-fraction, and regulatory links can be introduced using the 'add-a-node' button.
4. A random entity an be removed using the 'remove-a-node' button.
5. A random link can be removed using the 'remove-a-link' button.
6. A reporter for any entity can be added using the 'add-a-reporter' button. The 'node-to-report' and the 'reporter-inactive-fraction' can be specified. The 'perfect?' on/off switch can be toggled to make a perfect or imperfect reporter. A perfect reporter of an entity receives the same regulatory input as the entity/sensor of interest does. An imperfect reporter of an entity recieves input from the same sensors as the entity/sensor of interest does. but the direction and strength of the input can vary. Regulatory outputs of the entity are not recreated for any reporter.

Monitors reporting generation number, total molecules, stable generations since last instability, the number of generations of stability for considering a regulatory architecture stable (stability gen), the value of the perturbed entity (perturb-value), the duration of each perturbation (perturb-time), the frequency of the perturbations (perturb-freq) and the identity of the perturbed entity (perturbed node) are included.

For simulating change over time, this model uses a combination of deterministic and stochastic functions. The values of each entity/sensor changes at each tick using a deterministic equation:

value @ t+1 = value @ t + sum of inputs from all sensors.

The change in value contributed by each sensor for a given entity = round((unit of entity)*(value of sensor)*(1 - inactive-fraction of sensor) / (weight of regulatory link)). In other words, change = round(k*(value of sensor)), where k is a different constant for each sensor of each entity and round indicates rounding to the nearest integer. For positive regulators (link color grey), this change in value is added and for negative regulators (link color black), it is subtracted.

The order of operation on the entities/sensors varies with every tick based on the random-seed. After every two ticks, only about half the number of each entity/sensor, selected using a random number generator, is kept to simulate dilution and random partition upon cell division.

The value of each entity/sensor is plotted relative to the most abundant entity/sensor. This profile at each time point can be considered as the 'phenotype' of the system. These scaled values are also used to depict each entity/sensor in the regulatory architecture at each time point.

## COPYRIGHT AND LICENSE


Copyright 2008 Uri Wilensky for NetLogo and 2020 Antony Jose for this simulation of ESP systems.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

small-arrow-link
1.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
4
Line -1184463 true 150 150 105 195
Line -1184463 true 150 150 195 195
@#$#@#$#@
1
@#$#@#$#@
