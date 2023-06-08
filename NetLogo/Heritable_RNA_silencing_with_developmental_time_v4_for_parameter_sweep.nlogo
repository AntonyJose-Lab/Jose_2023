directed-link-breed [active-links active-link]  ; These are links that are visible.
; directed-link-breed [inactive-links inactive-link] ; These are links that are hidden. These can be used for further development of the model in response to induced changes.

turtles-own [
  val ; an entity/sensor's current 'activity' (e.g., concentration, % conformational change, sequence, etc.) represented visually as a scaled *relative* size of the node. The command 'ask turtles [ show val ]' will show the actual absolute numerical values.
  norm-val ; normalized value of each entity/sensor used for plotting that monitors the change across cell division cycles or equivalently across generations.
  property ; The steps of change by which activity changes if a positive or negative interaction crosses the threshold for change. Properties are characteristic of the entities/sensors themselves and do not change as the system evolves through interactions.
           ; Note that this specifies that every sensor sees the same property of a given entity/sensor.
  inactive-fraction ; Fraction of entity/sensor not available for regulatory interactions at each tick because of processes like protein folding, compartmentalization, diffusion, etc. This is a characteristic of each entity/sensor.
  turnover ; this is a term (0.0 to 1.0) for the fraction remaining after one tick because of turnover processes in the system that are not explicitly modeled.
  last-delay-time-values ; a convenient variable for remembering past n values (e.g., 10 needed for 5 germline zones and 2 generations) and using them later for calculations.
]

links-own [
  weight ; weight indicates the number of sensors needed to change one unit of property for each entity (characteristic of each 'interaction').
  delay ; delay indicates the time after increase in a sensor that the effect on the regulated entity/sensor can be seen.
]

globals [
  ;gene-is-sensor ; this variable holds the chance that the gene acts as a sensor instead of just an entity for the regulatory architecture. Slider on interface.
  ;cyc1-vs-cyc ; this variable stores the bias in regulatory output from each arbitrary component of core cycling stores, cyc1 and cyc2. Slider on interface.
  ;Additional globals on the interface are cycle-time, molecule-kinds, max-molecules, stasis-level, max-ever-molecules.

  delay-time ; the maximal delay with which a molecule can have an effect. This can be increased beyond a generation to allow for intergenerational effects and in fact can be increased indefinitely to allow for multigenerational effects.
  max-val ; maximum value among those of each entity/sensor in the system. Useful for scaling all entity/sensor sizes without getting them to be too big.
  sum-val ; total value of all entities/sensors in the system. This is equal to the total of the number of molecules in the system.
  min-val ; minimum value of all entities/sensors in the system. This is useful for setting the limits of the loss-of-function perturbations.
  generation-time ; the duration in ticks of one generation. Custom information can be used for this variable to indicate the progress of one cell through divisions and delay
  divide-vs-grow ; a list of 0s and 1s with 0 signifying growth and 1 signifying division. This is probably the simplest way to encode a germ cell life cycle. [e.g., 0 1 0 1 0 0 0 0 1 0 1 etc.]
                 ; for E. coli, this is just [ 0 1 ], with each tick = 10 min., but for  C. elegans, it is [ 0 1 0 1 0 1 0 1 0 1 ... including many 0 0 0 0 strings where there are no cell divisions.]
                 ; By looking at the literature and assembling a table of developmental time from fertilization to fertilization, I have made an approximate string for 365 ticks,
                 ; where every tick = 15 min and encodes cell divisions along the germ lineage from one generation ot the next. There are a total of 14 cell divisions to go from one
                 ; generation to the next. THE NEXT STEP IS TO INCLUDE THE IMPACT OF SOMATIC CELLS BY KEEPING TRACK OF EVERY CELL AT EVERY BRANCH POINT...THIS IS A PATH TOWARDS
                 ; DYNAMICLALY SIMULATING THE ENTIRE C ELEGANS CELL CODE GENERATION AFTER GENERATION!
]

to setup
  clear-all ; Note this command in the set up procedure clears all global variables. If you want to preserve any value through this process, first assign that variable to another local variable using let before the clear-all and then assign the local variable again to the global after the clear-all.

 ; SET RANDOM SEED FROM INPUT
  random-seed behaviorspace-run-number

  ; INITIALIZE PARAMETERS
  set-default-shape turtles "circle"
  set-default-shape links "small-arrow-link"
  set divide-vs-grow [ 1 0 1 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]
  ; the above long list encodes the minimal duration of one generation from zygote to zygote in C. elegans as estimated based on papers. This is a rough estimate of events in 15 min. intervals.
  set generation-time length divide-vs-grow ; This is time in ticks. For C. elegans, true time is 15 min. * generation-time.
  set delay-time (generation-time * ancestral-effect-generations) - 1; this is the maximal delayed action allowed in the simulation and accounts for maternal effects from each region of the germline on each region of the germline.

  ; COLOR WORLD
  ask patches [ set pcolor white ]
  ; HATCH all turtles, position them, and give them their basic attributes.

  ; UNKNOWN CYCLING STORES
  ask patch -8 4 [; this is prt1 of cycling stores that are maintaining form and function in each generation.
    sprout 1 [
      set color yellow
      ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;  ask patch -8 6 [ set plabel "cyc1" set plabel-color grey ]

      ask patch -8 -4 [; this is prt2 of the cycling stores that are maintaining form and function in each generation.
    sprout 1 [
      set color yellow
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
;    ask patch -8 -6 [ set plabel "cyc2" set plabel-color grey ]


    ; RELEVANT Gene Products
  ask patch 8 8 [ ; this is the gene for which a transcriptional reporter is being made.
    sprout 1 [
      set color black
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
;  ask patch 8 9.5 [ set plabel "gene" set plabel-color black ]

  ask patch 8 -8 [; this is a reporter for the gene of interest.
    sprout 1 [
      set color magenta
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;     ask patch 8 -9.5 [ set plabel "reporter" set plabel-color black ]

    ; RELEVANT REGULATORS cis regulators
  ask patch 0 2 [; this is positive small RNAs that increase expression of the gene.
    sprout 1 [
      set color cyan
            ; set label who
      set label-color cyan
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;   ask patch 0 1 [ set plabel "+s" set plabel-color black ]

    ask patch 0 -2 [; this is negative small RNAs that decrease expression of the gene.
    sprout 1 [
      set color blue
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;   ask patch 0 -3 [ set plabel "-s" set plabel-color black ]


  ; RELEVANT REGULATORS small RNAs

          ask patch 0 8 [; this is +ve small RNAs that regulate reporter mRNA.
    sprout 1 [
      set color lime
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
;   ask patch 0 7 [ set plabel "+p" set plabel-color black ]

          ask patch 0 5 [; this is -ve small RNAs that regulate reporter mRNA.
    sprout 1 [
      set color brown
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;  ask patch 0 4 [ set plabel "-p" set plabel-color black ]

          ask patch 0 -5 [; this is +ve small RNAs that regulate gene mRNA.
    sprout 1 [
      set color violet
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;  ask patch 0 -6 [ set plabel "+p" set plabel-color black ]


        ask patch 0 -8 [; this is -ve small RNAs that regulate gene mRNA.
    sprout 1 [
      set color gray
            ; set label who
      set label-color black
      set property 1 + random 9
      set val 1 + random 99
      set inactive-fraction random-float 1.0
      set turnover random-float 1.0
      set last-delay-time-values []
      set last-delay-time-values fput val last-delay-time-values
  ] ]
 ;  ask patch 0 -9 [ set plabel "-p" set plabel-color black ]

  ; Scale turtle sizes for displaying their relative values.
  ask turtles [ set max-val max [ val ] of turtles
  set size 0.1 + 2 * sqrt (val / max-val)
  ]

  ; BUILD ENTITY-SENSOR-PROPERTY SYSTEM
  if random-float 1.0 > gene-is-sensor [ ; a gene expressed in the germline can be either a sensor or an entity and be part of the cycling stores of heritable information
    ifelse random-float 1.0 > cyc1-vs-cyc2 [
    ask turtle 2 [ create-active-link-to turtle 0 ] ; this link makes a gene a sensor with regulatory input to cyc1
  ask active-link 2 0 [
    set color one-of [ grey black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
  ]]
    [
    ask turtle 2 [ create-active-link-to turtle 1 ] ; this link makes a gene a sensor with regulatory input to cyc2
  ask active-link 2 1 [
    set color one-of [ grey black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
  ]]
  ]

ifelse random-float 1.0 > cyc1-vs-cyc2 [ ; there can be general bias of regulatory inputs from the arbitrarily divided cyc1 vs. cyc2 to the different entities/sensors of interest.
  ask turtle 0 [ create-active-link-to turtle 4 ] ; from cyc1 or cyc2 to +s
  ask active-link 0 4 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]
  [
      ask turtle 1 [ create-active-link-to turtle 4 ] ; from cyc1 or cyc2 to +s
  ask active-link 1 4 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
    ]
  ]

ifelse random-float 1.0 > cyc1-vs-cyc2 [ ; there can be general bias of regulatory inputs from the arbitrarily divided cyc1 vs. cyc2 to the different entities/sensors of interest.
  ask turtle 0 [ create-active-link-to turtle 5 ] ; from cyc2 to -s
  ask active-link 0 5 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]
    [
      ask turtle 1 [ create-active-link-to turtle 5 ] ; from cyc2 to -s
  ask active-link 1 5 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
    ]

  ask turtle 4 [ create-active-link-to turtle 2 ] ; from +s to gene
  ask active-link 4 2 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
   ;; Reporter interactions with regulatory machinery that target the gene through UTR elements. Note that random-float 2.0 is used to simulate possibility of weaker (>1) or stronger (<1) regulation than that for the corresponding gene.
    ask turtle 4 [ create-active-link-to turtle 3 ]
  ask active-link 4 3 [ ; from +s to reporter
    set color [ color ] of active-link 4 2
    set weight round(([ weight ] of active-link 4 2) * (random-float 2.0))
    set thickness (0.5 - weight / 20)
    set delay round(([ delay ] of active-link 4 2) * (random-float 2.0))
]

  ask turtle 5 [ create-active-link-to turtle 2 ] ; from -s to gene
  ask active-link 5 2 [
    set color one-of [ black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ;; Reporter interactions with regulatory machinery that target the gene through UTR elements. Note that random-float 2.0 is used to simulate possibility of weaker (>1) or stronger (<1) regulation than that for the corresponding gene.
      ask turtle 5 [ create-active-link-to turtle 3 ]
    ask active-link 5 3 [ ; from -s to reporter
    set color [ color ] of active-link 5 2
    set weight round(([ weight ] of active-link 5 2) * (random-float 2.0))
    set thickness (0.5 - weight / 20)
    set delay round(([ delay ] of active-link 5 2) * (random-float 2.0))
]

  ;; INTERACTIONS THROUGH SMALL RNAs
  ;; interactions for + small RNAs and gene
ifelse random-float 1.0 > cyc1-vs-cyc2 [ ; bias of regulatory inputs from cyc1 vs. cyc2
  ask turtle 0 [ create-active-link-to turtle 6 ]
  ask active-link 0 6 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 6 [ create-active-link-to turtle 2 ]
  ask active-link 6 2 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]
  [
      ask turtle 1 [ create-active-link-to turtle 6 ]
  ask active-link 1 6 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 6 [ create-active-link-to turtle 2 ]
  ask active-link 6 2 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]

  ;; interactions for - small RNAs and gene
  ifelse random-float 1.0 > cyc1-vs-cyc2 [ ; bias of regulatory inputs from cyc1 vs. cyc2
  ask turtle 0 [ create-active-link-to turtle 7 ]
  ask active-link 0 7 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 7 [ create-active-link-to turtle 2 ]
  ask active-link 7 2 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]
  [
      ask turtle 1 [ create-active-link-to turtle 7 ]
  ask active-link 1 7 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 7 [ create-active-link-to turtle 2 ]
  ask active-link 7 2 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]

  ;; interactions for + small RNAs and reporter
ifelse random-float 1.0 > cyc1-vs-cyc2 [ ; bias of regulatory inputs from cyc1 vs. cyc2
  ask turtle 0 [ create-active-link-to turtle 8 ]
  ask active-link 0 8 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 8 [ create-active-link-to turtle 3 ]
  ask active-link 8 3 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]
  [
      ask turtle 1 [ create-active-link-to turtle 8 ]
  ask active-link 1 8 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 8 [ create-active-link-to turtle 3 ]
  ask active-link 8 3 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]

  ;; interactions for - small RNAs and reporter
  ifelse random-float 1.0 > cyc1-vs-cyc2 [ ; bias of regulatory inputs from cyc1 vs. cyc2
  ask turtle 0 [ create-active-link-to turtle 9 ]
  ask active-link 0 9 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 9 [ create-active-link-to turtle 3 ]
  ask active-link 9 3 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]
  [
      ask turtle 1 [ create-active-link-to turtle 9 ]
  ask active-link 1 9 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
     ask turtle 9 [ create-active-link-to turtle 3 ]
  ask active-link 9 3 [ ; from cyc2 to +small RNA targeting gene
    set color one-of [ black ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
  ]

  ;; Create the mutual postive regulation of cycling stores.
  ask turtle 0 [ create-active-link-to turtle 1 ] ; from cyc2 to cyc1
  ask active-link 0 1 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]
    ask turtle 1 [ create-active-link-to turtle 0 ] ; from cyc1 to cyc2
  ask active-link 1 0 [
    set color one-of [ grey ]
    set weight random-float 10
    set thickness (0.5 - weight / 20)
    set delay random delay-time
]

;; Ask links to show delay time in hours. Including the code below that is just supposed to label the links changes the behavior of systems such that different systems can become stable!
;; This is because the 'ask' command uses the random number generator, which changes subsequent particular random numbers used. But this is a convenient bit of code to know the delays of each link.
;  ask active-links [
;    set label-color blue
;    set label round(delay / 4)
;  ]
  reset-ticks
end

to go
  ifelse any? turtles [ ; Check in the beginning if there are turtles or not and stop if there are none becasue they have been lost through dilution or regulation.
    if (sum-val > max-ever-molecules) [ stop ] ; The maximum number of molecules there can ever be - puts an upper bound on the speed with which number of molecules can grow per tick. This simulates the cell death through explosive production of some molecules.
    ifelse (sum-val < stasis-level) [ ; Wait for molecules that may have accrued above a large number (i.e., carrying capacity of the cell) in previous cycle to drop down by dilution/regulation. This is equivalent to the system running out of 'food' molecules.
    ask turtles [
    let donors turtle-set [turtles-here] of in-active-link-neighbors ; Calculation of the input regulators for all entities/sensors, which automatically takes care of all outputs.
    let current-recipient who ; Capture name of current recipient to make all the changes that need to be made before moving to the next recipient.
    let change 0 ; Initializing the variable that will be used to change the value of each entity/sensor at each tick.

    ;; VALUE CHANGE CALCULATION
    ;; Overall logic implemented below:
    ;; number of each entity (cyc1, cyc2, gene, +s, -s, + small RNAs, - small RNAs or reporter) increased/decreased at one time = property (characteristic of an entity/sensor).
    ;; number of sensors needed to change one unit of property for each entity = weight of link (characteristic of each 'interaction').
    ;; only a fraction of the sensors are considered inactive to account for processes like protein folding, diffusion, etc.
    ;; that can limit the number of interactions between the sensor and the entity.
    ;; timing of change is after a delay of upto 2xgeneratoion-time ticks to account for duration of intermediate processes.

    if any? donors [ ; need to calculate property values to add or subtract from all donors to the recipient entity/sensor.
       ask donors [
            ask turtle who [
            let donor-link-weight [ weight ] of link-with turtle current-recipient
            let donor-delay [ delay ] of link-with turtle current-recipient
            if (donor-link-weight != 0) [  ; Check if a turtle that was already accepted as donor dies as a result of a previous operation, this avoids death of program.
            let donor-link-kind [ color ] of link-with turtle current-recipient

            ifelse length last-delay-time-values > donor-delay [
                 ;; Calculate change with delay
                    let i donor-delay
                    let delayed-val item i last-delay-time-values
                    let donor-val (val + delayed-val * ( 1 - inactive-fraction))

                  ifelse (donor-link-kind = grey) [
            set change change + round (donor-val / donor-link-weight)
            ]
            [
            set change change - round (donor-val / donor-link-weight)
            ]
              ]
                [
                  ;; Calculate change with whatever delay is possible. This only occurs in the beginning of the simulation.
                  let i (length last-delay-time-values) - 1
                    let delayed-val item i last-delay-time-values
                    let donor-val (val + delayed-val * ( 1 - inactive-fraction))

            ifelse (donor-link-kind = grey) [
            set change change + round (donor-val / donor-link-weight)
            ]
            [
            set change change - round (donor-val / donor-link-weight)
            ]
                ]
        ]
      ]
        ]
  set val (val + (change * property)) * turnover
          ; This adds the amount accrued. Note that only values are changing.
  if val < 0 [ set val 0 ]
    ]
    ]
    ]
    [
      ask turtles [ set val val * turnover]
    ]
  update-memory-of-past-values
  dilute-across-cell-divisions
  update-globals-and-visuals
  ]
  [
    stop ; To end prematurely if all turtles are lost.
  ]
  tick
end

to dilute-across-cell-divisions

     let time-now ticks mod generation-time
     if item time-now divide-vs-grow = 1 [; the procedure for dividing cells and continuing along one lineage. After every cycle-time ticks, only about 1/2 the number of each entity/sensor selected using a random number generator, is kept to simulate dilution across cell divisions.
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
end

to update-memory-of-past-values
  ask turtles [
      set last-delay-time-values fput val last-delay-time-values
      if length last-delay-time-values > delay-time [ set last-delay-time-values but-last last-delay-time-values ]
    ]
end

to update-globals-and-visuals
  set max-val max [ val ] of turtles
  set sum-val sum [ val ] of turtles

  if max-val != 0 [
      ask turtles with-min [ val ] [ set min-val val ]; Set the lowest non-zero value as the minimum value.
      ask turtles [ set size 0.1 + 2 * sqrt (val / max-val) ]; scale the size to be between 0.1 and 2.0
  ]

  ; Kill turtles with zero value.
  ask turtles with [ val = 0 ] [ die ] ; the associated links are automatically killed!

end

;;; Copyright 2020 Antony Jose for this model and 2008 Uri Wilensky for NetLogo.
;;; See Info tab for full copyright and license.
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

PLOT
890
10
1390
440
phenotype
time
relative number
0.0
3650.0
0.0
1.1
false
true
"ask turtle 2 [\n  create-temporary-plot-pen \"gene\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color black \n  ]\n  ask turtle 3 [\n  create-temporary-plot-pen \"reporter\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color magenta \n  ]\n  ask turtle 4 [\n  create-temporary-plot-pen \"+s\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color cyan \n  ]\n  ask turtle 5 [\n  create-temporary-plot-pen \"-s\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color blue \n  ]\n    ask turtle 6 [\n  create-temporary-plot-pen \"+p gene\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color lime \n  ]\n  ask turtle 7 [\n  create-temporary-plot-pen \"-p gene\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color brown \n  ]\n    ask turtle 8 [\n  create-temporary-plot-pen \"+p reporter\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color violet \n  ]\n      ask turtle 9 [\n  create-temporary-plot-pen \"-p reporter\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color gray \n  ]" " ask turtles [\n  set norm-val val / max-val\n] \nif turtle 2 != nobody [ ask turtle 2 [\n  create-temporary-plot-pen \"gene\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color black \n    plot norm-val\n  ]]\n  if turtle 3 != nobody [ ask turtle 3 [\n  create-temporary-plot-pen \"reporter\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color magenta \n    plot norm-val\n  ]]\n  if turtle 4 != nobody [ ask turtle 4 [\n  create-temporary-plot-pen \"+s\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color cyan  \n    plot norm-val\n  ]]\n  if turtle 5 != nobody [ ask turtle 5 [\n  create-temporary-plot-pen \"-s\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color blue \n    plot norm-val\n  ]]\n  if turtle 6 != nobody [ ask turtle 6 [\n  create-temporary-plot-pen \"+p gene\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color lime \n    plot norm-val\n  ]]\n  if turtle 7 != nobody [ ask turtle 7 [\n  create-temporary-plot-pen \"-p gene\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color brown \n    plot norm-val\n  ]]  \n  if turtle 8 != nobody [ ask turtle 8 [\n  create-temporary-plot-pen \"+p reporter\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color violet\n    plot norm-val\n  ]] \n    if turtle 9 != nobody [ ask turtle 9 [\n  create-temporary-plot-pen \"-p reporter\" ; this creates a temporary plot pen for all turtles\n  set-plot-pen-color gray\n    plot norm-val\n  ]] "
PENS

SLIDER
105
115
260
148
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
255
265
350
310
total molecules
sum-val
0
1
11

SLIDER
260
115
425
148
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
105
150
425
183
max-ever-molecules
max-ever-molecules
stasis-level
100 * stasis-level
50000.0
1
1
NIL
HORIZONTAL

MONITOR
180
265
257
310
generation
ticks / generation-time
17
1
11

BUTTON
205
225
260
258
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
150
225
205
258
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
315
225
370
258
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

BUTTON
205
40
325
73
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
260
225
315
258
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

SLIDER
105
185
260
218
gene-is-sensor
gene-is-sensor
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
260
185
425
218
cyc1-vs-cyc2
cyc1-vs-cyc2
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
145
80
377
113
ancestral-effect-generations
ancestral-effect-generations
0
100
2.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## ABOUT THE MODEL

This model simulates the activity of a gene expressed within the germline across generations in *C. elegans* to discover possible regulatory architectures that are compatible with observed RNA silencing phenomena.

It takes 14 cell divisions to go from one zygote to another. This duration frome one bottleneck stage to the next is characteristic for each organism and includes rapid cell divisions (during early development), extended quiescent periods (Z2 and Z3 until mid-L1 in *C. elegans*), proliferative divisions (to generate adult germline), and meiotic development (to generate oocytes), followed by fertilization. Therefore, it is best to decide on a time scale of simulation and encode cell division vs growth as 1 vs 0 in a list. For the *C. elegans* lineage, keeping each tick = 15 minutes resulted in 365 ticks for each generation, which is used in this simulation.

To account for intermediate steps between a sensor and the downstream sensor/entity, a delay is chosen randomly for each regulatory interaction. This delay can take up to two generations (i.e., 730 ticks) to account for parental effects.

To simulate cycling stores of regulatory information that can be stable across generations, two mutually activating sensors are used [cyc1 and cyc2]. These sensors can themselves be made of numerous entities and sensors in principle.

Along the cell lineage connecting two generations, the production of known and measurable regulators (e.g., mRNA and small RNAs) can be promoted by regulators that are part of these cycling stores and thereby be maintained across generations. When a reporter is created to express a fluorescent protein under the control of different cis regulatory sequences, regulatory inputs into promoter, 5’UTR and 3’UTR of the endogenous gene can be nearly recreated [summarized as +s for positive regulators and -s for negative regulators]. In contrast, the regulatory interactions of the coding regions of the mRNA or the translated protein product is likely not recreated. Therefore, in this agent-based model, all small RNAs are classified as either positive or negative regulators [+p and -p] that are produced through interactions with the cycling stores using distinct templates and that in turn regulate mRNAs of complementary sequence (i.e., the endogenous gene or its reporter).

There are a total of 10 entities/sensors in the model. Minimally 4 and maximally 8 of these can be measured experimentally. The measurable entities are:
1. gene mRNA - this entity/sesnor representing the mRNA of a gene of interest (e.g. gtbp-1)
2. reporter mRNA - this entity is a reporter (e.g., mCherry) for the gene of interest that is expressed using the same cis regulatory sequences (i.e., promoter and 3' UTR)
3. +s : this represents positive regulators that act through the cis sequences either at the DNA or RNA level.
4. -s : this represents negtive regulators that act through the cis sequences either at the DNA or RNA level.
5. +p : small RNAs that are used to increase the amount of mRNA (e.g., through stabilization, 'CSR-1 licensing', etc.). One for gene mRNA and one for reporter mRNA.
6. -p : small RNAs that are used to decrease the amount of 'active' mRNA (e.g., through decay, localization away from ribosomes, etc.). One for gene mRNA and one for reporter mRNA.
	Based on how much is known about a system and how many entities/sensors can be quantitatively measured, more can be added to the simulation.

There are two entities in the model that cannot be measured experimentally [cyc1 and cyc2], but are necessary for simulating replication: These are two sensors that activate each other and are the drivers of heredity across generations. These together form the minimal configuration of cycling stores of information (Jose, BioEssays, 2020) that need to allow replication of the regulatory architecture and need to be added to all models of heritable epigenetic information. Both sensors are made of numerous unknown entities and sensors - potentially the entire cell code except for the measurable entities explicitly simulated above have been arbitrarily separated into these two mutually enforssing sensors.

The gene(s) of interest could be part of this drive by acting as sensors that change the value of cyc1 or cyc2, or they can be entities that are changed by cyc1 or cyc2 without any reciprocal input.
When a reporter is added for a gene, the regulatory inputs into the gene are recreated for the reporter. However, the reporter is considered an entity because the gene products (mCherry protein, mRNA, pre-mRNA, etc.) are not expected to regulate anything within the system in the ideal case. This assumption can be relaxed to discover mechanisms through which extra copies of cis regulatory elements could effect a gene of interest.

In this model, some imperfection in recreation of the regulatory inputs is allowed to simulate reporters that do not quite perfectly capture endogenous gene regulation. Specifically, the regulatory direction and linked entities are preserved, but the sensitivity of regulation (weight of link) and speed of regulation (delay of link) are allowed to vary.

Each entity/sensor has the following attributes:
	val ; an entity/sensor's current 'activity' (e.g., concentration, % conformational change, sequence, etc.) represented visually as a scaled *relative* size of the node. The command 'ask turtles [ show val ]' will show the actual absolute numerical values.
  	norm-val ; normalized value of each entity/sensor used for plotting that monitors the change across cell division cycles or equivalently across generations.
  	property ; The steps of change by which activity changes if a positive or negative interaction crosses the threshold for change. Properties are characteristic of the entities/sensors themselves and do not change as the system evolves through interactions. Note that this specifies that every sensor sees the same property of a given entity/sensor.
  	inactive-fraction ; Fraction of entity/sensor not available for regulatory interactions at each tick because of processes like protein folding, compartmentalization, diffusion, etc. This is a characteristic of each entity/sensor.
  	turnover ; this is a term (0.0 to 1.0) for the fraction remaining after one tick because of turnover processes in the system that are not explicitly modeled.
  	last-delay-time-values ; a convenient variable for remembering past n values (e.g., 10 needed for 5 germline zones and 2 generations) and using them later for calculations.

Each regulatory link has the following attributes:
	weight ; weight indicates the number of sensors needed to change one unit of property for each entity (characteristic of each 'interaction').
  	delay ; delay indicates the time after increase in a sensor that the effect on the regulated entity/sensor can be seen.

Particular links are positive (grey) or negative (black) and have a weight that depicts the numbers of sensors needed to cause one unit of property change in the downstream regulated sensor. They also have a delay of 0 to [(generation-time * ancestral-effect-generations) - 1] ticks in transmitting the effect of change in the sensor to downstream sensors/entities.

The maximum number of molecules there can ever be [max-ever-molecules] puts an upper bound on the speed with which number of molecules can grow per tick. This simulates the cell death through explosive production of some molecules.
The stasis-level simulates waiting for molecules that may have accrued above a large number (i.e., carrying capacity of the cell) to drop down by turnover or dilution upon cell division before regulation resumes. This is equivalent to the system running out of 'food' molecules. To account for depletion of precursors, reactants, or building blocks, a maximal increase of 5000 molecules was allowed per tick. To account for the limited capacity of a cell, the simulation ends if total number of molecules increases beyond 50000.
    Experiments can guide changes to these values for total numbers of molecules.


## COPYRIGHT AND LICENSE


Copyright 2008 Uri Wilensky for NetLogo and 2021 Antony Jose for this simulation of regulatory architectures within the C. elegans germline.

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
<experiments>
  <experiment name="Heritable_RNA_silencing_exploration_with_365_tick_generation_100_generations" repetitions="100000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="36500"/>
    <exitCondition>count turtles &lt; 10</exitCondition>
    <metric>count turtles</metric>
    <metric>count active-links</metric>
    <metric>count active-links with [color = black]</metric>
  </experiment>
  <experiment name="Heritable_RNA_silencing_exploration_with_365_tick_generation_100_generations" repetitions="1000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="36500"/>
    <exitCondition>count turtles &lt; 6</exitCondition>
    <metric>count turtles</metric>
    <metric>count active-links</metric>
    <metric>count active-links with [color = black]</metric>
  </experiment>
  <experiment name="Heritable_RNA_silencing_exploration_with_365_tick_generation_100_generations_for_v2" repetitions="10000" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="36500"/>
    <exitCondition>count turtles &lt; 10</exitCondition>
    <metric>count turtles</metric>
    <metric>count active-links</metric>
    <metric>count active-links with [color = black]</metric>
  </experiment>
</experiments>
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
