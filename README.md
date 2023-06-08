# Heritable epigenetic changes are constrained by the dynamics of regulatory architectures

Antony M Jose<sup>*</sup>

<sup>*</sup>Corresponding author email: amjose@umd.edu

This paper establishes criteria for (1) heritable information in regulatory architectures, (2) the persistence of epigenetic changes across generations, (3) distinguishing between regulatory architectures using transient perturbations, (4) making unstable regulatory architectures in the germline heritable through interactions with somatic cells, and (5) generating epigenetic changes of defined magnitude and duration.

For this work, several programs were written to perform simulations and analyze data. Briefly, the possible weakly connected non-isomorphic digraphs that form regulatory architectures capable of indefinite persistence and the maximal information that can be stored using them (log<sub>2</sub>N) were calculated. Systems of ordinary differential equations were used to describe the rate of change for each interacting entity (x,y,z) in each of the 26 simplest heritable regulatory architectures (HRAs). Constraints based on steady state were used to obtain sets of parameters that include rates of promotion and/or inhibition as well as turnover for each entity/sensor. Particular parameter sets were then used to simulate steady states and illustrate the consequences of genetic or epigenetic changes to all 26 HRAs. The duration and extent of perturbation needed for heritable epigenetic changes were explored using numerical solutions of the equations describing the dynamics of each HRA. The analytical expressions for conditions that enable heritable epigenetic change were derived for the simplest heritable regulatory architecture where two entities (x and y) mutually promote each other’s production and these calculations were verified using numerical simulations.

Simple Entity-Sensor-Property (ESP) systems were simulated using custom programs in NetLogo (details are available within the ‘code’ and ‘info’ tabs of each program). Results from systematic explorations obtained by running NetLogo programs from the command line (i.e., ‘headless’) were analyzed using R. The developmental timing of cell division is *C. elegans* were curated manually from the literature and used to simulate ESP systems that incorporate this developmental timing and temporal delays in regulation.

Each program used is listed below grouped by the software needed. Details for each step are commented on within the programs. Please modify paths as required before use.


### Programs in NetLogo (v. 6.1.1):

(1) This program can be run 'headless' from the command line to identify sets of parameters for ESP systems that are heritable. Each system is identified by a `behaviorspace-run-number`.

 `ESP_systems_explorer_v1_for_parameter_sweep.nlogo`
  
(2) This program can be used to interactively explore ESP systems using `system-id` as the `behaviorspace-run-number`.

  `ESP_systems_single_system_explorer_v1.nlogo`

(3) This program can be run 'headless' from the command line to identify sets of parameters for ESP systems with HRAs that include developmental timing of cell divisions in *C. elegans* and parental regulation. Each system is identified by a `behaviorspace-run-number`.

  `Heritable_RNA_silencing_with_developmental_time_v4_for_parameter_sweep.nlogo`  
  
(4) This program can be used to interactively explore ESP systems with heritable regulatory architectures that include developmental timing of cell divisions in *C. elegans* and parental regulation using `system-id` as the `behaviorspace-run-number`.

  `Heritable_RNA_silencing_with_Celegans_developmental_time_single_system_explorer.nlogo`


### Programs in Python (v. 3.8.5):
  
(1) This program can be used to identify weakly-connected digraphs with regulation, which define regulatory architectures that can persist forever and thus are HRAs.
  
  `Heritable_Regulatory_Architectures_1-4_entities.py`

(2) This program can be used to identify parameter sets for each HRA that result in a steady state.

  `HRA_steady_state_params_d1.py`
  
(3) These programs can be used to examine the consequence of genetic or epigenetic changes of each entity/sensor of each HRA from steady state. The results were used in supplementary figures of the paper.
  
  `HRA_perturb_genetic_d3.py`
    
  `HRA_perturb_epigenetic_d3.py`
  
(4) These programs can be used to examine the consequence of genetic or epigenetic changes of each entity/sensor from steady state that were used to illustrate differences in the main figures of the paper.
  
  `HRA_perturb_genetic_vs_epigenetics_illustration_epigen_MainFig.py`
  
  `HRA_perturb_genetic_vs_epigenetics_illustration_gen_MainFig.py`
  
(5) This program was used to illustrate how the strength and duration of a perturbation alters the recovered steady state of the simplest HRA (type 'A' in the paper).
  
  `HRA_perturb_strength_duration_illustration.py`
  
(6) This program was used to verify the calculated analytical expression for the critical duration of perturbation needed for a perturbation of given strength to cause heritable epigenetic change by altering a component of a positive feedback loop.
  
  `HRA_A_tp_Analytical_Expression_check.py`
  

### Programs in R (v. 3.6.3):

(1) This program was used to analyze the parameters of stable/robust ESP systems obtained after running `ESP_systems_explorer_v1_for_parameter_sweep.nlogo` and results in several plots. 

  `ESP_parameter_extraction_for_paper.R`

