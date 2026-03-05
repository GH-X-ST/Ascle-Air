<p align="center">
  <sub>Student Group Design Project repository (Department of Aeronautics, Imperial College London).</sub><br>
  <sub>For educational and reference use only. Not an operational design manual, and definitely not certified for anything.</sub><br>
</p>

![Cover](General/01%20-%20Brand%20Guideline/Logo/HALO_PREVIEW.jpg)

## About
Helicopter Emergency Medical Services (HEMS) aircraft are currently adaptations of existing helicopter models, which, whilst adequate, are limited in their ability to serve the specific needs of the HEMS mission.

HALO explores a **3.5-tonne-class** conceptual design with **coaxial rotors**, targeting **~140 kt cruise** and **>600 km range**, powered by **twin Honeywell HTS9000** engines.  

The concept also looks at **active flow control** for cruise efficiency, an **H-tail** for rear loading and fast turnarounds, and onboard systems capable of supporting high auxiliary loads for medical equipment.

<p align="center">
  <img src="General/C%20-%20Final%20Submissions/Render2.png" width="800">
</p>

&nbsp;

## Contents
- [Deliverables](#deliverables)
- Catalogue
  - [Shared Resources](#shared-resources)
  - [HER01 — Project Coordinators](#her01--project-coordinators)
  - [HER02 — Aerodynamics](#her02--aerodynamics)
  - [HER03 — Structural Design & CAD](#her03--structural-design--cad)
  - [HER04 — Systems & Propulsion](#her04--systems--propulsion)
  - [HER05 — Flight Mechanics, Performance & Control](#her05--flight-mechanics-performance--control)

---

## Deliverables
Everything you’d typically want to open first:

- **Executive Summary:** [HER25_Executive_Summary.pdf](General/C%20-%20Final%20Submissions/HER25_Executive_Summary.pdf)
- **Final Presentation:** [final.pptx](General/C%20-%20Final%20Submissions/final.pptx)
- **Poster:** [PosterA1.pdf](General/C%20-%20Final%20Submissions/PosterA1.pdf)
&nbsp;

---

## Catalogue

### Shared Resources
Shared documents, templates, and reference libraries.

- **Briefing and admin documents:** [GDP guidance, kickoff slides, teamwork workshop material, etc.](General/00%20-%20Briefing/)
- **Brand assets & templates:** [Logo, fonts, report template zip, PowerPoint templates](General/01%20-%20Brand%20Guideline/)
- **Past reports library:** [“best” GDP reports](General/A%20-%20Past%20Best%20Reports/), [rotorcraft GDP reports](General/B%20-%20Past%20Helicopter%20Reports/GDP/), [VFS Student Design Competition](General/B%20-%20Past%20Helicopter%20Reports/VFS%20Student%20Design%20Competition/)
- **Final submissions:** [Final report, presentation, poster, etc.](General/C%20-%20Final%20Submissions/)

### HER01 — Project Coordinators
Project coordination, market/cost analysis, and cross-team integration.

- [Market Analysis](HER01/Market_analysis.m)
- [Cost](HER01/cost.m)
- [Duties](HER01/Duties%20GDP.docx)


### HER02 — Aerodynamics
Aerodynamic design tasks.

- Stability and force scripts: [Horizontal tailplane](HER02/Empennage/FHT.m), [Vertical tailplane](HER02/Empennage/FVT.m)
- Airfoil reference: [NACA0015](HER02/Empennage/NACA0015.m)


### HER03 — Structural Design & CAD
Structural sizing, mass estimation, and CAD-related calculations.

- Consolidated working script: [Integrated code](HER03/Big_code_of_all_mostly_working_sections.m)
- Structure weight calculation: [Structure weight estimation](HER03/afdd_structure_weight_estimation.m)
- Overall weight calculation: [Initial MTOW estimation](HER03/initial_mtow.m)


### HER04 — Systems & Propulsion
System architecture, weights, CG, propulsion integration, and equipment placement.

- **Weights and CG aggregation:** [Conceptual weights](HER04/conceptual_weights.mlx)
- **Propulsion:** [Propulsion calculation](HER04/Propulsion.mlx)
- **Mass input sheets:** [Configuration 1](HER04/mass_input_cfg_1.xlsx), [Configuration 2](HER04/mass_input_cfg_2.xlsx), [Configuration 3](HER04/mass_input_cfg_3.xlsx), [Configuration 4](HER04/mass_input_cfg_4.xlsx)
- **Placement tools:** [Algorithm](HER04/placement_algo_testing.mlx), [Input](HER04/placement_input.xlsx)


### HER05 — Flight Mechanics, Performance & Control
Mathematical model, trim, stability, control, and performance workflows.

- **01 — Mathematical Modelling:** [Coaxial rotor BEMT scripts, linearisation, Simulink models, and test notebooks](HER05/01%20-%20Mathematical%20Modeling/)
- **02 — Trim:** [Trim solving and plotting utilities](HER05/02%20-%20Trim/)
- **03 — Static Stability:** [Static Stability calculation](HER05/03%20-%20Static%20Stability/)
- **04 — Dynamic Stability:** [Linearisation tools and stability analysis scripts](HER05/04%20-%20Dynamic%20Stability/)
- **05 — Control Systems:** [PID](HER05/05%20-%20Control%20Systems/01%20-%20PID/), [H-infinity](HER05/05%20-%20Control%20Systems/02%20-%20H-infinity/)
- **06 — Performance:** [Performance and power calculation scripts](HER05/06%20-%20Performance/)
