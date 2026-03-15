<p align="center">
  <sub>A group design project presented to the Department of Aeronautics</sub><br>
  <sub>in partial fulfilment of the requirements for the degree of</sub><br>
  <sub>Master of Engineering (MEng) in Aeronautical Engineering</sub><br>
  <sub>at</sub><br>
  <sub>Imperial College London</sub><br>
</p>

![Cover light](General/01%20-%20Brand%20Guideline/Logo/HALO_PREVIEW.jpg#gh-light-mode-only)
![Cover dark](General/01%20-%20Brand%20Guideline/Logo/HALO_PREVIEW_dark.jpg#gh-dark-mode-only)

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

- **Briefing and Admin Documents:** [GDP guidance, kickoff slides, teamwork workshop material, etc.](General/00%20-%20Briefing/)
- **Brand Assets & Templates:** [Logo, fonts, report template zip, PowerPoint templates](General/01%20-%20Brand%20Guideline/)
- **Past Reports Library:** [“best” GDP reports](General/A%20-%20Past%20Best%20Reports/), [rotorcraft GDP reports](General/B%20-%20Past%20Helicopter%20Reports/GDP/), [VFS Student Design Competition](General/B%20-%20Past%20Helicopter%20Reports/VFS%20Student%20Design%20Competition/)
- **Final Submissions:** [Final report, presentation, poster, etc.](General/C%20-%20Final%20Submissions/)


### HER01 — Project Coordinators
Project coordination, market analysis, cost calculation, and cross-team integration.

- **Market Analysis:** [Market analysis script](HER01/Market_analysis.m)
- **Cost:** [Cost estimation](HER01/cost.m)
- **Duties:** [Duties allocation](HER01/Duties%20GDP.docx)


### HER02 — Aerodynamics
Aerodynamic design tasks.

- **Stability and Force Calculation:** [Horizontal tailplane](HER02/Empennage/FHT.m), [Vertical tailplane](HER02/Empennage/FVT.m)


### HER03 — Structural Design & CAD
Structural sizing, mass estimation, and CAD-related calculations.

- **Consolidated Working Script:** [Integrated code](HER03/Big_code_of_all_mostly_working_sections.m)
- **Weight Calculation:** [Structure weight estimation](HER03/afdd_structure_weight_estimation.m), [Initial MTOW estimation](HER03/initial_mtow.m)


### HER04 — Systems & Propulsion
System architecture, weights, CG, propulsion integration, and equipment placement.

- **Weights and CG Aggregation:** [Conceptual weights](HER04/conceptual_weights.mlx)
- **Propulsion:** [Propulsion calculation](HER04/Propulsion.mlx)
- **Placement Tools:** [Algorithm](HER04/placement_algo_testing.mlx), [Input](HER04/placement_input.xlsx)


### HER05 — Flight Mechanics, Performance & Control
Mathematical model, trim, stability, control, and performance workflows.

- **Mathematical Modelling:** [Coaxial rotor BEMT scripts, linearisation, Simulink models, and test notebooks](HER05/01%20-%20Mathematical%20Modeling/)
- **Trim:** [Trim solving and plotting utilities](HER05/02%20-%20Trim/)
- **Static Stability:** [Static Stability calculation](HER05/03%20-%20Static%20Stability/)
- **Dynamic Stability:** [Linearisation tools and stability analysis scripts](HER05/04%20-%20Dynamic%20Stability/)
- **Control Systems:** [PID](HER05/05%20-%20Control%20Systems/01%20-%20PID/), [H-infinity](HER05/05%20-%20Control%20Systems/02%20-%20H-infinity/)
- **Performance Analysis:** [Performance and power calculation scripts](HER05/06%20-%20Performance/)

---

## Stargazers Over Time
[![Stargazers over time light](https://starchart.cc/GH-X-ST/HALO.svg?background=%2300000000&axis=%23232333&line=%23708090)](https://starchart.cc/GH-X-ST/HALO#gh-light-mode-only)
[![Stargazers over time dark](https://starchart.cc/GH-X-ST/HALO.svg?background=%2300000000&axis=%23FFFFFF&line=%239eb3c9)](https://starchart.cc/GH-X-ST/HALO#gh-dark-mode-only)
