# Puma dispersal and ranging in Southeastern Brazil

This repository presents the code for the organization of data and the analyses 
of space use, dispersal, and movement of pumas in southeastern Brazil.
The code supplements the manuscript:

Niebuhr, B. B.; Cavalcanti, S. M. C.; Vilalba, E. A.; Alberico, V. V.; Gebin, J. C. Z.; 
Santos, D. C.; Barban, A. B.; de Oliveira, R.; Gurarie, E.; Morato, R. G. 2025. 
**Land use effects on the space use and dispersal of an apex predator in an ecotone between tropical biodiversity hotspots.** 
Diversity, 17.

Below we desribe the content of the repository.

## code

The folder **code/** presents the R scripts used to organize the data and perform the analyses presented in the paper.

0. [Organize raw GPS data](code/00_organize_raw_data.R) 
1. [Load and prepare GPS data](code/01_dispersal_load_data.md)
2. [Create R package to store the GPS data](code/02_create_data_package.md)
3. [Explore GPS data and make some basic plots](code/03_plot_explore.md)
4. [Finding dispersal timing: example for one individual](code/04_arima_fit_1ind.md)
5. [Finding dispersal timing: all pumas](code/05_arima_fit_multiple_animals_pkg.md)
6. [Understanding dispersal patterns](code/06_arima_understand_dispersal.md)
7. [Understanding movement patterns](code/07_understand_movement.md)
8. [Understanding residency patterns](code/08_understand_residency.md)
9. [Comparison with literature data](code/09_compare_literature.md)

## data

The GPS data collected and used in this study are stored on MoveBank (Movebank project ID 577226894 — Pumas from Tiete Project; 
Movebank project ID 594660300 — Onças do Legado), and access to it might be requested from the first author or through contact with 
Instituto Pró-Carnívoros. Data from the literature compilation is present in the **data/** folder in this repository. 

