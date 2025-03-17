# SCPF-Array-Modeling

This repository provides functions and codes to estimate the cross-spectral density and wavenumber-frequency spectrum for an array of area-averaging sensors with a finite surface area. 
The array of sensors was designed to measure the sub-convective wall pressure fluctuations underneath turbulent boundary layer flows.
The sample codes have been written assuming a sensor arrangement used by [Damani et al. 2024](https://arc.aiaa.org/doi/abs/10.2514/6.2024-3048). 
The sensors comprised 80 sub-resonant sensors arrange periodically along the flow direction and each sensor had a rectangular profile with a 5 mm length along the flow direction and 50 mm in the spanwise direction.
The mathematical approach taken towards the codes is described in section V of the conference publication.
Four wall pressure modeling codes have been provided for four well-known wall pressure models: [Corcos](https://www.sciencedirect.com/science/article/abs/pii/0022460X67901587), [Modified Corcos](https://apps.dtic.mil/sti/citations/ADA465806), [Smol'yakov](https://link.springer.com/article/10.1134/S1063771006030146) and [incompressible Chase](https://www.sciencedirect.com/science/article/abs/pii/0022460X80905532). 
Users are encouraged to use the codes to run different configurations of sensors and adapt these to flow-specific conditions.

## List of Functions and Scripts
+ AirDensityCalc - Calculates air density for a given temperature and pressure condition
+ AirSpeedOfSoundCalc - Calculates airspeed for a given temperature
+ AirViscosityCalc - Calculates air dynamic viscosity for a given temperature
+ ChaseArrayPrediction - Uses the [incompressible Chase](https://www.sciencedirect.com/science/article/abs/pii/0022460X80905532) model to evaluate wavenumber-frequency spectrum of the array
+ CorcosArrayPrediction - Uses the [Corcos](https://www.sciencedirect.com/science/article/abs/pii/0022460X67901587) model to evaluate wavenumber-frequency spectrum of the array
+ ModCorcosArrayPrediction - Uses the [Modified Corcos](https://apps.dtic.mil/sti/citations/ADA465806) model to evaluate wavenumber-frequency spectrum of the array
+ SmolyakovArrayPrediction - Uses the [Smol'yakov](https://link.springer.com/article/10.1134/S1063771006030146) model to evaluate wavenumber-frequency spectrum of the array

## References
+ [The Characteristics of Sub-Convective Wall Pressure Fluctuations in a Turbulent Boundary Layer](https://arc.aiaa.org/doi/abs/10.2514/6.2024-1904)
+ [Evaluating Models for Sub-Convective Pressure Fluctuations in Turbulent Boundary Layers](https://arc.aiaa.org/doi/abs/10.2514/6.2024-3048)
+ [Measurement and Analysis of Sub-Convective Pressure Fluctuations in Turbulent Boundary Layers: A Novel Methodology](https://vtechworks.lib.vt.edu/items/d3a198f2-1437-46fb-959f-310f25084c21)
+ For five boundary layer flow conditions and wind tunnel data: [DOI 10.7294/28395098](https://doi.org/10.7294/28395098.v1)
