# clusterCloud
The clusterCloud is an open-source library for OpenFOAM. It is compatible with OpenFOAM-9. 
The library provides an additional cloud called clusterCloud wich can be added via fvModels to an OpenFOAM solver. 
The available lagrangian particles are called ReactingClusterParcel and are based on the ReactingMultiphaseParcel. 
Additionally, a model to account for the turbulent clustering effect according to [1,2] is available. Besides the standard
oxidation model, a new reaction model based on kinetic-diffusion approach, adding variable CO/CO2 products was added. 
Similarly to the oxidation models, also gasification models (H2, CO2, H2O) are available. 

A similar approach has been implemented in Fluent, the results and the code (supplementary material) can be found in [3]

A diffusion function for the source terms of the particles has been optionally implemented as well and is based on the work from [4]. 

Two tutorial cases are available. A simple numerical example, similar to the one presented in [3] and an injection rig test case, which is 
based on the work from [5]

[1] Krüger, J., Haugen, N. E. L., & Løvås, T. (2017). Correlation effects between turbulence and the conversion rate of pulverized char particles. Combustion and Flame, 185, 160–172. https://doi.org/10.1016/j.combustflame.2017.07.008

[2] Haugen, N. E. L., Krüger, J., Mitra, D., & Løvås, T. (2018). The effect of turbulence on mass transfer rates of small inertial particles with surface reactions. Journal of Fluid Mechanics, 836, 932–951. https://doi.org/10.1017/jfm.2017.820

[3] Karchniwy, E., Haugen, N. E. L., Klimanek, A., Langørgen, Ø., & Sładek, S. (2021). The effect of turbulence on mass transfer in solid fuel combustion: RANS model. Combustion and Flame, 227, 65–78. https://doi.org/10.1016/j.combustflame.2020.12.040

[4] Zhang, J., Li, T., Ström, H., & Løvås, T. (2020). Grid-independent Eulerian-Lagrangian approaches for simulations of solid fuel particle combustion. Chemical Engineering Journal, 387(October 2019), 123964. https://doi.org/10.1016/j.cej.2019.123964

[5] Shen, Y. S., Guo, B. Y., Yu, A. B., & Zulli, P. (2009). A three-dimensional numerical study of the combustion of coal blends in blast furnace. Fuel, 88(2), 255–263. https://doi.org/10.1016/j.fuel.2008.08.013
