# HemeLB

A typical workflow with HemeLB consists of four steps:

1. A preprocessing step where you create a mesh or geometry file (we
   use extension .gmy) from a representation of the surface of your
   domain. For this, you need to use the application in the
   [geometry-tool](../geometry-tool)
   directory. [Documentation](user/geometry-tool.md)

2. Setting options in the configuration XML file, such as timestep and
   output data and
   frequency. [Documentation](user/XmlConfiguration.md)

3. Compiling and running the
   [main application](user/main-application.md).

4. Post processing the results, using the Python packages in the
   [python-tools](../python-tools) directory. [Documentation](user/python-tools.md)


[Developer documentation](dev)

# Publications

## Key code publications

- R.W. Nash, H.B. Carver, M.O. Bernabeu, J. Hetherington, D. Groen, T.
  Krüger, P.V. Coveney, "Choice of boundary condition for
  lattice-Boltzmann simulation of moderate-Reynolds-number flow in
  complex domains", Phys. Rev. E (2014).
  https://doi.org/10.1103/PhysRevE.89.023303

- D. Groen, J. Hetherington, H.B. Carver, R.W. Nash, M.O. Bernabeu,
  "Analysing and modelling the performance of the HemeLB
  lattice-Boltzmann simulation environment", J. Comput. Sci. (2013).
  https://doi.org/10.1016/j.jocs.2013.03.002

- M.D. Mazzeo & P.V. Coveney, "HemeLB: A high performance parallel
  lattice-Boltzmann code for large scale fluid flow in complex
  geometries", Comput. Phys. Commun. (2008)
  https://doi.org/10.1016/j.cpc.2008.02.013


## Papers using HemeLB

- Q. Zhou, K. Schirrmann, E. Doman, Q. Chen, N. Singh, P. Ravi Selvaganapathy,
  M.O. Bernabeu, O.E. Jensen, A. Juel, I.L. Chernyavsky, T. Krüger.
  Red blood cell dynamics in extravascular biological tissues modelled as
  canonical disordered porous media. Interface Focus 12, 20220037 (2022).
  https://dx.doi.org/10.1098/rsfs.2022.0037

- Q. Zhou, T. Perovic, I. Fechner, L.T. Edgar, P.R. Hoskins, H. Gerhardt,
  T. Krüger, M.O. Bernabeu. Association between erythrocyte dynamics and
  vessel remodelling in developmental vascular networks.
  J. R. Soc. Interface 18, 20210113 (2021).
  https://doi.org/10.1098/rsif.2021.0113

- R. Enjalbert, D. Hardman, T Krüger, M.O. Bernabeu. Compressed vessels
  bias red blood cell partitioning at bifurcations in a hematocrit-dependent
  manner: Implications in tumor blood flow. PNAS 118, e2025236118 (2021).
  https://doi.org/10.1073/pnas.2025236118

- Q. Zhou, J. Fidalgo, M.O. Bernabeu, M.S.N. Oliveira, T. Krüger.
  Emergent cell-free layer asymmetry and biased haematocrit partition
  in a biomimetic vascular network of successive bifurcations.
  Soft Matter 17, 3619-3633 (2021).
  https://doi.org/10.1039/D0SM01845G

- M.O. Bernabeu, J. Köry, J.A. Grogan, B. Markelc, A.B. Ricol, M. d’Avezac,
  R. Enjalbert, J. Kaeppler, N. Daly, J. Hetherington, T. Krüger, P.K. Maini,
  J.M. Pitt-Francis, R.J. Muschel, T. Alarcón, H.M. Byrne. Abnormal morphology
  biases haematocrit distribution in tumour vasculature and contributes to
  heterogeneity in tissue oxygenation. PNAS 117, 27811-27819 (2020).
  https://doi.org/10.1073/pnas.200777011

- Q. Zhou, J. Fidalgo, L. Calvi, M.O. Bernabeu, P.R. Hoskins, M.S.N. Oliveira,
  T. Krüger. Spatiotemporal Dynamics of Dilute Red Blood Cell Suspensions
  in Low-Inertia Microchannel Flow. Biophys. J. 118, 2561-2573 (2020).
  https://doi.org/10.1016/j.bpj.2020.03.019
