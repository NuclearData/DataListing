# ENDF/B Reaction Types {#sec:types}

This section lists some of the more useful reactions for use with the `FMn` input card and with the `MCNP`` `cross-section plotter. See the relevent section of the `MCNP`` `manual for more information. The complete ENDF/B list can be found in the ENDF/B manual \[MCL95\]. The tables and reaction types listed are:

| Table                                                                              | Description                              |
| ----------------------------------------------------------------------------------- | ----------------------------------------- |
| [1](#tab:NeutronReactions){reference-type="ref" reference="tab:NeutronReactions"}  | Neutron continuous-energy reactions      |
| [2](#tab:SaB){reference-type="ref" reference="tab:SaB"}                            | $S(\alpha, \beta)$ reactions             |
| [3](#tab:multigroup){reference-type="ref" reference="tab:multigroup"}              | Neutron and photon multigroup reactions  |
| [4](#tab:Photoatomic){reference-type="ref" reference="tab:Photoatomic"}            | Photoatomic reactions                    |
| [5](#tab:Photonuclear){reference-type="ref" reference="tab:Photonuclear"}          | Photonuclear reactions                   |
| [6](#tab:Electrons){reference-type="ref" reference="tab:Electrons"}                | Electron reactions                       |

In each of these tables, the MT column lists the ENDF/B reaction number and the FM column lists special `MCNP``6` reaction numbers that can be used with the `FM` card and cross-section plotter.

The nomenclature between `MCNP``6` and ENDF/B is inconsistent in that `MCNP``6` often refers to the number of the reaction type as R whereas ENDF/B uses MT, but they are the same. The problem arises because `MCNP``6` has an `MT` input card used for the $S(\alpha, \beta)$ thermal treatment. However, the nomenclature between Monte Carlo transport and deterministic transport techniques can be radically different. The reference \[FRA96c\] provides more information.

Generally only a subset of reactions is available for a particular nuclide. Some reaction data are eliminated by `MCNP``6` in cross-section processing if they are not required by the problem. Examples are photon production in a `MODE N` problem, or certain reaction cross sections not requested on an `FM` card. FM numbers should be used when available, rather than MT numbers. If an MT number is requested, the equivalent FM number will be displayed on the legend of cross-section plots.

::: {#tab:NeutronReactions}

| MT         | FM   | Microscopic Cross-Section Description                                                                                  |
| ---------- | ---- | ---------------------------------------------------------------------------------------------------------------------- |
| 1          | -1   | Total                                                                                                                  |
| 2          | -3   | Elastic                                                                                                                |
| 16         |      | $(n,2n)$                                                                                                               |
| 17         |      | $(n,3n)$                                                                                                               |
| 18         |      | Total fission $(n,fx)$ if and only if MT=18 is used to specify fission in the original evaluation.                     |
|            | -6   | Total fission cross section. (equal to MT=18 if MT=18 exists; otherwise equal to the sum of MTs 19, 20, 21, and 38.)   |
| 19         |      | $(n,f)$                                                                                                                |
| 20         |      | $(n,n'f)$                                                                                                              |
| 21         |      | $(n,2nf)$                                                                                                              |
| 22         |      | $(n,n'\alpha)$                                                                                                         |
| 28         |      | $(n,n'p)$                                                                                                              |
| 32         |      | $(n,n'd)$                                                                                                              |
| 33         |      | $(n,n't)$                                                                                                              |
| 38         |      | $(n,3nf)$                                                                                                              |
| 51         |      | $(n,n')$ to 1st excited state                                                                                          |
| 52         |      | $(n,n')$ to 2nd excited state                                                                                          |
| $\cdots$   |      | $\cdots$                                                                                                               |
| 90         |      | $(n,n')$ to 40th excited state                                                                                         |
| 91         |      | $(n,n')$ to continuum                                                                                                  |
| 101        | -2   | Absorption: sum of MT=102--117 (neutron disappearance; does not include fission)                                       |
| 102        |      | $(n,\gamma)$                                                                                                           |
| 103        |      | $(n,p)$                                                                                                                |
| 104        |      | $(n,d)$                                                                                                                |
| 105        |      | $(n,t)$                                                                                                                |
| 106        |      | $(n,3He)$                                                                                                              |
| 107        |      | $(n,\alpha)$                                                                                                           |
| 202        | -5   | Total photon production                                                                                                |
| 203        |      | Total proton production                                                                                                |
| 204        |      | Total deuterium production                                                                                             |
| 205        |      | Total tritium production                                                                                               |
| 206        |      | Total 3He production                                                                                                   |
| 207        |      | Total alpha production                                                                                                 |
| 301        | -4   | Average heating numbers (MeV/collision)                                                                                |
|            | -7   | Nubar (prompt or total)                                                                                                |
|            | -8   | Fission Q (in print table 98, but not plots)                                                                           |

  : Neutron continuous-energy and discrete reactions.
:::

At the time they are loaded, the total and elastic cross sections from
the data library are thermally adjusted by `MCNP``6` to the temperature
of the problem, if that temperature is different from the temperature at
which the cross-section set was processed. If different cells have
different temperatures, the cross sections are first adjusted to zero
degrees and adjusted again to the appropriate cell temperatures during
transport. The cross-section plot will *never* display the *transport*
adjustment. Therefore, for plotting, reactions 1 and -1 are equivalent
and reactions 2 and -3 are equivalent. However, for the `FM` card,
reactions -1 and -3 will use the zero-degree data and reactions 1 and 2
will use the transport-adjusted data. For example, if a library
evaluated at 300 K is used in a problem with cells at 400 K and 500 K,
the cross-section plotter and `MT=-1` and `MT=-3` options on the
`FM `card will use 0 K data. The `MT=1` and `MT=2` options on the `FM`
card will use 400 K and 500 K data.

The user looking for total production of p, d, t, , and should be warned
that in some evaluations, such processes are represented using reactions
with MT (or R) numbers other than the standard ones given in the above
list. This is of particular importance with the so-called "pseudolevel"
representation of certain reactions which take place in light isotopes.
For example, the ENDF/B-V evaluation of carbon includes cross sections
for the $(n,n’3\alpha)$ reaction in `MT=52` to `58`. The user interested
in particle production from light isotopes should check for the
existence of pseudolevels and thus possible deviations from the above
standard reaction list.

::: {#tab:SaB}
   MT   FM  Microscopic Cross-Section Description
  ---- ---- ---------------------------------------
   1        Total cross section
   2        Elastic scattering cross section
   4        Inelastic scattering cross section

  : $S(\alpha, \beta)$ reactions.
:::

::: {#tab:multigroup}
   MT    FM  Microscopic Cross-Section Description
  ----- ---- ------------------------------------------
    1    -1  Total cross section
   18    -2  Fission cross section
         -3  Nubar data
         -4  Fission chi data
   101   -5  Absorption cross section
         -6  Stopping powers
         -7  Momentum transfers
    n        Edit reaction n
   202       Photon production
   301       Heating number
   318       Fission Q
   401       Heating number times total cross section

  : Neutron and photon multigroup reactions.
:::

::: {#tab:Photoatomic}
   MT    FM  Microscopic Cross-Section Description
  ----- ---- ---------------------------------------
   501   -5  Total
   504   -1  Incoherent (Compton + Form Factor)
   502   -2  Coherent (Thomson + Form Factor)
   522   -3  Photoelectric with fluorescence
   516   -4  Pair production
   301   -6  Heating number

  : Photoatomic reactions.
:::

::: {#tab:Photonuclear}
   MT    FM   Microscopic Cross-Section Description
  ---- ------ ---------------------------------------
         1    Total
         2    Non-elastic
         3    Elastic
         4    Heating
         5    Other
        1005  Neutron production from reaction 5
        2005  Photon production from reaction 5
        9005  Proton production from reaction 5

  : Photonuclear reactions.
:::

::: {#tab:Electrons}
   MT   FM  Microscopic Cross-Section Description
  ---- ---- -----------------------------------------
        1   de/dx electron collision stopping power
        2   de/dx electron radiative stopping power
        3   de/dx total electron stopping power
        4   electron range
        5   electron radiation yield
        6   relativistic $\beta^2$
        7   stopping power density correction
        8   ratio of rad/col stopping powers
        9   drange
        10  dyield
        11  rng array values
        12  qav array values
        13  ear array values

  : Electon reactions.
:::

LANL maintains two electron-transport libraries, EL and EL03. The
electron transport algorithms and data in `MCNP``6` were adapted from
the ITS code \[HAL92\]. The EL library was developed and released in
1990 in conjunction with the addition of electron transport into
`MCNP``4`; the electron-transport algorithms and data correspond
(roughly) to that found in ITS version 1. The EL03 library \[ADA00\] was
developed and released in 2000 in conjunction with upgrades to the
electron physics package; these upgrades correspond (roughly) to that of
ITS version 3.The MT numbers for use in plotting the cross-section
values for these tables should be taken from print table 85 column
headings and are not from ENDF.


