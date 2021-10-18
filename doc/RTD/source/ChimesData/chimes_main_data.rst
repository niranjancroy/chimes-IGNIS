.. CHIMES main data 
   Alexander Richings, 19th February 2020 

.. _ChimesMainData_label: 

chimes_main_data
----------------

The ``chimes_main_data.hdf5`` file encodes all of the chemical reactions that go into the network, and tabulates the reaction rate coefficients. The reactions are organised into different groups according to how their rate coefficients are calculated (for example, whether they are constant, depend on temperature etc). The idea is that, in a given sub-step of the chemistry solver, the rate coefficients of all reactions in a given group are updated in the same way, which involves a single loop over all reactions in that group. 

The various reaction groups are described in more detail below. For each group, we also describe the various data arrays that are included in ``chimes_main_data.hdf5``. 

constant
^^^^^^^^

Reactions for which the rate coefficients are just constants. They do not need to be updated every sub-step in the chemistry solver, they can just be taken directly from the data tables. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. The first number gives the number of reactions that do not     |
|                                     | | involve molecules, while the second number gives the total number of       |
|                                     | | reactions, including molecules.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the         |
|                                     | | species indices of each reactant in each reaction. All reactions in this   |
|                                     | | group have exactly 2 reactants. Note that these indices correspond to the  |
|                                     | | position of the species in the full CHIMES network. When CHIMES is used    |
|                                     | | with a reduced network (i.e. if some metal elements have been switched     |
|                                     | | off), then CHIMES will internally re-map these indices to the              |
|                                     | | corresponding positions in the reduced network.                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 3) that gives the         |
|                                     | | species indices of each product in each reaction. If a reaction has fewer  |
|                                     | | than 3 product, the remaining indices are set to ``-1``. No reactions      |
|                                     | | have more than 3 products.  Note that these indices correspond to the      |
|                                     | | position of the species in the full CHIMES network. When CHIMES is used    |
|                                     | | with a reduced network (i.e. if some metal elements have been switched     |
|                                     | | off), then CHIMES will internally re-map these indices to the              |
|                                     | | corresponding positions in the reduced network.                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 1-dimensional array of length ``N_reactions`` that gives the constant    |
|                                     | | rate coefficients in cgs units (as linear, not log) for each reaction.     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. For each reaction we have an array |
|                                     | | of 9 integers that can be 0 or 1, indicating that the corresponding        |
|                                     | | element either is not or is involved in the given reaction, respectively.  |
|                                     | | The elements are encoded in the order C, N, O, Ne, Mg, Si, S, Ca, Fe. Then |
|                                     | | if one or more of the elements needed for a given reaction has been        |
|                                     | | switched off, that reaction will be skipped and will not be included in    |
|                                     | | the reduced network. Note that hydrogen and helium cannot be switched off, |
|                                     | | so we do not include flags for these.                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). CHIMES constructs two networks, one     |
|                                     | | including the molecular chemistry and one without it, and switches between |
|                                     | | them based on the temperature. This allows us to switch off molecules at   |
|                                     | | high temperatures. We can then use the ``molecular_flag`` from the tables  |
|                                     | | to determine which reactions can be included in the non-molecular network. |
|                                     | | Note that the reactions in each group are organised such that those that   |
|                                     | | include molecules are at the end of the group.                             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2_form_heating_``              | | Index giving the position of the H2 gas-phase formation reaction in this   |
| | ``reaction_index``                | | group. This is needed by the cooling routines, as this reaction            |
|                                     | | contributes to the heating rate, so we will need to access the             |
|                                     | | corresponding reaction rate.                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

T_dependent
^^^^^^^^^^^

Reactions for which the rate coefficients are functions of temperature only. If thermal evolution is switched on, they will need to be updated every sub-step of the integration. The rates are tabulated as a function of temperature, so to update the rate coefficients requires a 1-dimensional interpolation for each reaction. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 2-dimensional array of size (``N_reactions`` x 3) that gives the         |
|                                     | | species indices of each reactant in each reaction. If a reaction has fewer |
|                                     | | than 3 reactants, the remaining indices are set to ``-1``. No reactions    |
|                                     | | have more than 3 reactants. See the description under the **constant**     |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 3) that gives the         |
|                                     | | species indices of each product in each reaction. If a reaction has fewer  |
|                                     | | than 3 product, the remaining indices are set to ``-1``. No reactions      |
|                                     | | have more than 3 products. See the description under the **constant**      |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 2-dimensional array of size (``N_reactions`` x ``N_Temperatures``) that  |
|                                     | | tabulates the rate coefficients of each reaction as a function of          |
|                                     | | temperature. These are given as log10(rate) in cgs units.                  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2_collis_dissoc_``             | | Index giving the position of the H2 collisional dissociation reaction in   |
| | ``heating_reaction_``             | | this group. This is needed by the cooling routines, as this reaction       |
| | ``index``                         | | contributes to the heating rate, so we will need to access the             |
|                                     | | corresponding reaction rate.                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2_form_heating_``              | | Index giving the position of the H2 gas-phase formation reaction in this   |
| | ``reaction_index``                | | group. This is needed by the cooling routines, as this reaction            |
|                                     | | contributes to the heating rate, so we will need to access the             |
|                                     | | corresponding reaction rate.                                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

recombination_AB
^^^^^^^^^^^^^^^^

For the recombination of singly ionised hydrogen and helium, CHIMES switches between case A and case B recombination based on whether the gas is optically thin to the corresponding recombination radiation (as determined by the column densities of HI and HeI, respectively). In this reaction group, we have therefore tabulated the corresponding temperature-dependent case A and case B rate coefficients separately. The CHIMES solver will then select which one to use accordingly. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the         |
|                                     | | species indices of each reactant in each reaction. All reactions in this   |
|                                     | | group have exactly 2 reactants. See the description under the **constant** |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of the product in each reaction. All reactions in this group have  |
|                                     | | only 1 product, which is why we only use a 1-dimensional array here. See   |
|                                     | | the description under the **constant** group for further details.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates_caseA``                     | | A 2-dimensional array of size (``N_reactions`` x ``N_Temperatures``) that  |
|                                     | | tabulates the rate coefficients of the case A recombination reactions as a |
|                                     | | function of temperature. These are given as log10(rate) in cgs units.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates_caseB``                     | | A 2-dimensional array of size (``N_reactions`` x ``N_Temperatures``) that  |
|                                     | | tabulates the rate coefficients of the case B recombination reactions as a |
|                                     | | function of temperature. These are given as log10(rate) in cgs units.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

grain_recombination
^^^^^^^^^^^^^^^^^^^

Recombination of ions with electrons on the surface of dust grains. The rate coefficients depend on temperature and on the parameter ``Psi = G0 exp(-2.77 Av) T^0.5 / ne``, where ``G0`` is the strength of the radiation field in the 6-13.6 eV band in Habing units, ``Av`` is the dust extinction, ``T`` is the gas temperature, and ``ne`` is the electron density in cgs units. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the         |
|                                     | | species indices of each reactant in each reaction. All reactions in this   |
|                                     | | group have exactly 2 reactants. See the description under the **constant** |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of the product in each reaction. All reactions in this group have  |
|                                     | | only 1 product, which is why we only use a 1-dimensional array here. See   |
|                                     | | the description under the **constant** group for further details.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 3-dimensional array of size (``N_reactions`` x ``N_Temperatures`` x      |
|                                     | | ``N_Psi``) that tabulates the rate coefficients of each reaction as a      |
|                                     | | function of temperature and ``Psi``. These are given as log10(rate) in cgs |
|                                     | | units.                                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

H2_dust_formation
^^^^^^^^^^^^^^^^^

Formation of H2 on the surface of dust grains. The rate coefficient dependends on gas temperature and dust temperature. There is only one reaction in this group. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length 3 that gives the species indices of each   |
|                                     | | reactant in this reaction. Since there is only one reaction in this group, |
|                                     | | this array is only 1-dimensional. Also, this reaction only has 2           |
|                                     | | reactants, so the third index is -1 and is ignored by the code. See the    |
|                                     | | description under the **constant** group for further details.              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 1-dimensional array of length 3 that gives the species indices of the    |
|                                     | | product in this reaction. Since there is only one reaction in this group,  |
|                                     | | this array is only 1-dimensional. Also, this reaction only has 1 product,  |
|                                     | | so the second and third indices are -1 and are ignored by the code. See    |
|                                     | | the description under the **constant** group for further details.          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 2-dimensional array of size (``N_Temperatures`` x                        |
|                                     | | ``N_Dust_Temperatures``) that tabulates the rate coefficient of this       |
|                                     | | reaction as a function of gas temperature and dust temperature. These are  |
|                                     | | given as log10(rate) in cgs units. Note that, to get the reaction rate per |
|                                     | | unit volume, we need to multiply this rate coefficient by                  |
|                                     | | ``nHI * nHtot * dust_ratio``, and not by ``nHI^2``, where ``nHI`` and      |
|                                     | | ``nHtot`` are the HI and total hydrogen densities, respectively, and       |
|                                     | | ``dust_ratio`` is the dust-to-gas ratio relative to the Milky Way value.   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 1-dimensional array of length 9 that encodes which elements are involved |
|                                     | | in this reaction. See the description under the **constant** group for     |
|                                     | | details.                                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A single integer that can take the value of ``0`` (if no molecules are     |
|                                     | | involved in this reaction) or ``1`` (if molecules are involved). See the   |
|                                     | | description under the **constant** group for details.                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

H2_collis_dissoc
^^^^^^^^^^^^^^^^

Dissociation of H2 via collisions with electrons, HI and HeI. As well as temperature, the rate coefficients for these reactions also depend on the densities of HI, H2 and HeI, which determine whether the rate coefficient is in the low-density regime or LTE. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the         |
|                                     | | species indices of each reactant in each reaction. All reactions in this   |
|                                     | | group have exactly 2 reactants. See the description under the **constant** |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 3) that gives the         |
|                                     | | species indices of each product in each reaction. All reactions in this    |
|                                     | | group have exactly 3 products. See the description under the **constant**  |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``critical_density_H``              | | A 1-dimensional array of length ``N_Temperatures`` that tabulates the      |
|                                     | | critical density of H2 due to collisions with HI as a function of          |
|                                     | | temperature. These are given as log10(critical density) in cgs units.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``critical_density_H2``             | | A 1-dimensional array of length ``N_Temperatures`` that tabulates the      |
|                                     | | critical density of H2 due to collisions with H2 as a function of          |
|                                     | | temperature. These are given as log10(critical density) in cgs units.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``critical_density_He``             | | A 1-dimensional array of length ``N_Temperatures`` that tabulates the      |
|                                     | | critical density of H2 due to collisions with HI as a function of          |
|                                     | | temperature. These are given as log10(critical density) in cgs units.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``k0`` and ``kLTE``                 | | Two 2-dimensional arrays, each of size (``N_reactions`` x                  |
|                                     | | ``N_Temperatures``), that tabulate the low-density (``k0``) and LTE        |
|                                     | | (``kLTE``) rate coefficients for each reaction as a function of            |
|                                     | | temperature. These are given as log10(rate coefficient) in cgs units. The  |
|                                     | | overall rate coefficient for each reaction is then calculated as follows:  |
|                                     | | ``log10(rate coefficient) =``                                              |
|                                     | | ``(n_over_cr / (1 + n_over_cr)) * log10(kLTE) +``                          |
|                                     | | ``(1 / (1 + n_over_cr)) * log10(k0)``,                                     |
|                                     | | where:                                                                     |
|                                     | | ``n_over_cr = (nHI / critical_density_H) +``                               |
|                                     | | ``(2 * nH2 / critical_density_H2) + (nHeI / critical_density_He)``         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

cosmic_ray
^^^^^^^^^^

Ionisation and dissociation by cosmic rays. All rates are normalised relative to the cosmic ray ionisation rate of HI, which is given as a parameter in the gasVariables structure. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 3) that gives the         |
|                                     | | species indices of each product in each reaction. If a reaction has fewer  |
|                                     | | than 3 product, the remaining indices are set to ``-1``. No reactions      |
|                                     | | have more than 3 products. See the description under the **constant**      |
|                                     | | group for further details.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 1-dimensional array of length ``N_reactions``  that gives the rates of   |
|                                     | | each reaction relative to the HI cosmic ray ionisation rate. These are     |
|                                     | | stored as the linear, not the log, of this ratio. The rate of each         |
|                                     | | reaction in units of ``s^-1`` is then given by:                            |
|                                     | | ``rates * gasVariables.cr_rate``.                                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``secondary_ratio``                 | | A 2-dimensional array of size (2 x ``N_secondary_cosmic_ray_xHII``). This  |
|                                     | | gives the ratio of secondary to primary cosmic ray ionisation rate for HI  |
|                                     | | and HeI, tabulated as a function of the HII fraction, ``xHII``. The        |
|                                     | | corresponding cosmic ray ionisation rates are then multiplied by           |
|                                     | | ``(1 + secondary_ratio)``.                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``secondary_base_reaction``         | | A 1-dimensional array of length 2 that gives the positions of the HI and   |
|                                     | | HeI cosmic ray ionisation reactions in the **cosmic_ray** group. This      |
|                                     | | allows us to find those reactions so that we can update their rates to     |
|                                     | | include secondary ionisations, using the ``secondary_ratio`` array given   |
|                                     | | above.                                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

CO_cosmic_ray
^^^^^^^^^^^^^

Dissociation of CO by cosmic rays. This reaction has additional dependencies on the gas temperature and the H2 and CO abundances, and so is grouped separately from the rest of the cosmic ray reactions. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. In practice, there is only 1 reaction in this group.           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant. The only reaction in this group has 1 reactant.  |
|                                     | | See the description under the **constant** group for further details.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product. The only reaction in this group has 2 products.   |
|                                     | | See the description under the **constant** group for further details.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 2-dimensional array of size (``N_reactions`` x ``N_Temperatures``)  that |
|                                     | | gives the rate relative to the HI cosmic ray ionisation rate as a function |
|                                     | | of gas temperature, stored as log10(rate). Note that there are also extra  |
|                                     | | dependencies on the CO and H2 abundances (``xCO`` and ``xH2``,             |
|                                     | | respectively). The final rate per CO molecule, in units of ``s^-1``, is    |
|                                     | | then given by:                                                             |
|                                     | | ``rate * gasVariables.cr_rate * xH2 / sqrt(xCO)``.                         |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

photodissoc_group1
^^^^^^^^^^^^^^^^^^

Photodissociation of molecules (and also the negative ions ``C-`` and ``O-``) for which the rate is attenuated by dust extinction, ``Av``, simply as ``exp(-gamma * Av)``, where ``gamma`` is a parameter that we tabulate for  each reaction. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product in each reaction. All reactions in this group have |
|                                     | | exactly 2 products. See the description under the **constant** group for   |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gamma``                           | | A 1-dimensional array of length ``N_reactions`` that contains the gamma    |
|                                     | | parameter for each reaction, which is used for the dust attenuation (see   |
|                                     | | below).                                                                    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 1-dimensional array of length ``N_reactions`` that contains the          |
|                                     | | optically thin dissociation rates for each reaction, in units ``s^-1``.    |
|                                     | | These are stored as the linear rate, and NOT as the log. The reaction      |
|                                     | | rates then scale linearly with the strength of the radiation field in the  |
|                                     | | 6-13.6 eV band in Habing units, i.e. the ``G0`` parameter, and are         |
|                                     | | attenuated by dust extinction, ``Av``. The full rate is then given by:     |
|                                     | | ``rate * G0 * exp(-gamma Av)``.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

photodissoc_group2
^^^^^^^^^^^^^^^^^^

Photodissociation of molecules for which the rate is attenuated by dust extinction, ``Av``, with a more complex functional form as ``exp(-gamma_coeff[0] * Av)`` if ``Av > 15`` and ``exp(-gamma_coeff[1] * Av + gamma_coeff[2] * Av * Av)`` otherwise, where the ``gamma_coeff`` coefficients are the same for all reactions in this group (see below). 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product in each reaction. All reactions in this group have |
|                                     | | exactly 2 products. See the description under the **constant** group for   |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gamma_coeff``                     | | A 1-dimensional array of length 3 that contains the gamma coefficients     |
|                                     | | that are used to calculate the dust attenuation (see above). These are the |
|                                     | | same for all reactions in this group.                                      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 1-dimensional array of length ``N_reactions`` that contains the          |
|                                     | | optically thin dissociation rates for each reaction, in units ``s^-1``.    |
|                                     | | These are stored as the linear rate, and NOT as the log. The reaction      |
|                                     | | rates then scale linearly with the strength of the radiation field in the  |
|                                     | | 6-13.6 eV band in Habing units, i.e. the ``G0`` parameter, and are         |
|                                     | | attenuated by dust extinction, ``Av``. The full rate, ``R``, is then given |
|                                     | | by:                                                                        |
|                                     | | ``if (Av > 15)``                                                           |
|                                     | |   ``R = rate * G0 * exp(-gamma_coeff[0] * Av)``                            |
|                                     | | ``else``                                                                   |
|                                     | |   ``R = rate * G0 * exp(-gamma_coeff[1] * Av + gamma_coeff[2] * Av * Av)`` |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

H2_photodissoc
^^^^^^^^^^^^^^

Photodissociation of H2. Includes attenutation using the H2 self-shielding function of Richings et al. (2014b) and dust attenuation. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. In practice, there is only 1 reaction in this group.           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant. The only reaction in this group has 1 reactant.  |
|                                     | | See the description under the **constant** group for further details.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product. The only reaction in this group has 2 products.   |
|                                     | | See the description under the **constant** group for further details.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gamma``                           | | A 1-dimensional array of length ``N_reactions`` that contains the gamma    |
|                                     | | parameter, which is used for the dust attenuation (see below).             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``self_shielding``                  | | A 4-dimensional array of size ``(N_reactions x N_Temperatures x``          |
|                                     | | ``N_H2self_column_densities x N_b_turbulence)``, which tabulates the H2    |
|                                     | | self-shielding function with respect to gas temperature, H2 column density |
|                                     | | and turbulent broadening parameter (see Richings et al. 2014b).            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 1-dimensional array of length ``N_reactions`` that contains the          |
|                                     | | optically thin dissociation rate for this reaction, in units ``s^-1``.     |
|                                     | | This is stored as the linear rate, and NOT as the log. The reaction rate   |
|                                     | | then scales linearly with the strength of the radiation field in the       |
|                                     | | 12.24-13.51 eV band as parameterised by the ``H2_dissocJ`` parameter, and  |
|                                     | | is attenuated by dust extinction, ``Av``, and H2 self-shielding. The full  |
|                                     | | rate is then given by:                                                     |
|                                     | | ``rate * H2_dissocJ * isotropic_photon_density * speed_of_light``          |
|                                     | | ``* exp(-gamma Av) * S_H2``                                                |
|                                     | | where ``S_H2`` is calculated from the ``self_shielding`` array given above |
|                                     | | as a function of temperature, H2 column density and turbulent Doppler      |
|                                     | | broadening.                                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in this reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

CO_photodissoc
^^^^^^^^^^^^^^

Photodissociation of CO. Includes attenutation CO self-shielding, H2 cross-shielding, and dust shielding. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. In practice, there is only 1 reaction in this group.           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant. The only reaction in this group has 1 reactant.  |
|                                     | | See the description under the **constant** group for further details.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product. The only reaction in this group has 2 products.   |
|                                     | | See the description under the **constant** group for further details.      |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gamma``                           | | A 1-dimensional array of length ``N_reactions`` that contains the gamma    |
|                                     | | parameter, which is used for the dust attenuation (see below).             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``self_shielding``                  | | A 2-dimensional array of size ``(N_reactions x N_COself_column_densities`` |
|                                     | | ``x N_H2CO_column_densities``, which tabulates the shielding of CO by      |
|                                     | | itself and by H2 as a function of the CO and H2 column densities (see      |
|                                     | | Richings et al. 2014b for details).                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 1-dimensional array of length ``N_reactions`` that contains the          |
|                                     | | optically thin dissociation rate for this reaction, in units ``s^-1``.     |
|                                     | | This is stored as the linear rate, and NOT as the log. The reaction rate   |
|                                     | | then scales linearly with the strength of the radiation field in the       |
|                                     | | 6-13.6 eV  band in Habing units, i.e. the ``G0`` parameter, and is         |
|                                     | | attenuated by dust, CO and H2. The full rate is then given by:             |
|                                     | | ``rate * G0 * exp(-gamma Av) * S_CO``                                      |
|                                     | | where ``S_CO`` is calculated from the ``self_shielding`` array given above |
|                                     | | as a function of the CO and H2 column densities.                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in this reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+


photoion_fuv
^^^^^^^^^^^^

Photoionisation of species with an ionisation energy <13.6 eV. Additional information from the cross-sections tables for each UV spectrum will be needed to compute the rates of each reaction in this group; see the cross_sections section of this User Guide for details. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product in each reaction. All reactions in this group have |
|                                     | | exactly 2 products. See the description under the **constant** group for   |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gamma``                           | | A 1-dimensional array of length ``N_reactions`` that contains the gamma    |
|                                     | | parameter for each reaction. The attenuation of the photoionisation rate   |
|                                     | | due to dust extinction ``Av`` is then given as ``exp(-gamma * Av)``.       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``E_thresh``                        | | A 1-dimensional array of length ``N_reactions`` that contains the          |
|                                     | | ionisation energy of each reaction. Photons above this energy threshold    |
|                                     | | contribute to the given photoionisation reaction.                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+


photoion_euv
^^^^^^^^^^^^

Photoionisation of species with an ionisation energy >13.6 eV. Additional information from the cross-sections tables for each UV spectrum will be needed to compute the rates of each reaction in this group; see the cross_sections section of this User Guide for details. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product in each reaction. All reactions in this group have |
|                                     | | exactly 2 products. See the description under the **constant** group for   |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``E_thresh``                        | | A 1-dimensional array of length ``N_reactions`` that contains the          |
|                                     | | ionisation energy of each reaction. Photons above this energy threshold    |
|                                     | | contribute to the given photoionisation reaction.                          |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+


photoion_auger_fuv
^^^^^^^^^^^^^^^^^^

Auger photoionisation of species with an ionisation energy <13.6 eV, where a single photon is absorbed by an inner shell electron with enough energy to remove further electrons from the outer shells. Additional information from the cross-sections tables for each UV spectrum will be needed to compute the rates of each reaction in this group; see the cross_sections section of this User Guide for details. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product in each reaction. This array gives the products as |
|                                     | | the final ion and an electron. However, in practice multiple electrons are |
|                                     | | produced. The number of electrons from each reaction is calculated from    |
|                                     | | the difference between the initial and final ionisation states.            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``base_reaction``                   | | A 1-dimensional array of length ``N_reactions`` that contains the position |
|                                     | | of the base reaction (i.e. the photoionisation of the same reactant that   |
|                                     | | releases a single electron) in the **photoion_fuv** group.                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+


photoion_auger_euv
^^^^^^^^^^^^^^^^^^

Auger photoionisation of species with an ionisation energy >13.6 eV, where a single photon is absorbed by an inner shell electron with enough energy to remove further electrons from the outer shells. Additional information from the cross-sections tables for each UV spectrum will be needed to compute the rates of each reaction in this group; see the cross_sections section of this User Guide for details. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_reactions``                     | | A 1-dimensional array of length 2 containing the number of reactions in    |
|                                     | | this group. See the description under the **constant** group for details.  |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``reactants``                       | | A 1-dimensional array of length ``N_reactions`` that gives the species     |
|                                     | | indices of each reactant in each reaction. All reactions in this group     |
|                                     | | have only 1 reactant. See the description under the **constant** group for |
|                                     | | further details.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``products``                        | | A 2-dimensional array of size (``N_reactions`` x 2) that gives the species |
|                                     | | indices of each product in each reaction. This array gives the products as |
|                                     | | the final ion and an electron. However, in practice multiple electrons are |
|                                     | | produced. The number of electrons from each reaction is calculated from    |
|                                     | | the difference between the initial and final ionisation states.            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``base_reaction``                   | | A 1-dimensional array of length ``N_reactions`` that contains the position |
|                                     | | of the base reaction (i.e. the photoionisation of the same reactant that   |
|                                     | | releases a single electron) in the **photoion_euv** group.                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``element_incl``                    | | A 2-dimensional array of size (``N_reactions`` x 9) that encodes which     |
|                                     | | elements are involved in each reaction. See the description under the      |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``molecular_flag``                  | | A 1-dimensional array of integers of length ``N_reactions`` that can take  |
|                                     | | the value ``0`` (if no molecules are involved in the given reaction) or    |
|                                     | | ``1`` (if molecules are involved). See the description under the           |
|                                     | | **constant** group for details.                                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

cooling
^^^^^^^

The **cooling** group in ``chimes_main_data.hdf5`` contains the per-ion cooling and heating rates. These data arrays are described in detail below. 

+-------------------------------------+------------------------------------------------------------------------------+
| Data Array                          | Description                                                                  |
+=====================================+==============================================================================+
| ``N_coolants``                      | | The number of species for which the cooling rate per ion is simply a       |
|                                     | | function of temperature.                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``coolants``                        | | A 1-dimensional array of length ``N_coolants`` containing the species      |
|                                     | | indices for these coolants.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates``                           | | A 2-dimensional array of size (``N_coolants`` x ``N_Temperatures``)        |
|                                     | | containing the temperature-dependent cooling rates per ion for each        |
|                                     | | coolant. Stored as log10(rate) in cgs units.                               |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``N_coolants_2d``                   | | The number of species for which the cooling rate per ion is a function of  |
|                                     | | temperature and electron density at low temperatures. For these coolants,  |
|                                     | | the cooling rate does not simply scale linearly with the electron density, | 
|                                     | | so we cannot simply multiply by the electron density as we do for the      |
|                                     | | standard coolants. At high temperatures, we switch to purely temperature-  |
|                                     | | dependent rates.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``coolants_2d``                     | | A 1-dimensional array of length ``N_coolants_2d`` containing the species   |
|                                     | | indices for these coolants.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates_2d``                        | | A 3-dimensional array of size (``N_coolants_2d`` x                         |
|                                     | | ``N_cool_2d_Temperatures`` x ``N_cool_2d_ElectronDensities``) containing   |
|                                     | | the cooling rates per ion for these coolants as a function of temperature  |
|                                     | | and electron density. Stored as log10(rate) in cgs units. These are used   |
|                                     | | at low temperatures.                                                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates_hiT_2d``                    | | A 2-dimensional array of size (``N_coolants_2d`` x                         |
|                                     | | ``N_cool_2d_hiT_Temperatures``) containing the cooling rates per ion for   |
|                                     | | these coolants as a function of temperature only. Stored as log10(rate) in |
|                                     | | cgs units. These are used at high temperatures.                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``N_coolants_4d``                   | | The number of species for which the cooling rate per ion is a function of  |
|                                     | | temperature, HI density, electron density and HII density at low           |
|                                     | | temperatures. At high temperatures, we switch to purely temperature-       |
|                                     | | dependent rates.                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``coolants_4d``                     | | A 1-dimensional array of length ``N_coolants_4d`` containing the species   |
|                                     | | indices for these coolants.                                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates_4d``                        | | A 5-dimensional array of size (``N_coolants_4d`` x                         |
|                                     | | ``N_cool_4d_Temperatures`` x ``N_cool_4d_HIDensities`` x                   |
|                                     | | ``N_cool_4d_ElectronDensities`` x ``N_cool_4d_HIIDensities``) containing   | 
|                                     | |  the cooling rates per ion for these coolants as a function of temperature |
|                                     | | and the densities of HI, electrons and HII. Stored as log10(rate) in cgs   |
|                                     | | units. These are used at low temperatures.                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``rates_hiT_4d``                    | | A 2-dimensional array of size (``N_coolants_4d`` x                         |
|                                     | | ``N_cool_4d_hiT_Temperatures``) containing the cooling rates per ion for   |
|                                     | | these coolants as a function of temperature only. Stored as log10(rate) in |
|                                     | | cgs units. These are used at high temperatures.                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_cool_LTE``                     | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent H2 cooling rate in the LTE limit, used in the H2     |
|                                     | | cooling function in CHIMES (see Glover & Abel 2008; Richings et al 2014a). |
|                                     | | Stored as log10(rate) in cgs units.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_cool_lowDens_H2``              | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent H2 cooling rate in the low-density limit from        |
|                                     | | collisions with H2. Used in the H2 cooling function in CHIMES (see Glover  |
|                                     | | & Abel 2008; Richings et al 2014a). Stored as log10(rate) in cgs units.    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_cool_lowDens_HI``              | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent H2 cooling rate in the low-density limit from        |
|                                     | | collisions with HI. Used in the H2 cooling function in CHIMES (see Glover  |
|                                     | | & Abel 2008; Richings et al 2014a). Stored as log10(rate) in cgs units.    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_cool_lowDens_HII``             | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent H2 cooling rate in the low-density limit from        |
|                                     | | collisions with HII. Used in the H2 cooling function in CHIMES (see Glover |
|                                     | | & Abel 2008; Richings et al 2014a). Stored as log10(rate) in cgs units.    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_cool_lowDens_HeI``             | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent H2 cooling rate in the low-density limit from        |
|                                     | | collisions with HeI. Used in the H2 cooling function in CHIMES (see Glover |
|                                     | | & Abel 2008; Richings et al 2014a). Stored as log10(rate) in cgs units.    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``H2_cool_lowDens_elec``            | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent H2 cooling rate in the low-density limit from        |
|                                     | | collisions with electrons. Used in the H2 cooling function in CHIMES (see  |
|                                     | | Glover & Abel 2008; Richings et al 2014a). Stored as log10(rate) in cgs    |
|                                     | | units.                                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``CO_cool_rot_L0``                  | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent CO cooling rate from rotational transitions in the   |
|                                     | | low-density limit. Used in the CO cooling function in CHIMES (see Neufeld  |
|                                     | | & Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al.   |
|                                     | | 2014a). Stored as log10(rate) in cgs units.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``CO_cool_rot_Llte``                | | A 2-dimensional array of size (``N_mol_cool_Temperatures`` x               |
|                                     | | ``N_CO_cool_rot_ColumnDensities``)  containing the CO cooling rate from    |
|                                     | | rotational transitions in the LTE limit as a function of temperature and   |
|                                     | | effective CO column density. Used in the CO cooling function in CHIMES     |
|                                     | | (see Neufeld & Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010;      |
|                                     | | Richings et al. 2014a). Stored as log10(rate) in cgs units.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``CO_cool_rot_a``                   | | A 2-dimensional array of size (``N_mol_cool_Temperatures`` x               |
|                                     | | ``N_CO_cool_rot_ColumnDensities``)  containing the alpha parameter as a    |
|                                     | | function of temperature and effective CO column density. Used in the CO    |
|                                     | | cooling function in CHIMES (see Neufeld & Kaufman 1993; Neufeld et al.     |
|                                     | | 1995; Glover et al. 2010; Richings et al. 2014a). Stored as log10(rate) in |
|                                     | | cgs units.                                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``CO_cool_rot_nhalf``               | | A 2-dimensional array of size (``N_mol_cool_Temperatures`` x               |
|                                     | | ``N_CO_cool_rot_ColumnDensities``)  containing the nhalf parameter as a    |
|                                     | | function of temperature and effective CO column density. Used in the CO    |
|                                     | | cooling function in CHIMES (see Neufeld & Kaufman 1993; Neufeld et al.     |
|                                     | | 1995; Glover et al. 2010; Richings et al. 2014a). Stored as log10(rate) in |
|                                     | | cgs units.                                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``CO_cool_vib_L0``                  | | A 1-dimensional array of length ``N_mol_cool_Temperatures`` containing the |
|                                     | | temperature-dependent CO cooling rate from vibrational transitions in the  |
|                                     | | low-density limit. Used in the CO cooling function in CHIMES (see Neufeld  |
|                                     | | & Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al.   |
|                                     | | 2014a). Stored as log10(rate) in cgs units.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``CO_cool_vib_Llte``                | | A 2-dimensional array of size (``N_mol_cool_Temperatures`` x               |
|                                     | | ``N_CO_cool_vib_ColumnDensities``)  containing the CO cooling rate from    |
|                                     | | vibrational transitions in the LTE limit as a function of temperature and  |
|                                     | | effective CO column density. Used in the CO cooling function in CHIMES     |
|                                     | | (see Neufeld & Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010;      |
|                                     | | Richings et al. 2014a). Stored as log10(rate) in cgs units.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2O_cool_rot_``                 | | A 1-dimensional array of length ``N_H2O_cool_hiT_Temperatures`` containing |
| | ``hiT_L0``                        | | the temperature-dependent H2O cooling rate from rotational transitions in  |
|                                     | | the low-density limit, at high temperatures. Used in the H2O cooling       |
|                                     | | function in CHIMES (see Neufeld & Kaufman 1993; Neufeld et al. 1995;       |
|                                     | | Glover et al. 2010; Richings et al. 2014a). Stored as log10(rate) in cgs   |
|                                     | | units.                                                                     |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2O_cool_rot_``                 | | A 2-dimensional array of size (``N_H2O_cool_hiT_Temperatures`` x           |
| | ``hiT_Llte``                      | | ``N_H2O_cool_rot_ColumnDensities``)  containing the H2O cooling rate from  |
|                                     | | rotational transitions in the LTE limit as a function of temperature and   |
|                                     | | effective H2O column density, at high temperatures. Used in the H2O        |
|                                     | | cooling function in CHIMES (see Neufeld & Kaufman 1993; Neufeld et al.     |
|                                     | | 1995; Glover et al. 2010; Richings et al. 2014a). Stored as log10(rate) in |
|                                     | | cgs units.                                                                 |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2O_cool_rot_``                 | | A 2-dimensional array of size (``N_H2O_cool_hiT_Temperatures`` x           |
| | ``hiT_a``                         | | ``N_H2O_cool_rot_ColumnDensities``)  containing the alpha parameter as a   |
|                                     | | function of temperature and effective H2O column density, at high          |
|                                     | | temperatures. Used in the H2O cooling function in CHIMES (see Neufeld &    |
|                                     | | Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al.     |
|                                     | | 2014a). Stored as log10(rate) in cgs units.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2O_cool_rot_``                 | | A 2-dimensional array of size (``N_H2O_cool_hiT_Temperatures`` x           |
| | ``hiT_nhalf``                     | | ``N_H2O_cool_rot_ColumnDensities``)  containing the nhalf parameter as a   |
|                                     | | function of temperature and effective H2O column density, at high          |
|                                     | | temperatures. Used in the H2O cooling function in CHIMES (see Neufeld &    |
|                                     | | Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al.     |
|                                     | | 2014a). Stored as log10(rate) in cgs units.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Oortho_cool_``                | | A 1-dimensional array of length ``N_H2O_cool_lowT_Temperatures``           |
| | ``rot_lowT_L0``                   | | containing the temperature-dependent ortho-H2O cooling rate from           |
|                                     | | rotational transitions in the low-density limit, at low temperatures.      |
|                                     | | Used in the H2O cooling function in CHIMES (see Neufeld & Kaufman          |
|                                     | | 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al. 2014a).     |
|                                     | | Stored as log10(rate) in cgs units.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Oortho_cool_``                | | A 2-dimensional array of size (``N_H2O_cool_lowT_Temperatures`` x          |
| | ``rot_lowT_Llte``                 | | ``N_H2O_cool_rot_ColumnDensities``)  containing the ortho-H2O cooling      |
|                                     | | rate from rotational transitions in the LTE limit as a function of         |
|                                     | | temperature and effective H2O column density, at low temperatures.         |
|                                     | | Used in the H2O cooling function in CHIMES (see Neufeld & Kaufman          |
|                                     | | 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al. 2014a).     |
|                                     | | Stored as log10(rate) in cgs units.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Oortho_cool_``                | | A 2-dimensional array of size (``N_H2O_cool_lowT_Temperatures`` x          |
| | ``rot_lowT_a``                    | | ``N_H2O_cool_rot_ColumnDensities``)  containing the alpha parameter as a   |
|                                     | | function of temperature and effective H2O column density, at low           |
|                                     | | temperatures. Used in the H2O cooling function in CHIMES (see              |
|                                     | | Neufeld & Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010;           |
|                                     | | Richings et al. 2014a). Stored as log10(rate) in cgs units.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Oortho_cool_``                | | A 2-dimensional array of size (``N_H2O_cool_lowT_Temperatures`` x          |
| | ``rot_lowT_nhalf``                | | ``N_H2O_cool_rot_ColumnDensities``)  containing the nhalf parameter as a   |
|                                     | | function of temperature and effective H2O column density, at low           |
|                                     | | temperatures. Used in the H2O cooling function in CHIMES (see              |
|                                     | | Neufeld & Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010;           |
|                                     | | Richings et al. 2014a). Stored as log10(rate) in cgs units.                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Opara_cool_``                 | | A 1-dimensional array of length ``N_H2O_cool_lowT_Temperatures``           |
| | ``rot_lowT_L0``                   | | containing the temperature-dependent para-H2O cooling rate from            |
|                                     | | rotational transitions in the low-density limit, at low temperatures.      |
|                                     | | Used in the H2O cooling function in CHIMES (see Neufeld &                  |
|                                     | | Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et         |
|                                     | | al. 2014a). Stored as log10(rate) in cgs units.                            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Opara_cool_``                 | | A 2-dimensional array of size (``N_H2O_cool_lowT_Temperatures`` x          |
| | ``rot_lowT_Llte``                 | | ``N_H2O_cool_rot_ColumnDensities``)  containing the para-H2O cooling       |
|                                     | | rate from rotational transitions in the LTE limit as a function of         |
|                                     | | temperature and effective H2O column density, at low temperatures.         |
|                                     | | Used in the H2O cooling function in CHIMES (see Neufeld & Kaufman          |
|                                     | | 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al. 2014a).     |
|                                     | | Stored as log10(rate) in cgs units.                                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Opara_cool_``                 | | A 2-dimensional array of size (``N_H2O_cool_lowT_Temperatures`` x          |
| | ``rot_lowT_a``                    | | ``N_H2O_cool_rot_ColumnDensities``)  containing the alpha parameter as a   |
|                                     | | function of temperature and effective H2O column density, at low           |
|                                     | | temperatures. Used in the H2O cooling function in CHIMES (see Neufeld &    |
|                                     | | Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al.     |
|                                     | | 2014a). Stored as log10(rate) in cgs units.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``H2Opara_cool_``                 | | A 2-dimensional array of size (``N_H2O_cool_lowT_Temperatures`` x          |
| | ``rot_lowT_nhalf``                | | ``N_H2O_cool_rot_ColumnDensities``)  containing the nhalf parameter as a   |
|                                     | | function of temperature and effective H2O column density, at low           |
|                                     | | temperatures. Used in the H2O cooling function in CHIMES (see Neufeld &    |
|                                     | | Kaufman 1993; Neufeld et al. 1995; Glover et al. 2010; Richings et al.     |
|                                     | | 2014a). Stored as log10(rate) in cgs units.                                |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``gas_grain_transfer``              | | A 1-dimensional array of length ``N_Temperatures`` containing the cooling  |
|                                     | | rate due to energy exchange between the gas and the dust, as a function    |
|                                     | | of gas temperature, ``T_gas``. Stored as log10(rate) in cgs units. To get  |
|                                     | | the overall cooling rate, we also need the dust-to-gas ratio relative to   |
|                                     | | the local ISM (``dust_ratio``) and the temperature difference between      |
|                                     | | ``T_gas`` and the dust temperature ``T_dust``. The final rate is then:     |
|                                     | | ``rate * dust_ratio * (T_gas - T_dust)``.                                  |
|                                     | | See Richings et al. (2014a) for details.                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``grain_recombination``             | | A 2-dimensional array of size (``N_Temperatures`` x ``N_Psi``) containing  |
|                                     | | the grain recombination cooling rates as a function of gas temperature     |
|                                     | | ``T`` and the parameter ``Psi = G0 exp(-2.77 Av) T^0.5 / ne``, where       |
|                                     | | ``G0`` is the strength of the radiation field in the 6-13.6 eV band in     |
|                                     | | Habing units, ``Av`` is the dust extinction and ``ne`` is the electron     |
|                                     | | density. Stored as log10(rate) in cgs units. To get the overall cooling    |
|                                     | | rate, we also need the dust-to-gas ratio relative to the local ISM         |
|                                     | | (``dust_ratio``). The final rate is then:                                  |
|                                     | | ``rate * dust_ratio * ne``.                                                |
|                                     | | (See Glover & Jappsen 2007; Richings et al. 2014a).                        |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``photoelectric_heating``           | | A 2-dimensional array of size (``N_Temperatures`` x ``N_Psi``) containing  |
|                                     | | the photoelectric heating rates as a function of gas temperature ``T``     |
|                                     | | and the parameter ``Psi = G0 exp(-2.77 Av) T^0.5 / ne``, where ``G0`` is   |
|                                     | | the strength of the radiation field in the 6-13.6 eV band in Habing units, |
|                                     | | ``Av`` is the dust extinction and ``ne`` is the electron density. Stored   |
|                                     | | as log10(rate) in cgs units. To get the overall heating rate, we also need |
|                                     | | the dust-to-gas ratio relative to the local ISM (``dust_ratio``). The      |
|                                     | | final rate is then:                                                        |
|                                     | | ``rate * dust_ratio * G0 * exp(-2.77 Av) / nHtot``.                        |
|                                     | | (See Bakes & Tielens 1994; Wolfire et al. 2003; Richings et al. 2014a).    |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

TableBins
^^^^^^^^^

The **TableBins** group in ``chimes_main_data.hdf5`` contains all of the table bins that are used to tabulate and interpolate the various data arrays described above. 

References
^^^^^^^^^^

| `Bakes & Tielens (1994) <https://ui.adsabs.harvard.edu/abs/1994ApJ...427..822B>`_ 
| `Glover & Jappsen (2007) <https://ui.adsabs.harvard.edu/abs/2007ApJ...666....1G>`_ 
| `Glover & Abel (2008) <https://ui.adsabs.harvard.edu/abs/2008MNRAS.388.1627G>`_ 
| `Glover et al. (2010) <https://ui.adsabs.harvard.edu/abs/2010MNRAS.404....2G>`_
| `Neufeld & Kaufman (1993) <https://ui.adsabs.harvard.edu/abs/1993ApJ...418..263N>`_
| `Neufeld et al. (1995) <https://ui.adsabs.harvard.edu/abs/1995ApJS..100..132N>`_
| `Richings et al. (2014a) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.440.3349R>`_
| `Richings et al. (2014b) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.442.2780R>`_
