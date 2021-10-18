.. Read Chimes 
   Alexander Richings, 6th March 2020 

.. _ReadChimes_label: 

Read Chimes
-----------

The ``read_chimes.py`` python module contains various dictionaries and functions that can be used to parse the abundance arrays that are output from the CHIMES module. These can be used to help the User find the abundance of a particular species in the abundance array, whether they are using the full CHIMES network or a reduced network with some of the elements switched off. 

There are two dictionaries that define the names of the species in the full CHIMES network: 

+-------------------------------------+------------------------------------------------------------------------------+
| Dictionary                          | Description                                                                  |
+=====================================+==============================================================================+
| ``chimes_dict``                     | | This dictionary gives a series of ``{Key : Value}`` pairs where the        |
|                                     | | ``Key`` is the name of the species and the ``Value`` is the position of    |
|                                     | | that species in the full CHIMES abundance array. So the User can use       |
|                                     | | ``chimes_dict["<species_name>"]`` to get the position of the species       |
|                                     | | called ``"<species_name>"`` in the abundance array. Note that the name     |
|                                     | | given by the User must exactly match the name defined in this dictionary.  |
|                                     | | To see a list of all possible species names, you can do:                   |
|                                     | | ``for k, v in sorted(chimes_dict.items(), key = lambda x: x[1]):``         |
|                                     | |     ``print(k)``                                                           |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| ``chimes_dict_inv``                 | | This gives the inverse of the ``chimes-dict`` dictionary, i.e. where the   |
|                                     | | ``Key`` is the position in the full CHIMES abundance array and the         |
|                                     | | ``Value`` is the name of the species at that position. So the User can     |
|                                     | | use ``chimes_dict_inv[<position>]`` to find the name of the species that   |
|                                     | | is located at ``<position>`` in the abundance array.                       |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+

The following functions can be used to construct dictionaries mapping species names to positions in the abundance array for a reduced network, and to read in the abundances of a given species from an HDF5 snapshot file: 

+-------------------------------------+------------------------------------------------------------------------------+
| Dictionary                          | Description                                                                  |
+=====================================+==============================================================================+
| | ``create_reduced_``               | | Creates a dictionary mapping species names to their position in the        |
| | ``chimes_dictionary(``            | | CHIMES abundance array when using a reduced CHIMES network, i.e. when      |
| | ``element_flags =``               | | some of the elements have been switched off. It takes as an optional       |
| | ``np.zeros(9))``                  | | argument an array of integers of length 9, ``element_flags``, that can     |
|                                     | | be ``0`` or    ``1`` that encode whether each metal element is included    |
|                                     | | in the network. These are given in the order: C, N, O, Ne, Mg, Si, S,      |
|                                     | | Ca, Fe. Note that H and He are always included in the network. If no       |
|                                     | | ``element_flags`` array is given to the function, then it defaults to      |
|                                     | | the primordial network with only H and He, i.e. with all metal elements    |
|                                     | | switched off.                                                              |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``create_reduced_``               | | Creates the inverse dictionary, i.e. mapping the species position in the   |
| | ``inverse_chimes_``               | | CHIMES abundance array to the corresponding species name, when using a     |
| | ``dictionary(``                   | | reduced CHIMES network, i.e. when some of the elements have been switched  |
| | ``element_flags =``               | | off. It takes as an optional argument an array of integers of length 9,    |
| | ``np.zeros(9))``                  | | ``element_flags``, that encode which elements are included in the network  |
|                                     | | (see above for details).                                                   |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``read_chimes(``                  | | Takes the HDF5 snapshot file specified by ``filename`` and reads the       |
| | ``filename,``                     | | abundances of the species given by the string ``chimes_species``. The      |
| | ``chimes_species,``               | | optional argument ``chimes_dataset`` gives the name of the dataset in the  |
| | ``chimes_dataset =``              | | HDF5 snapshot file containing the 2-dimensional CHIMES abundance array.    |
| | ``"PartType0/``                   | | If this is not specified by the User, it defaults to                       |
| | ``Abundances")``                  | | ``PartType0/ChimesAbundances``. Returns a 1-dimensional array containing   |
|                                     | | the abundance of ``chimes_species`` for each gas particle in the snapshot. |
|                                     | | This routine assumes that we are using the full CHIMES network.            |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+
| | ``read_reduced_``                 | | As ``read_chimes()``, but for when we are using a reduced network, i.e.    |
| | ``chimes(filename,``              | | with some of the elements switched off. The ``filename``,                  |
| | ``chimes_species,``               | | ``chimes_species`` and ``chimes_dataset`` arguments are as described       |
| | ``chimes_dataset =``              | | above. The optional ``element_flags`` argument gives an array of integers  |
| | ``"PartType0/``                   | | that encode which elements are included in the reduced network. See the    |
| | ``Abundances",``                  | | ``create_reduced_chimes_dictionary()`` description above for details. If   |
| | ``element_flags =``               | | this is not specified by the User, it defaults to the primordial network   |
| | ``np.zeros(9))``                  | | with only H and He, i.e. with all metal elements switched off.             |
|                                     |                                                                              |
+-------------------------------------+------------------------------------------------------------------------------+


