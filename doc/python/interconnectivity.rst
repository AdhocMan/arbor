.. _pyinterconnectivity:

Interconnectivity
#################

.. currentmodule:: arbor

.. class:: connection

    Describes a connection between two cells, defined by source and destination end points (that is pre-synaptic and
    post-synaptic respectively), a connection weight and a delay time.

    The :attr:`dest` does not include the gid of a cell, this is because a :class:`arbor.connection` is bound to the
    destination cell which means that the gid is implicitly known.

    .. function:: connection(source, destination, weight, delay)

        Construct a connection between the :attr:`source` and the :attr:`dest` with a :attr:`weight` and :attr:`delay`.

    .. attribute:: source

        The source end point of the connection (type: :class:`arbor.cell_global_label`, which can be initialized with a
        (gid, label) or a (gid, (label, policy)) tuple. If the policy is not indicated, the default
        :attr:`arbor.selection_policy.univalent` is used).

    .. attribute:: dest

        The destination end point of the connection (type: :class:`arbor.cell_local_label` representing the label of the
        destination on the cell, which can be initialized with just a label, in which case the default
        :attr:`arbor.selection_policy.univalent` is used, or a (label, policy) tuple). The gid of the cell is
        implicitly known.

    .. attribute:: weight

        The weight delivered to the target synapse. It is up to the target mechanism to interpret this quantity.
        For Arbor-supplied point processes, such as the ``expsyn`` synapse, a weight of ``1`` corresponds to an
        increase in conductivity in the target mechanism of ``1`` μS (micro-Siemens).

    .. attribute:: delay

        The delay time of the connection [ms]. Must be positive.

    .. note::

        An minimal full example of a connection reads as follows:
        (see :ref:`network tutorial <tutorialnetworkring>` for a more comprehensive example):

        .. code-block:: python

            import arbor

            # Create two locset labels, describing the endpoints of the connection.
            labels = arbor.label_dict()
            labels['synapse_site'] = '(location 1 0.5)'
            labels['root'] = '(root)'

            # Place 'expsyn' mechanism on "synapse_site", and a threshold detector at "root"
            decor = arbor.decor()
            decor.place('"synapse_site"', 'expsyn', 'syn')
            decor.place('"root"', arbor.threshold_detector(-10), 'detector')

            # Implement the connections_on() function on a recipe as follows:
            def connections_on(gid):
               # construct a connection from the "detector" source label on cell 2
               # to the "syn" target label on cell gid with weight 0.01 and delay of 10 ms.
               src  = (2, "detector") # gid and locset label of the source
               dest = "syn" # gid of the destination is determined by the argument to `connections_on`.
               w    = 0.01  # weight of the connection. Correspondes to 0.01 μS on expsyn mechanisms
               d    = 10 # delay in ms
               return [arbor.connection(src, dest, w, d)]

.. class:: gap_junction_connection

    Describes a gap junction between two gap junction sites.

    The :attr:`local` site does not include the gid of a cell, this is because a :class:`arbor.gap_junction_connection`
    is bound to the destination cell which means that the gid is implicitly known.

    .. note::

       A bidirectional gap-junction between two cells ``c0`` and ``c1`` requires two
       :class:`gap_junction_connection` objects to be constructed: one where ``c0`` is the
       :attr:`local` site, and ``c1`` is the :attr:`peer` site; and another where ``c1`` is the
       :attr:`local` site, and ``c0`` is the :attr:`peer` site.

    .. function::gap_junction_connection(peer, local, weight)

        Construct a gap junction connection between :attr:`peer` and :attr:`local` with weight :attr:`weight`.

    .. attribute:: peer

        The gap junction site: the remote half of the gap junction connection (type: :class:`arbor.cell_global_label`,
        which can be initialized with a (gid, label) or a (gid, label, policy) tuple. If the policy is not indicated,
        the default :attr:`arbor.selection_policy.univalent` is used).

    .. attribute:: local

        The gap junction site: the local half of the gap junction connection (type: :class:`arbor.cell_local_label`
        representing the label of the destination on the cell, which can be initialized with just a label, in which case
        the default :attr:`arbor.selection_policy.univalent` is used, or a (label, policy) tuple). The gid of the
        cell is implicitly known.

    .. attribute:: weight

        The unit-less weight of the gap junction connection.

.. class:: threshold_detector

    A threshold detector, generates a spike when voltage crosses a threshold. Can be used as source endpoint for an
    :class:`arbor.connection`.

    .. attribute:: threshold

        Voltage threshold of threshold detector [mV]



High-Level Network Description
------------------------------


.. py:function:: unique(pop: list[cell_global_range_label]) -> list[cell_global_range_label]

    Removes any duplicate global labels contained by merging overlapping intervals for matching local labels.

.. py:class:: network_selection

    Selects or rejects a connection or gap junction when queried.

    .. py:method:: bernoulli_random(seed: unsigned int, p: float) -> network_selection
       :staticmethod:

       Random selection using the bernoulli random distribution with given probability.

    .. py:method:: custom(func) -> network_selection
       :staticmethod:

       A custom selection using the provided function ``func`` with signature ``(src: cell_global_label, dest: cell_global_label) -> bool``,
       which must accept a ``src`` and ``dest`` ``cell_global_label`` and return true if the connection / gap junction should be generated,
       and false otherwise. Repeated calls with the same arguments must give the same result. For gap junction, ``func`` must also be symmetric.

    .. py:method:: all() -> network_selection
       :staticmethod:

       Select all connections / gap junctions.

    .. py:method:: none() -> network_selection
       :staticmethod:

       Select no connections / gap junctions.

    .. py:method:: invert(s: network_selection) -> network_selection
       :staticmethod:

       Invert the given selection, by accepting any connection not accepted by ``s`` and vice versa.

    .. py:method:: inter_cell() -> network_selection
       :staticmethod:

       Select only connections / gap junctions between different cells.

    .. py:method:: not_equal() -> network_selection
       :staticmethod:

       Select only connections / gap junctions between different items.
       Items identified by a different local celll label on the same cell are selected.

    .. py:method:: __call__(src: cell_global_label, dest: cell_global_label) -> bool
       :noindex:

       Accept or reject connection or gap junction between ``src`` and ``dest``.

    .. py:method:: __and__(other: network_selection) -> network_selection

       Logical "and" operation between two selections.

    .. py:method:: __or__(other: network_selection) -> network_selection

       Logical "or" operation between two selections.

    .. py:method:: __xor__(other: network_selection) -> network_selection

       Logical "xor" operation between two selections.


.. py:class:: network_value

    A value used in network generation.

    .. py:method:: uniform_distribution(seed: unsigned int, range: [float, float]) -> network_value
       :staticmethod:

       Uniform random distribution defined on the half open interval (range[0], range[1]].

    .. py:method:: normal_distribution(seed: unsigned int, mean: float, std_deviation: float) -> network_value
       :staticmethod:

       Normal random distribution with given mean and standard deviation.

    .. py:method:: truncated_normal_distribution(seed: unsigned int, mean: float, std_deviation: float, range: [float, float]) -> network_value
       :staticmethod:

       Truncated normal random distribution with given mean and standard deviation, truncated to (range[0], range[1]] by accept-reject sampling.
       A low acceptance rate can leed to poor performance, for example with a very small range or a mean far outside the range.

    .. py:method:: custom(func) -> network_value
       :staticmethod:

       A custom value using the provided function ``func`` with signature ``(src: cell_global_label, dest: cell_global_label) -> float``,
       which must accept a ``src`` and ``dest`` ``cell_global_label`` and return a float value for the connection / gap junction.
       Repeated calls with the same arguments must give the same result. For gap junction, ``func`` must also be symmetric.

    .. py:method:: __call__(src: cell_global_label, dest: cell_global_label) -> float
       :noindex:

       Generate value for connection or gap junction between ``src`` and ``dest``.


.. py:class:: cell_connection_network

    A network of cell connections.

    .. py:method:: __init__(weight: float | network_value, delay: float | network_value, selection: network_selection, src_pop: list[cell_global_range_label], dest_pop: list[cell_global_range_label])

       Create a network between the source population ``src_pop`` and destination population ``dest_pop`` based on the given selection.

    .. py:method:: generate(gid: int) -> list[connection]

       Generate cell connections for cell with given ``gid``.

    .. property:: weight

       The weight used in generated cell connections.

    .. property:: delay

       The delay used in generated cell connections.

    .. property:: selection

       The network selection used to generate cell connections.

    .. property:: source_population

       The source population used to generate cell connections.

    .. property:: destination_population

       The destination population used to generate cell connections.


.. py:class:: gap_junction_network

    A network of gap junctions.

    .. py:method:: __init__(weight: float | network_value, selection: network_selection, src_pop: list[cell_global_range_label], dest_pop: list[cell_global_range_label])

       Create a network between the source population ``src_pop`` and destination population ``dest_pop`` based on the given selection.

    .. py:method:: generate(gid: int) -> list[gap_junction_connection]

       Generate cell connections for cell with given ``gid``.

    .. property:: weight

       The weight used in generated cell connections.

    .. property:: selection

       The network selection used to generate cell connections.

    .. property:: source_population

       The source population used to generate cell connections.

    .. property:: destination_population

       The destination population used to generate cell connections.
