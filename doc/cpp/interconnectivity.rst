.. _cppinterconnectivity:

.. cpp:namespace:: arb

Interconnectivity
#################

.. cpp:class:: cell_connection

    Describes a connection between two cells: a pre-synaptic source and a
    post-synaptic destination. The source is typically a threshold detector on
    a cell or a spike source. The destination is a synapse on the post-synaptic cell.

    The :cpp:member:`dest` does not include the gid of a cell, this is because a
    :cpp:class:`cell_connection` is bound to the destination cell which means that the gid
    is implicitly known.

    .. cpp:member:: cell_global_label_type source

        Source end point, represented by a :cpp:type:`cell_global_label_type` which packages
        a cell gid, label of a group of sources on the cell, and source selection policy.

    .. cpp:member:: cell_local_label_type dest

        Destination end point on the cell, represented by a :cpp:type:`cell_local_label_type`
        which packages a label of a group of targets on the cell and a selection policy.
        The target cell's gid is implicitly known.

    .. cpp:member:: float weight

        The weight delivered to the target synapse.
        The weight is dimensionless, and its interpretation is
        specific to the synapse type of the target. For example,
        the `expsyn` synapse interprets it as a conductance
        with units μS (micro-Siemens).

    .. cpp:member:: float delay

        Delay of the connection (milliseconds).

.. cpp:class:: gap_junction_connection

    Describes a gap junction connection between two gap junction sites. The :cpp:member:`local` site does
    not include the gid of a cell, this is because a :cpp:class:`gap_junction_connection` is bound to the local
    cell which means that the gid is implicitly known.

    .. note::

       A bidirectional gap-junction connection between two cells ``c0`` and ``c1`` requires two
       :cpp:class:`gap_junction_connection` objects to be constructed: one where ``c0`` is the
       :cpp:member:`local` site, and ``c1`` is the :cpp:member:`peer` site; and another where ``c1`` is the
       :cpp:member:`local` site, and ``c0`` is the :cpp:member:`peer` site.

    .. cpp:member:: cell_global_label_type peer

        Peer gap junction site, represented by a :cpp:type:`cell_local_label_type` which packages a cell gid,
        a label of a group of gap junction sites on the cell, and a site selection policy.

    .. cpp:member:: cell_local_label_type local

        Local gap junction site on the cell, represented by a :cpp:type:`cell_local_label_type`
        which packages a label of a group of gap junction sites on the cell and a selection policy.
        The gid of the local site's cell is implicitly known.

    .. cpp:member:: float weight

        unit-less gap junction connection weight.


High-Level Network Description
------------------------------


.. cpp:function:: network_population unique(const network_population& pop)

    Removes any duplicate global labels contained by merging overlapping intervals for matching local labels.

.. cpp:function:: join(const network_population& a, network_population b)

    Combine network populations. Does not remove duplicates.

.. cpp:function:: template <typename... ARGS> network_population join(const network_population& a, network_population b, ARGS&&... args)

    Combine more than two network populations. Does not remove duplicates.

.. cpp:class:: network_selection

    Selects or rejects a connection or gap junction when queried.

    .. cpp:function:: static network_selection bernoulli_random(unsigned seed, double p)

       Random selection using the bernoulli random distribution with given probability.

    .. cpp:function:: static network_selection custom(std::function<bool(const cell_global_label_type&, const cell_global_label_type&)> func)

       A custom selection using the provided function ``func`` ,
       which must accept a ``src`` and ``dest`` ``cell_global_label_type`` and return true if the connection / gap junction should be generated,
       and false otherwise. Repeated calls with the same arguments must give the same result. For gap junction, ``func`` must also be symmetric.

    .. cpp:function:: static network_selection all()

       Select all connections / gap junctions.

    .. cpp:function:: static network_selection none()

       Select no connections / gap junctions.

    .. cpp:function:: static network_selection invert(network_selection s)

       Invert the given selection, by accepting any connection not accepted by ``s`` and vice versa.

    .. cpp:function:: static network_selection inter_cell()

       Select only connections / gap junctions between different cells.

    .. cpp:function:: static network_selection not_equal()

       Select only connections / gap junctions between different items.
       Items identified by a different local cell label on the same cell are selected.

    .. cpp:function:: bool operator()(const cell_global_label_type& src, const cell_global_label_type& dest) const

       Accept or reject connection or gap junction between ``src`` and ``dest``.

    .. cpp:function:: network_selection operator&(network_selection right) const

       Logical "and" operation between two selections.

    .. cpp:function:: network_selection operator|(network_selection right) const

       Logical "or" operation between two selections.

    .. cpp:function:: network_selection operator^(network_selection right) const

       Logical "xor" operation between two selections.


.. cpp:class:: network_value

    A value used in network generation.

    .. cpp:function:: network_value(double value)

       Constant uniform value for all connections / gap junctions.

    .. cpp:function:: static network_value uniform(double value)

       Constant uniform value for all connections / gap junctions.

    .. cpp:function:: static network_value uniform_distribution(unsigned seed, const std::array<double, 2>& range)

       Uniform random distribution defined on the half open interval (range[0], range[1]].

    .. cpp:function:: static network_value normal_distribution(unsigned seed, double mean, double std_deviation)

       Normal random distribution with given mean and standard deviation.

    .. cpp:function:: static network_value truncated_normal_distribution(unsigned seed, double mean, double std_deviation, const std::array<double, 2>& range)

       Truncated normal random distribution with given mean and standard deviation, truncated to (range[0], range[1]] by accept-reject sampling.
       A low acceptance rate can leed to poor performance, for example with a very small range or a mean far outside the range.

    .. cpp:function:: static network_value custom(std::function<double(const cell_global_label_type&, const cell_global_label_type&)> func)

       A custom value using the provided function ``func``,
       which must accept a ``src`` and ``dest`` ``cell_global_label_type`` and return a float value for the connection / gap junction.
       Repeated calls with the same arguments must give the same result. For gap junction, ``func`` must also be symmetric.

    .. cpp:function:: double operator()(const cell_global_label_type& src, const cell_global_label_type& dest) const

       Generate value for connection or gap junction between ``src`` and ``dest``.


.. cpp:class:: cell_connection_network

    A network of cell connections.

    .. cpp:function:: cell_connection_network()

       Create an empty network.

    .. cpp:function:: cell_connection_network(network_value weight, network_value delay, network_selection selection, network_population src_pop, network_population dest_pop)

       Create a network between the source population ``src_pop`` and destination population ``dest_pop`` based on the given selection.

    .. cpp:function:: std::vector<cell_connection> generate(cell_gid_type gid) const

       Generate cell connections for cell with given ``gid``.

    .. cpp:function:: const network_value& weight() const

       The weight used in generated cell connections.

    .. cpp:function:: const network_value& delay() const

       The delay used in generated cell connections.

    .. cpp:function:: const network_selection& selection() const

       The network selection used to generate cell connections.

    .. cpp:function:: const network_population& source_population() const

       The source population used to generate cell connections.

    .. cpp:function:: const network_population& destination_population() const

       The destination population used to generate cell connections.


.. cpp:class:: gap_junction_network

    A network of gap junctions.

    .. cpp:function:: gap_junction_network()

       Create an empty network.

    .. cpp:function:: gap_junction_network(network_value weight, network_selection selection, network_population src_pop, network_population dest_pop)

       Create a network between the source population ``src_pop`` and destination population ``dest_pop`` based on the given selection.

    .. cpp:function:: std::vector<gap_junction_connection> generate(cell_gid_type gid) const

       Generate cell connections for cell with given ``gid``.

    .. property:: const network_value& weight() const

       The weight used in generated cell connections.

    .. property:: const network_selection& selection() const

       The network selection used to generate cell connections.

    .. property:: const network_population& source_population() const

       The source population used to generate cell connections.

    .. property:: const network_population& destination_population() const

       The destination population used to generate cell connections.
