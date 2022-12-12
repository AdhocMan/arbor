.. _interconnectivity:

Interconnectivity
=================

Networks can be regarded as graphs, where the nodes are locations on cells and the edges
describe the communications between them. In Arbor, two sorts of edges are modelled: a
:term:`connection` abstracts the propagation of action potentials (:term:`spikes <spike>`)
through the network, while a :term:`gap junction connection` is used to describe a direct
electrical connection between two locations on two cells.
Connections only capture the propagation delay and attenuation associated with spike
connectivity: the biophysical modelling of the chemical synapses themselves is the
responsibility of the target cell model.

Connection sites and gap junction sites are defined on locations on cells as
part of the :ref:`cell description <modelcelldesc>`.
These sites as such are not connected yet, however the :ref:`recipe <modelrecipe>`
exposes a number of callbacks to form connections and gap junctions between sites.
The recipe callbacks are interrogated during simulation creation.

In addition, simulations may update their connectivity by building a new
connection table outside calls to `run`, for example

.. code-block:: python

    rec = recipe()
    dec = arb.domain_decomposition(rec, ctx)
    sim = arb.simulation(rec, ctx, dec)

    # run simulation for 0.25ms with the basic connectivity
    sim.run(0.25, 0.025)

    # extend the recipe to more connections
    rec.add_connections()
    #  use `connections_on` to build a new connection table
    sim.update_connections(rec)

    # run simulation for 0.25ms with the extended connectivity
    sim.run(0.5, 0.025)

This will completely replace the old table, previous connections to be retained
must be explicitly included in the updated callback. This can also be used to
update connection weights and delays. Note, however, that there is currently no
way to introduce new sites to the simulation, nor any changes to gap junctions.

The ``update_connections`` method accepts either a full ``recipe`` (but will
**only** use the ``connections_on`` and ``events_generators`` callbacks) or a
``connectivity``, which is a reduced recipe exposing only the relevant callbacks.
Currently ``connectivity`` is only available in C++; Python users have to pass a
full recipe.

.. warning::

   The semantics of connection updates are subtle and might produce surprising
   results if handled carelessly. In particular, spikes in-flight over a
   connection will *always* be delivered, even if the connection has been
   deleted before the time of delivery has passed (`= t_emitted +
   connection_delay`). As Arbor's connection model joins processes on the axon,
   the synaptic cleft, and the receiving synapse into a simple pair `(weight,
   delay)` it is unclear 'where' the action potential is located at the time of
   deletion relative to the locus of disconnection. Thus, it was decided to
   deliver spike events regardless. This is will not cause issues when the
   transition is slow and smooth, ie weights decays over time towards a small
   value and then the connection is removed. However, drastic and/or frequent
   changes across busy synapses might cause unexpected behaviour.

.. note::

   Arbor uses a lazily constructed network (from the ``recipe`` callbacks) for
   good reason; storing the full connectivity (for all ``gids``) in the
   ``recipe`` can lead to prohibitively large memory footprints. Keep this in
   mind when designing your connectivity and heed the consequences of doing I/O
   in these callbacks. This is doubly important when using models with dynamic
   connectivity where the temptation to store all connections is even larger and
   each call to ``update`` will re-evaluate the corresponding callbacks.

.. _modelconnections:

.. glossary::

   connection
      Connections implement chemical synapses between **source** and **target** cells and are characterized
      by having a transmission delay.

      Connections in Arbor are defined in two steps:

      1. Create labeled **source** and **target** on two separate cells as part of their
         :ref:`cell descriptions <modelcelldesc>` in the :ref:`recipe <modelrecipe>`. Sources typically
         generate spikes. Targets are typically synapses with associated biophysical model descriptions.
         Each labeled group of sources or targets may contain multiple items on possibly multiple locations
         on the cell.
      2. Declare the connection in the recipe *on the target cell*:  from a source identified using
         a :gen:`global_label`; a target identified using a :gen:`local_label` (:gen:`gid` of target is
         the argument of the recipe method); a connection delay and a connection weight.

         .. code-block:: python

             def connections_on(self, gid):
                 if gid + 1 < self.num_cells():
                     return [arbor.connection((gid + 1, "spike-source"), "synapse", weight, delay)]
                 else:
                     return []

   spike
   action potential
      Spikes travel over :term:`connections <connection>`. In a synapse, they generate an event.

   threshold detector
      :ref:`Placed <cablecell-place>` on a cell. Possible source of a connection.
      Detects crossing of a fixed threshold and generates corresponding events.
      Also used to record spikes for analysis. See :ref:`here
      <cablecell-threshold-detectors>` for more information.

   spike source cell
      Artificial cell to generate spikes on a given schedule, see :ref:`spike cell <spikecell>`.

   recording
      By default, spikes are used for communication, but not stored for analysis,
      however, :ref:`simulation <modelsimulation>` objects can be instructed to record spikes.

   event
      In a synapse :term:`spikes <spike>` generate events, which constitute stimulation of the synapse
      mechanism and the transmission of a signal. A synapse may receive events directly from an
      :term:`event generator`.

   event generator
      Externally stimulate a synapse. Events can be delivered on a schedule.
      See :py:class:`arbor.event_generator` for details.

.. _modelgapjunctions:

.. glossary::

   gap junction connection
      Gap junctions represent electrical synapses where transmission between cells is bidirectional and direct.
      They are modelled as a conductance between two **gap junction sites** on two cells.

      Similarly to `Connections`, Gap Junctions in Arbor are defined in two steps:

      1. Create labeled **gap junction sites** on two separate cells as part of
         their :ref:`cell descriptions <modelcelldesc>` in the :ref:`recipe
         <modelrecipe>`. Each labeled group of gap junctions may contain multiple
         items on possibly multiple locations on the cell.
      2. Declare the Gap Junction connections in the recipe *on the local cell*:
         from a peer **gap junction site** identified using a
         :gen:`global_label`; to a local **gap junction site** identified using
         a :gen:`local_label` (:gen:`gid` of the site is implicitly known); and
         a unit-less connection weight. Two of these connections are needed, on
         each of the peer and local cells. The callback `gap_junctions_on`
         returns a list of these items, eg

         .. code-block:: python

             def gap_junctions_on(self, gid):
                 n = self.num_cells
                 if gid + 1 < n and gid > 0:
                     return [arbor.gap_junction_connection((gid + 1, "gj"), "gj", weight),
                             arbor.gap_junction_connection((gid - 1, "gj"), "gj", weight),]
                 elif gid + 1 < n:
                     return [arbor.gap_junction_connection((gid + 1, "gj"), "gj", weight),]
                 if gid > 0:
                     return [arbor.gap_junction_connection((gid - 1, "gj"), "gj", weight),]
                 else:
                     return []

         Note that gap junction connections are symmetrical and thus the above
         example generates two connections, one incoming and one outgoing.

   .. Note::
      Only cable cells support gap junctions as of now.


High-Level Network Specification
--------------------------------
A high-level network specification API is available to generate cell connections and gap junctions.
The API is to some extent based on the `connection-set algebra <https://pubmed.ncbi.nlm.nih.gov/22437992/>`_ proposed by Mikael Djurfeldt.
After specifying cell populations as source and destination, a selection criterion is independently applied to each possible connection identified by a global source and destination label.
This allows for lazy network generation, by only evaluating the selection criterion for connections relevant for a requested cell. Network selections can be combined by the logical operations `and`, `or` and `xor`.
For each generated connection, a network value provides a parameter, such as a weight or delay, which is similarly independently generated.
An example for a network selection is `inter_cell` for selecting connections between different cells and `bernoulli_random` for a random selection.
For network values, uniform and normal random distributions are possible examples. In more complex cases, a custom selection or value function can also be provided.
A network is intended to be created once and called to generate from within a cell recipe when required.

An example for combination of a custom and random network selection in Python:

.. code-block:: python

   # Custom network selection
   class ring_selection:
       def __init__(self, ncells):
           self.ncells = ncells

       # Select based on the cell_global_label src and dest
       def __call__(self, src, dest):
           if src.gid + 1 == dest.gid or (dest.gid == 0 and src.gid == ncells - 1):
               return True

           return False


   # (5) Create a recipe that generates a network of connected cells.
   class random_ring_recipe(arbor.recipe):
       def __init__(self, ncells):
           # The base C++ class constructor must be called first, to ensure that
           # all memory in the C++ class is initialized correctly.
           arbor.recipe.__init__(self)
           self.ncells = ncells
           self.props = arbor.neuron_cable_properties()

           w = arbor.network_value.uniform_distribution(
               84, [0.01, 0.02]
           )  # Random between 0.01 and 0.02 μS on expsyn
           d = 5  # ms delay

           # Select all cells by specifying the range [0, ncells).
           # Use different local labels for source and destination.
           src_pop = [arbor.cell_global_range_label(0, ncells, "detector")]
           dest_pop = [arbor.cell_global_range_label(0, ncells, "syn")]

           # Select only inter-cell connections, which are either randomly selected with a 8%
           # probability or part of a ring.
           selection = arbor.network_selection.inter_cell() & (
               arbor.network_selection.bernoulli_random(42, 0.08)
               | arbor.network_selection.custom(ring_selection(ncells))
           )

           self.connection_network = arbor.cell_connection_network(
               w, d, selection, src_pop, dest_pop
           )

       # (8) Generate network
       def connections_on(self, gid):
           return self.connection_network.generate(gid)



API
---

* :ref:`Python <pyinterconnectivity>`
* :ref:`C++ <cppinterconnectivity>`
