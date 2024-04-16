"""
Unit and regression test for the montecarlo_file module.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import montecarlo_jg as montecarlo

def get_IsingHamiltonian(G, mus=None):
    if mus == None:
        mus = np.zeros(len(G.nodes()))

    if len(G.nodes()) != len(mus):
        error("DimensionMismatch")

    if len(G.nodes()) != len(mus):
        error(" Dimension Mismatch")
    J = [[] for i in G.nodes()]
    for e in G.edges:
        J[e[0]].append((e[1], G.edges[e]['weight']))
        J[e[1]].append((e[0], G.edges[e]['weight']))
    return montecarlo.IsingHamiltonian(J,mus)

def test_BitString():
    my_bs = montecarlo.BitString(8)
    my_bs.flip_site(2)
    my_bs.flip_site(2)
    my_bs.flip_site(2)
    my_bs.flip_site(7)
    my_bs.flip_site(0)
    assert(len(my_bs) == 8)

    my_bs = montecarlo.BitString(13)
    my_bs.set_config([0,1,1,0,0,1,0,0,1,0,1,0,0])

    assert(my_bs.on() == 5)
    assert(my_bs.off() == 8)

    assert(my_bs.int() == 3220)

    my_bs = montecarlo.BitString(20)
    my_bs.set_int_config(3221)

    # Let's make sure this worked:
    tmp = np.array([0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1])
    assert((my_bs.config == tmp).all())

    # We can provide an even stronger test here:
    for i in range(1000):
        my_bs.set_int_config(i) # Converts from integer to binary
        assert(my_bs.int() == i) # Converts back from binary to integer and tests

    my_bs1 = montecarlo.BitString(13)
    my_bs1.set_config([0,1,1,0,0,1,0,1,1,0,1,0,0])

    my_bs2 = montecarlo.BitString(13)
    my_bs2.set_int_config(3252)

    assert(my_bs1 == my_bs2)

    my_bs2.flip_site(5)
    assert(my_bs1 != my_bs2)



def test_IsingHamiltonian():
    """Test that the BitString class does what we expect."""
    
    N = 6
    Jval = 2.0
    G = nx.Graph()
    G.add_nodes_from([i for i in range(N)])
    G.add_edges_from([(i,(i+1)% G.number_of_nodes() ) for i in range(N)])
    for e in G.edges:
        G.edges[e]['weight'] = Jval

    # Define a new configuration instance for a 6-site lattice
    conf = montecarlo.BitString(N)
    ham = get_IsingHamiltonian(G)

    # Compute the average values for Temperature = 1
    E, M, HC, MS = ham.compute_average_values(1)

    assert(np.isclose(E,  -11.95991923))
    assert(np.isclose(M,   -0.00000000))
    assert(np.isclose(HC,   0.31925472))
    assert(np.isclose(MS,   0.01202961))