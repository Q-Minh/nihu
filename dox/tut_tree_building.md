Tree building {#tut_tree_building}
=============

\page tut_tree_building How to build and work with cluster trees

[TOC]

Introduction {#tut_tree_intro}
============

Cluster trees are hierarchical oc-trees built up based on a sequence of nodal locations.
Cluster trees provide the interface for various tree traversal methods, needed to implement Fast Multipole Methods.


Building a Cluster Tree {#tut_tree_build}
========================

The Tree Type {#tut_tree_cluster_type}
--------------------------------------

Cluster trees are a hierarchical collections of clusters, therefore, the cluster 
type needs to be defined first.
In our simplest example, we work with a cluster tree containing 3-dimensional 
NiHu::fmm::empty_cluster clusters.

\snippet tree_building.cpp cluster type

The nodal locations {#tut_tree_nodal_locations}
-----------------------------------------------

We create an array of random nodal locations.

\snippet tree_building.cpp nodal locations


Building the tree {#tut_tree_building}
--------------------------------------

The cluster tree is built with its constructor.
In the simplest case the constructor needs three parameters: 
	- the two iterators defining the range of nodal locations, 
	- and a cluster division functor. The functor receives a cluster as argument,
	and returns if the cluster needs to be divided into children or not. 
	As an example, we apply the built_in functor NiHu::fmm::divide_num_nodes 
	that divides a cluster if the number of nodes in a cluster is larger than its parameter.

\snippet tree_building.cpp construct tree

Alternatively, NiHu contains built in dividing functors for  building balanced trees
of specified depth or prescribed leaf level cluster diameter:

\snippet tree_building.cpp construct tree alternatives


Source and Receiver nodes {#tut_src_rec_nodes}
----------------------------------------------

In the above trees, each node is both source and receiver.
In order to define source and receiver nodes separately, the alternative constructor
can be used that takes two ranges as input:

\snippet tree_building.cpp construct tree source receiver


Building the tree from NiHu meshes
----------------------------------

In most boundary element applications, the cluster tree is not built based on a 
sequence of nodal locations, rather based on elements of a boundary element mesh.
NiHu provides built-in iterators NiHu::fmm::elem_center_iterator that traverse 
mesh elements but return the element center when dereferenced:

\snippet tree_building.cpp construct tree mesh


Traversing Cluster Trees {#tut_tree_traverse}
=============================================

Simple traversing
-----------------

The NiHu::fmm::cluster_tree is a vector container of clusters, and clusters can be accessed by indexing.
The following snippet, for example, traverses each cluster of a tree:

\snippet tree_building.cpp vector traversing

Level by level traversing
-------------------------

The clusters are stored in a level-contiguous way that makes level-by-level traversing possible
by getting the begin and end index of clusters within a specified level:

\snippet tree_building.cpp level traversing

Leaf Traversing
---------------

In an unbalanced tree, the leaf clusters are not stored continuously. Therefore,
NiHu cluster trees provide an index vector to traverse leaf clusters:

\snippet tree_building.cpp leaf traversing

Similarly, source and receiver leaf clusters can be traveresed by getting their
index vectors:

\snippet tree_building.cpp source receiver traversing

Depth First Search traversing
-----------------------------

Depth first search traverse of the tree is provided by the cluster class that contains
indices of its children in the tree.

\snippet tree_building.cpp dfs traverse
