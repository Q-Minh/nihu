General programming techniques {#tech_other}
==============================

\page tech_other

[TOC]

Introduction {#tech_other_intro}
============


The store pattern {#tech_other_store}
=================

The extremely useful store pattern is simply implemented as
~~~~~~~~~~~~~~~~
template <class C>
struct store
{
	static const C m_data;
};

template <class C>
const C store<C>::m_data;
~~~~~~~~~~~~~~~~


The So-called DRTP-pattern {#tech_other_drtp}
==========================

