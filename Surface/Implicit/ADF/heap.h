#ifndef SLIMHEAP_INCLUDED // -*- C++ -*-
#define SLIMHEAP_INCLUDED

#ifdef _MSC_VER
#pragma once
#pragma warning(disable : 4786)
#endif


/************************************************************************

  Heap data structure for greedy algorithms with local update.

  This heap implementation uses the standard approach of implementing
  a priority queue using a heap-ordered contiguous array.  However, it
  differs from the typical implemtnation in that each item in the heap
  keeps track of its location in the heap.  This allows the client
  algorithm to selectively update the keys of heaped items and
  efficiently restore the heap-ordering property.

  Copyright (C) 1998,2001 Michael Garland.
  
  $Id: heap.h,v 1.1.1.1 2004/02/04 03:19:51 wensu Exp $

 ************************************************************************/

// #include <gfx/array.h>
#include <vector>

class MxHeapable
{
private:
    double import;
    int token;

public:
    MxHeapable() { not_in_heap(); heap_key(0.0f); }

    inline bool is_in_heap() { return token != -47; }
    inline void not_in_heap() { token = -47; }
    inline int get_heap_pos() { return token; }
    inline void set_heap_pos(int t) { token=t; }

    inline void  heap_key(double k) { import=k; }
    inline double heap_key() const  { return import; }
};


class MxHeap : public std::vector<MxHeapable *>
{
private:
    void place(MxHeapable *x, unsigned int i);
    void swap(unsigned int i, unsigned int j);

    unsigned int parent(unsigned int i) { return (i-1)/2; }
    unsigned int left(unsigned int i) { return 2*i+1; }
    unsigned int right(unsigned int i) { return 2*i+2; }

    void upheap(unsigned int i);
    void downheap(unsigned int i);

    

public:
    MxHeap() : std::vector<MxHeapable *>(8) { }
    MxHeap(unsigned int n) : std::vector<MxHeapable *>(n) { }

    void insert(MxHeapable *t) { insert(t, t->heap_key()); }
    void insert(MxHeapable *, double);
    void update(MxHeapable *t) { update(t, t->heap_key()); }
    void update(MxHeapable *, double);
	// varray reset was just setting the size to 0
    void reset() { clear(); }

    unsigned int size() const { return (begin() == end() ? 0 : end() - begin()); }
    MxHeapable       *item(int i)       { return (*this)[i]; }
    const MxHeapable *item(int i) const { return (*this)[i]; }
    MxHeapable *extract();
    MxHeapable *top() { return (size()<1 ? (MxHeapable *)NULL : item(0)); }
    MxHeapable *remove(MxHeapable *);
};

// SLIMHEAP_INCLUDED
#endif
