"""
-----------------------------------------------------------------------
fmech
Copyright (C) William Dickson, 2008.
  
wbd@caltech.edu
www.willdickson.com

Released under the LGPL Licence, Version 3

This file is part of fmech.

fmech is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
fmech is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with fmech.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------   
mitrich.pyx

Purpose: pyrex wrapper for Brian Mirtich's mass properties code. 

Currently this is a complete hack job. I plan on cleaning this up
later.

Author: William Dickson 
------------------------------------------------------------------------
"""
import math, stl_tools

cdef extern from "volInt.h":

    # Data strutures
    ctypedef enum:
        MAX_VERTS 
        MAX_FACES 
        MAX_POLYGON_SZ
        
    ctypedef struct FACE:
        int numVerts
        float norm[3]
        float w
        int verts[MAX_POLYGON_SZ]
        void *poly
    
    ctypedef struct POLYHEDRON:
        int numVerts
        int numFaces
        float verts[MAX_VERTS][3]
        FACE faces[MAX_FACES]

    # Global variables
    float **raw_verts
    int **raw_faces

    float T0   
    float T1[3]
    float T2[3]
    float TP[3]

    # Volume integral Functions 
    cdef void readPolyhedron(POLYHEDRON *p, int numVerts,int numFaces)
    cdef void compVolumeIntegrals(POLYHEDRON *p)

    # Utility functions - bounds checking 
    cdef int get_max_verts()
    cdef int get_max_faces()
    cdef int get_max_polygon_sz()


# Useful constants
cdef int X
cdef int Y
cdef int Z
X = 0
Y = 1
Z = 2

def get_mass_props(filename, density):
    """
    Computes the mass properties of the polyhedral surface in the
    given stl file assuming a uniform density. Note, units must be
    consistent.

    Inputs:

       filename = name of the stl file to read
       density = density of the body

    Outputs: (mass, cm, it)

       mass = mass of the body
       cm =  center of mass (x,y,z)
       it = inertia tensor of the body (I11,I22,I33,I12,I13,I23)
    
    """

    cdef POLYHEDRON p
    cdef int numVerts
    cdef int numFaces

    # Read file and get facet list
    facet_list = stl_tools.read_stl(filename)    
    vertex_dict = stl_tools.get_vertex_dict(facet_list)

    # Check number of vertices and number of face with max allowed
    if len(facet_list) > get_max_faces():
        raise RuntimeError, '# faces (%d) exceeds max # allowed (%d)'%(len(facet_list), get_max_faces())
    if len(vertex_dict) > get_max_verts():
        raise RuntimeError, '# vertices (%d) exceeds max # allowed (%d)'%(len(vertex_dict), get_max_verts())

    # Get list of raw vertices
    vertex_list = vertex_dict.items()
    vertex_list.sort(stl_tools.vertex_item_cmp)
    numVerts = len(vertex_list)
    for i in range(0,numVerts):
        for j in range(0,3):
            raw_verts[i][j] = vertex_list[i][0][j]
            
    # Get array of raw faces order ccw w.r.t. outward normal
    numFaces = len(facet_list)
    for i in range(0,numFaces):
        facet = facet_list[i]
        facet.verts2CCW()
        if len(facet.vertices) > get_max_polygon_sz():
            n1 = len(facet.vertices)
            n2 = get_max_polygon_sz()
            raise RuntimeError, 'face %d, # of vertices(%d) exceeds max # allowed (%d) '%(i,n1,n2) 
        raw_faces[i][0] = len(facet.vertices)
        for j in range(0,len(facet.vertices)):
            vertex_ind = vertex_dict[facet.vertices[j]]
            raw_faces[i][j+1] = vertex_ind

    # Read the polyhedron
    readPolyhedron(&p, numVerts, numFaces)
    
    # Compute the volume intergals
    compVolumeIntegrals(&p)

    #Compute mass and center of mass    
    mass = T0*density
    cm = T1[X]/T0, T1[Y]/T0, T1[Z]/T0

    # Compute inertia tensor
    it = (
        density * (T2[Y] + T2[Z]),
        density * (T2[Z] + T2[X]),
        density * (T2[X] + T2[Y]),
        -density * TP[X],
        -density * TP[Y],
        -density * TP[Z]
        )

    # Translate inertia tensor to center of mass
    it = (
        it[0] - mass*(cm[Y]*cm[Y] + cm[Z]*cm[Z]),
        it[1] - mass*(cm[Z]*cm[Z] + cm[X]*cm[X]),
        it[2] - mass*(cm[X]*cm[X] + cm[Y]*cm[Y]),
        it[3] + mass*cm[X]*cm[Y],
        it[4] + mass*cm[Y]*cm[Z],
        it[5] + mass*cm[Z]*cm[X]
        )
    
    return mass, cm, it




    

    
    

    
