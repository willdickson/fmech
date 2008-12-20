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
stl_tools.py 

Purpose: a set of tools for manipulating stl files.

Author: William Dickson 
------------------------------------------------------------------------
"""
import sys, string, math, copy
try:
    import cgkit.cgtypes as cgtypes # cgkit 2
except ImportError, err:
    import cgtypes # cgkit 1

def cross_prod(a,b):
    """
    Compute the cross product of two 3 vectors

    """
    c0 = a[1]*b[2] - b[1]*a[2]
    c1 = b[0]*a[2] - a[0]*b[2]
    c2 = a[0]*b[1] - b[0]*a[1]
    return c0, c1, c2

def dot_prod(a,b):
    """
    Compute the dot product of two 3 vectors
    """
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vect_sub(a,b):
    """
    Subtract vector b from vector a
    """
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]

def vect_len(a):
    """
    Get the length of a vector
    """
    return math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)

def vect2unit(v):
    """
    Computes the unit vector in v direction
    """
    mag = vect_len(v)
    if mag != 0:
        return tuple([x/mag for x in v])
    else:
        return tuple([0 for x in v])

class stl_facet:
    """
    A really simple class for representing facets of a triangulated
    surface.

    self.vertices = list of facet vertices
    self.ow_normal = outward facing normal vertor
    """
    def __init__(self, vertices, ow_normal):
        self.vertices = vertices
        self.ow_normal = ow_normal

    def verts2CCW(self):
        """
        Reorder vertices in counter-clockwise order (when looking
        at the face from outside the polyhedron)
        """
        v0 = vect_sub(self.vertices[0], self.vertices[1])
        v1 = vect_sub(self.vertices[0], self.vertices[-1])
        normal = cross_prod(v0,v1)
        test = dot_prod(normal, self.ow_normal)
        if test < 0:
            self.vertices.reverse()

def read_facet(infile):
    """
    Reads a single facet from the stl file
    """
    line = infile.readline()
    line = string.split(line)
    if line[0] == 'endsolid':
        return  None
    else:
        assert line[0] == 'facet', 'No facet'
        ow_normal = float(line[2]), float(line[3]), float(line[4])
        line = infile.readline()
        line = string.split(line)
        assert line[0] == 'outer', 'No outer'
        vertices = []
        for i in range(0,3):
            line = infile.readline()
            line = string.split(line)
            assert line[0] == 'vertex', 'No vertex'
            vertex = float(line[1]), float(line[2]), float(line[3])
            vertices.append(vertex)
        line = infile.readline()
        line = string.split(line)
        assert line[0] == 'endloop', 'No endloop'
        line = infile.readline()
        line = string.split(line)
        assert line[0] == 'endfacet', 'No endfacet'
        return stl_facet(vertices, ow_normal)

def read_stl(file_name):
    """
    Read contents of the stl file
    """
    infile = open(file_name,'r')
    line = infile.readline()
    line = string.split(line)
    assert line[0] == 'solid', 'This does not apper to be an ascii stl file'
    facet_list = []
    while 1:
        out_val = read_facet(infile)
        if out_val == None:
            break
        facet_list.append(out_val)
    infile.close()
    return facet_list

def write_stl(filename, facet_list):
    """
    Write stl file of the given facet list
    """
    outfile = open(filename,'w')
    outfile.write('solid ascii\n')
    for facet in facet_list:
        outfile.write(' '*1+'facet normal %f %f %f\n'%facet.ow_normal)
        outfile.write( ' '*2+'outer loop\n')
        for i in range(0,3):
            outfile.write(' '*3+'vertex %f %f %f\n'%facet.vertices[i])
        outfile.write(' '*2+'endloop\n')
        outfile.write(' '*1+'endfacet\n')
    outfile.write('endsolid')
    outfile.close()

def get_vertex_dict(facet_list):
    """
    Get list of unique vertices from facet list. Assigns a unique index
    to each vertex based on thier order of occurance in the list of facets.
    """
    cnt = 0
    vertex_dict = {}
    for facet in facet_list:
        for vertex in facet.vertices:
            if not vertex_dict.has_key(vertex):
                vertex_dict[vertex] = cnt
                cnt+=1
    return vertex_dict

def get_edge_dict(facet_list, vertex_dict):
    """
    Get list of unique edges (in terms of the vertex indices. Assigns a
    unique index to each edge.
    """
    cnt = 0
    edge_dict = {}
    for facet in facet_list:
        for i in range(0,len(facet.vertices)):
            j = (i+1)%len(facet.vertices)
            ind_0 = vertex_dict[facet.vertices[i]]
            ind_1 = vertex_dict[facet.vertices[j]]
            edge_0 = ind_0, ind_1
            edge_1 = ind_1, ind_0
            if not ( edge_dict.has_key(edge_0) or  edge_dict.has_key(edge_1) ):
                edge_dict[edge_0] = cnt
                cnt+=1
    return edge_dict

def write_mirtich(file_name, facet_list):
    """
    Write polygonal surface file which is compatible with Brian Mirtich's
    VolInt mass properties program.
    """
    vertex_dict = get_vertex_dict(facet_list)
    # Write output file
    outfile  = open(file_name, 'w')
    # Write vertices
    outfile.write('%d\n\n'%(len(vertex_dict),))
    vertex_list = vertex_dict.items()
    vertex_list.sort(vertex_item_cmp)
    for vertex in vertex_list:
        for i in range(0,3):
            outfile.write('%1f '%(vertex[0][i],))
        outfile.write('\n')
    # Write facets
    outfile.write('\n%d\n\n'%(len(facet_list,)))
    for facet in facet_list:
        outfile.write('%d '%(len(facet.vertices),))
        facet.verts2CCW()
        for vertex in facet.vertices:
            outfile.write('%d '%(vertex_dict[vertex],))
        outfile.write('\n')
    outfile.close()

def vertex_item_cmp(x,y):
    """
    Comparison function for sorting vertex_dict.items() w.r.t. to the
    vertex indices (values).
    """
    return x[1]-y[1]

def get_euler_number(facet_list, vertex_dict, edge_dict):
    """
    Compute the euler number of the surface.
    """
    return len(vertex_dict) - len(edge_dict) +len(facet_list)


def scale_facet(facet, scale):
    """
    Scale facet vertices by the given scaling factor
    """
    vertices = []
    for vertex in facet.vertices:
        vertices.append((scale*vertex[0], scale*vertex[1], scale*vertex[2]))
    return stl_facet(vertices, copy.copy(facet.ow_normal))

def scale_facet_list(facet_list, scale):
    """
    Scale list of facets by the given scaling factor
    """
    new_facet_list = []
    for facet in facet_list:
        new_facet_list.append(scale_facet(facet, scale))
    return new_facet_list

def shift_facet(facet, p):
    """
    Shift facet by position vector p
    """
    vertices = []
    for vertex in facet.vertices:
        vertices.append((vertex[0]+p[0], vertex[1]+p[1], vertex[2]+p[2]))
    return stl_facet(vertices, copy.copy(facet.ow_normal))

def shift_facet_list(facet_list, p):
    """
    Shift facets in facet list by position vector p
    """
    new_facet_list = []
    for facet in facet_list:
        new_facet_list.append(shift_facet(facet,p))
    return new_facet_list

def rotate_facet(facet, ax, ang):
    """
    Rotate facet using axis and angle
    """
    ow_normal = rotate_vec(facet.ow_normal, ax, ang)
    vertices = []
    for vertex in facet.vertices:
        new_vertex = rotate_vec(vertex, ax, ang)
        vertices.append(new_vertex)
    return stl_facet(vertices, ow_normal)

def rotate_facet_list(facet_list, ax, ang):
    """
    Rotate all facets in list using axis and angle
    """
    new_facet_list = []
    for facet in facet_list:
        new_facet = rotate_facet(facet, ax, ang)
        new_facet_list.append(new_facet)
    return new_facet_list

def rotate_vec(v, ax, ang):
    """
    Rotate a vector by a given angle about a given axis.

    v = input vector
    ax = rotation axis
    ang = rotation angle (radians)

    """
    # Get rotation quaterion and its inverse
    rot_q = cgtypes.quat(ang, ax)
    rot_q_inv = rot_q.inverse()
    # Rotate vector
    vq = cgtypes.quat(0.0,v[0],v[1],v[2])
    vq_new = rot_q*vq*rot_q_inv
    return vq_new.x, vq_new.y, vq_new.z

def deg2rad(ang_deg):
    """
    Converts degrees to radians
    """
    return math.pi*ang_deg/180.0

def get_extent(facet_list, n):
    """
    Get (max-min) for the given dimension over all vertices

    """
    vertex_dict = get_vertex_dict(facet_list)
    if n==0:
        val_list = [x for x,y,z in vertex_dict]
    elif n==1:
        val_list = [y for x,y,z in vertex_dict]
    elif n==2:
        val_list = [z for x,y,z in vertex_dict]
    else:
        raise RuntimeError, 'dimension n must be 0,1 or 2'
    return max(val_list) - min(val_list)
# ---------------------------------------------------------------

if __name__=='__main__':


    filename = 'body_ascii.stl'
    facet_list = read_stl(filename)
    facet_list = rotate_facet_list(facet_list, (0,0,1), deg2rad(-90.0))
    facet_list = shift_facet_list(facet_list, (4,0,0))
    vertex_dict = get_vertex_dict(facet_list)
    x_list = [x for x,y,z in vertex_dict.keys()]
    y_list = [y for x,y,z in vertex_dict.keys()]
    z_list = [z for x,y,z in vertex_dict.keys()]
    print 'Before scaling' + '-'*40
    print 'max(x): ', max(x_list), ', min(x): ', min(x_list)
    print 'max(y): ', max(y_list), ', min(y): ', min(y_list)
    print 'max(z): ', max(z_list), ', min(z): ', min(z_list)
    model_len = max(y_list)-min(y_list)
    body_len = 2.5
    scale = body_len/float(model_len)
    facet_list = scale_facet_list(facet_list, scale)
    vertex_dict = get_vertex_dict(facet_list)
    x_list = [x for x,y,z in vertex_dict.keys()]
    y_list = [y for x,y,z in vertex_dict.keys()]
    z_list = [z for x,y,z in vertex_dict.keys()]
    print 'After scaling' + '-'*40
    print 'max(x): ', max(x_list), ', min(x): ', min(x_list)
    print 'max(y): ', max(y_list), ', min(y): ', min(y_list)
    print 'max(z): ', max(z_list), ', min(z): ', min(z_list)
    write_stl('body_scaled.stl', facet_list)

