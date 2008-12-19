

	/*******************************************************
        *                                                      *
	*  volInt.h                                            *
	*                                                      *
	*  This code computes volume integrals needed for      *
	*  determining mass properties of polyhedral bodies.   *
	*                                                      *
	*  For more information, see the accompanying README   *
	*  file, and the paper                                 *
	*                                                      *
	*  Brian Mirtich, "Fast and Accurate Computation of    *
	*  Polyhedral Mass Properties," journal of graphics    *
	*  tools, volume 1, number 1, 1996.                    *
	*                                                      *
	*  This source code is public domain, and may be used  *
	*  in any way, shape or form, free of charge.          *
	*                                                      *
	*  Copyright 1995 by Brian Mirtich                     *
	*                                                      *
	*  mirtich@cs.berkeley.edu                             *
	*  http://www.cs.berkeley.edu/~mirtich                 *
        *                                                      *
	*******************************************************/

/*
	Revision history

	26 Jan 1996	Program creation.

	 3 Aug 1996	Corrected bug arising when polyhedron density
			is not 1.0.  Changes confined to function main().
			Thanks to Zoran Popovic for catching this one.

	27 May 1997     Corrected sign error in translation of inertia
	                product terms to center of mass frame.  Changes 
			confined to function main().  Thanks to 
			Chris Hecker.

        27 Dec 2004     (William Dickson) Modified program into a python 
                        extension
*/


#include <stdio.h>
#include <string.h>
#include <math.h>

/*
   ============================================================================
   constants
   ============================================================================
*/

#define MAX_VERTS 30000   /* maximum number of polyhedral vertices */
#define MAX_FACES 30000   /* maximum number of polyhedral faces */
#define MAX_POLYGON_SZ 20 /* maximum number of verts per polygonal face */

#define X 0
#define Y 1
#define Z 2

/*
   ============================================================================
   macros
   ============================================================================
*/

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

/*
   ============================================================================
   data structures
   ============================================================================
*/

typedef struct {
  int numVerts;
  double norm[3];
  double w;
  int verts[MAX_POLYGON_SZ];
  struct polyhedron *poly;
} FACE;

typedef struct polyhedron {
  int numVerts, numFaces;
  double verts[MAX_VERTS][3];
  FACE faces[MAX_FACES];
} POLYHEDRON;



/*
   ============================================================================
   globals
   ============================================================================
*/


float raw_verts[MAX_VERTS][3];
int raw_faces[MAX_FACES][MAX_POLYGON_SZ+1];

int A;   /* alpha */
int B;   /* beta */
int C;   /* gamma */

/* projection integrals */
double P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

/* face integrals */
double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;

/* volume integrals */
double T0, T1[3], T2[3], TP[3];

/* Volume integral functions */
void readPolyhedron(POLYHEDRON *p,int numVerts,int numFaces);
void compProjectionIntegrals(FACE *f);
void compFaceIntegrals(FACE *f);
void compVolumeIntegrals(POLYHEDRON *p);

/* Utility functions - to provde a it of error checking in the pyrex code */
int get_max_verts(void);
int get_max_faces(void);
int get_max_polygon_sz(void);


