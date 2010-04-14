#!/usr/bin/env python
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
view_sim.py

Purpose: a simple viewer for the insect flight simulations.

Author: William Dickson 
------------------------------------------------------------------------
"""
import sys
import vtk, string, cPickle, os
from fmech.sim_util import quat2mat
from fmech.sim_util import read_config
from time import sleep

# Rendering step
step = 1

# Get name of log file
log_file = sys.argv[1] 

# Camera position
cam_pos = (20,0,0)
cam_fpt = (0,0,0)

# Polyhedra data files
config = read_config()
body_stl_file = config['body_stl_file']
wing_stl_file = config['wing_stl_file']

def get_transform(p,q):
    """
    Gets the vtk transformation for the given position vector
    and orientation quaternion.
    """
    mat4x4 = vtk.vtkMatrix4x4()
    mat4x4.Identity()
    
    # Set rotation from quaternion
    r = quat2mat(q)
    cnt = 0
    for i in range(0,3):
        for j in range(0,3):
            mat4x4.SetElement(i,j, r[cnt])
            cnt += 1
    # Set translation
    for i in range(0,3):
        mat4x4.SetElement(i,3, p[i])

    # Set transform
    trans = vtk.vtkTransform()
    trans.SetMatrix(mat4x4)
    return trans

# Setup data sources
body_src = vtk.vtkSTLReader()
body_src.SetFileName(body_stl_file)
wing_src = vtk.vtkSTLReader()
wing_src.SetFileName(wing_stl_file)

# Try a normal filter on the body to get shading to work
body_filt = vtk.vtkPolyDataNormals()
body_filt.ConsistencyOn()
body_filt.SplittingOn()
body_filt.ComputePointNormalsOn()
body_filt.SetInput(body_src.GetOutput())

wing_filt = vtk.vtkPolyDataNormals()
wing_filt.ConsistencyOn()
wing_filt.SplittingOn()
wing_filt.ComputePointNormalsOn()
wing_filt.SetInput(wing_src.GetOutput())

# Set up mappers
body_mapper = vtk.vtkPolyDataMapper()
body_mapper.SetInput(body_filt.GetOutput())
wing_mapper = vtk.vtkPolyDataMapper()
wing_mapper.SetInput(wing_filt.GetOutput())

# Setup body actors
body_actor = vtk.vtkActor()
body_actor.GetProperty().SetColor((0.6, 0.6, 0.6))
body_actor.GetProperty().SetInterpolationToGouraud()
body_actor.GetProperty().SetSpecular(0.2)
#body_actor.GetProperty().SetAmbient(0.2)
body_actor.SetMapper(body_mapper)
# Setup left wing actor
l_wing_actor = vtk.vtkActor()
l_wing_actor.GetProperty().SetColor((0.6, 0.6, 0.6))
l_wing_actor.GetProperty().SetOpacity(0.7)
l_wing_actor.GetProperty().SetInterpolationToGouraud()
#l_wing_actor.GetProperty().SetAmbient(0.3)
l_wing_actor.SetMapper(wing_mapper)
# Setup right wing actor
r_wing_actor = vtk.vtkActor()
r_wing_actor.GetProperty().SetColor((0.6, 0.6, 0.6))
r_wing_actor.GetProperty().SetOpacity(0.7)
r_wing_actor.GetProperty().SetInterpolationToGouraud()
#r_wing_actor.GetProperty().SetAmbient(0.3)
r_wing_actor.SetMapper(wing_mapper)

# Read log file
log_fid = open(log_file, 'r')
log = cPickle.load(log_fid)
log_fid.close()

# Create the Renderer, RenderWindow, and RenderWindowInteractor
ren = vtk.vtkRenderer()
ren_win = vtk.vtkRenderWindow()
ren_win.AddRenderer(ren)
ren.SetBackground(0.0, 0.0, 0.0)
ren_win.FullScreenOn()
ren_win.SwapBuffersOn()

# Create balls
if 1:
    ball_grid = (4,3)
    ball_spacing = 100
    ball_z = -20.0
    ball = vtk.vtkSphereSource()
    ball.SetRadius(10.0)
    ball.SetThetaResolution(25)
    ball.SetPhiResolution(25)
    mapBall = vtk.vtkPolyDataMapper()
    mapBall.SetInput( ball.GetOutput())
    for i in range(0,ball_grid[0]):
        for j in range(0,ball_grid[1]):
            ball_x = i*ball_spacing - 50
            ball_y = -0.5*ball_spacing*ball_grid[1] +  j*ball_spacing
            ballActor = vtk.vtkActor()
            ballActor.GetProperty().SetDiffuseColor((1,0,0))
            ballActor.GetProperty().SetSpecular(.3)
            ballActor.GetProperty().SetSpecularPower(30)
            ballActor.SetMapper(mapBall)
            ballActor.SetPosition(ball_x, ball_y, ball_z)
            ren.AddActor(ballActor)        

# Add the actors to the render; set the background and size
ren.AddActor(body_actor)
ren.AddActor(l_wing_actor)
ren.AddActor(r_wing_actor)
cam = ren.GetActiveCamera()
cam.SetPosition(cam_pos)
cam.SetFocalPoint(cam_fpt)
cam.SetViewUp((0,0,1))

# Add text actor
text_mapper = vtk.vtkTextMapper()
text_mapper.GetTextProperty().SetFontSize(16)
text_mapper.SetInput('')
text = vtk.vtkTextActor()
text.SetMapper(text_mapper)
text.SetDisplayPosition(5,5)
text.GetProperty().SetColor(0.8,0.8,0.8)
ren.AddActor(text)
       
# Loop through log and render
n = len(log.body_p)
n_list = range(0,n,step)

for i in n_list:
    # Set body position and orientation
    p = log.body_p[i]
    q = log.body_q[i]
    trans = get_transform(p,q) 
    body_actor.SetUserTransform(trans)

    cam.SetFocalPoint(p)
    cam.SetPosition((p[0]+20, p[1] + 20, p[2] ))

    # Set left wing position and orientation
    p = log.l_wing_p[i]
    q = log.l_wing_q[i]
    trans = get_transform(p,q) 
    l_wing_actor.SetUserTransform(trans)
 
    # Set right wing position and orientation
    p = log.r_wing_p[i]
    q = log.r_wing_q[i]
    trans = get_transform(p,q) 
    r_wing_actor.SetUserTransform(trans)

    # Add time text
    text_mapper.SetInput('t %4.1f ms'%(log.t[i],))
    
    # Reset clipping range and render
    ren.ResetCameraClippingRange() 
    ren_win.Render()
