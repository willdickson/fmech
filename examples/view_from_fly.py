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

Author: William Dickson 
------------------------------------------------------------------------
"""
from __future__ import print_function
import os
import sys
import vtk 
import time
import string 
import cPickle 
from fmech.sim_util import *

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


# ----------------------------------------------------------------------
if __name__ == '__main__':

    # Rendering step
    step = 1
    fullscreen = False

    # Read log file
    log_file = sys.argv[1] 
    log_fid = open(log_file, 'r')
    log = cPickle.load(log_fid)
    log_fid.close()
    n = len(log.body_p)
    n_list = range(0,n,step)
    
    # Create the Renderer, RenderWindow, and RenderWindowInteractor
    ren = vtk.vtkRenderer()
    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren.SetBackground(0.0, 0.0, 0.0)
    if fullscreen:
        ren_win.FullScreenOn()
    else:
        ren_win.SetSize(1200,800)
    ren_win.SwapBuffersOn()

    # Setup light sources
    light = vtk.vtkLight()
    light.SetIntensity(1.5)
    light.SetPosition(50, 0, 50)
    light.SetDiffuseColor(1, 1, 1)
    ren.AddLight(light)

    light = vtk.vtkLight()
    light.SetIntensity(0.4)
    light.SetPosition(-50, 0, 50)
    light.SetDiffuseColor(1, 1, 1)
    ren.AddLight(light)

    # Create checkerboard floor
    floor_size = 600.0
    floor_zpos = -20.0
    floor_res = 40
    floor_src = vtk.vtkPlaneSource()
    floor_src.SetCenter(0.0, 0.0, floor_zpos)
    floor_src.SetOrigin(-floor_size, -floor_size, floor_zpos)
    floor_src.SetPoint1( floor_size, -floor_size, floor_zpos)
    floor_src.SetPoint2(-floor_size,  floor_size, floor_zpos)
    floor_src.SetResolution(floor_res, floor_res)
    floor_src.Update()
    grid = floor_src.GetOutput()
    floor_colors = vtk.vtkUnsignedCharArray()
    floor_colors.SetNumberOfComponents(3)
    floor_colors.SetNumberOfTuples(grid.GetNumberOfCells())
    floor_color_0 = (200,200,200)
    floor_color_1 = (10,100,10)
    for i in range(grid.GetNumberOfCells()):
        if (i//floor_res)%2 == 0:
            if i%floor_res%2 == 0:
                floor_colors.InsertTuple(i,floor_color_0)
            else:
                floor_colors.InsertTuple(i,floor_color_1)
        else:
            if i%floor_res%2 == 0:
                floor_colors.InsertTuple(i,floor_color_1)
            else:
                floor_colors.InsertTuple(i,floor_color_0)
    grid.GetCellData().SetScalars(floor_colors)
    floor_mapper = vtk.vtkPolyDataMapper()
    floor_mapper.SetInputConnection(floor_src.GetOutputPort())
    floor_actor = vtk.vtkActor()
    floor_actor.SetMapper(floor_mapper)
    ren.AddActor(floor_actor)

    # Looming ball
    ball = vtk.vtkSphereSource()
    ball.SetRadius(10.0)
    ball.SetThetaResolution(25)
    ball.SetPhiResolution(25)
    ball_map = vtk.vtkPolyDataMapper()
    ball_map.SetInputConnection( ball.GetOutputPort())
    ball_actor = vtk.vtkActor()
    ball_actor.GetProperty().SetDiffuseColor((1,0,0))
    ball_actor.GetProperty().SetSpecular(.3)
    ball_actor.GetProperty().SetSpecularPower(30)
    ball_actor.SetMapper(ball_map)
    ball_start = 3*floor_size, 0.0, 0.0
    ball_final = log.body_p[int(0.9*(n-1))]
    ball_actor.SetPosition(ball_start)
    ren.AddActor(ball_actor)        
    ball_vec = vec_sub(ball_final, ball_start)

    # Setup camera
    cam = ren.GetActiveCamera()
    cam.SetViewUp((0,0,1))
    cam.SetViewAngle(90.0)

    ## Add text actor
    text = vtk.vtkTextActor()
    text.GetTextProperty().SetFontSize(16)
    text.SetInput('')
    text.SetDisplayPosition(5,5)
    text.GetProperty().SetColor(0.8,0.8,0.8)
    ren.AddActor2D(text)
    
    # Loop through log and render
    for i in n_list:

        # Let transeint die out before starting visualization.
        if log.t[i] < 300.0:
            continue 

        # Set camera
        p = log.body_p[i]
        q = log.body_q[i]
        cam.SetPosition((p[0], p[1], p[2] ))

        # Fly's view
        focal_point = vec_add(qrotate_vec((10,0,0),q),p)
        cam.SetFocalPoint(focal_point)

        # Overhead view
        #cam.SetPosition((p[0]-500, p[1] - 500, p[2]+100 ))
        #cam.SetFocalPoint(p)

        ball_pos = vec_add(ball_start, vec_scalar_mul(float(i)/int(0.9*(n-1)),ball_vec))
        ball_actor.SetPosition(ball_pos)
            
        # Add time text
        text.SetInput('t %4.1f ms'%(log.t[i],))
        
        # Reset clipping range and render
        ren.ResetCameraClippingRange() 
        ren_win.Render()

        if False and i==0:
            raw_input('press enter to start')

        #time.sleep(0.01)
