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
from __future__ import print_function
import sys
import vtk, string, cPickle, os
from time import sleep
from fmech.sim_util import vec_mag 
from fmech.sim_util import vec2unit 
from fmech.sim_util import vec_cross
from fmech.sim_util import vec_angle
from fmech.sim_util import quat2mat
from fmech.sim_util import read_config
from fmech.sim_util import vec_scalar_mul 
from fmech.sim_util import axis_angle2quat
from write_image import write_image


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

    # Camera position
    cam_pos = (20,0,0)
    cam_fpt = (0,0,0)
    
    # Get name of log file
    log_file = sys.argv[1] 

    # Read log file
    log_fid = open(log_file, 'r')
    log = cPickle.load(log_fid)
    log_fid.close()
    
    # Polyhedra data files
    config = read_config()
    body_stl_file = config['body_stl_file']
    wing_stl_file = config['wing_stl_file']
    
    # Setup data sources for fly
    body_src = vtk.vtkSTLReader()
    body_src.SetFileName(body_stl_file)
    wing_src = vtk.vtkSTLReader()
    wing_src.SetFileName(wing_stl_file)
    
    # Try a normal filter on the body to get shading to work
    body_filt = vtk.vtkPolyDataNormals()
    body_filt.ConsistencyOn()
    body_filt.SplittingOn()
    body_filt.ComputePointNormalsOn()
    #body_filt.SetInput(body_src.GetOutput())
    body_filt.SetInputConnection(body_src.GetOutputPort())
    
    wing_filt = vtk.vtkPolyDataNormals()
    wing_filt.ConsistencyOn()
    wing_filt.SplittingOn()
    wing_filt.ComputePointNormalsOn()
    #wing_filt.SetInput(wing_src.GetOutput())
    wing_filt.SetInputConnection(wing_src.GetOutputPort())
    
    # Set up mappers
    body_mapper = vtk.vtkPolyDataMapper()
    #body_mapper.SetInput(body_filt.GetOutput())
    body_mapper.SetInputConnection(body_filt.GetOutputPort())
    wing_mapper = vtk.vtkPolyDataMapper()
    #wing_mapper.SetInput(wing_filt.GetOutput())
    wing_mapper.SetInputConnection(wing_filt.GetOutputPort())
    
    # Setup body actors
    body_actor = vtk.vtkActor()
    #body_actor.GetProperty().SetColor((0.6, 0.6, 0.6))
    body_actor.GetProperty().SetColor((0.4, 0.4, 0.4))
    body_actor.GetProperty().SetInterpolationToGouraud()
    body_actor.GetProperty().SetSpecular(0.2)
    #body_actor.GetProperty().SetAmbient(0.2)
    body_actor.SetMapper(body_mapper)
    
    # Setup left wing actor
    l_wing_actor = vtk.vtkActor()
    #l_wing_actor.GetProperty().SetColor((0.6, 0.6, 0.6))
    l_wing_actor.GetProperty().SetColor((0.4, 0.4, 0.4))
    l_wing_actor.GetProperty().SetOpacity(1.0)
    l_wing_actor.GetProperty().SetInterpolationToGouraud()
    #l_wing_actor.GetProperty().SetAmbient(0.3)
    l_wing_actor.SetMapper(wing_mapper)
    
    # Setup right wing actor
    r_wing_actor = vtk.vtkActor()
    r_wing_actor.GetProperty().SetColor((0.6, 0.6, 0.6))
    r_wing_actor.GetProperty().SetOpacity(1.0)
    r_wing_actor.GetProperty().SetInterpolationToGouraud()
    #r_wing_actor.GetProperty().SetAmbient(0.3)
    r_wing_actor.SetMapper(wing_mapper)
    
    # Set up force arrows 
    arrow_src = vtk.vtkArrowSource()
    arrow_src.SetTipResolution(8)
    arrow_src.SetTipLength(0.3)
    arrow_src.SetTipRadius(0.1)
    arrow_mapper = vtk.vtkPolyDataMapper()
    arrow_mapper.SetInputConnection(arrow_src.GetOutputPort())
    num_blade_elem = len(log.r_wing_forces[0])
    l_arrow_actor_list = []
    r_arrow_actor_list = []
    arrow_actor_dict = {'l': l_arrow_actor_list, 'r': r_arrow_actor_list}
    for key, arrow_actor_list in arrow_actor_dict.items():
        for i in range(num_blade_elem):
            arrow_actor = vtk.vtkActor()
            arrow_actor.GetProperty().SetColor((0.0, 0.0, 1.0))
            arrow_actor.GetProperty().SetOpacity(1.0)
            arrow_actor.GetProperty().SetInterpolationToGouraud()
            arrow_actor.SetMapper(arrow_mapper)
            arrow_actor_list.append(arrow_actor)
    #print(log.r_wing_forces[10][0])
    #assert 1==0
    
    # Create the Renderer, RenderWindow, and RenderWindowInteractor
    ren = vtk.vtkRenderer()
    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren.SetBackground(0.0, 0.0, 0.0)
    #ren_win.FullScreenOn()
    ren_win.SetSize(1200,800)
    ren_win.SwapBuffersOn()

    light = vtk.vtkLight()
    light.SetIntensity(1.5)
    light.SetPosition(20, 0, 20)
    light.SetDiffuseColor(1, 1, 1)
    ren.AddLight(light)

    # Create floor
    if 1:
        floor_size = 400.0
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

    
    # Create balls
    if 0:
        ball_grid = (4,3)
        ball_spacing = 100
        ball_z = -20.0
        ball = vtk.vtkSphereSource()
        ball.SetRadius(10.0)
        ball.SetThetaResolution(25)
        ball.SetPhiResolution(25)
        mapBall = vtk.vtkPolyDataMapper()
        mapBall.SetInputConnection( ball.GetOutputPort())
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

    # Add arrow actors
    for key, arrow_actor_list in arrow_actor_dict.items():
        for arrow_actor in arrow_actor_list:
            ren.AddActor(arrow_actor)

    # Setup camera
    cam = ren.GetActiveCamera()
    cam.SetPosition(cam_pos)
    cam.SetFocalPoint(cam_fpt)
    cam.SetViewUp((0,0,1))

    ## Add text actor
    text = vtk.vtkTextActor()
    text.GetTextProperty().SetFontSize(16)
    text.SetInput('')
    text.SetDisplayPosition(5,5)
    text.GetProperty().SetColor(0.8,0.8,0.8)
    ren.AddActor2D(text)
           
    # Loop through log and render
    n = len(log.body_p)
    n_list = range(0,n,step)

    
    for i in n_list:

        if log.t[i] < 300.0:
            continue 

        # Set body position and orientation
        p = log.body_p[i]
        q = log.body_q[i]
        trans = get_transform(p,q) 
        body_actor.SetUserTransform(trans)

        # Set camera
        cam.SetFocalPoint(p)
        cam.SetPosition((p[0]+20, p[1] + 20, p[2] ))
        cam.SetPosition((p[0]+20, p[1] + 20, p[2]+5 ))
        #cam.SetPosition((p[0]+40, p[1] + 0, p[2]+20 ))
    
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

        # Set wing forces
        force_disp_scale = 600.0
        wing_forces_dict = {'l': log.l_wing_forces[i], 'r': log.r_wing_forces[i]}
        for key, arrow_actor_list in arrow_actor_dict.items():
            for elem, arrow_actor in enumerate(arrow_actor_list):
                elem_pos, elem_force =  wing_forces_dict[key][elem]
                elem_force_mag = vec_mag(elem_force)
                elem_force_unit = vec2unit(elem_force)
                elem_rot_angle = vec_angle((1,0,0), elem_force_unit)
                elem_rot_axis = vec2unit(vec_cross((1,0,0), elem_force_unit))
                elem_rot_quat = axis_angle2quat(elem_rot_axis, elem_rot_angle)
                elem_transform = get_transform(elem_pos, elem_rot_quat)
                #arrow_actor.SetScale((force_disp_scale*elem_force_mag, 1.0,1.0))
                arrow_actor.SetScale(force_disp_scale*elem_force_mag)
                arrow_actor.SetUserTransform(elem_transform)
            
        # Add time text
        text.SetInput('t %4.1f ms'%(log.t[i],))
        
        # Reset clipping range and render
        ren.ResetCameraClippingRange() 
        ren_win.Render()

        image_name = './images/frame_{:06d}'.format(i)
        write_image(image_name, ren_win)

        if False and i==0:
            raw_input('press enter to start')

        #sleep(0.01)
