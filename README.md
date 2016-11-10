# isogeometricURSS
A repository for my code created as part of the URSS at Warwick in 2016

### Installation
In order to install: 

1) Cut the file geo_load.m from the downloaded code and replace the existing geoPDEs file with the same name, this should be located at '[your geoPDEs install]/inst/geometry/geo_load.m'.

2) Add the folder with the remaining isogeometricURSS code to your matlab path including files and subfolders and it should be ready to use.

### Usage
My library provides several common function calls in the geoPDEs library and nurbs library with versions for periodic splines. The most common of these are detailed below. There are other background functions of geoPDEs that have analogues in my library but I shall leave them from this description as they aren't called manually. 

Analogues of the nurbs library:
- perbspmak = an analogue of nrbmak
- perbspplot = an analogue of nrbplot

Analogues of the geoPDEs library:
- geo_load = replaces the geo_load of geoPDEs and can be given an object created by perbspmak or periodicCurveInterpolate
- sp_perbsp = an analogue of sp_nurbs forming a space of periodic basis splines instead of nurbs basis splines
- sp_scalar_perbsp = a periodic basis version of the sp_scalar class

New functionality:
- op_f_v_tp_param = a version of op_f_v_tp that takes values in the parametric domain instead of on the curve/mesh
- op_u_v_tp_param = a version of op_u_v_tp that takes values in the parametric domain instead of on the curve/mesh

From here operators in sp_scalar and msh operations can be called as usual.

### Examples
Many of the files have a [filename]Test.m file alongside them containing my own tests of that code which should be useful for examining functionality. 

The code in /problemAttempts/ has my own use of this library. It was used to solve the equations in the paper 'Modelling cell motility and chemotaxis with evolving surface finite elements' DOI:10.1098/rsif.2012.0276 by Elliot et al. as such many examples of my own attempts feature problems like curve shortening flow and versions thereof with lagrange multipliers designed to keep the area enclosed constant. 

It is also worth examining Listing 1 of the expository paper 'A new design for the implementation of isogeometric analysis in Octave and Matlab: GeoPDEs 3.0' by R. Vazquez as this code is what my library aimed to replicate the structure of, just swapping out the nurbs specific function calls with my own function calls. In particular the analogues are:
- call geo_load on a periodic spline object e.g. the output of periodicCurveInterpolate or perbspmak
- geometry.nurbs.knots becomes geometry.perbspline.knots
- sp nurbs(geometry.nurbs, msh) becomes sp_perbsp(geometry.perbspline, msh)

### Whats missing
I have yet to implement periodic basis splines of degree greater than 3, and there is also no analogue of the sp_vector space so those operations are not available in vector forms. 
I also have yet to implement any sophisticated forms of output so you CANT call things like sp_to_vtk for example. I do have a perbspplot function as an analogue of the nrbplot function.
I have not worked with any boundary conditions yet, other than the periodic conditions and others intrinsic to the problem I was writing this for. 
