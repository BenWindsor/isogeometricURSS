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

From here all operators in sp_scalar and msh operations can be called as usual.

### Whats missing
I have yet to implement periodic basis splines of degree greater than 3, and there is also no analogue of the sp_vector space so those operations are not available in vector forms. 
