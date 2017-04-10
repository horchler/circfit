circfit
========
##### Matlab function for least squares fitting of two-dimensional data to a circle.
###### Version 1.3, 4-6-16
##### Download Repository: [ZIP Archive](https://github.com/horchler/circfit/archive/master.zip)

--------

[```circfit(X,Y)```](https://github.com/horchler/circfit/blob/master/circfit.m) returns scalar radius ```R``` of a fitted circle. ```X``` and ```Y``` are 1-D arrays of position data in a rectilinear coordinate system. ```X``` and ```Y``` must be the same length and must contain at least three non-colinear points in order for a valid solution to be found. The function can also return position of the center of the fitted circle and the root mean squared error of the fit.  

Additional functions to fit curvature, check colinearity, calculate root mean squared error, and plot the fitted result are included.
&nbsp;  

--------

Andrew D. Horchler, *horchler @ gmail . com*, [biorobots.case.edu](http://biorobots.case.edu/)  
Created: 5-12-07, Revision: 1.3, 4-6-16  

This version tested with Matlab 9.0.0.341360 (R2016a)  
Mac OS X 10.11.4 (Build: 15E65), Java 1.7.0_75-b13  
Compatibility maintained back through Matlab 7.4 (R2007a)  
&nbsp;  

--------

Copyright &copy; 2007&ndash;2017, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.