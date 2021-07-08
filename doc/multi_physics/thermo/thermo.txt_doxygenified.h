/**

\mainpage Thermoelasticity: Combining the Heat equation and non-linear Solid Mechanics

\image html undef.gif "Undeformed configuration of an elastic block. " 
\image latex undef.eps "Undeformed configuration of an elastic block. " width=0.75\textwidth

 We consider the uniform steady thermal expansion of an elastic
 body that is differentially heated. The top surface is heated and the bottom
 surface is maintained at the reference temperature, which leads to
 a uniform temperature gradient throughout the material. The material
 expands more near the upper surface than near the lower surface, deforming
 the initially rectangular block into an curved configuration.

\image html def.gif "Deformed elastica when the top surface is heated and the lower surface is maintained at the reference temperature. The contours indicate the temperature of the body. " 
\image latex def.eps "Deformed elastica when the top surface is heated and the lower surface is maintained at the reference temperature. The contours indicate the temperature of the body. " width=0.75\textwidth

<HR>
<HR>
\section main The driver code
 
The driver code for this example is given below:

\include thermo.cc

<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/multi_physics/thermo/
">
demo_drivers/multi_physics/thermo/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/multi_physics/thermo/thermo.cc
">
demo_drivers/multi_physics/thermo/thermo.cc
</A>
</CENTER>
.


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

