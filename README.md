Application of FEM on beam element using the python codes

( I.)	The lengths of the beam elements are taken as input element wise.
( II. )	The value of Elastic Modulus and Moment of Inertia are also taken as input. (In case  the beam has variying Modulus of Rigidity , then the value of Modulus of Rigidity is taken element wise. A separate code has been provided)
( III.	) Now, the Loading type i.e. udl or concentrated load at mid point ,is selected and the value of load is entered. 
( IV.  )	If settlement is there. Then node numbers where the settlement is occuring is taken and the value of settlement occuring at those specified nodes are entered sequentially. Downward displacement is taken as negative for sign convention. For the settlement a separate code has also been provided.
( V.	) Now, the boundary conditions are applied in the code. The Degrees of freedom where displacement is Zero is taken as input while running the code. 
( VI.	) As, Displacement and rotation at each node has been obtained using the FEM for getting the Moments at the joints Slope Deflection Method has been used.
