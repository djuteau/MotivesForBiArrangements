#! @Chapter The Orlik-Solomon bicomplex of a bi-arrangement

# DeclareGlobalFunction( "ORLIK_SOLOMON_INTERNAL_NR_OF_FLAT" );

####################################
##
#! @Section GAP Categories
##
####################################

#! @Description
#! The GAP category of Orlik-Solomon bicomplex records.
#! Those are records containing a matroid, a coloring function,
#! a vector space for each flat and an arrow for each adjacent
#! pair of flats.
#! @Arguments record

DeclareCategory( "IsOrlikSolomonBicomplexRecord",
                 IsRecord );
                 
#! @Description
#! The GAP category of Orlik-Solomon bicomplex records of projective
#! bi-arrangements.
#! Those are records containing a matroid, a coloring function,
#! a vector space for each flat and an arrow for each adjacent
#! pair of flats.
#! They admit the properties <C>IsBlueExact</C>, <C>IsRedExact</C> and
#! <C>IsBiExact</C>.
#! @Arguments record

DeclareCategory( "IsProjectiveOrlikSolomonBicomplexRecord",
                 IsOrlikSolomonBicomplexRecord );

####################################
##
#! @Section Constructors
##
####################################


#! @Description
#!  Creates a matroid using the constructor from the alcove package.
#!  The input is a list <A>L</A> of vectors representing the linear forms
#!  defining the hyperplanes, and an <C>Homalg</C> ring <A>R</A>.
#! @Arguments L, R
#! @Returns a matroid
DeclareOperation( "Matroid",
        [ IsList, IsHomalgRing ] );
        
#! @Description
#!  The argument is an Orlik-Solomon bicomplex record <A>A</A>. The function computes the corresponding graded motive
#!  and adds it to the record, assuming the bicomplex is blue or red exact (according to the color of the zero flat).
#! @Arguments A
#! @Returns nothing.
DeclareOperation( "SemisimplifiedMotive",
	[ IsRecord ] );

#! @Description
#!  The arguments are a matroid <A>m</A>, and integer <A>k</A>, a boolean <A>default</A> representing
#!  the default color, and a boolean <A>last</A> representing the color of the zero flat.
#!  The output is a coloring function on the matroid <A>m</A> such that
#!  the first <A>k</A> hyperplanes are blue, the others are red, and a non-zero flat which can be expressed
#!  as an intersection of blue hyperplanes is blue, and similarly for red; a non-zero flat whose color is
#!  not determined by those conditions is set to the <A>default</A> color, and the zero flat is set to <A>last</A>.
#! @Arguments m, k, default, last
#! @Returns a coloring function
DeclareOperation( "Coloring",
	[ IsMatroid, IsInt, IsBool, IsBool ] );

#! @BeginGroup OrlikSolomonBicomplexRecord
#! @Label for IsMatroid, IsFunction[, IsCapCategory]
#! @Description
#!  In the first form, returns the Orlik-Solomon bicomplex of the bi-arrangement <A>m</A>
#!  with coloring function <A>chi</A>, as a bicomplex of objects in the
#! category <A>cat</A>. If no category is specified, then
#! by default <C>MatrixCategory( HomalgFieldOfRationals() )</C> is used.
#! Alternatively, the biarrangement may be specified by a list of equations for the blue
#! hyperplanes, a list of equations for the red hyperplanes, a default color
#! for the non-zero flats that are not obviously blue or red, and a color for the zero flat.
#! @Arguments m, chi, cat
#! @Returns a record
DeclareOperation( "OrlikSolomonBicomplexRecord",
        [ IsMatroid, IsFunction, IsCapCategory ] );

#! @Description
#! @Arguments m, chi
#! @Returns a record
DeclareOperation( "OrlikSolomonBicomplexRecord",
        [ IsMatroid, IsFunction ] );
        
#! @Description
#! @Arguments L, M, default
#! @Returns a record
DeclareOperation( "OrlikSolomonBicomplexRecord",
        [ IsList, IsList, IsBool, IsBool ] );
#! @EndGroup
        
# DeclareOperation( "RedMultizetaBiOS",
#        [ IsList, IsCapCategory ] );

# DeclareOperation( "RedMultizetaBiOS",
#        [ IsList ] );

# DeclareOperation( "BlueMultizetaBiOS",
#        [ IsList, IsCapCategory ] );

# DeclareOperation( "BlueMultizetaBiOS",
#        [ IsList ] );

# DeclareOperation( "MultizetaBiOS",
#        [ IsList, IsCapCategory ] );

# DeclareOperation( "MultizetaBiOS",
#        [ IsList ] );

#! @Description
#!  The arguments are an integer <A>n</A> and a permutation <A>w</A> on $n$ letters.
#!  The output is a matrix giving the equations of the hyperplanes of the corresponding
#!  cellular arrangement.
#! @Arguments n, w
#! Returns a matrix
DeclareOperation( "CellularArrangement",
		[ IsInt, IsPerm ] );

#! @Description
#!  The argument is a list of integers <A>a</A>. The output is a matrix giving the equations
#!  of the hyplerplanes of the corresponding iterated integral arrangement.
#! @Arguments a
#! @Returns a matrix
DeclareOperation( "IteratedIntegralArrangement",
		[ IsList ] );

#! @Description
#!  The argument is an integer <A>n</A>. The output is a matrix giving the equations
#!  of the hyplerplanes of the standard simplex in the $n$-dimensional projective space.
#! @Arguments n
#! @Returns a matrix
DeclareOperation( "SimplexArrangement",
		[ IsInt ] );

#! @Description
#!  The arguments are a list <A>w</A> representing a permutation, and two Booleans representing
#!  the default color and the color of the zero flat. The output is the corresponding
#!  cellular integral Orlik-Solomon bicomplex record.
#! @Arguments w, default, last
#! Returns an Orlik-Solomon bicomplex record
DeclareOperation( "CellularIntegralOrlikSolomonBicomplexRecord",
		[ IsList, IsBool, IsBool ] );

#! @Description
#!  The arguments are a list of integers <A>a</A>, and two Booleans representing
#!  the default color and the color of the zero flat. The output is the corresponding
#!  iterated integral Orlik-Solomon bicomplex record.
#! @Arguments a, default, last
#! Returns an Orlik-Solomon bicomplex record
DeclareOperation( "IteratedIntegralOrlikSolomonBicomplexRecord",
		[ IsList, IsBool, IsBool ] );


####################################
##
#! @Section Attributes
##
####################################


####################################
##
#! @Section Operations
##
####################################

#! @BeginGroup OrlikSolomonBicomplexObject
#! @Description
#!  The arguments are an Orlik-Solomon bicomplex record <A>A</A>, a list <A>S</A> representing a flat,
#!  and integers <A>i</A> and <A>j</A> denoting a position in the bicomplex.
#!  The output is the object at the position $(i,j)$ in the local Orlik-Solomon bicomplex
#!  of the flat $S$ of codimension $k$ and number $s$.
#!  Although it is not an attribute because there are several arguments, it has the role of an attribute,
#!  and the corresponding tester and setter are implemented. For the setter, the additional argument
#!  <A>V</A> is the vector space to be stored. On the other hand, in that case $k = i + j $ is omitted.
#! @Arguments A, FG, i, j, k, s
#! @Returns a vector space (or an object in the CAP abelian category attached to <A>A</A>).
DeclareOperation( "OrlikSolomonBicomplexObject",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ] );

#! @Arguments A, FG, i, j, k, s
DeclareOperation( "HasOrlikSolomonBicomplexObject",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ] );

#! @Arguments A, FG, i, j, s, V
DeclareOperation( "SetOrlikSolomonBicomplexObject",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsCapCategoryObject ] );

#! @EndGroup

#! @BeginGroup OrlikSolomonBicomplexDifferentialComponent_group
#! @Description
#!  The arguments are an Orlik-Solomon bicomplex record <A>A</A>, a string <A>FG</A>,
#!  integers <A>i</A> and <A>j</A> denoting a position in the bicomplex,
#!  and integers <A>s</A> and <A>t</A> denoting two flats $S$ and $T$ of codimension $i+j$,
#!  The output is the morphism starting at position $(i,j)$ between the local contributions
#!  of the flats $S$ and $T$.
#!  Although it is not an attribute because there are several arguments, it has the role of an attribute,
#!  and the corresponding tester and setter are implemented. For the setter, the additional argument
#!  <A>f</A> is the linear map to be stored.
#! @Arguments A, FG, i, j, s, t
#! @Returns a linear map (or a morphism in the CAP abelian category attached to <A>A</A>).
DeclareOperation( "OrlikSolomonBicomplexHorizontalDifferentialComponent",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ] ); 

DeclareOperation( "OrlikSolomonBicomplexVerticalDifferentialComponent",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ] ); 

DeclareOperation( "HasOrlikSolomonBicomplexHorizontalDifferentialComponent",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ] ); 

DeclareOperation( "HasOrlikSolomonBicomplexVerticalDifferentialComponent",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ] );

#! @Arguments A, FG, i, j, s, t, f
DeclareGlobalFunction( "SetOrlikSolomonBicomplexHorizontalDifferentialComponent" );

DeclareGlobalFunction( "SetOrlikSolomonBicomplexVerticalDifferentialComponent" );

#! @EndGroup

#! @BeginGroup Differentials
#! @Description
#!  The arguments are an Orlik-Solomon bicomplex record <A>A</A>,
#!  a string <A>FG</A> equal to "F" or "G",
#!  integers <A>i</A> and <A>j</A> denoting a position,
#!  an integer <A>k</A> denoting the codimension of the flat S,
#!  and an integer <A>s</A> denoting the number of S
#!  in the list of flats of codimension <A>k</A>.
#!  The output is the horizontal or vertical morphism starting at the position
#!  $(i,j)$ in the local Orlik-Solomon bicomplex for the flat <A>S</A>, corresponding to the
#!  the sheaf <A>FG</A>: for "F", all the black strata are set red, while for "G" all the black
#! strata are set blue.
#!
#!  A tester and setter were implemented, but it turned out to be more efficient not to store the result
#!  (if non-zero, it is a direct sum of stored components).
#! @Arguments A, FG, i, j, k, s
#! @Returns a linear map (or a morphism in the CAP abelian category attached to <A>A</A>).
DeclareOperation( "OrlikSolomonBicomplexHorizontalDifferential",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ]);

DeclareOperation( "OrlikSolomonBicomplexVerticalDifferential",
		[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ]);
#! @EndGroup

#! @BeginGroup Morphisms
#! @Description
#!  Returns the morphism from the "F" to the "G" version of the object at $(i,j)$
#!  for the local Orlik-Solomon bicomplex for the stratum S defined by its codimension <A>k</A>
#!  and its numnber <A>s</A>. In the setter, $k = i + j$ is not needed, but the linear map to be stored
#!  is required.
#! @Arguments A, FG, i, j, k, s
#! @Returns a linear map (or a morphism in the CAP abelian category attached to <A>A</A>).
DeclareOperation( "OrlikSolomonBicomplexMorphism",
		[ IsRecord, IsInt, IsInt, IsInt, IsInt ]);

#! @Description
#! @Arguments A, FG, i, j, s, f
DeclareOperation( "SetOrlikSolomonBicomplexMorphism",
		[ IsRecord, IsInt, IsInt, IsInt, IsCapCategoryMorphism ]);
#! @EndGroup

#! @Description
#!  The arguments are a matroid <A>mat</A>, a coloring function <A>chi</A>, and
#!  the codimension <A>k</A> and number <A>s</A> of a flat $S$. The function says
#!  whether the local bi-arrangement of the flat $S$ is tame.
#! @Returns a Boolean
#! @Arguments mat, chi, k, s
DeclareOperation( "IsTameFlat",
	[ IsMatroid, IsFunction, IsInt, IsInt ] );
	

#! @Description
#! 	The arguments are a matroid <A>mat</A> and a coloring function <A>chi</A>.
#! 	The output says whether the bi-matroid is tame. If not, the function also displays
#! 	the flats which are not tame, along with their codimension and number.
#! @Returns a Boolean
#! @Arguments mat, chi
DeclareOperation( "IsTame",
 	[ IsMatroid, IsFunction ] );

#! @Description
#!  The arguments are two matrices (lists of equations for the blue/red hyperplanes), and two
#!  Booleans representing the default color and color of the zero flat.
#!  The output says whether the bi-arrangement is tame. If not, the function also displays
#!  the flats which are not tame, along with their codimension and number.
#! @Arguments L, M, default, last
#! @Returns a Boolean
DeclareOperation( "IsTameBiarrangement",
	[ IsList, IsList, IsBool, IsBool ] );

#! @Description
#! @Arguments w, default, last
DeclareOperation( "IsTameCellularIntegralBiarrangement",
	[ IsList, IsBool, IsBool ] );

#! @Description
#! @Arguments a, default, last
DeclareOperation( "IsTameIteratedIntegralBiarrangement",
	[ IsList, IsBool, IsBool ] );

####################################
##
#! @Section Obsolete
##
####################################

#! @BeginGroup OrlikSolomonBicomplexHomologyObjects
#! @Description
#! The arguments are an Orlik-Solomon bicomplex record <A>A</A>, a list <A>S</A> representing
#! a flat, and integers <A>i</A> and <A>j</A>.
#! The output is the horizontal or vertical homology object at the position $(i,j)$ in the Orlik-Solomon bicomplex
#! of the local bi-arrangement of flats less than or equal to <A>S</A>.
#! @Arguments A, S, i, j
#! @Returns a vector space
DeclareOperation( "OrlikSolomonBicomplexHorizontalHomologyObject",
		[ IsRecord, IsList, IsInt, IsInt ] );

#! @Arguments A, S, i, j
DeclareOperation( "OrlikSolomonBicomplexVerticalHomologyObject",
		[ IsRecord, IsList, IsInt, IsInt ] );
#! @EndGroup

#! @BeginGroup IsBlueOrRedExact
#! @Description
#! The arguments are an Orlik-Solomon bicomplex record and a list representing a flat.
#! The output is a boolean telling us whether the flat is blue or red exact respectively.
#! @Arguments A, S
#! @Returns a boolean
DeclareOperation( "IsBlueExact",
		[ IsRecord, IsList ] );
		
#! @Arguments A, S
DeclareOperation( "IsRedExact",
		[ IsRecord, IsList ] );
#! @EndGroup

#! @BeginGroup BlueOrRedNonExactness
#! @Description
#! The arguments are an Orlik-Solomon bicomplex record and a list representing a flat.
#! The output is a matrix recording the dimensions of the horizontal or vertical homology of the
#! bicomplex not living on the diagonal (thus the defect of blue or red exactness).
#! @Arguments A, S
#! @Returns a matrix
DeclareOperation( "BlueNonExactness",
		[ IsRecord, IsList ] );
		
#! @Arguments A, S
DeclareOperation( "RedNonExactness",
		[ IsRecord, IsList ] );
#! @EndGroup

#! @Description
#! The arguments are an Orlik-Solomon bicomplex record <A>A</A> and a list <A>S</A> representing a flat.
#! The output is a matrix recording the dimensions of the components of the local
#! bi-arrangement of flats less than or equal to <A>S</A>.
#! @Arguments A, S
#! @Returns a matrix
DeclareOperation( "OrlikSolomonBicomplexDimensions",
		[ IsRecord, IsList ] );

#! @Description
#! The arguments is a matrix recording the dimensions of a projective Orlik-Solomon bicomplex.
#! Assuming bi-exactness, the output is the graded dimension of the corresponding motive
#! (each component is computed as the Euler characteristic of the total complex associated with
#! a rectangular subquotient of the bicomplex).
#! @Arguments mat
DeclareOperation( "Euler",
		[ IsList ] );
		
#! @Description
#! The arguments are an Orlik-Solomon bicomplex record and a list representing a flat.
#! The output is a bicomplex.
#! @Arguments A, S
#! @Returns a bicomplex 
DeclareOperation( "OrlikSolomonBicomplex",
        [ IsRecord, IsList ] );
		


