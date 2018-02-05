#! @Chapter The Orlik-Solomon bicomplex of a bi-arrangement

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
#!  The arguments are a matroid <A>m</A>, and integer <A>k</A> and a boolean <A>b</A> representing
#!  the default color. The output is a coloring function on the matroid <A>m</A> such that
#!  the first <A>k</A> hyperplanes are blue, the others are red, and a stratum which can be expressed
#!  as an intersection of blue hyperplanes is blue, and similarly for red; the strata whose color is
#!  not determined by those conditions are set to the default color.
#! @Arguments m, k default
#! @Returns a coloring function
DeclareOperation( "Coloring",
	[ IsMatroid, IsInt, IsBool ] );

#! @BeginGroup OrlikSolomonBicomplexRecord
#! @Label for IsMatroid, IsFunction[, IsCapCategory]
#! @Description
#!  In the first form, returns the Orlik-Solomon bicomplex of the bi-arrangement <A>m</A>
#!  with coloring function <A>chi</A>, as a bicomplex of objects in the
#! category <A>cat</A>. If no category is specified, then
#! by default <C>MatrixCategory( HomalgFieldOfRationals() )</C> is used.
#! Alternatively, the biarrangement may be specified by a list of equations for the blue
#! hyperplanes, a list of equations for the red hyperplanes, and a default color
#! for the strata that are not obviously blue or red.
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
        [ IsList, IsList, IsBool ] );
#! @EndGroup

#! @Description
#! The arguments are an Orlik-Solomon bicomplex record and a list representing a flat.
#! The output is a bicomplex.
#! @Arguments A, S
#! @Returns a bicomplex 
DeclareOperation( "OrlikSolomonBicomplex",
        [ IsRecord, IsList ] );
        
DeclareOperation( "RedMultizetaBiOS",
        [ IsList, IsCapCategory ] );

DeclareOperation( "RedMultizetaBiOS",
        [ IsList ] );

DeclareOperation( "BlueMultizetaBiOS",
        [ IsList, IsCapCategory ] );

DeclareOperation( "BlueMultizetaBiOS",
        [ IsList ] );

DeclareOperation( "MultizetaBiOS",
        [ IsList, IsCapCategory ] );

DeclareOperation( "MultizetaBiOS",
        [ IsList ] );
        
DeclareOperation( "CellularBiArrangement",
		[ IsInt, IsPerm, IsBool ] );

DeclareOperation( "CellularBiArrangement",
		[ IsInt, IsPerm ] );
		
DeclareOperation( "CellMotive",
		[ IsList, IsBool ] );

DeclareOperation( "CellMotive",
		[ IsList ] );


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

#! @Description
#! The arguments are a matroid <C>m</C> and an integer $k$.
#! The output is the list of flats of rank $k$ in <C>m</C>.
#! The function uses <C>FlatsOfRank</C>, but also returns
#! a result (the empty list) if $k$ is too large or negative.
#! @Arguments m, k
#! @Returns a list of flats
DeclareOperation( "FlatsOfRankExtended",
		[ IsMatroid, IsInt ] );

#! @Description
#! The argument is a boolean (<C>true</C>, <C>false</C> or <C>fail</C>).
#! The output is the color it represents (blue, red or black respectively).
#! In terms of sheaf operations, they correspond to $j_*$, $j_!$ and $j_{!*}$.
#! @Arguments b
#! @Returns a string 
DeclareOperation( "ColorBool",
		[ IsBool ] );

#! @Group DisplayColoring
#! @Label for IsMatroid, IsFunction | IsOrlikSolomonBicomplexRecord
#! @Description
#! The arguments are a matroid and a coloring function in the first version,
#! an Orlik-Solomon bicomplex record in the second version.
#! The function displays the colors of all the flats.
#! The colors are given by a boolean: <C>true</C> means <E>blue</E>,
#! <C>false</C> means <E>red</E>, <C>fail</C> means <E>black</E>.
#! @Arguments m, chi
#! @Returns displays all flats and their colors
DeclareOperation( "DisplayColoring",
		[ IsMatroid, IsFunction ] );

#! @Group DisplayColoring
#! @Arguments A
DeclareOperation( "DisplayColoring",
		[ IsRecord ] );
#		[ IsOrlikSolomonBicomplexRecord ] );


#! @BeginGroup OrlikSolomonBicomplexObject
#! @Description
#!  The arguments are an Orlik-Solomon bicomplex record <A>A</A>, a list <A>S</A> representing a flat,
#!  and integers <A>i</A> and <A>j</A> denoting a position in the bicomplex.
#!  The output is the object representing the local contribution of the flat at the position $(i,j)$.
#!  Although it is not an attribute because there are several arguments, it has the role of an attribute,
#!  and the corresponding tester and setter are implemented. For the setter, the additional argument
#!  <A>V</A> is the vector space to be stored.
#! @Arguments A, S, i, j
#! @Returns a vector space (or an object in the CAP abelian category attached to <A>A</A>).
DeclareOperation( "OrlikSolomonBicomplexObject",
		[ IsRecord, IsList, IsInt, IsInt ] );

#! @Arguments A, S, i, j
DeclareOperation( "HasOrlikSolomonBicomplexObject",
		[ IsRecord, IsList, IsInt, IsInt ] );

#! @Arguments A, S, i, j, V
DeclareOperation( "SetOrlikSolomonBicomplexObject",
		[ IsRecord, IsList, IsInt, IsInt, IsCapCategoryObject ] );

#! @EndGroup

#! @BeginGroup OrlikSolomonBicomplexDifferentialComponent_group
#! @Description
#!  The arguments are an Orlik-Solomon bicomplex record <A>A</A>, a list <A>S</A> (resp. <A>T</A>) representing a flat,
#!  and integers <A>i</A> and <A>j</A> (resp. <A>k</A> and <A>l</A>) denoting a position in the bicomplex.
#!  The output is the morphism between the local contributions of the flats <A>S</A> and <A>T</A>
#!  at the respective positions $(i,j)$ and $(k,l)$.
#!  Although it is not an attribute because there are several arguments, it has the role of an attribute,
#!  and the corresponding tester and setter are implemented. For the setter, the additional argument
#!  <A>f</A> is the linear map to be stored.
#! @Arguments A, S, i, j, T, k, l
#! @Returns a linear map (or a morphism in the CAP abelian category attached to <A>A</A>).
DeclareGlobalFunction( "OrlikSolomonBicomplexDifferentialComponent" );

#! @Arguments A, S, i, j, T, k, l
DeclareGlobalFunction( "HasOrlikSolomonBicomplexDifferentialComponent" );

#! @Arguments A, S, i, j, T, k, l, f
DeclareGlobalFunction( "SetOrlikSolomonBicomplexDifferentialComponent" );
#! @EndGroup

#! @Description
#!  The arguments are an Orlik-Solomon bicomplex record <A>A</A>, a list <A>S</A> representing a flat,
#!  and integers <A>i</A> and <A>j</A> (resp. <A>k</A> and <A>l</A>) denoting a position in the bicomplex.
#!  The output is the morphism between the the objects at the respective positions
#!  $(i,j)$ and $(k,l)$ in the global Orlik-Solomon bicomplex for the flat <A>S</A>.
#!
#!  A tester and setter were implemented, but it turned out to be more efficient not to store the result
#!  (if non-zero, it is a direct sum of stored components).
#! @Arguments A, S, i, j, k, l
#! @Returns a linear map (or a morphism in the CAP abelian category attached to <A>A</A>).
DeclareGlobalFunction( "OrlikSolomonBicomplexDifferential" );

#DeclareGlobalFunction( "HasOrlikSolomonBicomplexDifferential" );

#DeclareGlobalFunction( "SetOrlikSolomonBicomplexDifferential" );

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

