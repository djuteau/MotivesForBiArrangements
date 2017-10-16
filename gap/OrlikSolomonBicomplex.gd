#! @Chapter The Orlik-Solomon bicomplex of a bi-arrangement

#! @Section Constructors

#! @Description
#!  Return the Orlik-Solomon bicomplex of the bi-arrangement <A>m</A>
#!  with coloring function <A>chi</A>.
#! @Arguments m, chi
#! @Returns a record

DeclareOperation( "Matroid",
        [ IsList, IsHomalgRing ] );

DeclareOperation( "FlatsOfRankExtended",
		[ IsMatroid, IsInt ] );

DeclareOperation( "OrlikSolomonBicomplex",
        [ IsMatroid, IsFunction, IsCapCategory ] );

DeclareOperation( "OrlikSolomonBicomplex",
        [ IsMatroid, IsFunction ] );

DeclareOperation( "HasOrlikSolomonBicomplexObject",
		[ IsRecord, IsList, IsInt, IsInt ] );

DeclareOperation( "SetOrlikSolomonBicomplexObject",
		[ IsRecord, IsList, IsInt, IsInt, IsCapCategoryObject ] );

DeclareOperation( "OrlikSolomonBicomplexObject",
		[ IsRecord, IsList, IsInt, IsInt ] );

DeclareGlobalFunction( "HasOrlikSolomonBicomplexDifferentialComponent" );

DeclareGlobalFunction( "SetOrlikSolomonBicomplexDifferentialComponent" );

DeclareGlobalFunction( "OrlikSolomonBicomplexDifferentialComponent" );

DeclareGlobalFunction( "HasOrlikSolomonBicomplexDifferential" );

DeclareGlobalFunction( "SetOrlikSolomonBicomplexDifferential" );

DeclareGlobalFunction( "OrlikSolomonBicomplexDifferential" );

