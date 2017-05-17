#! @Chapter The Orlik-Solomon bicomplex of a bi-arrangement

#! @Section Constructors

#! @Description
#!  Return the Orlik-Solomon bicomplex of the bi-arrangement <A>m</A>
#!  with coloring function <A>chi</A>.
#! @Arguments m, chi
#! @Returns a record
DeclareOperation( "OrlikSolomonBicomplex",
        [ IsMatroid, IsFunction ] );

DeclareOperation( "Matroid",
        [ IsList, IsHomalgRing ] );
