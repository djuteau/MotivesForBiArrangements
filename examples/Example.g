LoadPackage( "MotivesForBiArrangements" );


A := BlueMultizetaBiOS( [ 5 ] );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
IsBlueExact( A, A.Smin );

B := RedMultizetaBiOS( [ 5 ] );;
PrintArray( OrlikSolomonBicomplexDimensions( B, B.Smin ) );
IsRedExact( B, B.Smin );

# Benchmarking (before removing all the zeros... now the run times are a bit better)

A := BlueMultizetaBiOS( [ 2 ], MatrixCategory( HomalgFieldOfRationals( ) ) );; time;				#  618
A := BlueMultizetaBiOS( [ 2 ], LeftPresentations( HomalgFieldOfRationals( ) ) );; time;				#  840
A := BlueMultizetaBiOS( [ 2 ], CategoryOfHomalgLeftModules( HomalgFieldOfRationals( ) ) );; time;	# 3129

A := BlueMultizetaBiOS( [ 3 ], MatrixCategory( HomalgFieldOfRationals( ) ) );; time;				#  4131
A := BlueMultizetaBiOS( [ 3 ], LeftPresentations( HomalgFieldOfRationals( ) ) );; time;				#  5318
A := BlueMultizetaBiOS( [ 3 ], CategoryOfHomalgLeftModules( HomalgFieldOfRationals( ) ) );; time;	# 19126

A := BlueMultizetaBiOS( [ 4 ], MatrixCategory( HomalgFieldOfRationals( ) ) );; time;				#  30947
A := BlueMultizetaBiOS( [ 4 ], LeftPresentations( HomalgFieldOfRationals( ) ) );; time;				#  45738
A := BlueMultizetaBiOS( [ 4 ], CategoryOfHomalgLeftModules( HomalgFieldOfRationals( ) ) );; time;	# 130685

A := BlueMultizetaBiOS( [ 5 ], MatrixCategory( HomalgFieldOfRationals( ) ) );; time;				#  323539
A := BlueMultizetaBiOS( [ 5 ], LeftPresentations( HomalgFieldOfRationals( ) ) );; time;				#  775797
A := BlueMultizetaBiOS( [ 5 ], CategoryOfHomalgLeftModules( HomalgFieldOfRationals( ) ) );; time;	# 


OSStatistics := function( A )
	local obj, mor, null_obj, null_mor, o, f;

	obj := 0;
	null_obj := 0;
	
	mor := 0;
	null_mor := 0;

	for s in RecNames( A ) do
		if IsCapCategoryObject( A.(s) ) and HasCapCategory( A.(s) ) then
			obj := obj + 1;
			if IsZero( A.(s) ) then
				null_obj := null_obj + 1;
			fi;
		fi;
		if IsCapCategoryMorphism( A.(s) ) and HasCapCategory( A.(s) ) then
			mor := mor + 1;
			if IsZero( A.(s) ) then
				null_mor := null_mor + 1;
			fi;
		fi;
	od;

	Print( null_obj, " out of ", obj, " objects are zero.\n" );
	Print( null_mor, " out of ", mor, " morphisms are zero.\n" );
	
end;


A := BlueMultizetaBiOS( [ 2 ] );; OSStatistics( A );
# 24 out of 45 objects are zero.
# 233 out of 300 morphisms are zero.

A := BlueMultizetaBiOS( [ 3 ] );; OSStatistics( A );
# 107 out of 170 objects are zero.
# 2612 out of 2908 morphisms are zero.

A := BlueMultizetaBiOS( [ 4 ] );; OSStatistics( A );
# 382 out of 565 objects are zero.
# 19569 out of 20678 morphisms are zero.

A := BlueMultizetaBiOS( [ 5 ] );; OSStatistics( A );
# 1219 out of 1739 objects are zero.
# 113421 out of 117215 morphisms are zero.
 



################

Q := HomalgFieldOfRationals( );

m := [
      [ 1 ],
      ];

m := Matroid( m, Q );

chi := function( flat )
    if flat = [ ] then
        Error( "\n" );
    fi;
	return false;
end;

A := OrlikSolomonBicomplex( m, chi );

OrlikSolomonBicomplexDimensions( A, [ 1 ] );

################

Q := HomalgFieldOfRationals( );

m := [
      [ 1, 0 ],
      [ 0, 1 ],
      [ 1, 1 ],
      [ 1, -1 ]
      ];

m := Matroid( m, Q );

chi := function( flat )
    if flat = [ ] then
        Error( "\n" );
    fi;
    if flat = [ 3 ] then
        return false;
    fi;
    return true;
end;

A := OrlikSolomonBicomplex( m, chi );

################

Q := HomalgFieldOfRationals( );

M := [
	[ 1, 0, 0 ],
	[ 0, 1, 0 ],
	[ 0, 0, 1 ],
	[ 1, 1, 1 ]
	];

M := Matroid( M, Q );

chi := function( flat )
    if flat = [ ] then
        Error( "\n" );
    fi;
    if flat = [ 4 ] then
        return false;
    fi;
    return true;
end;

A := OrlikSolomonBicomplex( M, chi );

Length( RecNames( A ) );

for s in RecNames( A ) do Print( s, "\n" ); Display( A.(s) ); Print( "\n" ); od;


################

Q := HomalgFieldOfRationals( );

M := [
	[ 1, 0, 0 ],
	[ 0, 1, 0 ],
	[ 0, 0, 1 ]
	];

M := Matroid( M, Q );

chi := function( flat )
    if flat = [ ] then
        Error( "\n" );
    fi;
    if flat = [ 4 ] then
        return false;
    fi;
    return true;
end;

A := OrlikSolomonBicomplex( M, chi );

names := RecNames( A );

l := Filtered( names, s -> StartsWith( s, "[ 1, 2, 3, 4 ],2,0" ) );

###################

u_1 = (-1,1,0,0)
u_2 = (-1,0,1,0)
u_3 = (0,0,1,0)
u_4 = (0,0,0,1)

vecteurs rouges :
u_5 = (0,1,0,0)
u_6 = (0,-1,1,0)
u_7 = (0,0,-1,1)
u_8 = (-1,0,0,1)

Strates de rang 2 dont la couleur est importante :
{1,2,6} bleu
{3,4,7} bleu
{1,3}, {1,4}, {2,3}, {2,4} bleu
{3,5,6} rouge
{2,7,8} rouge
{5,7}, {5,8}, {6,7}, {6,8} rouge

Strates de rang 3 dont la couleur est importante :
{1,2,3,5,6} bleu
{2,3,4,7,8} bleu
{1,2,4,6}, {1,3,4,7} bleu
{3,4,5,6,7} rouge
{1,2,6,7,8} rouge
{3,5,6,8}, {2,5,7,8} rouge

m := [
      [ -1,1,0,0 ],
      [ -1,0,1,0 ],
      [ 0,0,1,0 ],
      [ 0,0,0,1 ],
      
      [ 0,1,0,0 ],
      [ 0,-1,1,0 ],
      [ 0,0,-1,1 ],
      [ -1,0,0,1 ]
      ];

m := Matroid( m, HomalgFieldOfRationals( ) );

chi := function( flat )
	return flat in [ 
		[ 1 ], [ 2 ], [ 3 ], [ 4 ],
		[ 1,2,6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
		[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
		[ 1 .. 8 ]	
	];
end;

A := OrlikSolomonBicomplex( m, chi );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );

[ [   1,   4,   6,   4,   0 ],
  [   4,  12,  12,   4,   0 ],
  [   6,  11,   6,   0,   0 ],
  [   4,   3,   0,   0,   0 ],
  [   1,   0,   0,   0,   0 ] ]

psi := function( flat )
	return flat in [ 
		[ 1 ], [ 2 ], [ 3 ], [ 4 ],
		[ 1, 2, 6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
		[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
	];
end;

A := OrlikSolomonBicomplex( m, psi );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );

[ [   1,   4,   6,   4,   1 ],
  [   4,  12,  12,   4,   0 ],
  [   6,  11,   5,   0,   0 ],
  [   4,   4,   0,   0,   0 ],
  [   0,   0,   0,   0,   0 ] ]
  
  
# Now with [ 1, 4, 5, 8 ] blue

chi2 := function( flat )
	return flat in [ 
		[ 1 ], [ 2 ], [ 3 ], [ 4 ],
		[ 1,2,6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
		[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
		[ 1, 4, 5, 8 ],
		[ 1 .. 8 ]	
	];
end;

A := OrlikSolomonBicomplex( m, chi2 );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );

[ [   1,   4,   6,   4,   0 ],
  [   4,  12,  11,   4,   0 ],
  [   6,  12,   5,   0,   0 ],
  [   4,   4,   0,   0,   0 ],
  [   1,   0,   0,   0,   0 ] ]
  
psi2 := function( flat )
	return flat in [ 
		[ 1 ], [ 2 ], [ 3 ], [ 4 ],
		[ 1, 2, 6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
		[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
		[ 1, 4, 5, 8 ]
	];
end;

A := OrlikSolomonBicomplex( m, psi2 );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );

[ [   1,   4,   6,   4,   1 ],
  [   4,  12,  11,   3,   0 ],
  [   6,  12,   6,   0,   0 ],
  [   4,   4,   0,   0,   0 ],
  [   0,   0,   0,   0,   0 ] ]

# Now with [ 1, 4, 5, 8 ] black

chi3 := function( flat )
	if flat = [ 1, 4, 5, 8 ] then
		return fail;
	else
		return flat in [ 
			[ 1 ], [ 2 ], [ 3 ], [ 4 ],
			[ 1,2,6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
			[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
			[ 1 .. 8 ]	
		];
	fi;
end;

A := OrlikSolomonBicomplex( m, chi3 );;
IsBlueExact( A, A.Smin );
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );

[ [   1,   4,   6,   4,   0 ],
  [   4,  12,  11,   4,   0 ],
  [   6,  11,   5,   0,   0 ],
  [   4,   3,   0,   0,   0 ],
  [   1,   0,   0,   0,   0 ] ]
  
psi3 := function( flat )
	if flat = [ 1, 4, 5, 8 ] then
		return fail;
	else
		return flat in [ 
			[ 1 ], [ 2 ], [ 3 ], [ 4 ],
			[ 1,2,6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
			[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
		];
	fi;
end;
  
A := OrlikSolomonBicomplex( m, psi3 );;
IsRedExact( A, A.Smin );
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );

[ [   1,   4,   6,   4,   1 ],
  [   4,  12,  11,   3,   0 ],
  [   6,  11,   5,   0,   0 ],
  [   4,   4,   0,   0,   0 ],
  [   0,   0,   0,   0,   0 ] ]
  
# Now with 0 also black

phi := function( flat )
	if flat = [ 1, 4, 5, 8 ] or flat = [ 1 .. 8 ] then
		return fail;
	else
		return flat in [ 
			[ 1 ], [ 2 ], [ 3 ], [ 4 ],
			[ 1,2,6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
			[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
		];
	fi;
end;

A := OrlikSolomonBicomplex( m, phi );;
IsRedExact( A, A.Smin );
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );


# 
# chi3 := function( flat )
# 	if flat = [ 1, 4, 5, 8 ] then
# 		return fail;
# 	else
# 		return flat in [ 
# 			[ 1 ], [ 2 ], [ 3 ], [ 4 ],
# 			[ 1,2,6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
# 			[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
# 			[ 1 .. 8 ]	
# 		];
# 	fi;
# end;
# 
# A := OrlikSolomonBicomplex( m, chi3 );;
# IsBlueExact( A, A.Smin );
# PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
# 
# [ [   1,   4,   6,   4,   0 ],
#   [   4,  12,  11,   4,   0 ],
#   [   6,  11,   5,   0,   0 ],
#   [   4,   3,   0,   0,   0 ],
#   [   1,   0,   0,   0,   0 ] ]
#   
# psi2 := function( flat )
# 	return flat in [ 
# 		[ 1 ], [ 2 ], [ 3 ], [ 4 ],
# 		[ 1, 2, 6 ], [ 3, 4, 7 ], [ 1, 3 ], [ 1, 4 ], [ 2, 3 ], [ 2, 4 ],
# 		[ 1, 2, 3, 5, 6 ], [ 2, 3, 4, 7, 8 ], [ 1, 2, 4, 6 ], [ 1, 3, 4, 7 ],
# 		[ 1, 4, 5, 8 ]
# 	];
# end;
# 
# A := OrlikSolomonBicomplex( m, psi2 );;
# PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );


gap> A := BrownMotiveBlue( 5 );;
gap> IsBlueExact( A, A.Smin );
true
gap> PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
[ [  1,  3,  3,  0 ],
  [  3,  5,  3,  0 ],
  [  3,  2,  0,  0 ],
  [  1,  0,  0,  0 ] ]

gap> A := BrownMotiveBlue( 6 );;
gap> IsBlueExact( A, A.Smin );
true
gap> PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
[ [   1,   4,   6,   4,   0 ],
  [   6,  15,  12,   4,   0 ],
  [  14,  21,   6,   0,   0 ],
  [  13,  10,   0,   0,   0 ],
  [   4,   0,   0,   0,   0 ] ]

B := BrownMotive( 6 );;
gap> IsRedExact( B, B.Smin );
true
gap> PrintArray( OrlikSolomonBicomplexDimensions( B, B.Smin ) );
[ [   1,   4,   6,   4,   1 ],
  [   6,  15,  12,   3,   0 ],
  [  14,  21,   7,   0,   0 ],
  [  13,  13,   0,   0,   0 ],
  [   0,   0,   0,   0,   0 ] ]

gap> A := BrownMotiveBlue( 7 );;
gap> IsBlueExact( A, A.Smin );
true
gap> PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
[ [    1,    5,   10,   10,    5,    0 ],
  [   10,   34,   42,   22,    5,    0 ],
  [   41,   91,   63,   12,    0,    0 ],
  [   82,  112,   31,    0,    0,    0 ],
  [   72,   50,    0,    0,    0,    0 ],
  [   22,    0,    0,    0,    0,    0 ] ]
gap> B := BrownMotive( 7 );;
gap> PrintArray( OrlikSolomonBicomplexDimensions( B, B.Smin ) );
[ [    1,    5,   10,   10,    5,    1 ],
  [   10,   34,   42,   22,    4,    0 ],
  [   41,   91,   63,   13,    0,    0 ],
  [   82,  112,   31,    0,    0,    0 ],
  [   72,   72,    0,    0,    0,    0 ],
  [    0,    0,    0,    0,    0,    0 ] ]
gap> IsRedExact( B, B.Smin );
false

gap> A8 := BrownMotiveBlue( 8 );;
Error, reached the pre-set memory limit
(change it with the -o command line option) in
  res[i] := func( C[i] ); at /Users/imjprg/juteau/gap4r8/lib/coll.gi:746 called from 
List( [ 1 .. arg[1] ], function ( i )
      return [  ];
  end ) at /Users/imjprg/juteau/gap4r8/pkg/homalg_project/Gauss/gap/SparseMatrix.gi:422 called from
SparseZeroMatrix( NrRows( C ), NrColumns( C ), R!.ring ) at /Users/imjprg/juteau/gap4r8/pkg/homalg_project/GaussForHomalg/gap/GaussTools.gi:55 called from
RP!.ZeroMatrix( C ) at /Users/imjprg/juteau/gap4r8/pkg/homalg_project/MatricesForHomalg/gap/Tools.gi:1236 called from
Eval( B ) at /Users/imjprg/juteau/gap4r8/pkg/homalg_project/GaussForHomalg/gap/GaussTools.gi:191 called from
RP!.UnionOfColumns( A, B ) at /Users/imjprg/juteau/gap4r8/pkg/homalg_project/MatricesForHomalg/gap/Tools.gi:673 called from
...  at line 1426 of *stdin*
you can 'return;'
