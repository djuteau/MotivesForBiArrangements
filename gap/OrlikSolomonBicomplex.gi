####################################
#
# methods for operations:
#
####################################

##
InstallMethod( Matroid,
        [ IsList, IsHomalgRing ],
  
  function( L, R )
    local r, c;
    
    if L = [ ] then
        Error( "empty input\n" );
    fi;
    
    r := Length( L[1] );
    c := Length( L );
    
    return Matroid( HomalgMatrix( TransposedMat( L ), R ) );

end );

InstallMethod( FlatsOfRankExtended,
		[ IsMatroid, IsInt ],
		
	function( matroid, k )
		local rank;
		
		rank := RankOfMatroid( matroid );
		
		if k < 0 or k > rank then
			return [];
		fi;
		
		return FlatsOfRank( matroid, k );
end );

DisplayColoring := function( m, chi )
	local S;
	for S in Concatenation( Flats( m ) ) do
		Print( S, " ", chi( S ), "\n" );
	od;
end;

InstallMethod( HasOrlikSolomonBicomplexObject,
		[ IsRecord, IsList, IsInt, IsInt ],
		
	function( A, S, i, j )

		return IsBound( A.(JoinStringsWithSeparator( [ S, i, j ] )) );

end );

InstallMethod( SetOrlikSolomonBicomplexObject,
		[ IsRecord, IsList, IsInt, IsInt, IsCapCategoryObject ],
		
	function( A, S, i, j, V )
		
		A.(JoinStringsWithSeparator( [ S, i, j ] )) := V;
		
end );

InstallMethod( OrlikSolomonBicomplexObject,
		[ IsRecord, IsList, IsInt, IsInt ],
		
	function( A, Sigma, i, j )
		local V, SS;
		
		if HasOrlikSolomonBicomplexObject( A, Sigma, i, j ) then
			V := A.(JoinStringsWithSeparator( [ Sigma, i, j ] ));
		elif i + j < A.rank( Sigma ) then
			SS := Filtered( FlatsOfRankExtended( A.matroid, i + j ), S -> IsSubset( Sigma, S ) );
			V := DirectSum( A.cat, List( SS, S -> OrlikSolomonBicomplexObject( A, S, i, j ) ) );
			SetOrlikSolomonBicomplexObject( A, Sigma, i, j, V );
		else
			V := ZeroObject( A.cat );
			SetOrlikSolomonBicomplexObject( A, Sigma, i, j, V );
		fi;
		
		return V;
end );

InstallGlobalFunction( HasOrlikSolomonBicomplexDifferentialComponent,
#		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt ],

	function( A, S, i, j, T, k, l )

		return IsBound( A.(JoinStringsWithSeparator( [ S, i, j, T, k, l ] )) );

end );

InstallGlobalFunction( SetOrlikSolomonBicomplexDifferentialComponent,
#		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt, IsCapCategoryMorphism ],

	function( A, S, i, j, T, k, l, f )

		A.(JoinStringsWithSeparator( [ S, i, j, T, k, l ] )) := f;
		
end );

InstallGlobalFunction( OrlikSolomonBicomplexDifferentialComponent,
#		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt ],
		
	function( A, S, i, j, T, k, l )
		local V, W, f;
		
		if HasOrlikSolomonBicomplexDifferentialComponent( A, S, i, j, T, k, l ) then
			f := A.(JoinStringsWithSeparator( [ S, i, j, T, k, l ] ));
		else
			V := OrlikSolomonBicomplexObject( A, S, i, j );
			W := OrlikSolomonBicomplexObject( A, T, k, l );
			f := ZeroMorphism( V, W );
			SetOrlikSolomonBicomplexDifferentialComponent( A, S, i, j, T, k, l, f );
		fi;

		return f;
end );

InstallGlobalFunction( HasOrlikSolomonBicomplexDifferential,
#		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt ],
		
	function( A, Sigma, i, j, k, l )

		return IsBound( A.(JoinStringsWithSeparator( [ Sigma, i, j, k, l ] )) );
				
end );

InstallGlobalFunction( SetOrlikSolomonBicomplexDifferential,
#		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt, IsCapCategoryMorphism ],
		
	function( A, Sigma, i, j, k, l, f )

		A.(JoinStringsWithSeparator( [ Sigma, i, j, k, l ] )) := f;
				
end );

InstallGlobalFunction( OrlikSolomonBicomplexDifferential,
#		[ IsRecord, IsList, IsInt, IsInt, IsInt, IsInt ],
		
	function( A, Sigma, i, j, k, l )
		local V, W, f, SS, TT;
		
		if HasOrlikSolomonBicomplexDifferential( A, Sigma, i, j, k, l ) then
			f := A.(JoinStringsWithSeparator( [ Sigma, i, j, k, l ] ));
		else
            SS := Filtered( FlatsOfRankExtended( A.matroid, i + j ), S -> IsSubset( Sigma, S ) );
            TT := Filtered( FlatsOfRankExtended( A.matroid, k + l ), T -> IsSubset( Sigma, T ) );
			f := MorphismBetweenDirectSums(
            	DirectSum( A.cat, List( SS, S -> OrlikSolomonBicomplexObject( A, S, i, j ) ) ),
                List( SS, S -> List( TT, T -> OrlikSolomonBicomplexDifferentialComponent( A, S, i, j, T, k, l ) ) ),
                DirectSum( A.cat, List( TT, T -> OrlikSolomonBicomplexObject( A, T, k, l ) ) )
 			);
			SetOrlikSolomonBicomplexDifferential( A, Sigma, i, j, k, l, f );
		fi;

		return f;
end );

####################################
#
# methods for constructors:
#
####################################

##
InstallMethod( OrlikSolomonBicomplex,
        [ IsMatroid, IsFunction, IsCapCategory ],
        
  function( m, chi, cat )
    local A, k, i, Sigma, SS, TT, S, T, phi, d, D, s, t, psi, obj;

    A := rec(
    		cat := cat,
    		matroid := m,
    		rank := RankFunction( m ),
    		coloring := chi,
             (JoinStringsWithSeparator( [ [ ], 0, 0 ] )) := TensorUnit( cat ),
    		);
    
    for k in [ 1 .. RankOfMatroid( m ) ] do
        for Sigma in FlatsOfRankExtended( m, k ) do

            SS := Filtered( FlatsOfRankExtended( m, k - 1 ), S -> IsSubset( Sigma, S ) );
            TT := Filtered( FlatsOfRankExtended( m, k - 2 ), T -> IsSubset( Sigma, T ) );
            
            if chi( Sigma ) then
                for i in Reversed( [ 1 .. k ] ) do
                	phi := OrlikSolomonBicomplexDifferential( A, Sigma, i - 1, k - i, i - 2, k - i );
                    d := KernelEmbedding( phi );
                    SetOrlikSolomonBicomplexObject( A, Sigma, i, k - i, Source( d ) );
                    D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i - 1, k - i ) );
                    for s in [ 1 .. Length( SS ) ] do
                    	SetOrlikSolomonBicomplexDifferentialComponent( A, Sigma, i, k - i, SS[s], i - 1, k - i, PreCompose( d, ProjectionInFactorOfDirectSum( D, s ) ) );
                    od;
                    if i < k then
                        D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i, k - i - 1 ) );
                        psi := PreCompose(
                                       OrlikSolomonBicomplexDifferential( A, Sigma, i, k - i - 1, i - 1, k - i - 1 ),
                                       OrlikSolomonBicomplexDifferential( A, Sigma, i - 1, k - i - 1, i - 1, k - i )
                        );
                        psi := KernelLift( phi, psi );
                        for s in [ 1 .. Length( SS ) ] do
                            SetOrlikSolomonBicomplexDifferentialComponent(
                            	A, SS[s], i, k - i - 1, Sigma, i, k - i, PreCompose( InjectionOfCofactorOfDirectSum( D, s ), psi )
                            );
                        od;
                    fi;
                od;
            else
                for i in [ 0 .. k - 1] do
                    phi := OrlikSolomonBicomplexDifferential( A, Sigma, i, k - i - 2, i, k - i - 1 );
                    d := CokernelProjection( phi );
                    SetOrlikSolomonBicomplexObject( A, Sigma, i, k - i, Range( d ) ); #Error("");
                    D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i, k - i - 1 ) ); #Error("");
                    for s in [ 1 .. Length( SS ) ] do
                        SetOrlikSolomonBicomplexDifferentialComponent( A, SS[s], i, k - i - 1, Sigma, i, k - i, PreCompose( InjectionOfCofactorOfDirectSum( D, s ), d ) );
                    od; #Error("");
                    if i > 0 then
                        D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i - 1, k - i ) );
                        psi := PreCompose(
                                       OrlikSolomonBicomplexDifferential( A, Sigma, i, k - i - 1, i - 1, k - i - 1 ),
                                       OrlikSolomonBicomplexDifferential( A, Sigma, i - 1, k - i - 1, i - 1, k - i )
                        );
                        psi := CokernelColift( phi, psi );
                        for s in [ 1 .. Length( SS ) ] do
                            SetOrlikSolomonBicomplexDifferentialComponent(
                            	A, Sigma, i, k - i, SS[s], i - 1, k - i, PreCompose( psi, ProjectionInFactorOfDirectSum( D, s ) )
                            );
                        od;
                    fi;
                od;
            fi;
        od;
    od;
    
    return A;
    
end );

InstallMethod( OrlikSolomonBicomplex,
        [ IsMatroid, IsFunction ],
        
	function( m, chi )
	
		return OrlikSolomonBicomplex( m, chi, MatrixCategory( HomalgFieldOfRationals() ) );

end );

SimplexArrangement := function( n )
	local res, i, j;
	
	res := [ ];
	
	
	for i in [ 1 .. n + 1 ] do
		res[i] := [ ];
		for j in [ 1 .. n + 1 ] do
			res[i][j] := 0;
		od;
	od;

	res[1][2] := 1;

	for i in [ 2 .. n ] do
		res[i][i] := 1;
		res[i][i + 1] := -1;
	od;
	
	res[n + 1][1] := -1;
	res[n + 1][n + 1] := 1;
	
	return res;
end;

Multizeta01Word := function( ni_list )
	local res, n, j;
	res := [ ];
	for n in ni_list do
		Add( res, 1 );
		for j in [ 1 .. n - 1 ] do
			Add( res, 0 );
		od;
	od;
	return res;
end;

IteratedIntegralBlueArrangement := function( a ) 			# a is the list of a_i's
	local n, res, i, j;
	
	n := Length( a );
	
	res := [];
	
	for i in [ 1 .. n + 1 ] do
		res[i] := [];
		for j in [ 1 .. n + 1 ] do
			res[i][j] := 0;
		od;
	od;
	
	res[1][1] := 1;
	
	for i in [ 2 .. n + 1] do
		res[i][1] := a[i - 1];
		res[i][i] := -1;
	od;
	
	return res;
end;

MultizetaBlueArrangement := function( ni_list )
	return IteratedIntegralBlueArrangement( Multizeta01Word ( ni_list ) );
end;

RedMultizetaBiOS := function( ni_list )
	local Q, n, m, rk, chi;
	
	Q := HomalgFieldOfRationals( );
	n := Sum( ni_list );
	
	m := Concatenation( MultizetaBlueArrangement( ni_list ), SimplexArrangement( n ) );
	m := Matroid( m, Q );
	rk := RankFunction( m );
	
	chi := function( flat )
#		return not ForAll( flat, i -> i > n + 1 );
		return not rk( flat ) = Length( Filtered( flat, i -> i > n + 1 ) );
	end;
	
	return OrlikSolomonBicomplex( m, chi );
end;

BlueMultizetaBiOS := function( ni_list )
	local Q, n, m, rk, chi;
	
	Q := HomalgFieldOfRationals( );
	n := Sum( ni_list );
	
	m := Concatenation( MultizetaBlueArrangement( ni_list ), SimplexArrangement( n ) );
	m := Matroid( m, Q );
	rk := RankFunction( m );
	
	chi := function( flat )
		if rk( flat ) = n + 1 then
			return true;
		else
			return not rk( flat ) = Length( Filtered( flat, i -> i > n + 1 ) );
		fi;
	end;
	
	return OrlikSolomonBicomplex( m, chi );
end;

OrlikSolomonBicomplexDimensions := function( A, Sigma )
	local r, res, i, j;
	
	r := A.rank( Sigma );
	
	res := [];
	for i in [ 0 .. r ] do
		res[i + 1] := [];
		for j in [ 0 .. r - i ] do
			res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexObject( A, Sigma, i, j ) );
		od;
		for j in [ r - i + 1 .. r ] do
			res[i + 1][j + 1] := 0;
		od;
	od;
	return res;
end;

OrlikSolomonBicomplexDifferentials := function( A, Sigma )
	local r, hor, ver, i, j;
	
	r := A.rank( Sigma );

	hor := [];
	for i in [ 0 .. r - 1] do
		hor[i + 1] := [];
		for j in [ 0 .. r ] do
			hor[i + 1][j + 1] := OrlikSolomonBicomplexDifferential( A, Sigma, r - i, j, r - i - 1, j );
		od;
	od;
	
	ver := [];
	for i in [ 0 .. r ] do
		ver[i + 1] := [];
		for j in [ 0 .. r - 1 ] do
			ver[i + 1][j + 1] := OrlikSolomonBicomplexDifferential( A, Sigma, r - i, j, r - i, j + 1 );
			Print( "\n" );		od;
	od;
	
	return rec( hor := hor, ver := ver );
end;

# DisplayOrlikSolomonBicomplexDifferentials := function( A, Sigma )
# 	local record, hor, ver, i, j;
# 	
# 	record := OrlikSolomonBicomplexDifferentials( A, Sigma );
# 
# 	hor := [];	
# 	for i in [ 1 .. Length( record.hor ) ] do
# 		hor[i] := [];
# 		for j in [ 1 .. Length( record.hor[i] ) ] do
# 			hor[i][j] := UnderlyingMatrix( record.hor[i][j] );
# 		od; 
# 	od;
# 	
# 	ver := [];	
# 	for i in [ 1 .. Length( record.ver ) ] do
# 		ver[i] := [];
# 		for j in [ 1 .. Length( record.ver[i] ) ] do
# 			ver[i][j] := UnderlyingMatrix( record.ver[i][j] );
# 		od; 
# 	od;
# 	
# 	Print(hor);
# 	Print(ver);
# 	
# end;

DisplayOrlikSolomonBicomplexDifferentials := function( A, Sigma )
	local r, i, j;
	
	r := A.rank( Sigma );

	for j in [ 0 .. r - 1 ] do
		for i in [ 0 .. r - j - 1 ] do
			Print( JoinStringsWithSeparator( [ i + 1, j, i, j ] ), "\n");
			Display( OrlikSolomonBicomplexDifferential( A, Sigma, i + 1, j, i, j ) );
			Print( "\n" );
		od;
	od;
	
	for i in [ 0 .. r - 1 ] do
		for j in [ 0 .. r - i - 1 ] do
			Print( JoinStringsWithSeparator( [ i, j, i, j + 1 ] ), "\n");
			Display( OrlikSolomonBicomplexDifferential( A, Sigma, i, j, i, j + 1 ) );
			Print( "\n" );		od;
	od;
end;
