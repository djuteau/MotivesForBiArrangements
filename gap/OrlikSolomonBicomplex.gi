####################################
#
# methods for constructors:
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

##
InstallMethod( OrlikSolomonBicomplexRecord,
        [ IsMatroid, IsFunction, IsCapCategory ],

  function( m, chi, cat )
    local A, k, i, Sigma, SS, TT, S, T, phi, d, D, s, t, psi, obj, c, b, r;

    A := rec(
    		cat := cat,
    		matroid := m,
    		flats := Flats( m ),
    		rank := RankFunction( m ),
    		coloring := chi,
             (JoinStringsWithSeparator( [ [ ], 0, 0 ] )) := TensorUnit( cat ),
    		);

    A.Smin := A.flats[Length( A.flats )][1];

    for k in [ 1 .. RankOfMatroid( m ) ] do
    	Print( "\nk = ", k, "\n" );
        for Sigma in FlatsOfRankExtended( m, k ) do
        	
#        	Print( ".\c" );

            SS := Filtered( FlatsOfRankExtended( m, k - 1 ), S -> IsSubset( Sigma, S ) );
            TT := Filtered( FlatsOfRankExtended( m, k - 2 ), T -> IsSubset( Sigma, T ) );
            
            c := chi( Sigma );
            
            if c = fail then
             	b := BlueNonExactness( A, Sigma );
             	r := RedNonExactness( A, Sigma );
             	if IsZero( b ) and IsZero( r ) then
             		Print( Sigma, " is black exact.\n" );
             	else
             		if not IsZero( b ) then
             			Print( Sigma, " is black but not blue exact:\n" );
             			PrintArray( b );
             			Print( "\n" );
             		fi;
             		if not IsZero( r ) then
             			Print( Sigma, " is black but not red exact:\n" );
             			PrintArray( r );
             			Print( "\n" );
             		fi;
             	fi;
            elif c = true then
            	b := BlueNonExactness( A, Sigma );
            	if IsZero( b ) then
             		Print( Sigma, " is blue exact.\n" );
             	else
             		if not IsZero( b ) then
             			Print( Sigma, " is blue but not blue exact:\n" );
             			PrintArray( b );
             			Print( "\n" );
             		fi;
             	fi;
            else
             	r := RedNonExactness( A, Sigma );
            	if IsZero( r ) then
             		Print( Sigma, " is red exact.\n" );
             	else
             		if not IsZero( r ) then
             			Print( Sigma, " is red but not red exact:\n" );
             			PrintArray( r );
             			Print( "\n" );
             		fi;
             	fi;
            fi;
            
            if Sigma = A.Smin then
            	return A;
            fi;
             
            if c = fail then
            	for i in [ 1 .. k - 1 ] do
                    psi := PreCompose(
                                    OrlikSolomonBicomplexDifferential( A, Sigma, i, k - i - 1, i - 1, k - i - 1 ),
                                    OrlikSolomonBicomplexDifferential( A, Sigma, i - 1, k - i - 1, i - 1, k - i )
                    );
                	SetOrlikSolomonBicomplexObject( A, Sigma, i, k - i, ImageObject( psi ) );
                	
                	d := ImageEmbedding( psi );
                	D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i - 1, k - i ) );
                	for s in [ 1 .. Length( SS ) ] do
                		SetOrlikSolomonBicomplexDifferentialComponent( A, Sigma, i, k - i, SS[s], i - 1, k - i, PreCompose( d, ProjectionInFactorOfDirectSum( D, s ) ) );
                	od;
                	
                	d := CoastrictionToImage( psi );
                    D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i, k - i - 1 ) );
                    for s in [ 1 .. Length( SS ) ] do
                        SetOrlikSolomonBicomplexDifferentialComponent( A, SS[s], i, k - i - 1, Sigma, i, k - i, PreCompose( InjectionOfCofactorOfDirectSum( D, s ), d ) );
                    od;
             	od;
           
            elif c then
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
                    SetOrlikSolomonBicomplexObject( A, Sigma, i, k - i, Range( d ) );
                    D := List( [ 1 .. Length( SS ) ], s -> OrlikSolomonBicomplexObject( A, SS[s], i, k - i - 1 ) );
                    for s in [ 1 .. Length( SS ) ] do
                        SetOrlikSolomonBicomplexDifferentialComponent( A, SS[s], i, k - i - 1, Sigma, i, k - i, PreCompose( InjectionOfCofactorOfDirectSum( D, s ), d ) );
                    od;
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
    
    Print( Euler( OrlikSolomonBicomplexDimensions( A, A.Smin ) ), "\n" );

    return A;
    
end );



InstallMethod( OrlikSolomonBicomplexRecord,
        [ IsMatroid, IsFunction ],
        
	function( m, chi )
	
		return OrlikSolomonBicomplexRecord( m, chi, MatrixCategory( HomalgFieldOfRationals() ) );

end );

# InstallMethod( ProjectiveOrlikSolomonBicomplexRecord,
#         [ IsMatroid, IsFunction, IsCapCategory ],
# 
# 	function( m, chi, cat )
# 		local ;
# 		
# 		A := OrlikSolomonBicomplexRecord( m, chi, cat );
# 		
# 		
# 
# end );
# 
# InstallMethod( GradedCohomologyOfProjectiveBiarrangement,
# 		[ IsRecord ],
# 
# 	GradedCohomologyOfProjectiveBiarrangement := function( A )
# 		local n, betti, k, B, i, j;
# 		
# 		n := RankOfMatroid( A.matroid ) - 1;
# 		
# 		betti := [];
# 				
# 		for k in [ 0 .. n ] do
# 			B := DoubleCochainComplex(
# 				A.cat,
# 				function( i, j )
# 					if i < -k - 1 or i > 0 then
# 						return IdentityMorphism( ZeroObject( A.cat ) );
# 					elif i = -k - 1 then
# 						return UniversalMorphismFromZeroObject( OrlikSolomonBicomplexObject( A, A.Smin, k, j ) );
# 					elif i = 0 then
# 						return UniversalMorphismIntoZeroObject( OrlikSolomonBicomplexObject( A, A.Smin, 0, j ) );
# 					else
# 						return OrlikSolomonBicomplexDifferential( A, A.Smin, -i, j, -i - 1, j );
# 					fi;
# 				end,
# 				function( i, j )
# 					if j < -1 or j > n - k then
# 						return  ZeroMorphism( A.cat );
# 					elif j = -1 then
# 						return UniversalMorphismFromZeroObject( OrlikSolomonBicomplexObject( A, A.Smin, -i, 0 ) );
# 					elif j = n - k then
# 						return UniversalMorphismIntoZeroObject( OrlikSolomonBicomplexObject( A, A.Smin, -i, n - k ) );
# 					else
# 						return ( -1 )^j * OrlikSolomonBicomplexDifferential( A, A.Smin, -i, j, -i, j + 1 );
# 					fi;
# 				end
# 			);
# 			
# 		SetAboveBound( B, n - k + 1 );
# 		SetBelowBound( B, -1 );
# 		SetRightBound( B, 1 );
# 		SetLeftBound( B, -k - 1 );
# 		
# #		Error( "" );
# 
# 		B := TotalCochainComplex( B );
# 		
# 		betti[k] := List( [ -k, n - k ], i -> Dimension( DefectOfExactnessAt( B, i ) ) );
# 
# 		od;
# 		
# 		return betti;
# 		
# end );



InstallMethod( OrlikSolomonBicomplex,
		[ IsRecord, IsList ],

	function( A, S ) 
		local i, j, r;
		
		if IsBound( A.( String( S ) ) ) then
			return A.( String( S ) );
		fi;
		
		r := RankOfMatroid( A.matroid );

		A.( String( S ) ) := DoubleCochainComplex(
			A.cat,
			function( i, j ) return OrlikSolomonBicomplexDifferential( A, S, -i, j, -i - 1, j ); end,
			function( i, j ) return OrlikSolomonBicomplexDifferential( A, S, -i, j, -i, j + 1 ); end
		);
		SetAboveBound( A.( String( S ) ), r + 1 );
		SetBelowBound( A.( String( S ) ), -1 );
		SetRightBound( A.( String( S ) ), 1 );
		SetLeftBound( A.( String( S ) ), -r - 1 );

	return A.( String( S ) );
end ); 

InstallMethod( OrlikSolomonBicomplex,
		[ IsRecord, IsList ],

	function( A, S ) 
		local i, j, r, res;
		
		r := RankOfMatroid( A.matroid );

		res := DoubleCochainComplex(
			A.cat,
			function( i, j ) return OrlikSolomonBicomplexDifferential( A, S, -i, j, -i - 1, j ); end,
			function( i, j ) return OrlikSolomonBicomplexDifferential( A, S, -i, j, -i, j + 1 ); end
		);
		SetAboveBound( res, r + 1 );
		SetBelowBound( res, -1 );
		SetRightBound( res, 1 );
		SetLeftBound( res, -r - 1 );

	return res;
end ); 


####################################
#
# methods for operations:
#
####################################

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

InstallMethod( ColorBool,
		[ IsBool ],

	function ( b )
		if b = true then
			return "blue";
		elif b = false then
			return "red";
		else
			return "black";
		fi;
end );

InstallMethod( DisplayColoring,
		[ IsMatroid, IsFunction ],

	function( m, chi )
		local S;
		for S in Concatenation( Flats( m ) ) do
			Print( S, " ", chi( S ), "\n" );
		od;
end );

InstallMethod( DisplayColoring,
		[ IsRecord ],

	function( A )
		DisplayColoring( A.matroid, A.coloring );
end );

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
#			SetOrlikSolomonBicomplexObject( A, Sigma, i, j, V );
		else
			V := ZeroObject( A.cat );
#			SetOrlikSolomonBicomplexObject( A, Sigma, i, j, V );
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
#			SetOrlikSolomonBicomplexDifferentialComponent( A, S, i, j, T, k, l, f );
		fi;

		return f;
end );

# InstallGlobalFunction( HasOrlikSolomonBicomplexDifferential,
# #		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt ],
# 		
# 	function( A, Sigma, i, j, k, l )
# 
# 		return IsBound( A.(JoinStringsWithSeparator( [ Sigma, i, j, k, l ] )) );
# 				
# end );
# 
# InstallGlobalFunction( SetOrlikSolomonBicomplexDifferential,
# #		[ IsRecord, IsList, IsInt, IsInt, IsList, IsInt, IsInt, IsCapCategoryMorphism ],
# 		
# 	function( A, Sigma, i, j, k, l, f )
# 
# 		A.(JoinStringsWithSeparator( [ Sigma, i, j, k, l ] )) := f;
# 				
# end );

InstallGlobalFunction( OrlikSolomonBicomplexDifferential,
#		[ IsRecord, IsList, IsInt, IsInt, IsInt, IsInt ],
		
	function( A, Sigma, i, j, k, l )
		local V, W, f, SS, TT;
		
# 		if HasOrlikSolomonBicomplexDifferential( A, Sigma, i, j, k, l ) then
# 			f := A.(JoinStringsWithSeparator( [ Sigma, i, j, k, l ] ));
# 		else
            SS := Filtered( FlatsOfRankExtended( A.matroid, i + j ), S -> IsSubset( Sigma, S ) );
            TT := Filtered( FlatsOfRankExtended( A.matroid, k + l ), T -> IsSubset( Sigma, T ) );
#			f := 
			return MorphismBetweenDirectSums(
            	DirectSum( A.cat, List( SS, S -> OrlikSolomonBicomplexObject( A, S, i, j ) ) ),
                List( SS, S -> List( TT, T -> OrlikSolomonBicomplexDifferentialComponent( A, S, i, j, T, k, l ) ) ),
                DirectSum( A.cat, List( TT, T -> OrlikSolomonBicomplexObject( A, T, k, l ) ) )
 			);
#			SetOrlikSolomonBicomplexDifferential( A, Sigma, i, j, k, l, f );
#		fi;

#		return f;
end );

####################################
#
# methods for constructors:
#
####################################


InstallMethod( OrlikSolomonBicomplexHorizontalHomologyObject,
		[ IsRecord, IsList, IsInt, IsInt ],

	function( A, S, i, j )
  		local alpha, beta, iota, lambda;

		alpha := OrlikSolomonBicomplexDifferential( A, S, i + 1, j, i, j );
		beta := OrlikSolomonBicomplexDifferential( A, S, i, j, i - 1, j );

		if not IsZero( PreCompose( alpha, beta ) ) then
      
	      Error( "the composition of the given morphisms has to be zero" );
      
		fi;
  
		iota := ImageEmbedding( alpha );
  
		lambda := KernelLift( beta, iota );
  
		return CokernelObject( lambda );
  
end );

InstallMethod( IsBlueExact,
		[ IsRecord, IsList ],

	function( A, S )
		local r;
	
		r := A.rank( S );
		
		return ForAll(
			[ 0 .. r - 2 ],
			i -> ForAll(
				[ 0 .. r - 2 - i ],
				j -> IsZero( OrlikSolomonBicomplexHorizontalHomologyObject( A, S, i, j ) ) 
				)
		);

end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyObject,
		[ IsRecord, IsList, IsInt, IsInt ],

	function( A, S, i, j )
		local alpha, beta, iota, lambda;

		alpha := OrlikSolomonBicomplexDifferential( A, S, i, j - 1, i, j );
		beta := OrlikSolomonBicomplexDifferential( A, S, i, j, i, j + 1 );
  
		if not IsZero( PreCompose( alpha, beta ) ) then
      
      		Error( "the composition of the given morphisms has to be zero" );
      
  		fi;
  
  		iota := ImageEmbedding( alpha );
  
  		lambda := KernelLift( beta, iota );
  
  	return CokernelObject( lambda );
  
end );

InstallMethod( IsRedExact,
		[ IsRecord, IsList ],

	function( A, S )
		local r;
	
		r := A.rank( S );
		
		return ForAll(
			[ 0 .. r - 2 ],
			i -> ForAll(
				[ 0 .. r - 2 - i ],
				j -> IsZero( OrlikSolomonBicomplexVerticalHomologyObject( A, S, i, j ) ) 
				)
		);

end );

InstallMethod( BlueNonExactness,
		[ IsRecord, IsList ],

	function( A, S )
		local r, res, i, j;
		
		r := A.rank( S );
	
		res := NullMat( r - 1, r - 1, 0 );
	
		for i in [ 0 .. r - 2 ] do
			for j in [ 0 .. r - 2 - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexHorizontalHomologyObject( A, S, i, j ) );
			od;
		od;

	return res;
end );

InstallMethod( RedNonExactness,
		[ IsRecord, IsList ],

	function( A, S )
		local r, res, i, j;
		
		r := A.rank( S );
	
		res := NullMat( r - 1, r - 1 );
	
		for i in [ 0 .. r - 2 ] do
			for j in [ 0 .. r - 2 - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexVerticalHomologyObject( A, S, i, j ) );
			od;
		od;

	return res;
end );

InstallMethod( Euler,
		[ IsList ],

	function( mat )
		local M, n, i, j;
		M := StructuralCopy( mat );
		n := Length( M );
		for i in [ 1 .. n ] do
			for j in [ 1 .. n ] do
				M[i][j] := ( -1 )^( i + j ) * M[i][j];
			od;
		od;
		return ( - 1)^( n ) * List( [ 0 .. n - 2 ], k -> Sum( Sum( M{ [ 1 .. k + 1 ] } { [ 1 .. n - 1 - k ] } ) ) );
end);

InstallMethod( OrlikSolomonBicomplex,
		[ IsRecord, IsList ],

	function( A, S ) 
		local i, j;
		
		if IsBound( A.( String( S ) ) ) then
			return A.( String( S ) );
		fi;

		A.( String( S ) ) := DoubleCochainComplex(
			A.cat,
			function( i, j ) return OrlikSolomonBicomplexDifferential( A, S, -i, j, -i - 1, j ); end,
			function( i, j ) return OrlikSolomonBicomplexDifferential( A, S, -i, j, -i, j + 1 ); end
		);

	return A.( String( S ) );
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
	
	res[n + 1][1] := 1;
	res[n + 1][n + 1] := -1;
	
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

CellIntegralBlueArrangement := function( n, w )
	local res, i, j, a, b, c;
	
	res := [];
	
	for i in [ 1 .. n - 2 ] do
		res[i] := [];
		for j in [ 1 .. n - 2 ] do
			res[i][j] := 0;
		od;
	od;
	
	res[1][1] := 1;
	
	c := 2;
	
	for i in [ 1 .. n ] do
		a := i^w;
		if i = n then
			b := 1^w;
		else
			b := ( i + 1 )^w;
		fi;
		
		if a in [ 1 .. n - 3 ] then
			if b in [ 1 .. n - 3 ] then
				res[c][a + 1] := 1;
				res[c][b + 1] := -1;
				c := c + 1;
			elif b = n - 2 then
				res[c][a + 1] := 1;
				res[c][1] := -1;
				c := c + 1;
			elif b = n then
				res[c][a + 1] := 1;
				c := c + 1;
			fi;
		elif a = n - 2 then
			if b in [ 1 .. n - 3 ] then
				res[c][1] := 1;
				res[c][b + 1] := -1;
				c := c + 1;
			fi;
		elif a = n then
			if b in [ 1 .. n - 3 ] then
				res[c][b + 1] := -1;
				c := c + 1;
			fi;
		fi;
	od;
	
	return res;
end;

# BrownMotive := function( n )
# 	local L, M, matroid, chi, cat, c, rk, res;
# 
# 	L := List( Combinations( [ 1.. n - 2 ], 2 ), function( c ) local res; res := List( [ 1 .. n - 2 ], i -> 0 ); res{ c } := [ 1, - 1 ]; return res; end );
# 	Append( L, List( [ 1 .. n - 2 ], function( c ) local res; res := List( [ 1 .. n - 2 ], i -> 0 ); res[c] := 1; return res; end ) );
# 	L := Set( L );
# 	M := SimplexArrangement( n - 3 );
# 	L := Difference( L, M );
# 
# 	matroid := Matroid( Concatenation( L, M ), HomalgFieldOfRationals( ) );
# 	rk := RankFunction( matroid );
# 	
# 	chi := function( flat )
# 		return not rk( flat ) = Length( Filtered( flat, i -> i > Length( L ) ) );
# 	end;
# 
# 	cat := MatrixCategory ( HomalgFieldOfRationals( ) );
# 
# #	Print( "L = "); PrintArray( L ); Print( "\n" );
# #	Print( "M = "); PrintArray( M ); Print( "\n" );
# 
# #	return Concatenation( L, M );
# 
# 	res := OrlikSolomonBicomplexRecord( matroid, chi, cat );
# 	Error( "" );
# 	return res;
# end;

BrownMotive := function( n )
	local L, M, matroid, chi, c, rk, A;

	L := List( Combinations( [ 1.. n - 2 ], 2 ), function( c ) local res; res := List( [ 1 .. n - 2 ], i -> 0 ); res{ c } := [ 1, - 1 ]; return res; end );
	Append( L, List( [ 1 .. n - 2 ], function( c ) local res; res := List( [ 1 .. n - 2 ], i -> 0 ); res[c] := 1; return res; end ) );
	L := Set( L );
	M := SimplexArrangement( n - 3 );
	L := Difference( L, M );

	matroid := Matroid( Concatenation( L, M ), HomalgFieldOfRationals( ) );
	rk := RankFunction( matroid );
	
	chi := function( flat )
		if flat = [ 1 .. Length( Flats( matroid )[2] ) ] then
			return fail;
		else
			return not rk( flat ) = Length( Filtered( flat, i -> i > Length( L ) ) );
		fi;
	end;
	
	A := OrlikSolomonBicomplexRecord( matroid, chi );	
	A.L := L;
	A.M := M;

	return A;
#	return OrlikSolomonBicomplexRecord( matroid, chi );
end;

CellIntegralRedBiOS := function( n, w )
		local m, rk, i, chi, cat;
		
		cat := MatrixCategory ( HomalgFieldOfRationals( ) );
	
		m := Concatenation( CellIntegralBlueArrangement( n, w ), SimplexArrangement( n - 3 ) );
		m := Matroid( m, HomalgFieldOfRationals( ) );
		rk := RankFunction( m );
	
		chi := function( flat )
			return not rk( flat ) = Length( Filtered( flat, i -> i > n - 2 ) );
		end;
	
	return OrlikSolomonBicomplexRecord( m, chi, cat );
	
end;

CellIntegralBlueBiOS := function( n, w )
		local m, rk, i, chi, cat;
		
		cat := MatrixCategory ( HomalgFieldOfRationals( ) );
	
		m := Concatenation( CellIntegralBlueArrangement( n, w ), SimplexArrangement( n - 3 ) );
		m := Matroid( m, HomalgFieldOfRationals( ) );
		rk := RankFunction( m );

		chi := function( flat )
			if rk( flat ) = n - 2 then
				return true;
			else
				return not rk( flat ) = Length( Filtered( flat, i -> i > n - 2 ) );
			fi;
		end;
	
		return OrlikSolomonBicomplexRecord( m, chi, cat );
	
end;

CellIntegralBiOS := function( n, w )
		local m, rk, i, chi, cat;
		
		cat := MatrixCategory ( HomalgFieldOfRationals( ) );
	
		m := Concatenation( CellIntegralBlueArrangement( n, w ), SimplexArrangement( n - 3 ) );
		m := Matroid( m, HomalgFieldOfRationals( ) );
		rk := RankFunction( m );

		chi := function( flat )
			if rk( flat ) = n - 2 then
				return fail;
			else
				return not rk( flat ) = Length( Filtered( flat, i -> i > n - 2 ) );
			fi;
		end;
	
		return OrlikSolomonBicomplexRecord( m, chi, cat );
	
end;

MultizetaBlueArrangement := function( ni_list )
	return IteratedIntegralBlueArrangement( Multizeta01Word ( ni_list ) );
end;

# RedMultizetaBiOS := function( ni_list, cat )
# 	local Q, n, m, rk, chi;
# 	
# 	Q := HomalgFieldOfRationals( );
# 	n := Sum( ni_list );
# 	
# 	m := Concatenation( MultizetaBlueArrangement( ni_list ), SimplexArrangement( n ) );
# 	m := Matroid( m, Q );
# 	rk := RankFunction( m );
# 	
# 	chi := function( flat )
# #		return not ForAll( flat, i -> i > n + 1 );
# 		return not rk( flat ) = Length( Filtered( flat, i -> i > n + 1 ) );
# 	end;
# 	
# 	return OrlikSolomonBicomplexRecord( m, chi, cat );
# end;

InstallMethod( RedMultizetaBiOS,
        [ IsList, IsCapCategory ],
        
	function( ni_list, cat )

		local Q, n, m, rk, chi;
	
		Q := HomalgFieldOfRationals( );
		n := Sum( ni_list );
	
		m := Concatenation( MultizetaBlueArrangement( ni_list ), SimplexArrangement( n ) );
		m := Matroid( m, Q );
		rk := RankFunction( m );
	
		chi := function( flat )
			return not rk( flat ) = Length( Filtered( flat, i -> i > n + 1 ) );
		end;
	
	return OrlikSolomonBicomplexRecord( m, chi, cat );
	
end );

InstallMethod( RedMultizetaBiOS,
        [ IsList ],
        
	function( ni_list )
	
		return RedMultizetaBiOS( ni_list, MatrixCategory( HomalgFieldOfRationals() ) );

end );

InstallMethod( BlueMultizetaBiOS,
		[ IsList, IsCapCategory ],

	function( ni_list, cat )

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
		
		return OrlikSolomonBicomplexRecord( m, chi, cat );

end );

InstallMethod( BlueMultizetaBiOS,
        [ IsList ],
        
	function( ni_list )
	
		return BlueMultizetaBiOS( ni_list, MatrixCategory( HomalgFieldOfRationals() ) );

end );

InstallMethod( MultizetaBiOS,
		[ IsList, IsCapCategory ],

	function( ni_list, cat )

		local Q, n, m, rk, chi;
		
		Q := HomalgFieldOfRationals( );
		n := Sum( ni_list );
		
		m := Concatenation( MultizetaBlueArrangement( ni_list ), SimplexArrangement( n ) );
		m := Matroid( m, Q );
		rk := RankFunction( m );
		
		chi := function( flat )
			if rk( flat ) = n + 1 then
				return fail;
			else
				return not rk( flat ) = Length( Filtered( flat, i -> i > n + 1 ) );
			fi;
		end;
		
		return OrlikSolomonBicomplexRecord( m, chi, cat );

end );

InstallMethod( MultizetaBiOS,
        [ IsList ],
        
	function( ni_list )
	
		return MultizetaBiOS( ni_list, MatrixCategory( HomalgFieldOfRationals() ) );

end );

InstallMethod( OrlikSolomonBicomplexDimensions,
		[ IsRecord, IsList ],

	function( A, Sigma )
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
end );

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
		od;
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




CellularArrangement := function( n, w ) 
	# choice of affine coordinates t_i with (z_1, ..., z_n) = (t_1, t_2, t_3, ..., t_{n-3}, 1, infty, 0)
	# choice of projective coordinates (t_0, t_1, ..., t_n) with {t_0 = 0} the hyperplane at infinity (for consistency with SimplexArrangement)
	local res, i, j, a, b, c;
	
	res := [];
	
	for i in [ 1 .. n - 2 ] do
		res[i] := [];
		for j in [ 1 .. n - 2 ] do
			res[i][j] := 0;
		od;
	od;
	
	c := 1;
	
	for i in [ 1 .. n ] do
		a := i^w;
		if i = n then
			b := 1^w;
		else
			b := ( i + 1 )^w;
		fi;
		
		if a in [ 1 .. n - 3 ] then
			if b in [ 1 .. n - 3 ] then
				res[c][a + 1] := 1;
				res[c][b + 1] := -1;
				c := c + 1;
			elif b = n - 2 then
				res[c][a + 1] := 1;
				res[c][1] := -1;
				c := c + 1;
			elif b = n then
				res[c][a + 1] := 1;
				c := c + 1;
			fi;
		elif a = n - 2 then
			if b in [ 1 .. n - 3 ] then
				res[c][1] := 1;
				res[c][b + 1] := -1;
				c := c + 1;
			fi;
		elif a = n then
			if b in [ 1 .. n - 3 ] then
				res[c][b + 1] := -1;
				c := c + 1;
			fi;
		fi;
	od;
	
	if c = n-2 then
		res[n-2][1] := 1; # adding the hyperplane at infinity
	fi;
	
	return res;
end;

CellularBiArrangementBlue := function( n, w )
		local l, m, rk, i, chi, cat;
		
		cat := MatrixCategory ( HomalgFieldOfRationals( ) );
	
		l := Concatenation( CellularArrangement( n, w ), SimplexArrangement( n - 3 ) );
		m := Matroid( l, HomalgFieldOfRationals( ) );
		rk := RankFunction( m );

		chi := function( flat )
			if rk( flat ) = n - 2 then # flat is the maximal stratum {0}
				return true;
			elif rk( flat ) = Length( Filtered( flat, i -> i <= n - 2 ) ) then # flat is an intersection of blue hyperplanes
				return true;
			elif rk( flat ) = Length( Filtered( flat, i -> i > n-2 ) ) then # flat is an intersection of red hyperplanes
				return false;
			else 
				return fail;
			fi;
		end;
	
		return OrlikSolomonBicomplexRecord( m, chi, cat ); # or replace with the relevant command
end;

CellularBiArrangementRed := function( n, w )
		local l, m, rk, i, chi, cat;
		
		cat := MatrixCategory ( HomalgFieldOfRationals( ) );
	
		l := Concatenation( CellularArrangement( n, w ), SimplexArrangement( n - 3 ) );
		m := Matroid( l, HomalgFieldOfRationals( ) );
		rk := RankFunction( m );

		chi := function( flat )
			if rk( flat ) = n - 2 then # flat is the maximal stratum {0}
				return false;
			elif rk( flat ) = Length( Filtered( flat, i -> i <= n - 2 ) ) then # flat is an intersection of blue hyperplanes
				return true;
			elif rk( flat ) = Length( Filtered( flat, i -> i > n-2 ) ) then # flat is an intersection of red hyperplanes
				return false;
			else 
				return fail;
			fi;
		end;
	
		return OrlikSolomonBicomplexRecord( m, chi, cat ); # or replace with the relevant command
end;

InstallMethod( CellularBiArrangement,
	[ IsInt, IsPerm, IsBool ],

	function( n, w, default )
		local l, m, rk, i, chi, cat;
		
		l := Concatenation( CellularArrangement( n, w ), SimplexArrangement( n - 3 ) );
		m := Matroid( l, HomalgFieldOfRationals( ) );
		rk := RankFunction( m );

		chi := function( flat )
			if rk( flat ) = n - 2 then # flat is the maximal stratum {0}
				return fail;
			elif rk( flat ) = Length( Filtered( flat, i -> i <= n - 2 ) ) then # flat is an intersection of blue hyperplanes
				return true;
			elif rk( flat ) = Length( Filtered( flat, i -> i > n-2 ) ) then # flat is an intersection of red hyperplanes
				return false;
			else 
				return default;
			fi;
		end;
	
		return OrlikSolomonBicomplexRecord( m, chi ); # or replace with the relevant command
end );

InstallMethod( CellularBiArrangement,
	[ IsInt, IsPerm ],

	function( n, w )
		return CellularBiArrangement( n, w, fail );

end );

#l:=List(Elements(DoubleCosets(S5,D5,D5)[4]), w -> ListPerm(w,5));

DihedralDoubleCoset := function( w )
	local n, Sym, Dih, v;
	n := Length( w );
	Sym := SymmetricGroup( n );
	Dih := AsSubgroup( Sym, DihedralGroup( IsPermGroup, 2 * n ) );
	w := PermList( w );
	return List( Elements( DoubleCoset( Dih, w, Dih ) ), v -> ListPerm( v, n ) );
end;

InstallMethod( CellMotive,
	[ IsList, IsBool ],

	function( pi, default )
		local w, n, A, M, i, j, k, motive;
		w := PermList( pi );
		n := Length( pi );
		
		A := CellularBiArrangement( n, w, default );;
		M := OrlikSolomonBicomplexDimensions( A, A.Smin );
		Print( "\n\n" );
		PrintArray( M );
		Print( "\n" );
	
		for i in [ 1 .. n - 1 ] do
			for j in [ 1 .. n - 1 ] do
				M[i][j] := ( -1 )^( i + j ) * M[i][j];
			od;
		od;

		motive := ( - 1)^( n + 1 ) * List( [ 0 .. n - 3 ], k -> Sum( Sum( M{ [ 1 .. k + 1 ] } { [ 1 .. n - 2 - k ] } ) ) );
	
		Print( motive, "\n" );
	
		return A;
end );

InstallMethod( CellMotive,
	[ IsList ],

	function( pi )
		return CellMotive( pi, fail );
		
end );

TestBlueRedStrataExactness := function( A )
	local i, S, chi;
	chi := A.coloring;
	for i in [ 1 .. Length( A.flats ) ] do
		for S in A.flats[i] do
			if chi( S ) = true then
            	Print( S, " blue exact: ", IsBlueExact( A, S ), "\n" );
			elif chi( S ) = false then
                Print( S, " red exact: ", IsRedExact( A, S ), "\n" );
            fi;
		od;
	od;
end;

AllBlackToBlue := function( chi )
	return function( S )
		if chi( S ) = fail then return true; fi;
		return chi( S );
	end;
end;

AllBlackToRed := function( chi )
	return function( S )
		if chi( S ) = fail then return false; fi;
		return chi( S );
	end;
end;

BlackToBlue := function( chi, Sigma )
	return function( S )
		if S <> Sigma then return chi( S ); fi;
		if chi( S ) = fail then return true; fi;
		return chi( S );
	end;
end;

BlackToRed := function( chi, Sigma )
	return function( S )
		if S <> Sigma then return chi( S ); fi;
		if chi( S ) = fail then return false; fi;
		return chi( S );
	end;
end;