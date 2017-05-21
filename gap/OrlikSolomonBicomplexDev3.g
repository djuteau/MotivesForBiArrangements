LoadPackage("M2");
LoadPackage("complex");
LoadPackage("alcove");
LoadPackage("MotivesForBiArrangements");

InstallOtherMethod( MakeInfList,
	[ IsDenseList, IsInt, IsObject ],
function( l, i, x)
	local tails;
	tails := RepeatListN( [ x ] );
	return Concatenate( tails, i, l, tails );
end );

InstallOtherMethod( CochainMorphism,
	[ IsCochainComplex, IsCochainComplex, IsList, IsInt ],
function( C, D, l, n )
	local zeroId, minC, minD, min, maxC, maxD, max, ll, middle;
	zeroId := IdentityMorphism( ZeroObject( CatOfComplex( C ) ) );
	minC := ActiveLowerBound( C ) + 1;
	minD := ActiveLowerBound( D ) + 1;
	min := Minimum( minC, minD );
	maxC := ActiveUpperBound( C ) - 1;
	maxD := ActiveUpperBound( D ) - 1;
	max := Maximum( maxC, maxD );
	middle := List(
		[ min, max ], 
		function(i)
			if i in [ n .. n + Length( l ) - 1 ] then
				return l[i - n + 1];
			else
				return ZeroMorphism( ObjectAt( C, i ), ObjectAt( D, i ) );
			fi;
		end
	);
	ll := MakeInfList( ll, min, zeroId );
	return CochainMorphism( C, D, ll );
end );

InstallOtherMethod( CochainComplex,
	[ IsCapCategory, IsList, IsInt ],
function( A, c, n )
	if c = [] then
		return ZeroObject( CochainComplexCategory ( A ) );
	fi;
	return CochainComplex( c, n );
end );

## StalkBicomplex( IsCapCategoryObject, IsInt, IsInt )

StalkDoubleChainComplex := function( M, i, j )
	return DoubleChainComplex( StalkChainComplex( StalkChainComplex( M, j ), i ) );
end;

StalkDoubleCochainComplex := function( M, i, j )
	return DoubleCochainComplex( StalkCochainComplex( StalkCochainComplex( M, j ), i ) );
end;

##################################################################################################################

BiOS := function( m, chi )
    local flats, Q, zero, zeroId, zeroComplex, one, oneId, oneComplex, InfListOfZeroMorphisms, Vect, A, k, i, j, l, SIGMA, SS, S, T, rows, cols, rows_j, cols_i, summands, s, phi, psi;

    flats := MakeInfList( Flats( m ), 0, [] );

	Q := HomalgFieldOfRationals();
    zero := 0 * Q;
 	zeroId := IdentityMorphism( zero );       
    zeroComplex := CochainComplex( [ zeroId ] );
    one := 1 * Q;
    oneId := IdentityMorphism( one );
    oneComplex := CochainComplex( [ oneId ] );
    InfListOfZeroMorphisms := MakeInfList( [], 0, zeroId );
    Vect := CapCategory( zero );

    A := rec(
             (String( [ ] )) := StalkDoubleCochainComplex( 1 * Q, 0, 0 ),
             );

    for k in [ 1 .. RankOfMatroid( m ) ] do
		for SIGMA in flats[k] do
			SS := MakeInfList( List( Sublist( flats, 0, k - 1 ), l -> Filtered( l, S -> IsSubset( SIGMA, S ) ) ), 0, [] );
			
            rows := MakeInfList(
            	List(
            		[ 0 .. k - 1 ],
            		j -> MakeInfList(
            			List(
            				[ j - k .. 0 ],
            				i -> MorphismBetweenDirectSums(
            					DirectSum( Vect, List( SS[i + j], S -> ObjectAt( A.( String( S ) ), i, j ) ) ),
            					List(
            						SS[i + j],
            						S -> List(
            							SS[i + j + 1],
            							function(T)
            								if IsSubset( S, T ) then
            									return A.( JoinStringsWithSeparator( [ S, T, i, j ] ) );
            								fi;
            								return ZeroMorphism( ObjectAt( A.( String( S ) ), i, j ), ObjectAt( A.( String( T ) ), i + 1, j ) );
            							end
            						)
            					),
            					DirectSum( Vect, List( SS[i + j + 1], T -> ObjectAt( A.( String( T ) ), i + 1, j ) ) )
            				)
            			),
            			j - k,
            			zeroId
            		)
            	),
            	0,
            	InfListOfZeroMorphisms
            );
            
			cols := MakeInfList(
            	List(
            		[ -k + 1 .. 0 ],
            		i -> MakeInfList(
            			List(
	          	  			[ -1 .. k + i - 2 ],
	            			j -> MorphismBetweenDirectSums(
    	        				DirectSum( Vect, List( SS[i + j], T -> ObjectAt( A.( String( T ) ), i, j ) ) ),
    	        				List(
    	        					SS[i + j], 
    	        					T -> List(
    	        						SS[i + j + 1],  
	    	        					function(S)
	    	        						if IsSubset( S, T ) then 
	    	        							return A.( JoinStringsWithSeparator( [ T, S, i, j ] ) );
	    	    							fi;
	    	    							return ZeroMorphism( ObjectAt( A.( String( T ) ), i, j ), ObjectAt( A.( String( S ) ), i, j + 1 ) );
	            						end
		        					)
	        					),
    	        				DirectSum( Vect, List( SS[i + j + 1], S -> ObjectAt( A.( String( S ) ), i, j + 1) ) )
    	        			)
            			),
            			-1,
            			zeroId
            		)
            	),
            	-k + 1,
            	InfListOfZeroMorphisms
            );
            
            for j in [ 0 .. k - 1 ] do
            	i := j - k + 1;
				summands := List( SS[k - 1], S -> ObjectAt( A.( String( S ) ), i, j ) );
				phi := KernelEmbedding( rows[j][i] );
				psi := CokernelProjection( cols[i][j - 1] );
            	for s in [ 1 .. Length( SS[k - 1] ) ] do
            		S := SS[k - 1][s];

	            	A.( JoinStringsWithSeparator( [ SIGMA, S, i - 1, j ] ) ) := PreCompose(
    	        		phi,
    	        		ProjectionInFactorOfDirectSum( summands, s )
    	        	);
					rows_j := Replace( rows[j], i - 1, [ UniversalMorphismFromZeroObject( Source( phi ) ), phi ] );
					rows := Replace( rows, j, [ rows_j ] );
					
        	    	A.( JoinStringsWithSeparator( [ S, SIGMA, i, j ] ) ) := PreCompose(
            			InjectionOfCofactorOfDirectSum( summands, s ),
            			psi
            		);
            		cols_i := Replace( cols[i], j, [ psi, UniversalMorphismIntoZeroObject( Range( psi ) ) ] );
            		cols := Replace( cols, i, [ cols_i ] );
            	od;
            od;
		od;
		A.( String( SIGMA ) ) := DoubleCochainComplex( Vect, rows, cols );
	od;
	return A;
end;



