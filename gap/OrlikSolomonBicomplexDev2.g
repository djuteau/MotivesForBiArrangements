
LoadPackage("M2");
LoadPackage("complex");
LoadPackage("alcove");
LoadPackage("MotivesForBiArrangements");

## StalkBicomplex( IsCapCategoryObject, IsInt, IsInt )

StalkDoubleChainComplex := function( M, i, j )
	return DoubleChainComplex( StalkChainComplex( StalkChainComplex( M, j ), i ) );
end;

StalkDoubleCochainComplex := function( M, i, j )
	return DoubleCochainComplex( StalkCochainComplex( StalkCochainComplex( M, j ), i ) );
end;


##### CHANGE i TO -i !!!!

## C is a bicomplex, concentrated in 0 <= - i, j and - i + j < k
## HorizontalKillDiagonalCohomology(C) returns a bicomplex which coincides with A up to the diagonal
## - i + j = k - 1, with a new diagonal - i + j = k mapping to the diagonal - i + j = k - 1
## by the kernels of the horizontal differentials from that diagonal.
##
## The new vertical differential differentials, from - i + j = k - 1 to - i + j = k,
## are constructed by the universal property of the kernel (lift of the diagonal morphism
## = composition of horizontal and vertical differentials)

HorizontalKillDiagonalCohomology := function( A, k )
	local cat, ZeroEndoOfZeroObj, rows, cols, i, j, rows_j, cols_i, phi, res;
	cat := CatOfDoubleComplex( A );
	ZeroEndoOfZeroObj := IdentityMorphism( ZeroObject( cat ) );
	rows := Rows( A );
	cols := Columns( A );
	for j in [ 0 .. k - 1 ] do
		i := j - k;
	    phi := KernelEmbedding( rows[j][i + 1] );	
		rows_j := Replace( rows[j], i - 1, [ UniversalMorphismFromZeroObject( Source( phi ) ), phi ] );
		rows := Replace( rows, j, [ rows_j ] );
	od;
	for i in [ -k .. 0 ] do
		j := k + i - 1;
		phi := KernelLift( rows[j + 1][i + 1], PreCompose( HorizontalDifferentialAt( A, i, j ), VerticalDifferentialAt( A, i + 1, j ) ) );
		cols_i := Replace( cols[i], j, [ phi, UniversalMorphismIntoZeroObject( Range( phi ) ) ] );
		cols := Replace( cols, i, [ cols_i ] );
	od;
	res := DoubleChainComplex( cat, rows, cols );
	SetAboveBound( res, k );
	SetBelowBound( res, 0 );
	SetLeftBound( res, - k );
	SetRightBound( res, 0 );
	return res;
end;

VerticalKillDiagonalCohomology := function( A, k )
	local cat, ZeroEndoOfZeroObj, rows, cols, i, j, rows_j, cols_i, phi, res;
	cat := CatOfDoubleComplex( A );
	ZeroEndoOfZeroObj := IdentityMorphism( ZeroObject( cat ) );
	rows := Rows( A );
	cols := Columns( A );
	for i in [ -k .. 0 ] do
		j := k + i - 1;
	    phi := CokernelProjection( cols[i][j - 1] );	
		cols_i := Replace( cols[i], j, [ phi, UniversalMorphismIntoZeroObject( Range( phi ) ) ] );
		cols := Replace( cols, j, [ cols_i ] );
	od;
	for j in [ 1 .. k ] do
		i := j - k;
		phi := CokernelColift( cols[i][j - 2], PreCompose( HorizontalDifferentialAt( A, i, j - 1), VerticalDifferentialAt( A, i + 1, j - 1) ) );
		rows_j := Replace( rows[j], i - 1, [ UniversalMorphismFromZeroObject( Source( phi ) ), phi ] );
#		Error();
		rows := Replace( rows, j, [ rows_j ] );
	od;
	res := DoubleChainComplex( cat, rows, cols );
	SetAboveBound( res, k );
	SetBelowBound( res, 0 );
	SetLeftBound( res, - k );
	SetRightBound( res, 0 );
	return res;
end;
	

KillDiagonalCohomology := function( A, k, color )
	if color then
		return HorizontalKillDiagonalCohomology( A, k );
	else
		return VerticalKillDiagonalCohomology( A, k );
	fi;
end;

## cat := CategoryOfHomalgLeftModules(Q);

Q := HomalgFieldOfRationals();
m := Matroid([[1,0,0],[0,1,0],[0,0,1],[1,1,0]], Q);
flats := Flats(m);

# BlueArrangement := function( m )
# 	return [ m, x -> true ];
# end;

InstallOtherMethod( MakeInfList,
	[ IsDenseList, IsInt, IsObject ],
function( l, i, x)
	local tails;
	tails := RepeatListN( [ x ] );
	return Concatenate( tails, i, l, tails );
end );


InstallOtherMethod( CochainMorphism,
	[ IsCochainComplex, IsCochainComplex, IsList, IsInt ],
function( C, D, l , n )
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
	return ChainMorphism( C, D, ll );
end );


# DeclareMethod( "MorphismBetweenDirectSums", [ IsList, IsCapCategory ] );
# 
# InstallMethod( MorphismBetweenDirectSums,
# 	[ IsList, IsInt, IsInt, IsCapCategory ],
# function( l, m, n, cat )
# 	if m > 0 and n > 0 then
# 		return MorphismBetweenDirectSums( l );
# 	elif m > 0 and n = 0 then
# 		return UniversalMorphismIntoZeroObject( 
# 	if n = 0 then
# 	
# 	fi;
# end );

BiOS := function( m, chi )
    local flats, zero, zeroID, zeroComplex, A, k, i, j, l, SIGMA, SS, S, T, cols, summands, s;

    flats := MakeInfList( Flats( m ), 0, [] );
        
    zero := 0 * HomalgFieldOfRationals();
 	zeroId := IdentityMorphism( zero );       
    zeroComplex := CochainComplex( [ zeroId ] );

    A := rec(
             (String( [ ] )) := StalkDoubleCochainComplex( 1 * HomalgFieldOfRationals(), 0, 0 ),
             );

    for k in [ 1 .. RankOfMatroid( m ) ] do
    	for SIGMA in flats[k] do
            SS := MakeInfList( List( Sublist( flats, 0, k ), l -> Filtered( l, S -> IsSubset( SIGMA, S ) ) ), 0, [] );
            cols := MakeInfList(
            	List(
            		[ -k + 1 .. 0 ],
            		i -> CochainComplex(
	            		List(
	            			[ 0 .. k + i - 2 ],
	            			j -> MorphismBetweenDirectSums(
    	        				DirectSum( List( SS[i + j], T -> ObjectAt( A.( String( T ) ), i, j ) ) ),
    	        				List(
    	        					List( SS[i + j], T ->
    	        						List(
    	        							SS[i + j + 1], S -> 
	    	        						function(S)
	    	        							if IsSubset( S, T ) then 
	    	        								return A.( JoinStringsWithSeparator( [ String( T ), String(S), i, j ] ) );
	    	        							fi;
	    	        							return ZeroMorphism( ObjectAt( A.( String( T ) ), i, j ), ObjectAt( A.( String( S ) ), i, j + 1) );
	    	        						end
	    	        					)
    	        					)
    	        				),
    	        				DirectSum( List( SS[i + j + 1], S -> ObjectAt( A.( String( S ) ), i, j + 1) ) )
    	        			)
        	 		   	),
        	 		   	0
            		),
            		-k + 1,
            		zeroComplex
            	)
            );
            A.( String(SIGMA) ) := KillDiagonalCohomology(
            	DoubleCochainComplex(
            		List(
            			[ -k .. -1 ],
            			j -> CochainMorphism(
            				cols[i],
            				cols[i + 1],
            				List(
            					[ 0 .. k + i - 2 ],
	            				j -> MorphismBetweenDirectSums(
    	        					DirectSum( List( SS[i + j + 1], S -> ObjectAt( A.( String( S ) ), i, j ) ) ),
    	        					List(
    	        						List( SS[i + j + 1], S ->
    	        							List(
    	        								SS[i + j], T -> 
	    	        							function(S)
	    	        								if IsSubset( S, T ) then 
	    	        									return A.( JoinStringsWithSeparator( [ String( S ), String(T), i, j ] ) );
	    	        								fi;
	    	        								return ZeroMorphism( ObjectAt( A.( String( S ) ), i, j ), ObjectAt( A.( String( T ) ), i + 1, j ) );
	    	        							end
	    	        						)
    	        						)
    	        					),
    	        					DirectSum( List( SS[i + j], T -> ObjectAt( A.( String( T ) ), i + 1, j ) ) ) 	        					
    	        				)
            				),
            				0
            			)
            		)
            	),
            	k,
            	chi( SIGMA )
            );
            for j in [ 0 .. k ] do
            	i := j - k + 1;
				summands := List( SS[k - 1], S -> ObjectAt( A.( String( S ) ), i, j ) );
            	for s in [ 1 .. Length( SS[k - 1] ) ] do
            		S := SS[k - 1][s];
	            	A.( JoinStringsWithSeparator( [ SIGMA, S, i - 1, j ] ) ) :=
    	        		PreCompose(
    	        			HorizontalDifferentialAt( A.( String( SIGMA ) ), i - 1, j),
    	        			ProjectionInFactorOfDirectSum( summands, s)
    	        		);
        	    	A.( JoinStringsWithSeparator( [ S, SIGMA, i, j ] ) ) :=
            			PreCompose(
            				InjectionOfCofactorOfDirectSum( summands, s ),
            				VerticalDifferentialAt( A.( String( SIGMA ) ), i, j )
            			);
            	od;
            od;
    	od;
    od;
	return A;
end;











    for k in [ 2 .. RankOfMatroid( m ) ] do
        for Sigma in flats[ k ] do
            SS := MakeInfList( List( Sublist( flats, 0, k ), l -> Filtered( l, S -> IsSubset( Sigma, S ) ) ), 0, [] );
            cols := List( [ -k + 1 .. 0 ], i -> CochainComplex( List( [ 0 .. i + k - 1 ], function(j)
            		if j - i < 0 then
            			return UniversalMorphismFromZeroObject( DirectSum( List( SS[j - i + 1], S -> ObjectAt( A.( String( S ) ), i, j + 1 ) ) ) );
            		fi;
            		return MorphismBetweenDirectSums(
            		List(
            			SS[j - i],
            				T -> List(
            					SS[j - i + 1],
            					 	function( S ) 
            					 	   if not IsList( S ) then
            					 	   		Error( "not a list\n" );
            					 	   fi;
            							if IsSubset( S, T ) then
            								return VerticalDifferentialAt( A.( JoinStringsWithSeparator( [ T, S, i, j ] ) ) );
            							else
            								return ZeroMorphism( ObjectAt( A.( String( T ) ), i, j ), ObjectAt( A.( String( S ) ), i, j + 1 ) );
            							fi;
            						end )));
            			end ), 0 ) );
		od;
	od;
	
	return cols;
end;





            H := function( i , j)
            	if i < 0 or j < 0 or i + j > k then 
            		return ZeroEndoOfZeroObj;
            	elif
            	else return MorphismBetweenDirectSums(
            		List(
            			SS,
            			S - > List(
            						TT,
            					 	function( T ) 
            							if IsSubset( S, T ) then
            								return HorizontalDifferentialAt( A.( JoinStringsWithSeparator( [ S, T, i, j ] ) ) );
            							else
            								return ZeroMorphism( ObjectAt( A.( String( S ), i, j ) ), A.( String( T ), i - 1, j ) );
            							fi;
            						end;
            					)
            			)
            		)
            end;
            V := function( i, j)
            
            end;
			A.( String( Sigma ) ) := Bicomplex( Qev, H, V )
		od;
		
	od;            
            
            
