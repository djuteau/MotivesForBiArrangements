## StalkBicomplex( IsCapCategoryObject, IsInt, IsInt )

StalkDoubleChainComplex := function( M, i, j )
	return StalkChainComplex( StalkChainComplex( M , j ), i );
end;

StalkDoubleChainComplex := function( M, i, j )
	return DoubleChainComplex( StalkChainComplex( StalkChainComplex( M, j ), i ) );
end;


##### CHANGE i TO -i !!!!

## C is a bicomplex of objects of A, concentrated in 0 <= i, j and i + j < k
## HorizontalKill returns a bicomplex which coincides with A up to the diagonal
## i + j = k - 1, with a new diagonal i + j = m mapping to the diagonal i + j = k - 1
## by the kernels of the horizontal differentials from that diagonal.
##
## The new vertical differential differentials, from i + j = k - 1 to i + j = k,
## are constructed by the universal property of the kernel (lift of the diagonal morphism
## = composition of horizontal and vertical differentials)
##
## Possible improvement: Use and update bounds

HorizontalKillDiagonalCohomology := function( A, C, k )
	local ZeroEndoOfZeroObj, H, V, phi, res;
	ZeroEndoOfZeroObj := IdentityMorphism( ZeroObject( A ) );
	phi := List( [ 0 .. k - 1 ], j -> KernelEmbedding( HorizontalDifferentialAt( C, k - 1 - j, j ) ) );
	H := function( i, j )
		if i < 0 or j < 0 or i + j > k + 1 then 
			return ZeroEndoOfZeroObj;
		elif i = 0 then 
			return UniversalMorphismIntoZeroObject( ObjectAt( C, 0, j) );
		elif i + j = k + 1 then 
			return UniversalMorphismFromZeroObject( Source ( phi[j] ) );
		elif i + j < k then
			return HorizontalDifferentialAt( C, i, j);
		else return phi[j];
		fi;
	end;
	V := function( i, j )
		if i < 0 or j < - 1 or i + j > k then
			return ZeroEndoOfZeroObj;
		elif j = - 1 then
			return UniversalMorphismFromZeroObject( ObjectAt( C, i, 0 ) );
		elif i + j = k then
			return UniversalMorphismIntoZeroObject( ObjectAt( C, i, j ) );
		elif i + j < k - 1 then
			return VerticalDifferentialAt( C, i, j - 1 );
		else return LiftAlongMonomorphism(
						phi[j],
						PreCompose( HorizontalDifferentialAt( C, i, j - 1 ), VerticalDifferentialAt( C, i - 1, j - 1 ) )
						);
		fi;
	end;
	res := DoubleChainComplex( A, H, V );
	SetAboveBound( res, k );
	SetBelowBound( res, 0 );
	SetLeftBound( res, 0 );
	SetRightBound( res, k );
	return res;
end;

KillDiagonalCohomology := function( A, C, k, color )
	if color then
		return HorizontalKillDiagonalCohomology( A, C, k );
	else return VerticalKillDiagonalCohomology( A, C, k );
end;

## A := CategoryOfHomalgLeftModules(Q);



BiOS := function( m, chi )
    local flats, Q, zero, F, idF, A, k, i, Sigma, SS, TT, phi, d, D, s, t, psi;
    
    flats := Flats( m );
    
    Q := HomalgRing( m );
    
    zero := 0 * Q;
    
    F := 1 * Q;
    idF := IdentityMorphism( F );
    
    Qev := HomalgCategory(zero);
    
    A := rec(
             (String( [ ] )) := StalkBicomplex( F, 0, 0 ),
             );
    
    for k in [ 1 .. RankOfMatroid( m ) ] do
        for Sigma in flats[ k + 1 ] do
            SS := List( flats{ [ O .. k ] }, l -> Filtered( l, S -> IsSubset( Sigma, S ) ) );
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
            
            
