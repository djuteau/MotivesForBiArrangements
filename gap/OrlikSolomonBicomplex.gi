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

####################################
#
# methods for constructors:
#
####################################

##
InstallMethod( OrlikSolomonBicomplex,
        [ IsMatroid, IsFunction ],
        
  function( m, chi )
    local Q, zero, zeroId, one, oneId, A, k, i, Sigma, SS, TT, S, T, phi, d, D, s, t, psi, obj, cat;
    
    Q := HomalgFieldOfRationals();
    cat := LeftPresentations( Q );
    
    zero := ZeroObject( cat );
    zeroId := IdentityMorphism( zero );
    one := TensorUnit( cat );
    oneId := IdentityMorphism( one );


    A := rec(
             (JoinStringsWithSeparator( [ [ ], 0, 0 ] )) := one,
             );

    for Sigma in FlatsOfRank( m, 1 ) do
        if chi( Sigma ) then
            A.(JoinStringsWithSeparator( [ Sigma, 1, 0 ] )) := one;
            A.(JoinStringsWithSeparator( [ Sigma, 0, 1 ] )) := zero;
            A.(JoinStringsWithSeparator( [ [ ], 1, -1, Sigma, 1, 0 ] )) := ZeroMorphism( zero, one );
            A.(JoinStringsWithSeparator( [ Sigma, 1, 0, [ ], 0, 0 ] )) := oneId;
            A.(JoinStringsWithSeparator( [ [ ], 0, 0, Sigma, 0, 1 ] )) := ZeroMorphism( one, zero );
            A.(JoinStringsWithSeparator( [ Sigma, 0, 1, [ ], -1, 1 ] )) := zeroId;
        else
            A.(JoinStringsWithSeparator( [ Sigma, 1, 0 ] )) := zero;
            A.(JoinStringsWithSeparator( [ Sigma, 0, 1 ] )) := one;
            A.(JoinStringsWithSeparator( [ [ ], 1, -1, Sigma, 1, 0 ] )) := zeroId;
            A.(JoinStringsWithSeparator( [ Sigma, 1, 0, [ ], 0, 0 ] )) := ZeroMorphism( zero, one );
            A.(JoinStringsWithSeparator( [ [ ], 0, 0, Sigma, 0, 1 ] )) := oneId;
            A.(JoinStringsWithSeparator( [ Sigma, 0, 1, [ ], -1, 1 ] )) := ZeroMorphism( one, zero );
        fi;
    od;

    
    for k in [ 2 .. RankOfMatroid( m ) ] do
        for Sigma in FlatsOfRank( m, k ) do
            SS := Filtered( FlatsOfRank( m, k - 1 ), S -> IsSubset( Sigma, S ) );
            TT := Filtered( FlatsOfRank( m, k - 2 ), T -> IsSubset( Sigma, T ) );
            
            for i in [ 1 .. k ] do
                phi := List( SS, S -> List( TT,
                               function( T )
                               	 if not IsBound( A.(JoinStringsWithSeparator( [ S, i - 1, k - i ] )) ) then
                               	 	A.(JoinStringsWithSeparator( [ S, i - 1, k - i ] )) := zero;
                               	 fi;
                               	 if not IsBound( A.(JoinStringsWithSeparator( [ T, i - 2, k - i ] )) ) then
                               	 	A.(JoinStringsWithSeparator( [ T, i - 2, k - i ] )) := zero;
                               	 fi;
                                 if not IsBound( A.(JoinStringsWithSeparator( [ S, i - 1, k - i, T, i - 2, k - i ] )) ) then
                                     A.(JoinStringsWithSeparator( [ S, i - 1, k - i, T, i - 2, k - i ] )) :=
                                       ZeroMorphism( A.(JoinStringsWithSeparator( [ S, i - 1, k - i ] )),  A.(JoinStringsWithSeparator( [ T, i - 2, k - i ] )) );
                                 fi;
                                 return A.(JoinStringsWithSeparator( [ S, i - 1, k - i, T, i - 2, k - i ] ));
                             end ) );
                phi := MorphismBetweenDirectSums( phi );
                A.(JoinStringsWithSeparator( [ Sigma, i - 1, k - i, i - 2, k - i ] )) := phi;
            od;
            
            for i in [ 0 .. k - 1 ] do
                phi := List( TT, T -> List( SS,
                               function( S )
                               	 if not IsBound( A.(JoinStringsWithSeparator( [ T, i, k - i - 2 ] )) ) then
                               	 	A.(JoinStringsWithSeparator( [ T, i, k - i - 2 ] )) := zero;
                               	 fi;
                               	 if not IsBound( A.(JoinStringsWithSeparator( [ S, i, k - i - 1 ] )) ) then
                               	 	A.(JoinStringsWithSeparator( [ S, i, k - i - 1 ] )) := zero;
                               	 fi; 
                               	 
                                 if not IsBound( A.(JoinStringsWithSeparator( [ T, i, k - i - 2, S, i, k - i - 1 ] )) ) then
                                     A.(JoinStringsWithSeparator( [ T, i, k - i - 2, S, i, k - i - 1 ] )) :=
                                       ZeroMorphism( A.(JoinStringsWithSeparator( [ T, i, k - i - 2 ] )), A.(JoinStringsWithSeparator( [ S, i, k - i - 1 ] )) );
                                 fi;
                                 return A.(JoinStringsWithSeparator( [ T, i, k - i - 2, S, i, k - i - 1 ] ));
                             end ) );
                
                phi := MorphismBetweenDirectSums( phi );
                A.(JoinStringsWithSeparator( [ Sigma, i, k - i - 2, i, k - i - 1 ] )) := phi;
            od;
            
            if chi( Sigma ) then
                for i in Reversed( [ 1 .. k ] ) do
                    phi := A.(JoinStringsWithSeparator( [ Sigma, i - 1, k - i, i - 2, k - i ] ));
                    d := KernelEmbedding( phi );
                    A.(JoinStringsWithSeparator( [ Sigma, i, k - i ] )) := Source( d );
                    D := List( [ 1 .. Length( SS ) ], s -> A.(JoinStringsWithSeparator( [ SS[s], i - 1, k - i ] )) );
                    for s in [ 1 .. Length( SS ) ] do
                        A.(JoinStringsWithSeparator( [ Sigma, i, k - i, SS[s], i - 1, k - i ] )) := PreCompose( d, ProjectionInFactorOfDirectSum( D, s ) );
                    od;
                    if i < k then
                        D := List( [ 1 .. Length( SS ) ], s -> A.(JoinStringsWithSeparator( [ SS[s], i, k - i - 1 ] )) );
                        psi := PreCompose(
                                       A.(JoinStringsWithSeparator( [ Sigma, i, k - i - 1, i - 1, k - i - 1 ] )),
                                       A.(JoinStringsWithSeparator( [ Sigma, i - 1, k - i - 1, i - 1, k - i ] )) );
                        
                        psi := KernelLift( phi, psi );
                        for s in [ 1 .. Length( SS ) ] do
                            A.(JoinStringsWithSeparator( [ SS[s], i, k - i - 1, Sigma, i, k - i ] )) :=
                              PreCompose( InjectionOfCofactorOfDirectSum( D, s ), psi );
                        od;
                    fi;
                od;
                A.(JoinStringsWithSeparator( [ Sigma, 0, k ] )) := zero;
                obj := A.(JoinStringsWithSeparator( [ Sigma, k, 0 ] )); ###
                for S in SS do
                	A.(JoinStringsWithSeparator( [ Sigma, 0, k, S, -1, k ] )) := zeroId;
                	phi := UniversalMorphismFromZeroObject( obj );
                	A.(JoinStringsWithSeparator( [ S, k, -1, Sigma, k, 0 ] )) := phi;
                od;
            else
                for i in [ 1 .. k ] do
                    phi := A.(JoinStringsWithSeparator( [ Sigma, i, k - i - 2, i, k - i - 1 ] ));
                    d := CokernelProjection( phi );
                    A.(JoinStringsWithSeparator( [ Sigma, i, k - i ] )) := Range( d );
                    D := List( [ 1 .. Length( SS ) ], s -> A.(JoinStringsWithSeparator( [ SS[s], i, k - i - 1 ] )) );
                    for s in [ 1 .. Length( SS ) ] do
                        A.(JoinStringsWithSeparator( [ SS[s], i, k - i - 1, Sigma, i, k - i ] )) := PreCompose( InjectionOfCofactorOfDirectSum( D, s ), d );
                    od;
                    if i > 1 then
                        D := List( [ 1 .. Length( SS ) ], s -> A.(JoinStringsWithSeparator( [ SS[s], i - 1, k - i ] )) );
                        psi := PreCompose(
                                       A.(JoinStringsWithSeparator( [ Sigma, i, k - i - 1, i - 1, k - i - 1 ] )),
                                       A.(JoinStringsWithSeparator( [ Sigma, i - 1, k - i - 1, i - 1, k - i ] )) );
                        
                        psi := CokernelColift( phi, psi );
                        for s in [ 1 .. Length( SS ) ] do
                            A.(JoinStringsWithSeparator( [ Sigma, i, k - i, SS[s], i - 1, k - i ] )) :=
                              PreCompose( psi, ProjectionInFactorOfDirectSum( D, s ) );
                        od;
                    fi;
                od;
                A.(JoinStringsWithSeparator( [ Sigma, k, 0 ] )) := zero;
                obj := A.(JoinStringsWithSeparator( [ Sigma, 0, k ] )); ###
#                obj := A.(JoinStringsWithSeparator( [ Sigma, 0, k ] )); ###
                for S in SS do
                	A.(JoinStringsWithSeparator( [ S, k, -1, Sigma, k, 0 ] )) := zeroId;
                	phi := UniversalMorphismIntoZeroObject( obj );
                	A.(JoinStringsWithSeparator( [ Sigma, 0, k, S, -1, k ] )) := phi;
                od;

            fi;
        od;
    od;
    
    return A;
    
end );
