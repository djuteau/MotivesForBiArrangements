# LoadPackage( "alcove");
# LoadPackage( "LinearAlgebraForCAP");
# LoadPackage( "complex" );
# LoadPackage( "Bicomplex" );

# cat := MatrixCategory( HomalgFieldOfRationals( ) );

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

InstallMethod( Coloring,

    [ IsMatroid, IsInt, IsBool, IsBool ],
    
    function( m, k, default, last )
    
    local flats, i, j;
    
    flats := Flats( m );
    
    return Concatenation(
        List(
            [ 1 .. Rank( m ) - 1 ],
            i -> List(
                flats[i + 1],
                function( S )
                    if Length( Filtered( S, j -> j <= k ) ) = i then
                        return true;
                    elif Length( Filtered( S, j -> j > k ) ) = i then
                        return false;
                    else
                        return default;
                    fi;
                end 
                )
        ),
        [ [ last ] ]
        );
        
end );

InstallMethod( Coloring,

    [ IsMatroid, IsInt, IsBool, IsBool, IsBool ],
    
    function( m, k, conflict, default, last )
    
    local flats, i, j;
    
    flats := Flats( m );
    
    return Concatenation(
        List(
            [ 1 .. Rank( m ) - 1 ],
            i -> List(
                flats[i + 1],
                function( S )
                    if Length( Filtered( S, j -> j <= k ) ) = i and Length( Filtered( S, j -> j > k ) ) = i then
                        return conflict;
                    elif Length( Filtered( S, j -> j <= k ) ) = i then
                        return true;
                    elif Length( Filtered( S, j -> j > k ) ) = i then
                        return false;
                    else
                        return default;
                    fi;
                end 
                )
        ),
        [ [ last ] ]
        );
        
end );

# InstallMethod( Coloring,
# 
#     [ IsMatroid, IsInt, IsBool, IsBool, IsBool ],
#     
#     function( m, k, conflict, default, last )
#     
#     local flats, blue, red, i, j;
#     
#     flats := Flats( m );
#     blue := [ 1 .. k ];
#     red := [ k + 1 .. Length( flats[2] ) ];
#     
#     return Concatenation(
#         List(
#             [ 1 .. Rank( m ) - 1 ],
#             i -> List(
#                 flats[i + 1],
#                 function( S )
#                     if S = cl( Intersection( S, red ) ) then
#                         return false;
#                     elif S = cl( Intersection( S, blue ) ) then
#                         return true;
#                     else
#                         return default;
#                     fi;
#                 end 
#                 )
#         ),
#         [ [ last ] ]
#         );
#         
# end );
# 
# S = cl( Intersection( S, blue ) )  
# S = cl( Intersection( S, red ) )  

InstallMethod( ColoringFunction,
	[ IsMatroid, IsInt, IsBool, IsBool ],
	
	function( m, k, default, last )
	 	local rk;
 	
	 	rk := RankFunction( m );
 	
	 	return function( flat )
				if rk( flat ) = Rank( m ) then # flat is the maximal stratum {0}
					return last;
				elif rk( flat ) = Length( Filtered( flat, i -> i <= k ) ) then # flat is an intersection of blue hyperplanes
					return true;
				elif rk( flat ) = Length( Filtered( flat, i -> i > k ) ) then # flat is an intersection of red hyperplanes
					return false;
				else 
					return default;
				fi;
			end;
 	
 end );


InstallMethod( ArrangementFromGraph,

	[ IsInt, IsList ],

	function( n, L )
	
	local res, m, i, j;
	
	res := [ ];
	
	m := Length( L );
	
	for i in [ 1 .. m ] do
		res[i] := [ ];
		for j in [ 1 .. n + 1 ] do
			res[i][j] := 0;
		od;
		if L[i][1] > 1 then
			res[i][L[i][1] - 1] := -1;
		fi;
		if L[i][2] > 1 then
			res[i][L[i][2] - 1] := 1;
		fi;
	od;
	
	return res;

end );

InstallMethod( ArrangementFromGraph,

	[ IsList ],

	function( L )
	
	return ArrangementFromGraph( MaximumList( Union( L ) ) - 2, L );

end );

InstallMethod( ProjectiveSpaceOrlikSolomonBicomplexRecord,

    [ IsInt, IsList, IsList ],

	function ( n, L, M )
	
	return OrlikSolomonBicomplexRecord(
		ArrangementFromGraph( n, L ),
		ArrangementFromGraph ( n, M ),
		fail,
		true
	);
	
end );

InstallMethod( ProjectiveSpaceOrlikSolomonBicomplexRecord,

    [ IsList, IsList ],

	function ( L, M )
	
	return OrlikSolomonBicomplexRecord(
		ArrangementFromGraph( L ),
		ArrangementFromGraph ( M ),
		fail,
		false
	);
	
end );

InstallMethod( OrlikSolomonBicomplexRecord,

	[ IsList, IsList ],
 	
 	function( L, M )
 		local m, chi;
 		
 		m := Matroid ( Concatenation( L, M ), HomalgFieldOfRationals( ) );
 		chi := Coloring( m, Length( L ), fail, true );
 		
 		return OrlikSolomonBicomplexRecord( m, chi );
 
 end );

InstallMethod( OrlikSolomonBicomplexRecord,

	[ IsList, IsList, IsBool, IsBool ],
 	
 	function( L, M, default, last )
 		local m, chi;
 		
 		m := Matroid ( Concatenation( L, M ), HomalgFieldOfRationals( ) );
 		chi := Coloring( m, Length( L ), default, last );
 		
 		return OrlikSolomonBicomplexRecord( m, chi );
 
 end );
 
 InstallMethod( OrlikSolomonBicomplexRecord,
        [ IsMatroid, IsList ],
        
	function( m, chi )
		local matrixcat;
		
		matrixcat := MatrixCategory( HomalgFieldOfRationals() );
		
		CapCategorySwitchLogicOff( matrixcat );
		
		DeactivateCachingOfCategory( matrixcat );
		
		DisableInputSanityChecks( matrixcat );
		
		SetCachingToCrisp( matrixcat, "ZeroMorphism" );
		
		return OrlikSolomonBicomplexRecord( m, chi, matrixcat );

end );

InstallMethod( OrlikSolomonBicomplexRecord,
	[ IsMatroid, IsList, IsCapCategory ],

	function( m, chi, cat )
    	local A, k, i, j, phi, d, D, s, t, u, tt, uu, a, psi, obj, c, b, r, res, phiF, phiG, psiF, psiG, FG;

	    A := rec(
    		cat := cat,
    		matroid := m,
    		flats := Flats( m ),
    		rank := RankOfMatroid( m ),
    		);

		A.nr_flats := List( A.flats, Length );
		A.coloring := chi;
		A.incidence := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank ], j ->
						List( A.flats[i + 1], s -> List( A.flats[j + 1], t -> IsSubset( s, t ) ) ) ) );


		A.F := rec(
			objects := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank - i ], j -> [ ] ) ),
			differentials_horizontal := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) ),
			differentials_vertical := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) )
			);
	
    	A.F.objects[1][1][1] := TensorUnit( cat );

		
		A.G := rec(
			objects := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank - i ], j -> [ ] ) ),
			differentials_horizontal := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) ),
			differentials_vertical := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) )
			);
	
    	A.G.objects[1][1][1] := TensorUnit( cat );
    	
    	A.morphisms := List(
    		[ 0 .. A.rank ],
    		i -> List(
    			[ 0 .. A.rank - i ],
    			j -> [ ]
    			)
    		);
    
    	A.morphisms[1][1][1] := IdentityMorphism( TensorUnit( cat ) );
	
	    for k in [ 1 .. A.rank ] do
			if k = A.rank then
#				Print( "\nCodimension ", k, " (1 flat)\n");
			else
#				Print( "\nCodimension ", k, " (", A.nr_flats[k + 1], " flats)\n");
			fi;
		
			for s in [ 1 .. A.nr_flats[k + 1] ] do
			
#				Print( ".\c" );

				c := A.coloring[k][s];
				
				tt := Positions( A.incidence[k + 1][k][s], true );

				if c = true then

					for i in Reversed( [ 0 .. k ] ) do
						for FG in [ "F", "G" ] do
							j := k - i;
							phi := OrlikSolomonBicomplexHorizontalDifferential( A, FG, i - 1, j, k, s );
							d := KernelEmbedding( phi );
							SetOrlikSolomonBicomplexObject( A, FG, i, j, s, Source( d ) );
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i - 1, j, k - 1, t ) );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
									A, FG, i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( d, D, a )
								);
							od;
							if i < k then
								D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j - 1, k - 1, t ) );
								psi := PreCompose(
											   OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, j - 1, k, s ),
											   OrlikSolomonBicomplexVerticalDifferential( A, FG, i - 1, j - 1, k, s )
								);
								psi := KernelLift( phi, psi );
								for a in [ 1 .. Length( tt ) ] do
									SetOrlikSolomonBicomplexVerticalDifferentialComponent(
										A, FG, i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( psi, D, a )
									);
								od;
							fi;
						od;

						SetOrlikSolomonBicomplexMorphism(
							A, i, j, s, KernelObjectFunctorial(
								OrlikSolomonBicomplexHorizontalDifferential( A, "F", i - 1, j, k, s ),
								OrlikSolomonBicomplexMorphism( A, i - 1, j, k, s ),
								OrlikSolomonBicomplexHorizontalDifferential( A, "G", i - 1, j, k, s )
								)
						);
					od;

				elif c = false then

					for i in [ 0 .. k ] do
						for FG in [ "F", "G" ] do
							j := k - i;
							phi := OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j - 2, k, s );
							d := CokernelProjection( phi );
							SetOrlikSolomonBicomplexObject( A, FG, i, j, s, Range( d ) );
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j - 1, k - 1, t ) );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexVerticalDifferentialComponent(
									A, FG, i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( d, D, a )
								);
							od;
							if i > 0 then
								D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i - 1, j, k - 1, t ) );
								psi := PreCompose(
											   OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, j - 1, k, s ),
											   OrlikSolomonBicomplexVerticalDifferential( A, FG, i - 1, j - 1, k, s )
								);
								psi := CokernelColift( phi, psi );
								for a in [ 1 .. Length( tt ) ] do
									SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
										A, FG, i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( psi, D, a )
									);
								od;
							fi;
						od;
						
						SetOrlikSolomonBicomplexMorphism(
							A, i, j, s, CokernelObjectFunctorial(
								OrlikSolomonBicomplexVerticalDifferential( A, "F", i, j - 2, k, s ),
								OrlikSolomonBicomplexMorphism( A, i, j - 1, k, s ),
								OrlikSolomonBicomplexVerticalDifferential( A, "G", i, j - 2, k, s )
								)
						);

					od;


				else
				
					for i in [ 0 .. k ] do
						j := k - i;
						phiF := OrlikSolomonBicomplexVerticalDifferential( A, "F", i, j - 2, k, s );
						d := CokernelProjection( phiF );
						SetOrlikSolomonBicomplexObject( A, "F", i, j, s, Range( d ) );
						D := List( tt, t -> OrlikSolomonBicomplexObject( A, "F", i, j - 1, k - 1, t ) );
						for a in [ 1 .. Length( tt ) ] do
							SetOrlikSolomonBicomplexVerticalDifferentialComponent(
								A, "F", i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( d, D, a )
							);
						od;
#						if i > 0 then
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, "F", i - 1, j, k - 1, t ) );
							psiF := PreCompose(
										   OrlikSolomonBicomplexHorizontalDifferential( A, "F", i, j - 1, k, s ),
										   OrlikSolomonBicomplexVerticalDifferential( A, "F", i - 1, j - 1, k, s )
							);
							psiF := CokernelColift( phiF, psiF );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
									A, "F", i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( psiF, D, a )
								);
							od;
#						fi;
						
						phiG := OrlikSolomonBicomplexHorizontalDifferential( A, "G", i - 1, j, k, s );
						d := KernelEmbedding( phiG );
						SetOrlikSolomonBicomplexObject( A, "G", i, j, s, Source( d ) );
						D := List( tt, t -> OrlikSolomonBicomplexObject( A, "G", i - 1, j, k - 1, t ) );
						for a in [ 1 .. Length( tt ) ] do
							SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
								A, "G", i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( d, D, a )
							);
						od;
#						if i < k then
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, "G", i, j - 1, k - 1, t ) );
							psiG := PreCompose(
										   OrlikSolomonBicomplexHorizontalDifferential( A, "G", i, j - 1, k, s ),
										   OrlikSolomonBicomplexVerticalDifferential( A, "G", i - 1, j - 1, k, s )
							);
							psiG := KernelLift( phiG, psiG );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexVerticalDifferentialComponent(
									A, "G", i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( psiG, D, a )
								);
							od;
#						fi;

						SetOrlikSolomonBicomplexMorphism(
							A, i, j, s, KernelLift(
								phiG,
								PreCompose( psiF, OrlikSolomonBicomplexMorphism( A, i - 1, j, k, s ) )
							)
						);
						
					od;


				fi;
			od;
		od;
		
#		Print( "\n\n" );
		SemisimplifiedMotive( A );
 
 		return A;
 
end );

####################

InstallMethod( OrlikSolomonBicomplexRecord,
	[ IsMatroid, IsList, IsCapCategory ],

	function( m, chi, cat )
    	local A, k, i, j, phi, d, D, s, t, u, tt, uu, a, psi, obj, c, b, r, res, phiF, phiG, psiF, psiG, FG;

	    A := rec(
    		cat := cat,
    		matroid := m,
    		flats := Flats( m ),
    		rank := RankOfMatroid( m ),
    		);

		A.nr_flats := List( A.flats, Length );
		A.coloring := chi;
		A.incidence := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank ], j ->
						List( A.flats[i + 1], s -> List( A.flats[j + 1], t -> IsSubset( s, t ) ) ) ) );


		A.F := rec(
			objects := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank - i ], j -> [ ] ) ),
			differentials_horizontal := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) ),
			differentials_vertical := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) ),
            homology_objects := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank - i - 1 ], j -> [ ] ) )
			);
	
    	A.F.objects[1][1][1] := TensorUnit( cat );

		
		A.G := rec(
			objects := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank - i ], j -> [ ] ) ),
			differentials_horizontal := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) ),
			differentials_vertical := List( [ 0 .. A.rank ],
                                    i -> List( [ 0 .. A.rank - i ],
                                      j -> List( [ 1 .. A.nr_flats[i + j + 1] ], k -> [ ] ) ) ),
            homology_objects := List( [ 0 .. A.rank ], i -> List( [ 0 .. A.rank - i - 1 ], j -> [ ] ) )
			);
	
    	A.G.objects[1][1][1] := TensorUnit( cat );
    	
    	A.morphisms := List(
    		[ 0 .. A.rank ],
    		i -> List(
    			[ 0 .. A.rank - i ],
    			j -> [ ]
    			)
    		);
    
    	A.morphisms[1][1][1] := IdentityMorphism( TensorUnit( cat ) );
	
	    for k in [ 1 .. A.rank ] do
			if k = A.rank then
#				Print( "\nCodimension ", k, " (1 flat)\n");
			else
#				Print( "\nCodimension ", k, " (", A.nr_flats[k + 1], " flats)\n");
			fi;
		
			for s in [ 1 .. A.nr_flats[k + 1] ] do
			
#				Print( ".\c" );

				c := A.coloring[k][s];
				
				tt := Positions( A.incidence[k + 1][k][s], true );

				if c = true then

					for i in Reversed( [ 0 .. k ] ) do
						for FG in [ "F", "G" ] do
							j := k - i;
							phi := OrlikSolomonBicomplexHorizontalDifferential( A, FG, i - 1, j, k, s );
							d := KernelEmbedding( phi );
							SetOrlikSolomonBicomplexObject( A, FG, i, j, s, Source( d ) );
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i - 1, j, k - 1, t ) );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
									A, FG, i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( d, D, a )
								);
							od;
							if i < k then
								D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j - 1, k - 1, t ) );
								psi := PreCompose(
											   OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, j - 1, k, s ),
											   OrlikSolomonBicomplexVerticalDifferential( A, FG, i - 1, j - 1, k, s )
								);
								psi := KernelLift( phi, psi );
								for a in [ 1 .. Length( tt ) ] do
									SetOrlikSolomonBicomplexVerticalDifferentialComponent(
										A, FG, i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( psi, D, a )
									);
								od;
							fi;
						od;

						SetOrlikSolomonBicomplexMorphism(
							A, i, j, s, KernelObjectFunctorial(
								OrlikSolomonBicomplexHorizontalDifferential( A, "F", i - 1, j, k, s ),
								OrlikSolomonBicomplexMorphism( A, i - 1, j, k, s ),
								OrlikSolomonBicomplexHorizontalDifferential( A, "G", i - 1, j, k, s )
								)
						);
					od;

				elif c = false then

					for i in [ 0 .. k ] do
						for FG in [ "F", "G" ] do
							j := k - i;
							phi := OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j - 2, k, s );
							d := CokernelProjection( phi );
							SetOrlikSolomonBicomplexObject( A, FG, i, j, s, Range( d ) );
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j - 1, k - 1, t ) );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexVerticalDifferentialComponent(
									A, FG, i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( d, D, a )
								);
							od;
							if i > 0 then
								D := List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i - 1, j, k - 1, t ) );
								psi := PreCompose(
											   OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, j - 1, k, s ),
											   OrlikSolomonBicomplexVerticalDifferential( A, FG, i - 1, j - 1, k, s )
								);
								psi := CokernelColift( phi, psi );
								for a in [ 1 .. Length( tt ) ] do
									SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
										A, FG, i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( psi, D, a )
									);
								od;
							fi;
						od;
						
						SetOrlikSolomonBicomplexMorphism(
							A, i, j, s, CokernelObjectFunctorial(
								OrlikSolomonBicomplexVerticalDifferential( A, "F", i, j - 2, k, s ),
								OrlikSolomonBicomplexMorphism( A, i, j - 1, k, s ),
								OrlikSolomonBicomplexVerticalDifferential( A, "G", i, j - 2, k, s )
								)
						);

					od;


				else
				
					for i in [ 0 .. k ] do
						j := k - i;
						phiF := OrlikSolomonBicomplexVerticalDifferential( A, "F", i, j - 2, k, s );
						d := CokernelProjection( phiF );
						SetOrlikSolomonBicomplexObject( A, "F", i, j, s, Range( d ) );
						D := List( tt, t -> OrlikSolomonBicomplexObject( A, "F", i, j - 1, k - 1, t ) );
						for a in [ 1 .. Length( tt ) ] do
							SetOrlikSolomonBicomplexVerticalDifferentialComponent(
								A, "F", i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( d, D, a )
							);
						od;
#						if i > 0 then
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, "F", i - 1, j, k - 1, t ) );
							psiF := PreCompose(
										   OrlikSolomonBicomplexHorizontalDifferential( A, "F", i, j - 1, k, s ),
										   OrlikSolomonBicomplexVerticalDifferential( A, "F", i - 1, j - 1, k, s )
							);
							psiF := CokernelColift( phiF, psiF );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
									A, "F", i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( psiF, D, a )
								);
							od;
#						fi;
						
						phiG := OrlikSolomonBicomplexHorizontalDifferential( A, "G", i - 1, j, k, s );
						d := KernelEmbedding( phiG );
						SetOrlikSolomonBicomplexObject( A, "G", i, j, s, Source( d ) );
						D := List( tt, t -> OrlikSolomonBicomplexObject( A, "G", i - 1, j, k - 1, t ) );
						for a in [ 1 .. Length( tt ) ] do
							SetOrlikSolomonBicomplexHorizontalDifferentialComponent(
								A, "G", i, j, s, tt[a], ComponentOfMorphismIntoDirectSum( d, D, a )
							);
						od;
#						if i < k then
							D := List( tt, t -> OrlikSolomonBicomplexObject( A, "G", i, j - 1, k - 1, t ) );
							psiG := PreCompose(
										   OrlikSolomonBicomplexHorizontalDifferential( A, "G", i, j - 1, k, s ),
										   OrlikSolomonBicomplexVerticalDifferential( A, "G", i - 1, j - 1, k, s )
							);
							psiG := KernelLift( phiG, psiG );
							for a in [ 1 .. Length( tt ) ] do
								SetOrlikSolomonBicomplexVerticalDifferentialComponent(
									A, "G", i, j - 1, tt[a], s, ComponentOfMorphismFromDirectSum( psiG, D, a )
								);
							od;
#						fi;

						SetOrlikSolomonBicomplexMorphism(
							A, i, j, s, KernelLift(
								phiG,
								PreCompose( psiF, OrlikSolomonBicomplexMorphism( A, i - 1, j, k, s ) )
							)
						);
						
					od;


				fi;
			od;
		od;
		
#		Print( "\n\n" );
		SemisimplifiedMotive( A );
 
 		return A;
 
end );


 
InstallMethod( OrlikSolomonBicomplexMorphism,

	[ IsRecord, IsInt, IsInt, IsInt, IsInt ],

	function( A, i, j, k, s )
	
		local tt, t;
		
		if i < 0 or j < 0 or i + j > k then
			return ZeroObjectFunctorial( A.cat );
		fi;
		
		tt := Positions( A.incidence[k + 1][i + j + 1][s], true );
		
		return MorphismBetweenDirectSums(
			DirectSum( A.cat, List( tt, t -> OrlikSolomonBicomplexObject( A, "F", i, j, i + j, t ) ) ),
			List( tt, t1 -> List( tt, function( t2 )
				if t1 = t2 then
					return A.morphisms[i + 1][j + 1][t1];
				else
					return ZeroMorphism(
						OrlikSolomonBicomplexObject( A, "F", i, j, i + j, t1 ),
						OrlikSolomonBicomplexObject( A, "G", i, j, i + j, t2 )
					);
				fi;
			end ) ),
			DirectSum( A.cat, List( tt, t -> OrlikSolomonBicomplexObject( A, "G", i, j, i + j, t ) ) )	
		);
	
end );

InstallMethod( SetOrlikSolomonBicomplexMorphism,

	[ IsRecord, IsInt, IsInt, IsInt, IsCapCategoryMorphism ],

	function( A, i, j, s, f )
		
		A.morphisms[i + 1][j + 1][s] := f;

end );

InstallMethod( OrlikSolomonBicomplexObject,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, k, s )
		local tt;
	
		if i < 0 or j < 0 or i + j > k then
			return ZeroObject( A.cat );		
#		if HasOrlikSolomonBicomplexObject( A, FG, i, j, s ) then
		elif i + j = k then
			return A.(FG).objects[i + 1][j + 1][s];
		else
 			tt := Positions( A.incidence[k + 1][i + j + 1][s], true );
 			return DirectSum( A.cat, List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j, i + j, t ) ) );
		fi;
end );

InstallMethod( SetOrlikSolomonBicomplexObject,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsCapCategoryObject ],

	function( A, FG, i, j, s, V )

		A.(FG).objects[i + 1][j + 1][s] := V;

end );

InstallMethod( OrlikSolomonBicomplexHorizontalDifferential,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, k, s )
		local tt, uu;
		
		if i <= 0 or j < 0 then
			return UniversalMorphismIntoZeroObject( OrlikSolomonBicomplexObject( A, FG, i, j, k, s ) );
		fi;
		
		if i + j > k then
			return UniversalMorphismFromZeroObject( OrlikSolomonBicomplexObject( A, FG, i - 1, j, k, s ) );
		fi;
		
		tt := Positions( A.incidence[k + 1][i + j + 1][s], true ); 
		uu := Positions( A.incidence[k + 1][i + j][s], true );
		return MorphismBetweenDirectSums(
			DirectSum( A.cat, List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j, i + j, t ) ) ),
#			OrlikSolomonBicomplexObject( A, FG, i, j, k, s ),
			List( tt, t -> List( uu, u -> OrlikSolomonBicomplexHorizontalDifferentialComponent( A, FG, i, j, t, u ) ) ),
#			OrlikSolomonBicomplexObject( A, FG, i - 1, j, k, s )
			DirectSum( A.cat, List( uu, u -> OrlikSolomonBicomplexObject( A, FG, i - 1, j, i + j - 1, u ) ) )
 		);

end );

InstallMethod( OrlikSolomonBicomplexVerticalDifferential,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, k, s )
		local tt, uu;
		
		if i < 0 or j < 0 then
			return UniversalMorphismFromZeroObject( OrlikSolomonBicomplexObject( A, FG, i, j + 1, k, s ) );
		fi;

		if i + j >= k then
			return UniversalMorphismIntoZeroObject( OrlikSolomonBicomplexObject( A, FG, i, j, k, s ) );
		fi;
		
		tt := Positions( A.incidence[k + 1][i + j + 2][s], true ); 
		uu := Positions( A.incidence[k + 1][i + j + 1][s], true );
		return MorphismBetweenDirectSums(
			DirectSum( A.cat, List( uu, u -> OrlikSolomonBicomplexObject( A, FG, i, j, i + j, u ) ) ),
			List( uu, u -> List( tt, t -> OrlikSolomonBicomplexVerticalDifferentialComponent( A, FG, i, j, u, t ) ) ),
			DirectSum( A.cat, List( tt, t -> OrlikSolomonBicomplexObject( A, FG, i, j + 1, i + j + 1, t ) ) )
 		);

end );

InstallMethod( HasOrlikSolomonBicomplexHorizontalDifferentialComponent,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, t, u )
	
		return IsBound( A.(FG).differentials_horizontal[i + 1][j + 1][t][u] );

end );

InstallMethod( OrlikSolomonBicomplexHorizontalDifferentialComponent,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, t, u )
		
		if HasOrlikSolomonBicomplexHorizontalDifferentialComponent( A, FG, i, j, t, u ) then
			return A.(FG).differentials_horizontal[i + 1][j + 1][t][u];
		else
			return ZeroMorphism(
				OrlikSolomonBicomplexObject( A, FG, i, j, i + j, t ),
				OrlikSolomonBicomplexObject( A, FG, i - 1, j, i + j - 1, u )
			);
		fi;
end );

InstallGlobalFunction( SetOrlikSolomonBicomplexHorizontalDifferentialComponent,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt, IsCapCategoryMorphism ],

	function( A, FG, i, j, t, u, f )
	
		A.(FG).differentials_horizontal[i + 1][j + 1][t][u] := f;

end );

InstallMethod( HasOrlikSolomonBicomplexVerticalDifferentialComponent,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, t, u )
	
		return IsBound( A.(FG).differentials_vertical[i + 1][j + 1][t][u] );

end );


InstallMethod( OrlikSolomonBicomplexVerticalDifferentialComponent,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],

	function( A, FG, i, j, t, u )
		
		if HasOrlikSolomonBicomplexVerticalDifferentialComponent( A, FG, i, j, t, u ) then
			return A.(FG).differentials_vertical[i + 1][j + 1][t][u];
		else
			return ZeroMorphism(
				OrlikSolomonBicomplexObject( A, FG, i, j, i + j, t ),
				OrlikSolomonBicomplexObject( A, FG, i, j + 1, i + j + 1, u )
			);
		fi;
end );

InstallGlobalFunction( SetOrlikSolomonBicomplexVerticalDifferentialComponent,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt, IsCapCategoryMorphism ],

	function( A, FG, i, j, t, u, f )
	
		if j >= 0 then
			A.(FG).differentials_vertical[i + 1][j + 1][t][u] := f;
		fi;

end );
 
InstallMethod( OrlikSolomonBicomplexHorizontalHomologyObject,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],
	
	function ( A, FG, i, j, k, s )
	
	return HomologyObject(
	    OrlikSolomonBicomplexHorizontalDifferential( A, FG, i + 1, j, k, s ),
	    OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, j, k, s )
	);
	
end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyObject,

	[ IsRecord, IsString, IsInt, IsInt, IsInt, IsInt ],
	
	function ( A, FG, i, j, k, s )
	
	return HomologyObject(
	    OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j - 1, k, s ),
	    OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j, k, s )
	);
	
end );

InstallMethod( OrlikSolomonBicomplexHorizontalHomologyMorphism,

	[ IsRecord, IsInt, IsInt, IsInt, IsInt ],
	
	function ( A, i, j, k, s )
	
	return HomologyObjectFunctorial(
	    OrlikSolomonBicomplexHorizontalDifferential( A, "F", i + 1, j, k, s ),
	    OrlikSolomonBicomplexHorizontalDifferential( A, "F", i, j, k, s ),
	    OrlikSolomonBicomplexMorphism( A, i, j, k, s ),
	    OrlikSolomonBicomplexHorizontalDifferential( A, "G", i + 1, j, k, s ),
	    OrlikSolomonBicomplexHorizontalDifferential( A, "G", i, j, k, s )
	);
	
end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyMorphism,

	[ IsRecord, IsInt, IsInt, IsInt, IsInt ],
	
	function ( A, i, j, k, s )
	
	return HomologyObjectFunctorial(
	    OrlikSolomonBicomplexVerticalDifferential( A, "F", i, j - 1, k, s ),
	    OrlikSolomonBicomplexVerticalDifferential( A, "F", i, j, k, s ),
	    OrlikSolomonBicomplexMorphism( A, i, j, k, s ),
	    OrlikSolomonBicomplexVerticalDifferential( A, "G", i, j - 1, k, s ),
	    OrlikSolomonBicomplexVerticalDifferential( A, "G", i, j, k, s )
	);
	
end );

InstallMethod( SemisimplifiedMotive,

    [ IsRecord ],
    
    function( A )
    
    local n, c;
    
    n := A.rank;
    c := A.coloring[n][1];
    
    if c then
        return List( [ 0 .. n - 1 ], k -> Dimension( ImageObject( OrlikSolomonBicomplexVerticalHomologyMorphism( A, k + 1, n - 1 - k, n, 1 ) ) ) );
    else
        return List( [ 0 .. n - 1 ], k -> Dimension( ImageObject( OrlikSolomonBicomplexHorizontalHomologyMorphism( A, k, n - k, n, 1 ) ) ) );
    fi;
    
end );
    
# InstallMethod( SemisimplifiedMotive,
# 
# 	[ IsRecord, IsBool ],
# 
# 	function( A, last )
# 	
# 		local n, phi, k;
# 						
# 		A.dimgrBlue := [ ];
# 		A.dimgrRed := [ ];
# 		A.dimgrImage := [ ];
# 
# 		
# 		n := A.rank - 1;
# 
# 		if last then
# 		
# 			for k in [ 0 .. n ] do
# 				phi := CokernelObjectFunctorial(
# 					OrlikSolomonBicomplexVerticalDifferential( A, "F", k + 1, n - k - 1, n + 1, 1 ),
# 					OrlikSolomonBicomplexMorphism( A, k + 1, n - k, n + 1, 1 ),
# 					OrlikSolomonBicomplexVerticalDifferential( A, "G", k + 1, n - k - 1, n + 1, 1 )
# 				);
# 				Add( A.dimgrRed, Dimension( Source( phi ) ) );
# 				Add( A.dimgrBlue, Dimension( Range( phi ) ) );
# 				Add( A.dimgrImage, Dimension( ImageObject( phi ) ) );
# 			od;
# 		
# 		else
# 		
# 			for k in [ 0 .. n ] do
# 				phi := KernelObjectFunctorial(
# 					OrlikSolomonBicomplexHorizontalDifferential( A, "F", k, n - k + 1, n + 1, 1 ),
# 					OrlikSolomonBicomplexMorphism( A, k, n - k + 1, n + 1, 1 ),
# 					OrlikSolomonBicomplexHorizontalDifferential( A, "G", k, n - k + 1, n + 1, 1 )
# 				);
# 				Add( A.dimgrRed, Dimension( Source( phi ) ) );
# 				Add( A.dimgrBlue, Dimension( Range( phi ) ) );
# 				Add( A.dimgrImage, Dimension( ImageObject( phi ) ) );
# 			od;
# 		
# 		fi;
# 
# #	PrintArray( [ A.dimgrBlue, A.dimgrImage, A.dimgrRed ] );
# #	Print( "\n" );
# 
# end );


InstallMethod( OrlikSolomonBicomplexDimensions,

	[ IsRecord, IsString ],

	function( A, FG )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
			for j in [ 0 .. r - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexObject( A, FG, i, j, r, 1 ) );
			od;
			for j in [ r - i + 1 .. r ] do
				res[i + 1][j + 1] := 0;
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexDimensions,

	[ IsRecord ],

	function( A )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
			for j in [ 0 .. r - i ] do
				res[i + 1][j + 1] := Dimension( ImageObject( OrlikSolomonBicomplexMorphism( A, i, j, r, 1 ) ) );
			od;
			for j in [ r - i + 1 .. r ] do
				res[i + 1][j + 1] := 0;
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexHorizontalHomologyDimensions,

	[ IsRecord, IsString ],

	function( A, FG )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
			for j in [ 0 .. r - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexHorizontalHomologyObject( A, FG, i, j, r, 1 ) );
			od;
			for j in [ r - i + 1 .. r ] do
				res[i + 1][j + 1] := 0;
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexHorizontalHomologyDimensions,

	[ IsRecord ],

	function( A )
		local r, res, i, j;
	
		r := A.rank;

		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
		    for j in [ 0 .. r ] do
		        res[i + 1][j + 1] := 0;
		    od;
		od;

		for i in [ 0 .. r ] do
			for j in [ 0 .. r - i ] do
				res[i + 1][j + 1] := Dimension( ImageObject( OrlikSolomonBicomplexHorizontalHomologyMorphism( A, i, j, r, 1 ) ) );
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexHorizontalHomologyDimensions,

	[ IsRecord, IsString, IsInt, IsInt ],

	function( A, FG, k, s )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
		    for j in [ 0 .. r ] do
		        res[i + 1][j + 1] := 0;
		    od;
		od;
		for i in [ 0 .. k ] do
			for j in [ 0 .. k - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexHorizontalHomologyObject( A, FG, i, j, k, s ) );
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexHorizontalHomologyDimensions,

	[ IsRecord, IsInt, IsInt ],

	function( A, k, s )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
		    for j in [ 0 .. r ] do
		        res[i + 1][j + 1] := 0;
		    od;
		od;
		for i in [ 0 .. k ] do
			for j in [ 0 .. k - i ] do
				res[i + 1][j + 1] := Dimension( ImageObject( OrlikSolomonBicomplexHorizontalHomologyMorphism( A, i, j, k, s ) ) );
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyDimensions,

	[ IsRecord, IsString ],

	function( A, FG )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
			for j in [ 0 .. r - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexVerticalHomologyObject( A, FG, i, j, r, 1 ) );
			od;
			for j in [ r - i + 1 .. r ] do
				res[i + 1][j + 1] := 0;
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyDimensions,

	[ IsRecord ],

	function( A )
		local r, res, i, j;
	
		r := A.rank;

		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
		    for j in [ 0 .. r ] do
		        res[i + 1][j + 1] := 0;
		    od;
		od;

		for i in [ 0 .. r ] do
			for j in [ 0 .. r - i ] do
				res[i + 1][j + 1] := Dimension( ImageObject( OrlikSolomonBicomplexVerticalHomologyMorphism( A, i, j, r, 1 ) ) );
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyDimensions,

	[ IsRecord, IsString, IsInt, IsInt ],

	function( A, FG, k, s )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
		    for j in [ 0 .. r ] do
		        res[i + 1][j + 1] := 0;
		    od;
		od;
		for i in [ 0 .. k ] do
			for j in [ 0 .. k - i ] do
				res[i + 1][j + 1] := Dimension( OrlikSolomonBicomplexVerticalHomologyObject( A, FG, i, j, k, s ) );
			od;
		od;
		return res;
end );

InstallMethod( OrlikSolomonBicomplexVerticalHomologyDimensions,

	[ IsRecord, IsInt, IsInt ],

	function( A, k, s )
		local r, res, i, j;
	
		r := A.rank;
	
		res := [];
		for i in [ 0 .. r ] do
			res[i + 1] := [];
		    for j in [ 0 .. r ] do
		        res[i + 1][j + 1] := 0;
		    od;
		od;
		for i in [ 0 .. k ] do
			for j in [ 0 .. k - i ] do
				res[i + 1][j + 1] := Dimension( ImageObject( OrlikSolomonBicomplexVerticalHomologyMorphism( A, i, j, k, s ) ) );
			od;
		od;
		return res;
end );

# InstallMethod( IsExactOrlikSolomonBicomplex,
# 
#     [ IsRecord ],
#     
#     function( A )
#     
#     return ForAll(
#             [ 2 .. A.rank ],
#             k -> ForAll(
#                 [ 1 .. A.nr_flats[k + 1] ],
#                 s -> B.coloring[k][s] = true or 
#                     ( IsZero( OrlikSolomonBicomplexVerticalHomologyDimensions( B, "F", k, s ) ) 
#                         and IsZero( OrlikSolomonBicomplexVerticalHomologyDimensions( B, "G", k, s ) ) )
#                 )
#                 and
#                 ForAll(
#                 [ 1 .. A.nr_flats[k + 1] ],
#                 s -> A.coloring[k][s] = false or 
#                     ( IsZero( OrlikSolomonBicomplexHorizontalHomologyDimensions( B, "F", k, s ) ) 
#                         and IsZero( OrlikSolomonBicomplexHorizontalHomologyDimensions( B, "G", k, s ) ) )
#                 )        
#             );
# end );
# 
# InstallMethod( DefectOfExactness,
# 
#     [ IsRecord ],
#     
#     function( A )
#     
#     return List(
#             [ 2 .. B.rank ],
#             k -> [ Filtered(
#                     [ 1 .. B.nr_flats[k + 1] ],
#                     s -> B.coloring[k][s] = false and 
#                         ( not IsZero( OrlikSolomonBicomplexVerticalHomologyDimensions( B, "F", k, s ) ) 
#                             or not IsZero( OrlikSolomonBicomplexVerticalHomologyDimensions( B, "G", k, s ) ) )
#                 ),
#                 Filtered(
#                     [ 1 .. B.nr_flats[k + 1] ],
#                     s -> B.coloring[k][s] = true and 
#                         ( not IsZero( OrlikSolomonBicomplexHorizontalHomologyDimensions( B, "F", k, s ) ) 
#                             or not IsZero( OrlikSolomonBicomplexHorizontalHomologyDimensions( B, "G", k, s ) ) )
#                 ),
#                 Filtered(
#                     [ 1 .. B.nr_flats[k + 1] ],
#                     s -> B.coloring[k][s] = fail and 
#                         ( not IsZero( OrlikSolomonBicomplexVerticalHomologyDimensions( B, "F", k, s ) ) 
#                             or not IsZero( OrlikSolomonBicomplexHorizontalHomologyDimensions( B, "G", k, s ) ) )
#                 )
#                 ]
#             );
# end );



# OrlikSolomonBicomplexDimensionsImages :=
# #		[ IsRecord, IsList ],
# 
# 	function( A )
# 		local r, res, i, j;
# 	
# 		r := A.rank;
# 	
# 		res := [];
# 		for i in [ 0 .. r ] do
# 			res[i + 1] := [];
# 			for j in [ 0 .. r - i ] do
# 				res[i + 1][j + 1] := Dimension( ImageObject( OrlikSolomonBicomplexMorphism( A, i, j, r, 1 ) ) );
# 			od;
# 			for j in [ r - i + 1 .. r ] do
# 				res[i + 1][j + 1] := 0;
# 			od;
# 		od;
# 		return res;
# end;

InstallMethod( IsTameFlat,

	[ IsMatroid, IsFunction, IsInt, IsInt ],

	function( mat, chi, k, s )
	
		local S, C, l, cl, c;
	
		S := Flats( mat )[k + 1][s];
	
		cl := ClosureOperator( mat );
	
		c := chi( S );
	
		if c = true then
			C := [ true ];
		elif c = false then
			C := [ false ];
		else
			C := [ true, false ];
		fi;
	
		l := List( C, c -> Filtered( S, h -> chi( [ h ] ) = c ) );
	
		return ForAll( [ 1 .. Length( C ) ], i -> [ ] <> Difference( l[i], Union( Filtered( Circuits( mat ), T -> IsSubset( S, T ) and chi( cl( T ) ) = ( not C[i] ) ) ) ) );
	
end );

# InstallMethod( IsTame,
# 
# 	[ IsMatroid, IsFunction ],
# 
# 	function( mat, chi )
# 
# 	local flats;
# 	
# 	flats := Flats( mat );
# 	
# 	return ForAll( [ 1 .. Length( flats ) - 1 ], k -> ForAll( [ 1 .. Length( flats[k + 1] ) ], s -> IsTameFlat( mat, chi, k, s ) ) );
# 	
# end );

InstallMethod( IsTame,

 	[ IsMatroid, IsFunction ],

	function( mat, chi )

		local flats, k, s, res;
	
		res := true;
	
		flats := Flats( mat );
	
		for k in [ 1 .. Length( flats ) - 1 ] do
			for s in [ 1 .. Length( flats[k + 1] ) ] do
				if not IsTameFlat( mat, chi, k, s ) then
					Print( k, ", ", s, ", ", flats[k + 1][s], "\n" );
					res := false;
				fi;
			od;
		od;
	
		return res;
	
end );

InstallMethod( IsTameBiarrangement,

	[ IsList, IsFunction, IsBool, IsBool ],
 	
 	function( L, M, default, last )
 		local m, chi;
 		
 		m := Matroid ( Concatenation( L, M ), HomalgFieldOfRationals( ) );
 		chi := ColoringFunction( m, Length( L ), default, last );
 		
 		return IsTame( m, chi );
 
 end );


SemiSimplifiedMotiveByRectangles :=
	function( A )
	
	local n, rows, k, j, i, phi, res, morphism, L, FG;
	
	n := A.rank - 1;
	
	res := [ ];
	
	for FG in [ "F", "G" ] do

		rows := List(
			[ 0 .. n ],
			j -> ChainComplex(
				List(
					[ 1 .. n + 1 ],
					i -> OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, j, n + 1, 1 )
				),
				1
			)
		);
		
		phi := List(
			[ 1 .. n ],
			j -> ChainMorphism(
				rows[j],
				rows[j + 1],
				List( [ 0 .. n ], i -> OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j - 1, n + 1, 1 ) ),
				0
				)
		);
		
		A.(FG).complex := ChainComplex( Reversed( phi ), -1 );
		A.(FG).bicomplex := HomologicalBicomplex( ChainComplex( Reversed( phi ), -1 ) );
		
	od;

 L := List( [ -n - 1 .. 0 ], i -> ChainMorphism( A.G.complex[i], A.F.complex[i], List( [ 0 .. n + 1 - i ], j -> OrlikSolomonBicomplexMorphism( A, -j, i, n + 1, 1 ) ), 0 ) );
morphism := ChainMorphism( A.G.complex, A.F.complex, L, -n - 1 );
morphism := BicomplexMorphism( morphism );
	
	return morphism;
end;

SemiSimplifiedMotiveByRectangles :=
	function( A )
	
	local n, rows, cols, k, j, i, phi, morphism, L, FG;
	
	n := A.rank - 1;
	
	for FG in [ "F", "G" ] do

		cols := List(
			[ 0, -1 .. -n ],
			i -> ChainComplex(
				List(
					[ 0 .. n ],
					j -> OrlikSolomonBicomplexVerticalDifferential( A, FG, -i, j, n + 1, 1 )
				),
				1
			)
		);
		
		phi := List(
			[ 1 .. n ],
			j -> ChainMorphism(
				rows[j],
				rows[j + 1],
				List( [ 0 .. n ], i -> OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j - 1, n + 1, 1 ) ),
				0
				)
		);
		
		A.(FG).complex := ChainComplex( Reversed( phi ), -1 );
		A.(FG).bicomplex := HomologicalBicomplex( ChainComplex( Reversed( phi ), -1 ) );
		
	od;

 L := List( [ -n - 1 .. 0 ], i -> ChainMorphism( A.G.complex[i], A.F.complex[i], List( [ 0 .. n + 1 - i ], j -> OrlikSolomonBicomplexMorphism( A, -j, i, n + 1, 1 ) ), 0 ) );
morphism := ChainMorphism( A.G.complex, A.F.complex, L, -n - 1 );
morphism := BicomplexMorphism( morphism );
	
	return morphism;
end;
		
# 		for k in [ 0 .. n ] do
# 			BrutalTruncationAboveFunctor( CapCategory( complex[1] ), k + 1 );
# 			BrutalTruncation( BruntalTruncationFunctorial
# 		
# 		A.(FG).total := TotalComplex( A.(FG).bicomplex );
# 		
# 		res[k] := List( [ 0 .. 2 * n ], r -> Dimension( DefectOfExactnessAt( tot, 2 * k - r ) ) );
# 	od;
# 	
# 	res := Concatenation(
# 		[ List(
# 			[ 0 .. 2 * n ],
# 			r -> Dimension( DefectOfExactnessAt(
# 				ChainComplex( List( Reversed( [ 0 .. n - 1 ] ), j -> OrlikSolomonBicomplexVerticalDifferential( A, FG, 0, j, n + 1, 1 ) ), -n + 1 ),
# 				-r
# 			) )
# 		) ],
# 		res,
# 		[ List(
# 			[ 0 .. 2 * n ],
# 			r -> Dimension( DefectOfExactnessAt(
# 				ChainComplex( List( [ 1 .. n ], i -> OrlikSolomonBicomplexHorizontalDifferential( A, FG, i, 0, n + 1, 1 ) ), 1 ),
# 				r
# 			) )
# 		) ]
# 	);
# 	
# 	od;
# 	
# 	return res;
# 	
# end;

DisplayOrlikSolomonBicomplexDifferentials := function( A, FG )
	local n, i, j;
	
	n := A.rank - 1;

	for j in [ 0 .. n ] do
		for i in [ 0 .. n - j ] do
			Print( JoinStringsWithSeparator( [ i + 1, j, i, j ] ), "\n");
			Display( OrlikSolomonBicomplexHorizontalDifferential( A, FG, i + 1, j, n + 1, 1 ) );
			Print( "\n" );
		od;
	od;
	
	for i in [ 0 .. n ] do
		for j in [ 0 .. n - i ] do
			Print( JoinStringsWithSeparator( [ i, j, i, j + 1 ] ), "\n");
			Display( OrlikSolomonBicomplexVerticalDifferential( A, FG, i, j, n + 1, 1 ) );
			Print( "\n" );		od;
	od;
end;

CyclicIntervalFromList := function( list, k, l )
 local n;
 n := Length( list );
 if k + l - 1 <= n then
  return list{[ k .. k + l - 1 ]};
 else
  return list{Concatenation( [ k .. n ], [ 1 .. l - n + k - 1 ] )};
 fi;
end;

IsCyclicInterval := function ( list, n )
    local  p, k;
    p := Set( list );
    k := Length( p );
    if p[1] > 1 then
        return p = [ p[1] .. p[1] + k - 1 ];
    else
        return IsCyclicInterval( Difference( [ 1 .. n ], p ), n );
    fi;
end;


DihedralGroupPerm := function( n )
 return Group( PermList( Concatenation( [ 2 .. n ], [ 1 ] ) ), PermList( [ n, n - 1 .. 1 ] ) );
end;

IsConvergentConfiguration := function( pi )
 local n, std;
 n := Length( pi );
 std := [ 1 .. n];
 return ForAll( [ 2 .. n - 2 ], k -> ForAll( [ 1 .. n ], i -> not IsCyclicInterval( CyclicIntervalFromList( pi, i, k ), n ) ) );
end;


ConvergentConfigurations := function( n )
 local sym, dih, reps;
 sym := SymmetricGroup( n );
 dih := DihedralGroupPerm( n );
 reps := List( DoubleCosetRepsAndSizes( sym, dih, dih ), i -> i[1] );
 return Filtered( List( reps, w -> ListPerm( w, n ) ), IsConvergentConfiguration );
end;

CycleToTranspositions := function( c )
    return List( [ 1 .. Length(c) - 1 ], i -> c{[i, i + 1]} );
end;

InstallMethod( CellularIntegralOrlikSolomonBicomplexRecord,

	[ IsList ],

	function ( w )

	return OrlikSolomonBicomplexRecord(
		ArrangementFromGraph( CycleToTranspositions( w ) ),
		ArrangementFromGraph( CycleToTranspositions( [ 1 .. Length( w ) ] ) ),
		fail,
		true
	);
	
end );

BrownMotive := function( n )

    M := CycleToTranspositions( [ 1 .. n + 2 ] );
    L := Difference( Combinations( [ 1 .. n + 2 ], 2 ), M ); 

    return OrlikSolomonBicomplexRecord(
		ArrangementFromGraph( L ),
		ArrangementFromGraph( M ),
		false
		true,
		true
	);

end;


# IsTameCellularIntegralBiarrangement( [9, 2, 4, 1, 8, 6, 3, 5, 7], fail, true );
