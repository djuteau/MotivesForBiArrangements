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

InstallMethod( SemisimplifiedMotive,

	[ IsRecord ],

	function( A )
	
		local n, phi, last, k;
						
		A.dimgrBlue := [ ];
		A.dimgrRed := [ ];
		A.dimgrImage := [ ];

		
		n := A.rank - 1;

		last := A.coloring[ n + 1 ][1];

		if last then
		
			for k in [ 0 .. n ] do
				phi := CokernelObjectFunctorial(
					OrlikSolomonBicomplexVerticalDifferential( A, "F", k + 1, n - k - 1, n + 1, 1 ),
					OrlikSolomonBicomplexMorphism( A, k + 1, n - k, n + 1, 1 ),
					OrlikSolomonBicomplexVerticalDifferential( A, "G", k + 1, n - k - 1, n + 1, 1 )
				);
				Add( A.dimgrBlue, Dimension( Source( phi ) ) );
				Add( A.dimgrRed, Dimension( Range( phi ) ) );
				Add( A.dimgrImage, Dimension( ImageObject( phi ) ) );
			od;
		
		else
		
			for k in [ 0 .. n ] do
				phi := KernelObjectFunctorial(
					OrlikSolomonBicomplexHorizontalDifferential( A, "F", k, n - k + 1, n + 1, 1 ),
					OrlikSolomonBicomplexMorphism( A, k, n - k + 1, n + 1, 1 ),
					OrlikSolomonBicomplexHorizontalDifferential( A, "G", k, n - k + 1, n + 1, 1 )
				);
				Add( A.dimgrBlue, Dimension( Source( phi ) ) );
				Add( A.dimgrRed, Dimension( Range( phi ) ) );
				Add( A.dimgrImage, Dimension( ImageObject( phi ) ) );
			od;
		
		fi;

	PrintArray( [ A.dimgrBlue, A.dimgrImage, A.dimgrRed ] );
	Print( "\n" );

end );

InstallMethod( OrlikSolomonBicomplexRecord,
	[ IsMatroid, IsFunction, IsCapCategory ],

	function( m, chi, cat )
    	local A, k, i, j, phi, d, D, s, t, u, tt, uu, a, psi, obj, c, b, r, M, res, phiF, phiG, psiF, psiG, FG;

	    A := rec(
    		cat := cat,
    		matroid := m,
    		flats := Flats( m ),
    		rank := RankOfMatroid( m ),
    		);

		A.nr_flats := List( A.flats, Length );
		A.coloring := List( A.flats{ [ 2 .. Length( A.flats ) ] }, l -> List( l, chi ) );
		A.Smin := A.flats[Length( A.flats )][1];
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
				Print( "\nCodimension ", k, " (1 flat)\n");
			else
				Print( "\nCodimension ", k, " (", A.nr_flats[k + 1], " flats)\n");
			fi;
		
			for s in [ 1 .. A.nr_flats[k + 1] ] do
			
				Print( ".\c" );

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
				
					M := A;

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
		
		Print( "\n\n" );
		SemisimplifiedMotive( A );
 
 		return A;
 
end );

InstallMethod( OrlikSolomonBicomplexRecord,
        [ IsMatroid, IsFunction ],
        
	function( m, chi )
		local matrixcat;
		
		matrixcat := MatrixCategory( HomalgFieldOfRationals() );
		
		CapCategorySwitchLogicOff( matrixcat );
		
		DeactivateCachingOfCategory( matrixcat );
		
		DisableBasicOperationTypeCheck( matrixcat );
		
		SetCachingToCrisp( matrixcat, "ZeroMorphism" );
		
		return OrlikSolomonBicomplexRecord( m, chi, matrixcat );

end );


 
InstallMethod( OrlikSolomonBicomplexMorphism,

	[ IsRecord, IsInt, IsInt, IsInt, IsInt ],

	function( A, i, j, k, s )
	
		local tt, t;
		
		if i < 0 or j < 0 then
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

InstallMethod( IteratedIntegralArrangement,

	[ IsList ],

	function( a ) 			# a is the list of a_i's
	
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
end );

InstallMethod( CellularArrangement,

	[ IsInt, IsPerm ],

	function( n, w ) 

# choice of affine coordinates: t_i with (z_1, ..., z_n) = (t_1, t_2, t_3, ..., t_{n-3}, 1, infty, 0)
# choice of projective coordinates: (z_0, z_1, ..., z_n) with {z_0 = 0} the hyperplane at infinity,
# (for consistency with SimplexArrangement)

	local res, i, j, a, b, c;
	
	res := NullMat( n - 2, n - 2 );

# 	for i in [ 1 .. n - 2 ] do
# 		res[i] := [];
# 		for j in [ 1 .. n - 2 ] do
# 			res[i][j] := 0;
# 		od;
# 	od;
	
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
end );

InstallMethod( SimplexArrangement,

	[ IsInt ],

	function( n )
	
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
end );

InstallMethod( Coloring,

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

InstallMethod( OrlikSolomonBicomplexRecord,

	[ IsList, IsList, IsBool, IsBool ],
 	
 	function( L, M, default, last )
 		local m, chi;
 		
 		m := Matroid ( Concatenation( L, M ), HomalgFieldOfRationals( ) );
 		chi := Coloring( m, Length( L ), default, last );
 		
 		return OrlikSolomonBicomplexRecord( m, chi );
 
 end );


InstallMethod( CellularIntegralOrlikSolomonBicomplexRecord,

	[ IsList, IsBool, IsBool ],

	function ( w, default, last )

	return OrlikSolomonBicomplexRecord(
		CellularArrangement( Length( w ), PermList( w ) ),
		SimplexArrangement ( Length( w ) - 3 ),
		default,
		last
	);
	
end );

InstallMethod( IteratedIntegralOrlikSolomonBicomplexRecord,

	[ IsList, IsBool, IsBool ],

	function ( a, default, last )

	return OrlikSolomonBicomplexRecord(
		CellularArrangement( a ),
		SimplexArrangement ( Length( a ) ),
		default,
		last
	);
	
end );


InstallMethod( OrlikSolomonBicomplexDimensions,

	[ IsRecord, IsList ],

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

	[ IsList, IsList, IsBool, IsBool ],
 	
 	function( L, M, default, last )
 		local m, chi;
 		
 		m := Matroid ( Concatenation( L, M ), HomalgFieldOfRationals( ) );
 		chi := Coloring( m, Length( L ), default, last );
 		
 		return IsTame( m, chi );
 
 end );

InstallMethod( IsTameCellularIntegralBiarrangement,

	[ IsList, IsBool, IsBool ],

	function ( w, default, last )

		return IsTameBiarrangement(
			CellularArrangement( Length( w ), PermList( w ) ),
			SimplexArrangement ( Length( w ) - 3 ),
			default,
			last
		);
	
end );

InstallMethod( IsTameIteratedIntegralBiarrangement,

	[ IsList, IsBool, IsBool ],

	function ( a, default, last )

		return IsTameBiarrangement(
			IteratedIntegralArrangement( a ),
			SimplexArrangement ( Length( a ) ),
			default,
			last
		);
	
end );


# IsTameCellularIntegralBiarrangement( [9, 2, 4, 1, 8, 6, 3, 5, 7], fail, true );
