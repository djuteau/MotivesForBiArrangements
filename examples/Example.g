LoadPackage( "MotivesForBiArrangements" );


A := BlueMultizetaBiOS( [ 5 ] );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
IsBlueExact( A, A.Smin );

B := RedMultizetaBiOS( [ 5 ] );;
PrintArray( OrlikSolomonBicomplexDimensions( B, B.Smin ) );
IsRedExact( B, B.Smin );

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

m := [
      [ 1, 0 ],
      [ 0, 1 ],
      ];

m := Matroid( m, Q );




