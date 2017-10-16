LoadPackage( "MotivesForBiArrangements" );


A := BlueMultizetaBiOS( [ 5 ] );;
PrintArray( OrlikSolomonBicomplexDimensions( A, A.Smin ) );
IsBlueExact( A, A.Smin );

B := RedMultizetaBiOS( [ 5 ] );;
PrintArray( OrlikSolomonBicomplexDimensions( B, B.Smin ) );
IsRedExact( B, B.Smin );

UnderlyingMatrix( OrlikSolomonBicomplexDifferential( A, [ 1, 2, 3, 4, 5, 6 ], 0, 0, 0, 1 ) ) =
Involution( UnderlyingMatrix( OrlikSolomonBicomplexDifferential( B, [ 1, 2, 3, 4, 5, 6 ], 1, 0, 0, 0 ) ) );

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




