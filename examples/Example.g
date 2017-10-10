LoadPackage( "MotivesForBiArrangements" );

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
    
    if flat = [ 4 ] then
        return false;
    fi;
    
    return true;
    
end;

A := OrlikSolomonBicomplex( m, chi );

M := [
	[ 1, 0, 0 ],
	[ 0, 1, 0 ],
	[ 0, 0, 1 ],
	[ 1, 1, 1 ]
	];

M := Matroid( M, Q );

A := OrlikSolomonBicomplex( M, chi );
