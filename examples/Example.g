LoadPackage( "MotivesForBiArrangements" );

Q := HomalgFieldOfRationals( );

m := [
      [ 1, 0 ],
      [ 0, 1 ],
      [ 1, 1 ]
      ];

chi := function( flat )
    if flat = [ ] then
        Error( "\n" );
    fi;
    
    if flat = [ 3 ] then
        return false;
    fi;
    
    return true;
    
end;

m := Matroid( m, Q );

OS := OrlikSolomonBicomplex( m, chi );
