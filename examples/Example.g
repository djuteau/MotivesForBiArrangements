LoadPackage( "MotivesForBiArrangements" );

Q := HomalgFieldOfRationals( );

m := [
      [ 1, 0 ],
      [ 0, 1 ],
      [ 1, 1 ]
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

OS := OrlikSolomonBicomplex( m, chi );
