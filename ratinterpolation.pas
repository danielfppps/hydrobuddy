{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit ratinterpolation;
interface
uses Math, Sysutils, Ap, tsort;

procedure BuildFloaterHormannRationalInterpolant(X : TReal1DArray;
     N : AlglibInteger;
     D : AlglibInteger;
     var W : TReal1DArray);

implementation

(*************************************************************************
Rational barycentric interpolation without poles

The subroutine constructs the rational interpolating function without real
poles. It should be noted that the barycentric weights of the  interpolant
constructed are independent of the values of the given function.

Input parameters:
    X   -   interpolation nodes, array[0..N-1].
    N   -   number of nodes, N>0.
    D   -   order of the interpolation scheme, 0 <= D <= N-1.

Output parameters:
    W   -   array of the barycentric weights which  can  be  used  in  the
            BarycentricInterpolate subroutine. Array[0..N-1]

Note:
    this algorithm always succeeds and calculates the weights  with  close
    to machine precision.

  -- ALGLIB PROJECT --
     Copyright 17.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure BuildFloaterHormannRationalInterpolant(X : TReal1DArray;
     N : AlglibInteger;
     D : AlglibInteger;
     var W : TReal1DArray);
var
    S0 : Double;
    S : Double;
    V : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    Perm : TInteger1DArray;
    WTemp : TReal1DArray;
begin
    X := DynamicArrayCopy(X);
    Assert(N>0, 'BuildRationalInterpolantWithoutPoles: N<=0!');
    Assert((D>=0) and (D<=N), 'BuildRationalInterpolantWithoutPoles: incorrect D!');
    
    //
    // Prepare
    //
    SetLength(W, N-1+1);
    S0 := 1;
    K:=1;
    while K<=D do
    begin
        S0 := -S0;
        Inc(K);
    end;
    SetLength(Perm, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Perm[I] := I;
        Inc(I);
    end;
    TagSortFastI(X, Perm, N);
    
    //
    // Calculate Wk
    //
    K:=0;
    while K<=N-1 do
    begin
        
        //
        // Wk
        //
        S := 0;
        I:=Max(K-D, 0);
        while I<=Min(K, N-1-D) do
        begin
            V := 1;
            J:=I;
            while J<=I+D do
            begin
                if J<>K then
                begin
                    V := V/AbsReal(X[K]-X[J]);
                end;
                Inc(J);
            end;
            S := S+V;
            Inc(I);
        end;
        W[K] := S0*S;
        
        //
        // Next S0
        //
        S0 := -S0;
        Inc(K);
    end;
    
    //
    // Reorder W
    //
    SetLength(WTemp, N-1+1);
    APVMove(@WTemp[0], 0, N-1, @W[0], 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        W[Perm[I]] := WTemp[I];
        Inc(I);
    end;
end;


end.