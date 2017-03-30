{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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
unit inverseupdate;
interface
uses Math, Sysutils, Ap;

procedure RMatrixInvUpdateSimple(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     UpdColumn : AlglibInteger;
     UpdVal : Double);
procedure RMatrixInvUpdateRow(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     const V : TReal1DArray);
procedure RMatrixInvUpdateColumn(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdColumn : AlglibInteger;
     const U : TReal1DArray);
procedure RMatrixInvUpdateUV(var InvA : TReal2DArray;
     N : AlglibInteger;
     const U : TReal1DArray;
     const V : TReal1DArray);
procedure ShermanMorrisonSimpleUpdate(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     UpdColumn : AlglibInteger;
     UpdVal : Double);
procedure ShermanMorrisonUpdateRow(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     const V : TReal1DArray);
procedure ShermanMorrisonUpdateColumn(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdColumn : AlglibInteger;
     const U : TReal1DArray);
procedure ShermanMorrisonUpdateUV(var InvA : TReal2DArray;
     N : AlglibInteger;
     const U : TReal1DArray;
     const V : TReal1DArray);

implementation

(*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm updates matrix A^-1 when adding a number to an element
of matrix A.

Input parameters:
    InvA    -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    UpdRow  -   row where the element to be updated is stored.
    UpdColumn - column where the element to be updated is stored.
    UpdVal  -   a number to be added to the element.


Output parameters:
    InvA    -   inverse of modified matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixInvUpdateSimple(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     UpdColumn : AlglibInteger;
     UpdVal : Double);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    Lambda : Double;
    VT : Double;
    i_ : AlglibInteger;
begin
    Assert((UpdRow>=0) and (UpdRow<N), 'RMatrixInvUpdateSimple: incorrect UpdRow!');
    Assert((UpdColumn>=0) and (UpdColumn<N), 'RMatrixInvUpdateSimple: incorrect UpdColumn!');
    SetLength(T1, N-1+1);
    SetLength(T2, N-1+1);
    
    //
    // T1 = InvA * U
    //
    for i_ := 0 to N-1 do
    begin
        T1[i_] := InvA[i_,UpdRow];
    end;
    
    //
    // T2 = v*InvA
    //
    APVMove(@T2[0], 0, N-1, @InvA[UpdColumn][0], 0, N-1);
    
    //
    // Lambda = v * InvA * U
    //
    Lambda := UpdVal*InvA[UpdColumn,UpdRow];
    
    //
    // InvA = InvA - correction
    //
    I:=0;
    while I<=N-1 do
    begin
        VT := UpdVal*T1[I];
        VT := VT/(1+Lambda);
        APVSub(@InvA[I][0], 0, N-1, @T2[0], 0, N-1, VT);
        Inc(I);
    end;
end;


(*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm updates matrix A^-1 when adding a vector to a row
of matrix A.

Input parameters:
    InvA    -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    UpdRow  -   the row of A whose vector V was added.
                0 <= Row <= N-1
    V       -   the vector to be added to a row.
                Array whose index ranges within [0..N-1].

Output parameters:
    InvA    -   inverse of modified matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixInvUpdateRow(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     const V : TReal1DArray);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    Lambda : Double;
    VT : Double;
    i_ : AlglibInteger;
begin
    SetLength(T1, N-1+1);
    SetLength(T2, N-1+1);
    
    //
    // T1 = InvA * U
    //
    for i_ := 0 to N-1 do
    begin
        T1[i_] := InvA[i_,UpdRow];
    end;
    
    //
    // T2 = v*InvA
    // Lambda = v * InvA * U
    //
    J:=0;
    while J<=N-1 do
    begin
        VT := 0.0;
        for i_ := 0 to N-1 do
        begin
            VT := VT + V[i_]*InvA[i_,J];
        end;
        T2[J] := VT;
        Inc(J);
    end;
    Lambda := T2[UpdRow];
    
    //
    // InvA = InvA - correction
    //
    I:=0;
    while I<=N-1 do
    begin
        VT := T1[I]/(1+Lambda);
        APVSub(@InvA[I][0], 0, N-1, @T2[0], 0, N-1, VT);
        Inc(I);
    end;
end;


(*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm updates matrix A^-1 when adding a vector to a column
of matrix A.

Input parameters:
    InvA        -   inverse of matrix A.
                    Array whose indexes range within [0..N-1, 0..N-1].
    N           -   size of matrix A.
    UpdColumn   -   the column of A whose vector U was added.
                    0 <= UpdColumn <= N-1
    U           -   the vector to be added to a column.
                    Array whose index ranges within [0..N-1].

Output parameters:
    InvA        -   inverse of modified matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixInvUpdateColumn(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdColumn : AlglibInteger;
     const U : TReal1DArray);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    Lambda : Double;
    VT : Double;
begin
    SetLength(T1, N-1+1);
    SetLength(T2, N-1+1);
    
    //
    // T1 = InvA * U
    // Lambda = v * InvA * U
    //
    I:=0;
    while I<=N-1 do
    begin
        VT := APVDotProduct(@InvA[I][0], 0, N-1, @U[0], 0, N-1);
        T1[I] := VT;
        Inc(I);
    end;
    Lambda := T1[UpdColumn];
    
    //
    // T2 = v*InvA
    //
    APVMove(@T2[0], 0, N-1, @InvA[UpdColumn][0], 0, N-1);
    
    //
    // InvA = InvA - correction
    //
    I:=0;
    while I<=N-1 do
    begin
        VT := T1[I]/(1+Lambda);
        APVSub(@InvA[I][0], 0, N-1, @T2[0], 0, N-1, VT);
        Inc(I);
    end;
end;


(*************************************************************************
Inverse matrix update by the Sherman-Morrison formula

The algorithm computes the inverse of matrix A+u*v’ by using the given matrix
A^-1 and the vectors u and v.

Input parameters:
    InvA    -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    U       -   the vector modifying the matrix.
                Array whose index ranges within [0..N-1].
    V       -   the vector modifying the matrix.
                Array whose index ranges within [0..N-1].

Output parameters:
    InvA - inverse of matrix A + u*v'.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixInvUpdateUV(var InvA : TReal2DArray;
     N : AlglibInteger;
     const U : TReal1DArray;
     const V : TReal1DArray);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    Lambda : Double;
    VT : Double;
    i_ : AlglibInteger;
begin
    SetLength(T1, N-1+1);
    SetLength(T2, N-1+1);
    
    //
    // T1 = InvA * U
    // Lambda = v * T1
    //
    I:=0;
    while I<=N-1 do
    begin
        VT := APVDotProduct(@InvA[I][0], 0, N-1, @U[0], 0, N-1);
        T1[I] := VT;
        Inc(I);
    end;
    Lambda := APVDotProduct(@V[0], 0, N-1, @T1[0], 0, N-1);
    
    //
    // T2 = v*InvA
    //
    J:=0;
    while J<=N-1 do
    begin
        VT := 0.0;
        for i_ := 0 to N-1 do
        begin
            VT := VT + V[i_]*InvA[i_,J];
        end;
        T2[J] := VT;
        Inc(J);
    end;
    
    //
    // InvA = InvA - correction
    //
    I:=0;
    while I<=N-1 do
    begin
        VT := T1[I]/(1+Lambda);
        APVSub(@InvA[I][0], 0, N-1, @T2[0], 0, N-1, VT);
        Inc(I);
    end;
end;


procedure ShermanMorrisonSimpleUpdate(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     UpdColumn : AlglibInteger;
     UpdVal : Double);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    Lambda : Double;
    VT : Double;
    i_ : AlglibInteger;
begin
    SetLength(T1, N+1);
    SetLength(T2, N+1);
    
    //
    // T1 = InvA * U
    //
    for i_ := 1 to N do
    begin
        T1[i_] := InvA[i_,UpdRow];
    end;
    
    //
    // T2 = v*InvA
    //
    APVMove(@T2[0], 1, N, @InvA[UpdColumn][0], 1, N);
    
    //
    // Lambda = v * InvA * U
    //
    Lambda := UpdVal*InvA[UpdColumn,UpdRow];
    
    //
    // InvA = InvA - correction
    //
    I:=1;
    while I<=N do
    begin
        VT := UpdVal*T1[I];
        VT := VT/(1+Lambda);
        APVSub(@InvA[I][0], 1, N, @T2[0], 1, N, VT);
        Inc(I);
    end;
end;


procedure ShermanMorrisonUpdateRow(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdRow : AlglibInteger;
     const V : TReal1DArray);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    Lambda : Double;
    VT : Double;
    i_ : AlglibInteger;
begin
    SetLength(T1, N+1);
    SetLength(T2, N+1);
    
    //
    // T1 = InvA * U
    //
    for i_ := 1 to N do
    begin
        T1[i_] := InvA[i_,UpdRow];
    end;
    
    //
    // T2 = v*InvA
    // Lambda = v * InvA * U
    //
    J:=1;
    while J<=N do
    begin
        VT := 0.0;
        for i_ := 1 to N do
        begin
            VT := VT + V[i_]*InvA[i_,J];
        end;
        T2[J] := VT;
        Inc(J);
    end;
    Lambda := T2[UpdRow];
    
    //
    // InvA = InvA - correction
    //
    I:=1;
    while I<=N do
    begin
        VT := T1[I]/(1+Lambda);
        APVSub(@InvA[I][0], 1, N, @T2[0], 1, N, VT);
        Inc(I);
    end;
end;


procedure ShermanMorrisonUpdateColumn(var InvA : TReal2DArray;
     N : AlglibInteger;
     UpdColumn : AlglibInteger;
     const U : TReal1DArray);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    Lambda : Double;
    VT : Double;
begin
    SetLength(T1, N+1);
    SetLength(T2, N+1);
    
    //
    // T1 = InvA * U
    // Lambda = v * InvA * U
    //
    I:=1;
    while I<=N do
    begin
        VT := APVDotProduct(@InvA[I][0], 1, N, @U[0], 1, N);
        T1[I] := VT;
        Inc(I);
    end;
    Lambda := T1[UpdColumn];
    
    //
    // T2 = v*InvA
    //
    APVMove(@T2[0], 1, N, @InvA[UpdColumn][0], 1, N);
    
    //
    // InvA = InvA - correction
    //
    I:=1;
    while I<=N do
    begin
        VT := T1[I]/(1+Lambda);
        APVSub(@InvA[I][0], 1, N, @T2[0], 1, N, VT);
        Inc(I);
    end;
end;


procedure ShermanMorrisonUpdateUV(var InvA : TReal2DArray;
     N : AlglibInteger;
     const U : TReal1DArray;
     const V : TReal1DArray);
var
    T1 : TReal1DArray;
    T2 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    Lambda : Double;
    VT : Double;
    i_ : AlglibInteger;
begin
    SetLength(T1, N+1);
    SetLength(T2, N+1);
    
    //
    // T1 = InvA * U
    // Lambda = v * T1
    //
    I:=1;
    while I<=N do
    begin
        VT := APVDotProduct(@InvA[I][0], 1, N, @U[0], 1, N);
        T1[I] := VT;
        Inc(I);
    end;
    Lambda := APVDotProduct(@V[0], 1, N, @T1[0], 1, N);
    
    //
    // T2 = v*InvA
    //
    J:=1;
    while J<=N do
    begin
        VT := 0.0;
        for i_ := 1 to N do
        begin
            VT := VT + V[i_]*InvA[i_,J];
        end;
        T2[J] := VT;
        Inc(J);
    end;
    
    //
    // InvA = InvA - correction
    //
    I:=1;
    while I<=N do
    begin
        VT := T1[I]/(1+Lambda);
        APVSub(@InvA[I][0], 1, N, @T2[0], 1, N, VT);
        Inc(I);
    end;
end;


end.