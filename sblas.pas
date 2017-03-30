{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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
unit sblas;
interface
uses Math, Sysutils, Ap;

procedure SymmetricMatrixVectorMultiply(const A : TReal2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TReal1DArray;
     Alpha : Double;
     var Y : TReal1DArray);
procedure SymmetricRank2Update(var A : TReal2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TReal1DArray;
     const Y : TReal1DArray;
     var T : TReal1DArray;
     Alpha : Double);

implementation

procedure SymmetricMatrixVectorMultiply(const A : TReal2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TReal1DArray;
     Alpha : Double;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    BA1 : AlglibInteger;
    BA2 : AlglibInteger;
    BY1 : AlglibInteger;
    BY2 : AlglibInteger;
    BX1 : AlglibInteger;
    BX2 : AlglibInteger;
    N : AlglibInteger;
    V : Double;
begin
    N := I2-I1+1;
    if N<=0 then
    begin
        Exit;
    end;
    
    //
    // Let A = L + D + U, where
    //  L is strictly lower triangular (main diagonal is zero)
    //  D is diagonal
    //  U is strictly upper triangular (main diagonal is zero)
    //
    // A*x = L*x + D*x + U*x
    //
    // Calculate D*x first
    //
    I:=I1;
    while I<=I2 do
    begin
        Y[I-I1+1] := A[I,I]*X[I-I1+1];
        Inc(I);
    end;
    
    //
    // Add L*x + U*x
    //
    if IsUpper then
    begin
        I:=I1;
        while I<=I2-1 do
        begin
            
            //
            // Add L*x to the result
            //
            V := X[I-I1+1];
            BY1 := I-I1+2;
            BY2 := N;
            BA1 := I+1;
            BA2 := I2;
            APVAdd(@Y[0], BY1, BY2, @A[I][0], BA1, BA2, V);
            
            //
            // Add U*x to the result
            //
            BX1 := I-I1+2;
            BX2 := N;
            BA1 := I+1;
            BA2 := I2;
            V := APVDotProduct(@X[0], BX1, BX2, @A[I][0], BA1, BA2);
            Y[I-I1+1] := Y[I-I1+1]+V;
            Inc(I);
        end;
    end
    else
    begin
        I:=I1+1;
        while I<=I2 do
        begin
            
            //
            // Add L*x to the result
            //
            BX1 := 1;
            BX2 := I-I1;
            BA1 := I1;
            BA2 := I-1;
            V := APVDotProduct(@X[0], BX1, BX2, @A[I][0], BA1, BA2);
            Y[I-I1+1] := Y[I-I1+1]+V;
            
            //
            // Add U*x to the result
            //
            V := X[I-I1+1];
            BY1 := 1;
            BY2 := I-I1;
            BA1 := I1;
            BA2 := I-1;
            APVAdd(@Y[0], BY1, BY2, @A[I][0], BA1, BA2, V);
            Inc(I);
        end;
    end;
    APVMul(@Y[0], 1, N, Alpha);
end;


procedure SymmetricRank2Update(var A : TReal2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TReal1DArray;
     const Y : TReal1DArray;
     var T : TReal1DArray;
     Alpha : Double);
var
    I : AlglibInteger;
    TP1 : AlglibInteger;
    TP2 : AlglibInteger;
    V : Double;
begin
    if IsUpper then
    begin
        I:=I1;
        while I<=I2 do
        begin
            TP1 := I+1-I1;
            TP2 := I2-I1+1;
            V := X[I+1-I1];
            APVMove(@T[0], TP1, TP2, @Y[0], TP1, TP2, V);
            V := Y[I+1-I1];
            APVAdd(@T[0], TP1, TP2, @X[0], TP1, TP2, V);
            APVMul(@T[0], TP1, TP2, Alpha);
            APVAdd(@A[I][0], I, I2, @T[0], TP1, TP2);
            Inc(I);
        end;
    end
    else
    begin
        I:=I1;
        while I<=I2 do
        begin
            TP1 := 1;
            TP2 := I+1-I1;
            V := X[I+1-I1];
            APVMove(@T[0], TP1, TP2, @Y[0], TP1, TP2, V);
            V := Y[I+1-I1];
            APVAdd(@T[0], TP1, TP2, @X[0], TP1, TP2, V);
            APVMul(@T[0], TP1, TP2, Alpha);
            APVAdd(@A[I][0], I1, I, @T[0], TP1, TP2);
            Inc(I);
        end;
    end;
end;


end.