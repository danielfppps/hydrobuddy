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
unit hblas;
interface
uses Math, Sysutils, Ap;

procedure HermitianMatrixVectorMultiply(const A : TComplex2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TComplex1DArray;
     Alpha : Complex;
     var Y : TComplex1DArray);
procedure HermitianRank2Update(var A : TComplex2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TComplex1DArray;
     const Y : TComplex1DArray;
     var T : TComplex1DArray;
     Alpha : Complex);

implementation

procedure HermitianMatrixVectorMultiply(const A : TComplex2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TComplex1DArray;
     Alpha : Complex;
     var Y : TComplex1DArray);
var
    I : AlglibInteger;
    BA1 : AlglibInteger;
    BA2 : AlglibInteger;
    BY1 : AlglibInteger;
    BY2 : AlglibInteger;
    BX1 : AlglibInteger;
    BX2 : AlglibInteger;
    N : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
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
        Y[I-I1+1] := C_Mul(A[I,I],X[I-I1+1]);
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
            i1_ := (BA1) - (BY1);
            for i_ := BY1 to BY2 do
            begin
                Y[i_] := C_Add(Y[i_], C_Mul(V, Conj(A[I,i_+i1_])));
            end;
            
            //
            // Add U*x to the result
            //
            BX1 := I-I1+2;
            BX2 := N;
            BA1 := I+1;
            BA2 := I2;
            i1_ := (BA1)-(BX1);
            V := C_Complex(0.0);
            for i_ := BX1 to BX2 do
            begin
                V := C_Add(V,C_Mul(X[i_],A[I,i_+i1_]));
            end;
            Y[I-I1+1] := C_Add(Y[I-I1+1],V);
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
            i1_ := (BA1)-(BX1);
            V := C_Complex(0.0);
            for i_ := BX1 to BX2 do
            begin
                V := C_Add(V,C_Mul(X[i_],A[I,i_+i1_]));
            end;
            Y[I-I1+1] := C_Add(Y[I-I1+1],V);
            
            //
            // Add U*x to the result
            //
            V := X[I-I1+1];
            BY1 := 1;
            BY2 := I-I1;
            BA1 := I1;
            BA2 := I-1;
            i1_ := (BA1) - (BY1);
            for i_ := BY1 to BY2 do
            begin
                Y[i_] := C_Add(Y[i_], C_Mul(V, Conj(A[I,i_+i1_])));
            end;
            Inc(I);
        end;
    end;
    for i_ := 1 to N do
    begin
        Y[i_] := C_Mul(Alpha, Y[i_]);
    end;
end;


procedure HermitianRank2Update(var A : TComplex2DArray;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     const X : TComplex1DArray;
     const Y : TComplex1DArray;
     var T : TComplex1DArray;
     Alpha : Complex);
var
    I : AlglibInteger;
    TP1 : AlglibInteger;
    TP2 : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if IsUpper then
    begin
        I:=I1;
        while I<=I2 do
        begin
            TP1 := I+1-I1;
            TP2 := I2-I1+1;
            V := C_Mul(Alpha,X[I+1-I1]);
            for i_ := TP1 to TP2 do
            begin
                T[i_] := C_Mul(V, Conj(Y[i_]));
            end;
            V := C_Mul(Conj(Alpha),Y[I+1-I1]);
            for i_ := TP1 to TP2 do
            begin
                T[i_] := C_Add(T[i_], C_Mul(V, Conj(X[i_])));
            end;
            i1_ := (TP1) - (I);
            for i_ := I to I2 do
            begin
                A[I,i_] := C_Add(A[I,i_], T[i_+i1_]);
            end;
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
            V := C_Mul(Alpha,X[I+1-I1]);
            for i_ := TP1 to TP2 do
            begin
                T[i_] := C_Mul(V, Conj(Y[i_]));
            end;
            V := C_Mul(Conj(Alpha),Y[I+1-I1]);
            for i_ := TP1 to TP2 do
            begin
                T[i_] := C_Add(T[i_], C_Mul(V, Conj(X[i_])));
            end;
            i1_ := (TP1) - (I1);
            for i_ := I1 to I do
            begin
                A[I,i_] := C_Add(A[I,i_], T[i_+i1_]);
            end;
            Inc(I);
        end;
    end;
end;


end.