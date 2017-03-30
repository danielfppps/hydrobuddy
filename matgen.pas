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
unit matgen;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd;

procedure RMatrixRndOrthogonal(N : AlglibInteger; var A : TReal2DArray);
procedure RMatrixRndCond(N : AlglibInteger; C : Double; var A : TReal2DArray);
procedure CMatrixRndOrthogonal(N : AlglibInteger; var A : TComplex2DArray);
procedure CMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TComplex2DArray);
procedure SMatrixRndCond(N : AlglibInteger; C : Double; var A : TReal2DArray);
procedure SPDMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TReal2DArray);
procedure HMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TComplex2DArray);
procedure HPDMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TComplex2DArray);
procedure RMatrixRndOrthogonalFromTheRight(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
procedure RMatrixRndOrthogonalFromTheLeft(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
procedure CMatrixRndOrthogonalFromTheRight(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
procedure CMatrixRndOrthogonalFromTheLeft(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
procedure SMatrixRndMultiply(var A : TReal2DArray; N : AlglibInteger);
procedure HMatrixRndMultiply(var A : TComplex2DArray; N : AlglibInteger);

implementation

(*************************************************************************
Generation of a random uniformly distributed (Haar) orthogonal matrix

INPUT PARAMETERS:
    N   -   matrix size, N>=1
    
OUTPUT PARAMETERS:
    A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixRndOrthogonal(N : AlglibInteger; var A : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(N>=1, 'RMatrixRndOrthogonal: N<1!');
    SetLength(A, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                A[I,J] := 1;
            end
            else
            begin
                A[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    RMatrixRndOrthogonalFromTheRight(A, N, N);
end;


(*************************************************************************
Generation of random NxN matrix with given condition number and norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixRndCond(N : AlglibInteger; C : Double; var A : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : Double;
    L2 : Double;
begin
    Assert((N>=1) and AP_FP_Greater_Eq(C,1), 'RMatrixRndCond: N<1 or C<1!');
    SetLength(A, N-1+1, N-1+1);
    if N=1 then
    begin
        
        //
        // special case
        //
        A[0,0] := 2*RandomInteger(2)-1;
        Exit;
    end;
    L1 := 0;
    L2 := Ln(1/C);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    A[0,0] := Exp(L1);
    I:=1;
    while I<=N-2 do
    begin
        A[I,I] := Exp(RandomReal*(L2-L1)+L1);
        Inc(I);
    end;
    A[N-1,N-1] := Exp(L2);
    RMatrixRndOrthogonalFromTheLeft(A, N, N);
    RMatrixRndOrthogonalFromTheRight(A, N, N);
end;


(*************************************************************************
Generation of a random Haar distributed orthogonal complex matrix

INPUT PARAMETERS:
    N   -   matrix size, N>=1

OUTPUT PARAMETERS:
    A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixRndOrthogonal(N : AlglibInteger; var A : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(N>=1, 'CMatrixRndOrthogonal: N<1!');
    SetLength(A, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                A[I,J] := C_Complex(1);
            end
            else
            begin
                A[I,J] := C_Complex(0);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    CMatrixRndOrthogonalFromTheRight(A, N, N);
end;


(*************************************************************************
Generation of random NxN complex matrix with given condition number C and
norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : Double;
    L2 : Double;
    State : HQRNDState;
    V : Complex;
begin
    Assert((N>=1) and AP_FP_Greater_Eq(C,1), 'CMatrixRndCond: N<1 or C<1!');
    SetLength(A, N-1+1, N-1+1);
    if N=1 then
    begin
        
        //
        // special case
        //
        HQRNDRandomize(State);
        HQRNDUnit2(State, V.X, V.Y);
        A[0,0] := V;
        Exit;
    end;
    L1 := 0;
    L2 := Ln(1/C);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    A[0,0] := C_Complex(Exp(L1));
    I:=1;
    while I<=N-2 do
    begin
        A[I,I] := C_Complex(Exp(RandomReal*(L2-L1)+L1));
        Inc(I);
    end;
    A[N-1,N-1] := C_Complex(Exp(L2));
    CMatrixRndOrthogonalFromTheLeft(A, N, N);
    CMatrixRndOrthogonalFromTheRight(A, N, N);
end;


(*************************************************************************
Generation of random NxN symmetric matrix with given condition number  and
norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure SMatrixRndCond(N : AlglibInteger; C : Double; var A : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : Double;
    L2 : Double;
begin
    Assert((N>=1) and AP_FP_Greater_Eq(C,1), 'SMatrixRndCond: N<1 or C<1!');
    SetLength(A, N-1+1, N-1+1);
    if N=1 then
    begin
        
        //
        // special case
        //
        A[0,0] := 2*RandomInteger(2)-1;
        Exit;
    end;
    
    //
    // Prepare matrix
    //
    L1 := 0;
    L2 := Ln(1/C);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    A[0,0] := Exp(L1);
    I:=1;
    while I<=N-2 do
    begin
        A[I,I] := (2*RandomInteger(2)-1)*Exp(RandomReal*(L2-L1)+L1);
        Inc(I);
    end;
    A[N-1,N-1] := Exp(L2);
    
    //
    // Multiply
    //
    SMatrixRndMultiply(A, N);
end;


(*************************************************************************
Generation of random NxN symmetric positive definite matrix with given
condition number and norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random SPD matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure SPDMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : Double;
    L2 : Double;
begin
    
    //
    // Special cases
    //
    if (N<=0) or AP_FP_Less(C,1) then
    begin
        Exit;
    end;
    SetLength(A, N-1+1, N-1+1);
    if N=1 then
    begin
        A[0,0] := 1;
        Exit;
    end;
    
    //
    // Prepare matrix
    //
    L1 := 0;
    L2 := Ln(1/C);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    A[0,0] := Exp(L1);
    I:=1;
    while I<=N-2 do
    begin
        A[I,I] := Exp(RandomReal*(L2-L1)+L1);
        Inc(I);
    end;
    A[N-1,N-1] := Exp(L2);
    
    //
    // Multiply
    //
    SMatrixRndMultiply(A, N);
end;


(*************************************************************************
Generation of random NxN Hermitian matrix with given condition number  and
norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure HMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : Double;
    L2 : Double;
begin
    Assert((N>=1) and AP_FP_Greater_Eq(C,1), 'HMatrixRndCond: N<1 or C<1!');
    SetLength(A, N-1+1, N-1+1);
    if N=1 then
    begin
        
        //
        // special case
        //
        A[0,0] := C_Complex(2*RandomInteger(2)-1);
        Exit;
    end;
    
    //
    // Prepare matrix
    //
    L1 := 0;
    L2 := Ln(1/C);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    A[0,0] := C_Complex(Exp(L1));
    I:=1;
    while I<=N-2 do
    begin
        A[I,I] := C_Complex((2*RandomInteger(2)-1)*Exp(RandomReal*(L2-L1)+L1));
        Inc(I);
    end;
    A[N-1,N-1] := C_Complex(Exp(L2));
    
    //
    // Multiply
    //
    HMatrixRndMultiply(A, N);
    
    //
    // post-process to ensure that matrix diagonal is real
    //
    I:=0;
    while I<=N-1 do
    begin
        A[I,I].Y := 0;
        Inc(I);
    end;
end;


(*************************************************************************
Generation of random NxN Hermitian positive definite matrix with given
condition number and norm2(A)=1

INPUT PARAMETERS:
    N   -   matrix size
    C   -   condition number (in 2-norm)

OUTPUT PARAMETERS:
    A   -   random HPD matrix with norm2(A)=1 and cond(A)=C

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure HPDMatrixRndCond(N : AlglibInteger;
     C : Double;
     var A : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L1 : Double;
    L2 : Double;
begin
    
    //
    // Special cases
    //
    if (N<=0) or AP_FP_Less(C,1) then
    begin
        Exit;
    end;
    SetLength(A, N-1+1, N-1+1);
    if N=1 then
    begin
        A[0,0] := C_Complex(1);
        Exit;
    end;
    
    //
    // Prepare matrix
    //
    L1 := 0;
    L2 := Ln(1/C);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    A[0,0] := C_Complex(Exp(L1));
    I:=1;
    while I<=N-2 do
    begin
        A[I,I] := C_Complex(Exp(RandomReal*(L2-L1)+L1));
        Inc(I);
    end;
    A[N-1,N-1] := C_Complex(Exp(L2));
    
    //
    // Multiply
    //
    HMatrixRndMultiply(A, N);
    
    //
    // post-process to ensure that matrix diagonal is real
    //
    I:=0;
    while I<=N-1 do
    begin
        A[I,I].Y := 0;
        Inc(I);
    end;
end;


(*************************************************************************
Multiplication of MxN matrix by NxN random Haar distributed orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixRndOrthogonalFromTheRight(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
var
    Tau : Double;
    Lambda : Double;
    S : AlglibInteger;
    I : AlglibInteger;
    U1 : Double;
    U2 : Double;
    W : TReal1DArray;
    V : TReal1DArray;
    State : HQRNDState;
    i_ : AlglibInteger;
begin
    Assert((N>=1) and (M>=1), 'RMatrixRndOrthogonalFromTheRight: N<1 or M<1!');
    if N=1 then
    begin
        
        //
        // Special case
        //
        Tau := 2*RandomInteger(2)-1;
        I:=0;
        while I<=M-1 do
        begin
            A[I,0] := A[I,0]*Tau;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // General case.
    // First pass.
    //
    SetLength(W, M-1+1);
    SetLength(V, N+1);
    HQRNDRandomize(State);
    S:=2;
    while S<=N do
    begin
        
        //
        // Prepare random normal v
        //
        repeat
            I := 1;
            while I<=S do
            begin
                HQRNDNormal2(State, U1, U2);
                V[I] := U1;
                if I+1<=S then
                begin
                    V[I+1] := U2;
                end;
                I := I+2;
            end;
            Lambda := APVDotProduct(@V[0], 1, S, @V[0], 1, S);
        until AP_FP_Neq(Lambda,0);
        
        //
        // Prepare and apply reflection
        //
        GenerateReflection(V, S, Tau);
        V[1] := 1;
        ApplyReflectionFromTheRight(A, Tau, V, 0, M-1, N-S, N-1, W);
        Inc(S);
    end;
    
    //
    // Second pass.
    //
    I:=0;
    while I<=N-1 do
    begin
        Tau := 2*RandomInteger(2)-1;
        for i_ := 0 to M-1 do
        begin
            A[i_,I] := Tau*A[i_,I];
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Multiplication of MxN matrix by MxM random Haar distributed orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   Q*A, where Q is random MxM orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RMatrixRndOrthogonalFromTheLeft(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
var
    Tau : Double;
    Lambda : Double;
    S : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    U1 : Double;
    U2 : Double;
    W : TReal1DArray;
    V : TReal1DArray;
    State : HQRNDState;
begin
    Assert((N>=1) and (M>=1), 'RMatrixRndOrthogonalFromTheRight: N<1 or M<1!');
    if M=1 then
    begin
        
        //
        // special case
        //
        Tau := 2*RandomInteger(2)-1;
        J:=0;
        while J<=N-1 do
        begin
            A[0,J] := A[0,J]*Tau;
            Inc(J);
        end;
        Exit;
    end;
    
    //
    // General case.
    // First pass.
    //
    SetLength(W, N-1+1);
    SetLength(V, M+1);
    HQRNDRandomize(State);
    S:=2;
    while S<=M do
    begin
        
        //
        // Prepare random normal v
        //
        repeat
            I := 1;
            while I<=S do
            begin
                HQRNDNormal2(State, U1, U2);
                V[I] := U1;
                if I+1<=S then
                begin
                    V[I+1] := U2;
                end;
                I := I+2;
            end;
            Lambda := APVDotProduct(@V[0], 1, S, @V[0], 1, S);
        until AP_FP_Neq(Lambda,0);
        
        //
        // Prepare and apply reflection
        //
        GenerateReflection(V, S, Tau);
        V[1] := 1;
        ApplyReflectionFromTheLeft(A, Tau, V, M-S, M-1, 0, N-1, W);
        Inc(S);
    end;
    
    //
    // Second pass.
    //
    I:=0;
    while I<=M-1 do
    begin
        Tau := 2*RandomInteger(2)-1;
        APVMul(@A[I][0], 0, N-1, Tau);
        Inc(I);
    end;
end;


(*************************************************************************
Multiplication of MxN complex matrix by NxN random Haar distributed
complex orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixRndOrthogonalFromTheRight(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
var
    Lambda : Complex;
    Tau : Complex;
    S : AlglibInteger;
    I : AlglibInteger;
    W : TComplex1DArray;
    V : TComplex1DArray;
    State : HQRNDState;
    i_ : AlglibInteger;
begin
    Assert((N>=1) and (M>=1), 'CMatrixRndOrthogonalFromTheRight: N<1 or M<1!');
    if N=1 then
    begin
        
        //
        // Special case
        //
        HQRNDRandomize(State);
        HQRNDUnit2(State, Tau.X, Tau.Y);
        I:=0;
        while I<=M-1 do
        begin
            A[I,0] := C_Mul(A[I,0],Tau);
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // General case.
    // First pass.
    //
    SetLength(W, M-1+1);
    SetLength(V, N+1);
    HQRNDRandomize(State);
    S:=2;
    while S<=N do
    begin
        
        //
        // Prepare random normal v
        //
        repeat
            I:=1;
            while I<=S do
            begin
                HQRNDNormal2(State, Tau.X, Tau.Y);
                V[I] := Tau;
                Inc(I);
            end;
            Lambda := C_Complex(0.0);
            for i_ := 1 to S do
            begin
                Lambda := C_Add(Lambda,C_Mul(V[i_],Conj(V[i_])));
            end;
        until C_NotEqualR(Lambda,0);
        
        //
        // Prepare and apply reflection
        //
        ComplexGenerateReflection(V, S, Tau);
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheRight(A, Tau, V, 0, M-1, N-S, N-1, W);
        Inc(S);
    end;
    
    //
    // Second pass.
    //
    I:=0;
    while I<=N-1 do
    begin
        HQRNDUnit2(State, Tau.X, Tau.Y);
        for i_ := 0 to M-1 do
        begin
            A[i_,I] := C_Mul(Tau, A[i_,I]);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Multiplication of MxN complex matrix by MxM random Haar distributed
complex orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..M-1, 0..N-1]
    M, N-   matrix size

OUTPUT PARAMETERS:
    A   -   Q*A, where Q is random MxM orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure CMatrixRndOrthogonalFromTheLeft(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
var
    Tau : Complex;
    Lambda : Complex;
    S : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    W : TComplex1DArray;
    V : TComplex1DArray;
    State : HQRNDState;
    i_ : AlglibInteger;
begin
    Assert((N>=1) and (M>=1), 'CMatrixRndOrthogonalFromTheRight: N<1 or M<1!');
    if M=1 then
    begin
        
        //
        // special case
        //
        HQRNDRandomize(State);
        HQRNDUnit2(State, Tau.X, Tau.Y);
        J:=0;
        while J<=N-1 do
        begin
            A[0,J] := C_Mul(A[0,J],Tau);
            Inc(J);
        end;
        Exit;
    end;
    
    //
    // General case.
    // First pass.
    //
    SetLength(W, N-1+1);
    SetLength(V, M+1);
    HQRNDRandomize(State);
    S:=2;
    while S<=M do
    begin
        
        //
        // Prepare random normal v
        //
        repeat
            I:=1;
            while I<=S do
            begin
                HQRNDNormal2(State, Tau.X, Tau.Y);
                V[I] := Tau;
                Inc(I);
            end;
            Lambda := C_Complex(0.0);
            for i_ := 1 to S do
            begin
                Lambda := C_Add(Lambda,C_Mul(V[i_],Conj(V[i_])));
            end;
        until C_NotEqualR(Lambda,0);
        
        //
        // Prepare and apply reflection
        //
        ComplexGenerateReflection(V, S, Tau);
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheLeft(A, Tau, V, M-S, M-1, 0, N-1, W);
        Inc(S);
    end;
    
    //
    // Second pass.
    //
    I:=0;
    while I<=M-1 do
    begin
        HQRNDUnit2(State, Tau.X, Tau.Y);
        for i_ := 0 to N-1 do
        begin
            A[I,i_] := C_Mul(Tau, A[I,i_]);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Symmetric multiplication of NxN matrix by random Haar distributed
orthogonal  matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..N-1, 0..N-1]
    N   -   matrix size

OUTPUT PARAMETERS:
    A   -   Q'*A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure SMatrixRndMultiply(var A : TReal2DArray; N : AlglibInteger);
var
    Tau : Double;
    Lambda : Double;
    S : AlglibInteger;
    I : AlglibInteger;
    U1 : Double;
    U2 : Double;
    W : TReal1DArray;
    V : TReal1DArray;
    State : HQRNDState;
    i_ : AlglibInteger;
begin
    
    //
    // General case.
    //
    SetLength(W, N-1+1);
    SetLength(V, N+1);
    HQRNDRandomize(State);
    S:=2;
    while S<=N do
    begin
        
        //
        // Prepare random normal v
        //
        repeat
            I := 1;
            while I<=S do
            begin
                HQRNDNormal2(State, U1, U2);
                V[I] := U1;
                if I+1<=S then
                begin
                    V[I+1] := U2;
                end;
                I := I+2;
            end;
            Lambda := APVDotProduct(@V[0], 1, S, @V[0], 1, S);
        until AP_FP_Neq(Lambda,0);
        
        //
        // Prepare and apply reflection
        //
        GenerateReflection(V, S, Tau);
        V[1] := 1;
        ApplyReflectionFromTheRight(A, Tau, V, 0, N-1, N-S, N-1, W);
        ApplyReflectionFromTheLeft(A, Tau, V, N-S, N-1, 0, N-1, W);
        Inc(S);
    end;
    
    //
    // Second pass.
    //
    I:=0;
    while I<=N-1 do
    begin
        Tau := 2*RandomInteger(2)-1;
        for i_ := 0 to N-1 do
        begin
            A[i_,I] := Tau*A[i_,I];
        end;
        APVMul(@A[I][0], 0, N-1, Tau);
        Inc(I);
    end;
end;


(*************************************************************************
Hermitian multiplication of NxN matrix by random Haar distributed
complex orthogonal matrix

INPUT PARAMETERS:
    A   -   matrix, array[0..N-1, 0..N-1]
    N   -   matrix size

OUTPUT PARAMETERS:
    A   -   Q^H*A*Q, where Q is random NxN orthogonal matrix

  -- ALGLIB routine --
     04.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure HMatrixRndMultiply(var A : TComplex2DArray; N : AlglibInteger);
var
    Tau : Complex;
    Lambda : Complex;
    S : AlglibInteger;
    I : AlglibInteger;
    W : TComplex1DArray;
    V : TComplex1DArray;
    State : HQRNDState;
    i_ : AlglibInteger;
begin
    
    //
    // General case.
    //
    SetLength(W, N-1+1);
    SetLength(V, N+1);
    HQRNDRandomize(State);
    S:=2;
    while S<=N do
    begin
        
        //
        // Prepare random normal v
        //
        repeat
            I:=1;
            while I<=S do
            begin
                HQRNDNormal2(State, Tau.X, Tau.Y);
                V[I] := Tau;
                Inc(I);
            end;
            Lambda := C_Complex(0.0);
            for i_ := 1 to S do
            begin
                Lambda := C_Add(Lambda,C_Mul(V[i_],Conj(V[i_])));
            end;
        until C_NotEqualR(Lambda,0);
        
        //
        // Prepare and apply reflection
        //
        ComplexGenerateReflection(V, S, Tau);
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheRight(A, Tau, V, 0, N-1, N-S, N-1, W);
        ComplexApplyReflectionFromTheLeft(A, Conj(Tau), V, N-S, N-1, 0, N-1, W);
        Inc(S);
    end;
    
    //
    // Second pass.
    //
    I:=0;
    while I<=N-1 do
    begin
        HQRNDUnit2(State, Tau.X, Tau.Y);
        for i_ := 0 to N-1 do
        begin
            A[i_,I] := C_Mul(Tau, A[i_,I]);
        end;
        Tau := Conj(Tau);
        for i_ := 0 to N-1 do
        begin
            A[I,i_] := C_Mul(Tau, A[I,i_]);
        end;
        Inc(I);
    end;
end;


end.