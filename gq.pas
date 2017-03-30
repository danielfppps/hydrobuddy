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
unit gq;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, evd, gammafunc;

procedure GQGenerateRec(const Alpha : TReal1DArray;
     const Beta : TReal1DArray;
     Mu0 : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
procedure GQGenerateGaussLobattoRec(Alpha : TReal1DArray;
     Beta : TReal1DArray;
     Mu0 : Double;
     A : Double;
     B : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
procedure GQGenerateGaussRadauRec(Alpha : TReal1DArray;
     Beta : TReal1DArray;
     Mu0 : Double;
     A : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
procedure GQGenerateGaussLegendre(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
procedure GQGenerateGaussJacobi(N : AlglibInteger;
     Alpha : Double;
     Beta : Double;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
procedure GQGenerateGaussLaguerre(N : AlglibInteger;
     Alpha : Double;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
procedure GQGenerateGaussHermite(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);

implementation

(*************************************************************************
Computation of nodes and weights for a Gauss quadrature formula

The algorithm generates the N-point Gauss quadrature formula  with  weight
function given by coefficients alpha and beta  of  a  recurrence  relation
which generates a system of orthogonal polynomials:

P-1(x)   =  0
P0(x)    =  1
Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zeroth moment Mu0

Mu0 = integral(W(x)dx,a,b)

INPUT PARAMETERS:
    Alpha   –   array[0..N-1], alpha coefficients
    Beta    –   array[0..N-1], beta coefficients
                Zero-indexed element is not used and may be arbitrary.
                Beta[I]>0.
    Mu0     –   zeroth moment of the weight function.
    N       –   number of nodes of the quadrature formula, N>=1

OUTPUT PARAMETERS:
    Info    -   error code:
                * -3    internal eigenproblem solver hasn't converged
                * -2    Beta[i]<=0
                * -1    incorrect N was passed
                *  1    OK
    X       -   array[0..N-1] - array of quadrature nodes,
                in ascending order.
    W       -   array[0..N-1] - array of quadrature weights.

  -- ALGLIB --
     Copyright 2005-2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateRec(const Alpha : TReal1DArray;
     const Beta : TReal1DArray;
     Mu0 : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    I : AlglibInteger;
    D : TReal1DArray;
    E : TReal1DArray;
    Z : TReal2DArray;
begin
    if N<1 then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    
    //
    // Initialize
    //
    SetLength(D, N);
    SetLength(E, N);
    I:=1;
    while I<=N-1 do
    begin
        D[I-1] := Alpha[I-1];
        if AP_FP_Less_Eq(Beta[I],0) then
        begin
            Info := -2;
            Exit;
        end;
        E[I-1] := Sqrt(Beta[I]);
        Inc(I);
    end;
    D[N-1] := Alpha[N-1];
    
    //
    // EVD
    //
    if  not SMatrixTDEVD(D, E, N, 3, Z) then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // Generate
    //
    SetLength(X, N);
    SetLength(W, N);
    I:=1;
    while I<=N do
    begin
        X[I-1] := D[I-1];
        W[I-1] := Mu0*AP_Sqr(Z[0,I-1]);
        Inc(I);
    end;
end;


(*************************************************************************
Computation of nodes and weights for a Gauss-Lobatto quadrature formula

The algorithm generates the N-point Gauss-Lobatto quadrature formula  with
weight function given by coefficients alpha and beta of a recurrence which
generates a system of orthogonal polynomials.

P-1(x)   =  0
P0(x)    =  1
Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zeroth moment Mu0

Mu0 = integral(W(x)dx,a,b)

INPUT PARAMETERS:
    Alpha   –   array[0..N-2], alpha coefficients
    Beta    –   array[0..N-2], beta coefficients.
                Zero-indexed element is not used, may be arbitrary.
                Beta[I]>0
    Mu0     –   zeroth moment of the weighting function.
    A       –   left boundary of the integration interval.
    B       –   right boundary of the integration interval.
    N       –   number of nodes of the quadrature formula, N>=3
                (including the left and right boundary nodes).

OUTPUT PARAMETERS:
    Info    -   error code:
                * -3    internal eigenproblem solver hasn't converged
                * -2    Beta[i]<=0
                * -1    incorrect N was passed
                *  1    OK
    X       -   array[0..N-1] - array of quadrature nodes,
                in ascending order.
    W       -   array[0..N-1] - array of quadrature weights.

  -- ALGLIB --
     Copyright 2005-2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateGaussLobattoRec(Alpha : TReal1DArray;
     Beta : TReal1DArray;
     Mu0 : Double;
     A : Double;
     B : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    I : AlglibInteger;
    D : TReal1DArray;
    E : TReal1DArray;
    Z : TReal2DArray;
    PIM1A : Double;
    PIA : Double;
    PIM1B : Double;
    PIB : Double;
    T : Double;
    A11 : Double;
    A12 : Double;
    A21 : Double;
    A22 : Double;
    B1 : Double;
    B2 : Double;
    Alph : Double;
    Bet : Double;
begin
    Alpha := DynamicArrayCopy(Alpha);
    Beta := DynamicArrayCopy(Beta);
    if N<=2 then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    
    //
    // Initialize, D[1:N+1], E[1:N]
    //
    N := N-2;
    SetLength(D, N+2);
    SetLength(E, N+1);
    I:=1;
    while I<=N+1 do
    begin
        D[I-1] := Alpha[I-1];
        Inc(I);
    end;
    I:=1;
    while I<=N do
    begin
        if AP_FP_Less_Eq(Beta[I],0) then
        begin
            Info := -2;
            Exit;
        end;
        E[I-1] := Sqrt(Beta[I]);
        Inc(I);
    end;
    
    //
    // Caclulate Pn(a), Pn+1(a), Pn(b), Pn+1(b)
    //
    Beta[0] := 0;
    PIM1A := 0;
    PIA := 1;
    PIM1B := 0;
    PIB := 1;
    I:=1;
    while I<=N+1 do
    begin
        
        //
        // Pi(a)
        //
        T := (A-Alpha[I-1])*PIA-Beta[I-1]*PIM1A;
        PIM1A := PIA;
        PIA := T;
        
        //
        // Pi(b)
        //
        T := (B-Alpha[I-1])*PIB-Beta[I-1]*PIM1B;
        PIM1B := PIB;
        PIB := T;
        Inc(I);
    end;
    
    //
    // Calculate alpha'(n+1), beta'(n+1)
    //
    A11 := PIA;
    A12 := PIM1A;
    A21 := PIB;
    A22 := PIM1B;
    B1 := A*PIA;
    B2 := B*PIB;
    if AP_FP_Greater(AbsReal(A11),AbsReal(A21)) then
    begin
        A22 := A22-A12*A21/A11;
        B2 := B2-B1*A21/A11;
        Bet := B2/A22;
        Alph := (B1-Bet*A12)/A11;
    end
    else
    begin
        A12 := A12-A22*A11/A21;
        B1 := B1-B2*A11/A21;
        Bet := B1/A12;
        Alph := (B2-Bet*A22)/A21;
    end;
    if AP_FP_Less(Bet,0) then
    begin
        Info := -3;
        Exit;
    end;
    D[N+1] := Alph;
    E[N] := Sqrt(Bet);
    
    //
    // EVD
    //
    if  not SMatrixTDEVD(D, E, N+2, 3, Z) then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // Generate
    //
    SetLength(X, N+2);
    SetLength(W, N+2);
    I:=1;
    while I<=N+2 do
    begin
        X[I-1] := D[I-1];
        W[I-1] := Mu0*AP_Sqr(Z[0,I-1]);
        Inc(I);
    end;
end;


(*************************************************************************
Computation of nodes and weights for a Gauss-Radau quadrature formula

The algorithm generates the N-point Gauss-Radau  quadrature  formula  with
weight function given by the coefficients alpha and  beta  of a recurrence
which generates a system of orthogonal polynomials.

P-1(x)   =  0
P0(x)    =  1
Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zeroth moment Mu0

Mu0 = integral(W(x)dx,a,b)

INPUT PARAMETERS:
    Alpha   –   array[0..N-2], alpha coefficients.
    Beta    –   array[0..N-1], beta coefficients
                Zero-indexed element is not used.
                Beta[I]>0
    Mu0     –   zeroth moment of the weighting function.
    A       –   left boundary of the integration interval.
    N       –   number of nodes of the quadrature formula, N>=2
                (including the left boundary node).

OUTPUT PARAMETERS:
    Info    -   error code:
                * -3    internal eigenproblem solver hasn't converged
                * -2    Beta[i]<=0
                * -1    incorrect N was passed
                *  1    OK
    X       -   array[0..N-1] - array of quadrature nodes,
                in ascending order.
    W       -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 2005-2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateGaussRadauRec(Alpha : TReal1DArray;
     Beta : TReal1DArray;
     Mu0 : Double;
     A : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    I : AlglibInteger;
    D : TReal1DArray;
    E : TReal1DArray;
    Z : TReal2DArray;
    PolIM1 : Double;
    PolI : Double;
    T : Double;
begin
    Alpha := DynamicArrayCopy(Alpha);
    Beta := DynamicArrayCopy(Beta);
    if N<2 then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    
    //
    // Initialize, D[1:N], E[1:N]
    //
    N := N-1;
    SetLength(D, N+1);
    SetLength(E, N);
    I:=1;
    while I<=N do
    begin
        D[I-1] := Alpha[I-1];
        if AP_FP_Less_Eq(Beta[I],0) then
        begin
            Info := -2;
            Exit;
        end;
        E[I-1] := Sqrt(Beta[I]);
        Inc(I);
    end;
    
    //
    // Caclulate Pn(a), Pn-1(a), and D[N+1]
    //
    Beta[0] := 0;
    PolIM1 := 0;
    PolI := 1;
    I:=1;
    while I<=N do
    begin
        T := (A-Alpha[I-1])*PolI-Beta[I-1]*PolIM1;
        PolIM1 := PolI;
        PolI := T;
        Inc(I);
    end;
    D[N] := A-Beta[N]*PolIM1/PolI;
    
    //
    // EVD
    //
    if  not SMatrixTDEVD(D, E, N+1, 3, Z) then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // Generate
    //
    SetLength(X, N+1);
    SetLength(W, N+1);
    I:=1;
    while I<=N+1 do
    begin
        X[I-1] := D[I-1];
        W[I-1] := Mu0*AP_Sqr(Z[0,I-1]);
        Inc(I);
    end;
end;


(*************************************************************************
Returns nodes/weights for Gauss-Legendre quadrature on [-1,1] with N
nodes.

INPUT PARAMETERS:
    N           -   number of nodes, >=1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error   was   detected   when  calculating
                            weights/nodes.  N  is  too  large   to  obtain
                            weights/nodes  with  high   enough   accuracy.
                            Try  to   use   multiple   precision  version.
                    * -3    internal eigenproblem solver hasn't  converged
                    * -1    incorrect N was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateGaussLegendre(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    Alpha : TReal1DArray;
    Beta : TReal1DArray;
    I : AlglibInteger;
begin
    if N<1 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(Alpha, N);
    SetLength(Beta, N);
    I:=0;
    while I<=N-1 do
    begin
        Alpha[I] := 0;
        Inc(I);
    end;
    Beta[0] := 2;
    I:=1;
    while I<=N-1 do
    begin
        Beta[I] := 1/(4-1/AP_Sqr(I));
        Inc(I);
    end;
    GQGenerateRec(Alpha, Beta, Beta[0], N, Info, X, W);
    
    //
    // test basic properties to detect errors
    //
    if Info>0 then
    begin
        if AP_FP_Less(X[0],-1) or AP_FP_Greater(X[N-1],+1) then
        begin
            Info := -4;
        end;
        I:=0;
        while I<=N-2 do
        begin
            if AP_FP_Greater_Eq(X[I],X[I+1]) then
            begin
                Info := -4;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Returns  nodes/weights  for  Gauss-Jacobi quadrature on [-1,1] with weight
function W(x)=Power(1-x,Alpha)*Power(1+x,Beta).

INPUT PARAMETERS:
    N           -   number of nodes, >=1
    Alpha       -   power-law coefficient, Alpha>-1
    Beta        -   power-law coefficient, Beta>-1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes. Alpha or  Beta  are  too  close
                            to -1 to obtain weights/nodes with high enough
                            accuracy, or, may be, N is too large.  Try  to
                            use multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N/Alpha/Beta was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateGaussJacobi(N : AlglibInteger;
     Alpha : Double;
     Beta : Double;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    A : TReal1DArray;
    B : TReal1DArray;
    Alpha2 : Double;
    Beta2 : Double;
    APB : Double;
    T : Double;
    I : AlglibInteger;
    S : Double;
begin
    if (N<1) or AP_FP_Less_Eq(Alpha,-1) or AP_FP_Less_Eq(Beta,-1) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(A, N);
    SetLength(B, N);
    APB := Alpha+Beta;
    A[0] := (Beta-Alpha)/(APB+2);
    T := (APB+1)*Ln(2)+LnGamma(Alpha+1, S)+LnGamma(Beta+1, S)-LnGamma(APB+2, S);
    if AP_FP_Greater(t,Ln(MaxRealNumber)) then
    begin
        Info := -4;
        Exit;
    end;
    B[0] := Exp(t);
    if N>1 then
    begin
        Alpha2 := AP_Sqr(Alpha);
        Beta2 := AP_Sqr(Beta);
        A[1] := (Beta2-Alpha2)/((APB+2)*(APB+4));
        B[1] := 4*(Alpha+1)*(Beta+1)/((APB+3)*AP_Sqr(APB+2));
        I:=2;
        while I<=N-1 do
        begin
            A[I] := Double(0.25)*(Beta2-Alpha2)/(I*I*(1+Double(0.5)*APB/I)*(1+Double(0.5)*(APB+2)/I));
            B[I] := Double(0.25)*(1+Alpha/I)*(1+Beta/I)*(1+APB/I)/((1+Double(0.5)*(APB+1)/I)*(1+Double(0.5)*(APB-1)/I)*AP_Sqr(1+Double(0.5)*APB/I));
            Inc(I);
        end;
    end;
    GQGenerateRec(A, B, B[0], N, Info, X, W);
    
    //
    // test basic properties to detect errors
    //
    if Info>0 then
    begin
        if AP_FP_Less(X[0],-1) or AP_FP_Greater(X[N-1],+1) then
        begin
            Info := -4;
        end;
        I:=0;
        while I<=N-2 do
        begin
            if AP_FP_Greater_Eq(X[I],X[I+1]) then
            begin
                Info := -4;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Returns  nodes/weights  for  Gauss-Laguerre  quadrature  on  [0,+inf) with
weight function W(x)=Power(x,Alpha)*Exp(-x)

INPUT PARAMETERS:
    N           -   number of nodes, >=1
    Alpha       -   power-law coefficient, Alpha>-1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes. Alpha is too  close  to  -1  to
                            obtain weights/nodes with high enough accuracy
                            or, may  be,  N  is  too  large.  Try  to  use
                            multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N/Alpha was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateGaussLaguerre(N : AlglibInteger;
     Alpha : Double;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    A : TReal1DArray;
    B : TReal1DArray;
    T : Double;
    I : AlglibInteger;
    S : Double;
begin
    if (N<1) or AP_FP_Less_Eq(Alpha,-1) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(A, N);
    SetLength(B, N);
    A[0] := Alpha+1;
    T := LnGamma(Alpha+1, S);
    if AP_FP_Greater_Eq(T,Ln(MaxRealNumber)) then
    begin
        Info := -4;
        Exit;
    end;
    B[0] := Exp(T);
    if N>1 then
    begin
        I:=1;
        while I<=N-1 do
        begin
            A[I] := 2*I+Alpha+1;
            B[I] := I*(I+Alpha);
            Inc(I);
        end;
    end;
    GQGenerateRec(A, B, B[0], N, Info, X, W);
    
    //
    // test basic properties to detect errors
    //
    if Info>0 then
    begin
        if AP_FP_Less(X[0],0) then
        begin
            Info := -4;
        end;
        I:=0;
        while I<=N-2 do
        begin
            if AP_FP_Greater_Eq(X[I],X[I+1]) then
            begin
                Info := -4;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Returns  nodes/weights  for  Gauss-Hermite  quadrature on (-inf,+inf) with
weight function W(x)=Exp(-x*x)

INPUT PARAMETERS:
    N           -   number of nodes, >=1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes.  May be, N is too large. Try to
                            use multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N/Alpha was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    W           -   array[0..N-1] - array of quadrature weights.


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GQGenerateGaussHermite(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var W : TReal1DArray);
var
    A : TReal1DArray;
    B : TReal1DArray;
    I : AlglibInteger;
begin
    if N<1 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(A, N);
    SetLength(B, N);
    I:=0;
    while I<=N-1 do
    begin
        A[I] := 0;
        Inc(I);
    end;
    B[0] := Sqrt(4*ArcTan(1));
    if N>1 then
    begin
        I:=1;
        while I<=N-1 do
        begin
            B[I] := Double(0.5)*I;
            Inc(I);
        end;
    end;
    GQGenerateRec(A, B, B[0], N, Info, X, W);
    
    //
    // test basic properties to detect errors
    //
    if Info>0 then
    begin
        I:=0;
        while I<=N-2 do
        begin
            if AP_FP_Greater_Eq(X[I],X[I+1]) then
            begin
                Info := -4;
            end;
            Inc(I);
        end;
    end;
end;


end.