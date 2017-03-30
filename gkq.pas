{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2005-2009, Sergey Bochkanov (ALGLIB project).

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
unit gkq;
interface
uses Math, Sysutils, Ap, tsort, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, evd, gammafunc, gq;

procedure GKQGenerateRec(Alpha : TReal1DArray;
     Beta : TReal1DArray;
     Mu0 : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
procedure GKQGenerateGaussLegendre(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
procedure GKQGenerateGaussJacobi(N : AlglibInteger;
     Alpha : Double;
     Beta : Double;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
procedure GKQLegendreCalc(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
procedure GKQLegendreTbl(N : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray;
     var Eps : Double);

implementation

(*************************************************************************
Computation of nodes and weights of a Gauss-Kronrod quadrature formula

The algorithm generates the N-point Gauss-Kronrod quadrature formula  with
weight  function  given  by  coefficients  alpha  and beta of a recurrence
relation which generates a system of orthogonal polynomials:

    P-1(x)   =  0
    P0(x)    =  1
    Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

and zero moment Mu0

    Mu0 = integral(W(x)dx,a,b)


INPUT PARAMETERS:
    Alpha       –   alpha coefficients, array[0..floor(3*K/2)].
    Beta        –   beta coefficients,  array[0..ceil(3*K/2)].
                    Beta[0] is not used and may be arbitrary.
                    Beta[I]>0.
    Mu0         –   zeroth moment of the weight function.
    N           –   number of nodes of the Gauss-Kronrod quadrature formula,
                    N >= 3,
                    N =  2*K+1.

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -5    no real and positive Gauss-Kronrod formula can
                            be created for such a weight function  with  a
                            given number of nodes.
                    * -4    N is too large, task may be ill  conditioned -
                            x[i]=x[i+1] found.
                    * -3    internal eigenproblem solver hasn't converged
                    * -2    Beta[i]<=0
                    * -1    incorrect N was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
    WKronrod    -   array[0..N-1] - Kronrod weights
    WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
                    corresponding to extended Kronrod nodes).

  -- ALGLIB --
     Copyright 08.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GKQGenerateRec(Alpha : TReal1DArray;
     Beta : TReal1DArray;
     Mu0 : Double;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
var
    TA : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    T : TReal1DArray;
    S : TReal1DArray;
    WLen : AlglibInteger;
    WOffs : AlglibInteger;
    U : Double;
    M : AlglibInteger;
    L : AlglibInteger;
    K : AlglibInteger;
    XGTmp : TReal1DArray;
    WGTmp : TReal1DArray;
begin
    Alpha := DynamicArrayCopy(Alpha);
    Beta := DynamicArrayCopy(Beta);
    if (N mod 2<>1) or (N<3) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=Ceil(AP_Double(3*(n div 2))/2) do
    begin
        if AP_FP_Less_Eq(Beta[I],0) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // from external conventions about N/Beta/Mu0 to internal
    //
    N := N div 2;
    Beta[0] := Mu0;
    
    //
    // Calculate Gauss nodes/weights, save them for later processing
    //
    GQGenerateRec(Alpha, Beta, Mu0, N, Info, XGTmp, WGTmp);
    if Info<0 then
    begin
        Exit;
    end;
    
    //
    // Resize:
    // * A from 0..floor(3*n/2) to 0..2*n
    // * B from 0..ceil(3*n/2)  to 0..2*n
    //
    SetLength(TA, Floor(AP_Double(3*N)/2)+1);
    APVMove(@TA[0], 0, Floor(AP_Double(3*N)/2), @Alpha[0], 0, Floor(AP_Double(3*N)/2));
    SetLength(Alpha, 2*N+1);
    APVMove(@Alpha[0], 0, Floor(AP_Double(3*N)/2), @TA[0], 0, Floor(AP_Double(3*N)/2));
    I:=Floor(AP_Double(3*N)/2)+1;
    while I<=2*N do
    begin
        Alpha[I] := 0;
        Inc(I);
    end;
    SetLength(TA, Ceil(AP_Double(3*N)/2)+1);
    APVMove(@TA[0], 0, Ceil(AP_Double(3*N)/2), @Beta[0], 0, Ceil(AP_Double(3*N)/2));
    SetLength(Beta, 2*N+1);
    APVMove(@Beta[0], 0, Ceil(AP_Double(3*N)/2), @TA[0], 0, Ceil(AP_Double(3*N)/2));
    I:=Ceil(AP_Double(3*N)/2)+1;
    while I<=2*N do
    begin
        Beta[I] := 0;
        Inc(I);
    end;
    
    //
    // Initialize T, S
    //
    WLen := 2+N div 2;
    SetLength(T, WLen);
    SetLength(S, WLen);
    SetLength(TA, WLen);
    WOffs := 1;
    I:=0;
    while I<=WLen-1 do
    begin
        T[I] := 0;
        S[I] := 0;
        Inc(I);
    end;
    
    //
    // Algorithm from Dirk P. Laurie, "Calculation of Gauss-Kronrod quadrature rules", 1997.
    //
    T[WOffs+0] := Beta[N+1];
    M:=0;
    while M<=N-2 do
    begin
        U := 0;
        K:=(M+1) div 2;
        while K>=0 do
        begin
            L := M-K;
            U := U+(Alpha[K+N+1]-Alpha[L])*T[WOffs+K]+Beta[K+N+1]*S[WOffs+K-1]-Beta[L]*S[WOffs+K];
            S[WOffs+K] := U;
            Dec(K);
        end;
        APVMove(@TA[0], 0, WLen-1, @T[0], 0, WLen-1);
        APVMove(@T[0], 0, WLen-1, @S[0], 0, WLen-1);
        APVMove(@S[0], 0, WLen-1, @TA[0], 0, WLen-1);
        Inc(M);
    end;
    J:=N div 2;
    while J>=0 do
    begin
        S[WOffs+J] := S[WOffs+J-1];
        Dec(J);
    end;
    M:=N-1;
    while M<=2*N-3 do
    begin
        U := 0;
        K:=M+1-N;
        while K<=(M-1) div 2 do
        begin
            L := M-K;
            J := N-1-L;
            U := U-(Alpha[K+N+1]-Alpha[L])*T[WOffs+J]-Beta[K+N+1]*S[WOffs+J]+Beta[L]*S[WOffs+J+1];
            S[WOffs+J] := U;
            Inc(K);
        end;
        if M mod 2=0 then
        begin
            K := M div 2;
            Alpha[K+N+1] := Alpha[K]+(S[WOffs+J]-Beta[K+N+1]*S[WOffs+J+1])/T[WOffs+J+1];
        end
        else
        begin
            K := (M+1) div 2;
            Beta[K+N+1] := S[WOffs+J]/S[WOffs+J+1];
        end;
        APVMove(@TA[0], 0, WLen-1, @T[0], 0, WLen-1);
        APVMove(@T[0], 0, WLen-1, @S[0], 0, WLen-1);
        APVMove(@S[0], 0, WLen-1, @TA[0], 0, WLen-1);
        Inc(M);
    end;
    Alpha[2*N] := Alpha[N-1]-Beta[2*N]*S[WOffs+0]/T[WOffs+0];
    
    //
    // calculation of Kronrod nodes and weights, unpacking of Gauss weights
    //
    GQGenerateRec(Alpha, Beta, Mu0, 2*N+1, Info, X, WKronrod);
    if Info=-2 then
    begin
        Info := -5;
    end;
    if Info<0 then
    begin
        Exit;
    end;
    I:=0;
    while I<=2*N-1 do
    begin
        if AP_FP_Greater_Eq(X[I],X[I+1]) then
        begin
            Info := -4;
        end;
        Inc(I);
    end;
    if Info<0 then
    begin
        Exit;
    end;
    SetLength(WGauss, 2*N+1);
    I:=0;
    while I<=2*N do
    begin
        WGauss[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        WGauss[2*I+1] := WGTmp[I];
        Inc(I);
    end;
end;


(*************************************************************************
Returns   Gauss   and   Gauss-Kronrod   nodes/weights  for  Gauss-Legendre
quadrature with N points.

GKQLegendreCalc (calculation) or  GKQLegendreTbl  (precomputed  table)  is
used depending on machine precision and number of nodes.

INPUT PARAMETERS:
    N           -   number of Kronrod nodes, must be odd number, >=3.

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error   was   detected   when  calculating
                            weights/nodes.  N  is  too  large   to  obtain
                            weights/nodes  with  high   enough   accuracy.
                            Try  to   use   multiple   precision  version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes, ordered in
                    ascending order.
    WKronrod    -   array[0..N-1] - Kronrod weights
    WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
                    corresponding to extended Kronrod nodes).


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GKQGenerateGaussLegendre(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
var
    Eps : Double;
begin
    if AP_FP_Greater(MachineEpsilon,Double(1.0E-32)) and ((N=15) or (N=21) or (N=31) or (N=41) or (N=51) or (N=61)) then
    begin
        Info := 1;
        GKQLegendreTbl(N, X, WKronrod, WGauss, Eps);
    end
    else
    begin
        GKQLegendreCalc(N, Info, X, WKronrod, WGauss);
    end;
end;


(*************************************************************************
Returns   Gauss   and   Gauss-Kronrod   nodes/weights   for   Gauss-Jacobi
quadrature on [-1,1] with weight function

    W(x)=Power(1-x,Alpha)*Power(1+x,Beta).

INPUT PARAMETERS:
    N           -   number of Kronrod nodes, must be odd number, >=3.
    Alpha       -   power-law coefficient, Alpha>-1
    Beta        -   power-law coefficient, Beta>-1

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -5    no real and positive Gauss-Kronrod formula can
                            be created for such a weight function  with  a
                            given number of nodes.
                    * -4    an  error  was   detected   when   calculating
                            weights/nodes. Alpha or  Beta  are  too  close
                            to -1 to obtain weights/nodes with high enough
                            accuracy, or, may be, N is too large.  Try  to
                            use multiple precision version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N was passed
                    * +1    OK
                    * +2    OK, but quadrature rule have exterior  nodes,
                            x[0]<-1 or x[n-1]>+1
    X           -   array[0..N-1] - array of quadrature nodes, ordered in
                    ascending order.
    WKronrod    -   array[0..N-1] - Kronrod weights
    WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
                    corresponding to extended Kronrod nodes).


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GKQGenerateGaussJacobi(N : AlglibInteger;
     Alpha : Double;
     Beta : Double;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
var
    CLen : AlglibInteger;
    A : TReal1DArray;
    B : TReal1DArray;
    Alpha2 : Double;
    Beta2 : Double;
    APB : Double;
    T : Double;
    I : AlglibInteger;
    S : Double;
begin
    if (N mod 2<>1) or (N<3) then
    begin
        Info := -1;
        Exit;
    end;
    if AP_FP_Less_Eq(Alpha,-1) or AP_FP_Less_Eq(Beta,-1) then
    begin
        Info := -1;
        Exit;
    end;
    CLen := Ceil(AP_Double(3*(N div 2))/2)+1;
    SetLength(A, CLen);
    SetLength(B, CLen);
    I:=0;
    while I<=CLen-1 do
    begin
        A[I] := 0;
        Inc(I);
    end;
    APB := Alpha+Beta;
    A[0] := (Beta-Alpha)/(APB+2);
    T := (APB+1)*Ln(2)+LnGamma(Alpha+1, S)+LnGamma(Beta+1, S)-LnGamma(APB+2, S);
    if AP_FP_Greater(t,Ln(MaxRealNumber)) then
    begin
        Info := -4;
        Exit;
    end;
    B[0] := Exp(t);
    if CLen>1 then
    begin
        Alpha2 := AP_Sqr(Alpha);
        Beta2 := AP_Sqr(Beta);
        A[1] := (Beta2-Alpha2)/((APB+2)*(APB+4));
        B[1] := 4*(Alpha+1)*(Beta+1)/((APB+3)*AP_Sqr(APB+2));
        I:=2;
        while I<=CLen-1 do
        begin
            A[I] := Double(0.25)*(Beta2-Alpha2)/(I*I*(1+Double(0.5)*APB/I)*(1+Double(0.5)*(APB+2)/I));
            B[I] := Double(0.25)*(1+Alpha/I)*(1+Beta/I)*(1+APB/I)/((1+Double(0.5)*(APB+1)/I)*(1+Double(0.5)*(APB-1)/I)*AP_Sqr(1+Double(0.5)*APB/I));
            Inc(I);
        end;
    end;
    GKQGenerateRec(A, B, B[0], N, Info, X, WKronrod, WGauss);
    
    //
    // test basic properties to detect errors
    //
    if Info>0 then
    begin
        if AP_FP_Less(X[0],-1) or AP_FP_Greater(X[N-1],+1) then
        begin
            Info := 2;
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
Returns Gauss and Gauss-Kronrod nodes for quadrature with N points.

Reduction to tridiagonal eigenproblem is used.

INPUT PARAMETERS:
    N           -   number of Kronrod nodes, must be odd number, >=3.

OUTPUT PARAMETERS:
    Info        -   error code:
                    * -4    an  error   was   detected   when  calculating
                            weights/nodes.  N  is  too  large   to  obtain
                            weights/nodes  with  high   enough   accuracy.
                            Try  to   use   multiple   precision  version.
                    * -3    internal eigenproblem solver hasn't converged
                    * -1    incorrect N was passed
                    * +1    OK
    X           -   array[0..N-1] - array of quadrature nodes, ordered in
                    ascending order.
    WKronrod    -   array[0..N-1] - Kronrod weights
    WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
                    corresponding to extended Kronrod nodes).

  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GKQLegendreCalc(N : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray);
var
    Alpha : TReal1DArray;
    Beta : TReal1DArray;
    ALen : AlglibInteger;
    BLen : AlglibInteger;
    Mu0 : Double;
    K : AlglibInteger;
    I : AlglibInteger;
begin
    if (N mod 2<>1) or (N<3) then
    begin
        Info := -1;
        Exit;
    end;
    Mu0 := 2;
    ALen := Floor(AP_Double(3*(n div 2))/2)+1;
    BLen := Ceil(AP_Double(3*(n div 2))/2)+1;
    SetLength(Alpha, ALen);
    SetLength(Beta, BLen);
    K:=0;
    while K<=ALen-1 do
    begin
        Alpha[K] := 0;
        Inc(K);
    end;
    Beta[0] := 2;
    K:=1;
    while K<=BLen-1 do
    begin
        Beta[K] := 1/(4-1/AP_Sqr(K));
        Inc(K);
    end;
    GKQGenerateRec(Alpha, Beta, Mu0, N, Info, X, WKronrod, WGauss);
    
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
Returns Gauss and Gauss-Kronrod nodes for quadrature with N  points  using
pre-calculated table. Nodes/weights were  computed  with  accuracy  up  to
1.0E-32 (if MPFR version of ALGLIB is used). In standard double  precision
accuracy reduces to something about 2.0E-16 (depending  on your compiler's
handling of long floating point constants).

INPUT PARAMETERS:
    N           -   number of Kronrod nodes.
                    N can be 15, 21, 31, 41, 51, 61.

OUTPUT PARAMETERS:
    X           -   array[0..N-1] - array of quadrature nodes, ordered in
                    ascending order.
    WKronrod    -   array[0..N-1] - Kronrod weights
    WGauss      -   array[0..N-1] - Gauss weights (interleaved with zeros
                    corresponding to extended Kronrod nodes).


  -- ALGLIB --
     Copyright 12.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure GKQLegendreTbl(N : AlglibInteger;
     var X : TReal1DArray;
     var WKronrod : TReal1DArray;
     var WGauss : TReal1DArray;
     var Eps : Double);
var
    I : AlglibInteger;
    NG : AlglibInteger;
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
    Tmp : Double;
begin
    Assert((N=15) or (N=21) or (N=31) or (N=41) or (N=51) or (N=61), 'GKQNodesTbl: incorrect N!');
    SetLength(X, N-1+1);
    SetLength(WKronrod, N-1+1);
    SetLength(WGauss, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := 0;
        WKronrod[I] := 0;
        WGauss[I] := 0;
        Inc(I);
    end;
    Eps := Max(MachineEpsilon, Double(1.0E-32));
    if N=15 then
    begin
        NG := 4;
        WGauss[0] := Double(0.129484966168869693270611432679082);
        WGauss[1] := Double(0.279705391489276667901467771423780);
        WGauss[2] := Double(0.381830050505118944950369775488975);
        WGauss[3] := Double(0.417959183673469387755102040816327);
        X[0] := Double(0.991455371120812639206854697526329);
        X[1] := Double(0.949107912342758524526189684047851);
        X[2] := Double(0.864864423359769072789712788640926);
        X[3] := Double(0.741531185599394439863864773280788);
        X[4] := Double(0.586087235467691130294144838258730);
        X[5] := Double(0.405845151377397166906606412076961);
        X[6] := Double(0.207784955007898467600689403773245);
        X[7] := Double(0.000000000000000000000000000000000);
        WKronrod[0] := Double(0.022935322010529224963732008058970);
        WKronrod[1] := Double(0.063092092629978553290700663189204);
        WKronrod[2] := Double(0.104790010322250183839876322541518);
        WKronrod[3] := Double(0.140653259715525918745189590510238);
        WKronrod[4] := Double(0.169004726639267902826583426598550);
        WKronrod[5] := Double(0.190350578064785409913256402421014);
        WKronrod[6] := Double(0.204432940075298892414161999234649);
        WKronrod[7] := Double(0.209482141084727828012999174891714);
    end;
    if N=21 then
    begin
        NG := 5;
        WGauss[0] := Double(0.066671344308688137593568809893332);
        WGauss[1] := Double(0.149451349150580593145776339657697);
        WGauss[2] := Double(0.219086362515982043995534934228163);
        WGauss[3] := Double(0.269266719309996355091226921569469);
        WGauss[4] := Double(0.295524224714752870173892994651338);
        X[0] := Double(0.995657163025808080735527280689003);
        X[1] := Double(0.973906528517171720077964012084452);
        X[2] := Double(0.930157491355708226001207180059508);
        X[3] := Double(0.865063366688984510732096688423493);
        X[4] := Double(0.780817726586416897063717578345042);
        X[5] := Double(0.679409568299024406234327365114874);
        X[6] := Double(0.562757134668604683339000099272694);
        X[7] := Double(0.433395394129247190799265943165784);
        X[8] := Double(0.294392862701460198131126603103866);
        X[9] := Double(0.148874338981631210884826001129720);
        X[10] := Double(0.000000000000000000000000000000000);
        WKronrod[0] := Double(0.011694638867371874278064396062192);
        WKronrod[1] := Double(0.032558162307964727478818972459390);
        WKronrod[2] := Double(0.054755896574351996031381300244580);
        WKronrod[3] := Double(0.075039674810919952767043140916190);
        WKronrod[4] := Double(0.093125454583697605535065465083366);
        WKronrod[5] := Double(0.109387158802297641899210590325805);
        WKronrod[6] := Double(0.123491976262065851077958109831074);
        WKronrod[7] := Double(0.134709217311473325928054001771707);
        WKronrod[8] := Double(0.142775938577060080797094273138717);
        WKronrod[9] := Double(0.147739104901338491374841515972068);
        WKronrod[10] := Double(0.149445554002916905664936468389821);
    end;
    if N=31 then
    begin
        NG := 8;
        WGauss[0] := Double(0.030753241996117268354628393577204);
        WGauss[1] := Double(0.070366047488108124709267416450667);
        WGauss[2] := Double(0.107159220467171935011869546685869);
        WGauss[3] := Double(0.139570677926154314447804794511028);
        WGauss[4] := Double(0.166269205816993933553200860481209);
        WGauss[5] := Double(0.186161000015562211026800561866423);
        WGauss[6] := Double(0.198431485327111576456118326443839);
        WGauss[7] := Double(0.202578241925561272880620199967519);
        X[0] := Double(0.998002298693397060285172840152271);
        X[1] := Double(0.987992518020485428489565718586613);
        X[2] := Double(0.967739075679139134257347978784337);
        X[3] := Double(0.937273392400705904307758947710209);
        X[4] := Double(0.897264532344081900882509656454496);
        X[5] := Double(0.848206583410427216200648320774217);
        X[6] := Double(0.790418501442465932967649294817947);
        X[7] := Double(0.724417731360170047416186054613938);
        X[8] := Double(0.650996741297416970533735895313275);
        X[9] := Double(0.570972172608538847537226737253911);
        X[10] := Double(0.485081863640239680693655740232351);
        X[11] := Double(0.394151347077563369897207370981045);
        X[12] := Double(0.299180007153168812166780024266389);
        X[13] := Double(0.201194093997434522300628303394596);
        X[14] := Double(0.101142066918717499027074231447392);
        X[15] := Double(0.000000000000000000000000000000000);
        WKronrod[0] := Double(0.005377479872923348987792051430128);
        WKronrod[1] := Double(0.015007947329316122538374763075807);
        WKronrod[2] := Double(0.025460847326715320186874001019653);
        WKronrod[3] := Double(0.035346360791375846222037948478360);
        WKronrod[4] := Double(0.044589751324764876608227299373280);
        WKronrod[5] := Double(0.053481524690928087265343147239430);
        WKronrod[6] := Double(0.062009567800670640285139230960803);
        WKronrod[7] := Double(0.069854121318728258709520077099147);
        WKronrod[8] := Double(0.076849680757720378894432777482659);
        WKronrod[9] := Double(0.083080502823133021038289247286104);
        WKronrod[10] := Double(0.088564443056211770647275443693774);
        WKronrod[11] := Double(0.093126598170825321225486872747346);
        WKronrod[12] := Double(0.096642726983623678505179907627589);
        WKronrod[13] := Double(0.099173598721791959332393173484603);
        WKronrod[14] := Double(0.100769845523875595044946662617570);
        WKronrod[15] := Double(0.101330007014791549017374792767493);
    end;
    if N=41 then
    begin
        NG := 10;
        WGauss[0] := Double(0.017614007139152118311861962351853);
        WGauss[1] := Double(0.040601429800386941331039952274932);
        WGauss[2] := Double(0.062672048334109063569506535187042);
        WGauss[3] := Double(0.083276741576704748724758143222046);
        WGauss[4] := Double(0.101930119817240435036750135480350);
        WGauss[5] := Double(0.118194531961518417312377377711382);
        WGauss[6] := Double(0.131688638449176626898494499748163);
        WGauss[7] := Double(0.142096109318382051329298325067165);
        WGauss[8] := Double(0.149172986472603746787828737001969);
        WGauss[9] := Double(0.152753387130725850698084331955098);
        X[0] := Double(0.998859031588277663838315576545863);
        X[1] := Double(0.993128599185094924786122388471320);
        X[2] := Double(0.981507877450250259193342994720217);
        X[3] := Double(0.963971927277913791267666131197277);
        X[4] := Double(0.940822633831754753519982722212443);
        X[5] := Double(0.912234428251325905867752441203298);
        X[6] := Double(0.878276811252281976077442995113078);
        X[7] := Double(0.839116971822218823394529061701521);
        X[8] := Double(0.795041428837551198350638833272788);
        X[9] := Double(0.746331906460150792614305070355642);
        X[10] := Double(0.693237656334751384805490711845932);
        X[11] := Double(0.636053680726515025452836696226286);
        X[12] := Double(0.575140446819710315342946036586425);
        X[13] := Double(0.510867001950827098004364050955251);
        X[14] := Double(0.443593175238725103199992213492640);
        X[15] := Double(0.373706088715419560672548177024927);
        X[16] := Double(0.301627868114913004320555356858592);
        X[17] := Double(0.227785851141645078080496195368575);
        X[18] := Double(0.152605465240922675505220241022678);
        X[19] := Double(0.076526521133497333754640409398838);
        X[20] := Double(0.000000000000000000000000000000000);
        WKronrod[0] := Double(0.003073583718520531501218293246031);
        WKronrod[1] := Double(0.008600269855642942198661787950102);
        WKronrod[2] := Double(0.014626169256971252983787960308868);
        WKronrod[3] := Double(0.020388373461266523598010231432755);
        WKronrod[4] := Double(0.025882133604951158834505067096153);
        WKronrod[5] := Double(0.031287306777032798958543119323801);
        WKronrod[6] := Double(0.036600169758200798030557240707211);
        WKronrod[7] := Double(0.041668873327973686263788305936895);
        WKronrod[8] := Double(0.046434821867497674720231880926108);
        WKronrod[9] := Double(0.050944573923728691932707670050345);
        WKronrod[10] := Double(0.055195105348285994744832372419777);
        WKronrod[11] := Double(0.059111400880639572374967220648594);
        WKronrod[12] := Double(0.062653237554781168025870122174255);
        WKronrod[13] := Double(0.065834597133618422111563556969398);
        WKronrod[14] := Double(0.068648672928521619345623411885368);
        WKronrod[15] := Double(0.071054423553444068305790361723210);
        WKronrod[16] := Double(0.073030690332786667495189417658913);
        WKronrod[17] := Double(0.074582875400499188986581418362488);
        WKronrod[18] := Double(0.075704497684556674659542775376617);
        WKronrod[19] := Double(0.076377867672080736705502835038061);
        WKronrod[20] := Double(0.076600711917999656445049901530102);
    end;
    if N=51 then
    begin
        NG := 13;
        WGauss[0] := Double(0.011393798501026287947902964113235);
        WGauss[1] := Double(0.026354986615032137261901815295299);
        WGauss[2] := Double(0.040939156701306312655623487711646);
        WGauss[3] := Double(0.054904695975835191925936891540473);
        WGauss[4] := Double(0.068038333812356917207187185656708);
        WGauss[5] := Double(0.080140700335001018013234959669111);
        WGauss[6] := Double(0.091028261982963649811497220702892);
        WGauss[7] := Double(0.100535949067050644202206890392686);
        WGauss[8] := Double(0.108519624474263653116093957050117);
        WGauss[9] := Double(0.114858259145711648339325545869556);
        WGauss[10] := Double(0.119455763535784772228178126512901);
        WGauss[11] := Double(0.122242442990310041688959518945852);
        WGauss[12] := Double(0.123176053726715451203902873079050);
        X[0] := Double(0.999262104992609834193457486540341);
        X[1] := Double(0.995556969790498097908784946893902);
        X[2] := Double(0.988035794534077247637331014577406);
        X[3] := Double(0.976663921459517511498315386479594);
        X[4] := Double(0.961614986425842512418130033660167);
        X[5] := Double(0.942974571228974339414011169658471);
        X[6] := Double(0.920747115281701561746346084546331);
        X[7] := Double(0.894991997878275368851042006782805);
        X[8] := Double(0.865847065293275595448996969588340);
        X[9] := Double(0.833442628760834001421021108693570);
        X[10] := Double(0.797873797998500059410410904994307);
        X[11] := Double(0.759259263037357630577282865204361);
        X[12] := Double(0.717766406813084388186654079773298);
        X[13] := Double(0.673566368473468364485120633247622);
        X[14] := Double(0.626810099010317412788122681624518);
        X[15] := Double(0.577662930241222967723689841612654);
        X[16] := Double(0.526325284334719182599623778158010);
        X[17] := Double(0.473002731445714960522182115009192);
        X[18] := Double(0.417885382193037748851814394594572);
        X[19] := Double(0.361172305809387837735821730127641);
        X[20] := Double(0.303089538931107830167478909980339);
        X[21] := Double(0.243866883720988432045190362797452);
        X[22] := Double(0.183718939421048892015969888759528);
        X[23] := Double(0.122864692610710396387359818808037);
        X[24] := Double(0.061544483005685078886546392366797);
        X[25] := Double(0.000000000000000000000000000000000);
        WKronrod[0] := Double(0.001987383892330315926507851882843);
        WKronrod[1] := Double(0.005561932135356713758040236901066);
        WKronrod[2] := Double(0.009473973386174151607207710523655);
        WKronrod[3] := Double(0.013236229195571674813656405846976);
        WKronrod[4] := Double(0.016847817709128298231516667536336);
        WKronrod[5] := Double(0.020435371145882835456568292235939);
        WKronrod[6] := Double(0.024009945606953216220092489164881);
        WKronrod[7] := Double(0.027475317587851737802948455517811);
        WKronrod[8] := Double(0.030792300167387488891109020215229);
        WKronrod[9] := Double(0.034002130274329337836748795229551);
        WKronrod[10] := Double(0.037116271483415543560330625367620);
        WKronrod[11] := Double(0.040083825504032382074839284467076);
        WKronrod[12] := Double(0.042872845020170049476895792439495);
        WKronrod[13] := Double(0.045502913049921788909870584752660);
        WKronrod[14] := Double(0.047982537138836713906392255756915);
        WKronrod[15] := Double(0.050277679080715671963325259433440);
        WKronrod[16] := Double(0.052362885806407475864366712137873);
        WKronrod[17] := Double(0.054251129888545490144543370459876);
        WKronrod[18] := Double(0.055950811220412317308240686382747);
        WKronrod[19] := Double(0.057437116361567832853582693939506);
        WKronrod[20] := Double(0.058689680022394207961974175856788);
        WKronrod[21] := Double(0.059720340324174059979099291932562);
        WKronrod[22] := Double(0.060539455376045862945360267517565);
        WKronrod[23] := Double(0.061128509717053048305859030416293);
        WKronrod[24] := Double(0.061471189871425316661544131965264);
        WKronrod[25] := Double(0.061580818067832935078759824240055);
    end;
    if N=61 then
    begin
        NG := 15;
        WGauss[0] := Double(0.007968192496166605615465883474674);
        WGauss[1] := Double(0.018466468311090959142302131912047);
        WGauss[2] := Double(0.028784707883323369349719179611292);
        WGauss[3] := Double(0.038799192569627049596801936446348);
        WGauss[4] := Double(0.048402672830594052902938140422808);
        WGauss[5] := Double(0.057493156217619066481721689402056);
        WGauss[6] := Double(0.065974229882180495128128515115962);
        WGauss[7] := Double(0.073755974737705206268243850022191);
        WGauss[8] := Double(0.080755895229420215354694938460530);
        WGauss[9] := Double(0.086899787201082979802387530715126);
        WGauss[10] := Double(0.092122522237786128717632707087619);
        WGauss[11] := Double(0.096368737174644259639468626351810);
        WGauss[12] := Double(0.099593420586795267062780282103569);
        WGauss[13] := Double(0.101762389748405504596428952168554);
        WGauss[14] := Double(0.102852652893558840341285636705415);
        X[0] := Double(0.999484410050490637571325895705811);
        X[1] := Double(0.996893484074649540271630050918695);
        X[2] := Double(0.991630996870404594858628366109486);
        X[3] := Double(0.983668123279747209970032581605663);
        X[4] := Double(0.973116322501126268374693868423707);
        X[5] := Double(0.960021864968307512216871025581798);
        X[6] := Double(0.944374444748559979415831324037439);
        X[7] := Double(0.926200047429274325879324277080474);
        X[8] := Double(0.905573307699907798546522558925958);
        X[9] := Double(0.882560535792052681543116462530226);
        X[10] := Double(0.857205233546061098958658510658944);
        X[11] := Double(0.829565762382768397442898119732502);
        X[12] := Double(0.799727835821839083013668942322683);
        X[13] := Double(0.767777432104826194917977340974503);
        X[14] := Double(0.733790062453226804726171131369528);
        X[15] := Double(0.697850494793315796932292388026640);
        X[16] := Double(0.660061064126626961370053668149271);
        X[17] := Double(0.620526182989242861140477556431189);
        X[18] := Double(0.579345235826361691756024932172540);
        X[19] := Double(0.536624148142019899264169793311073);
        X[20] := Double(0.492480467861778574993693061207709);
        X[21] := Double(0.447033769538089176780609900322854);
        X[22] := Double(0.400401254830394392535476211542661);
        X[23] := Double(0.352704725530878113471037207089374);
        X[24] := Double(0.304073202273625077372677107199257);
        X[25] := Double(0.254636926167889846439805129817805);
        X[26] := Double(0.204525116682309891438957671002025);
        X[27] := Double(0.153869913608583546963794672743256);
        X[28] := Double(0.102806937966737030147096751318001);
        X[29] := Double(0.051471842555317695833025213166723);
        X[30] := Double(0.000000000000000000000000000000000);
        WKronrod[0] := Double(0.001389013698677007624551591226760);
        WKronrod[1] := Double(0.003890461127099884051267201844516);
        WKronrod[2] := Double(0.006630703915931292173319826369750);
        WKronrod[3] := Double(0.009273279659517763428441146892024);
        WKronrod[4] := Double(0.011823015253496341742232898853251);
        WKronrod[5] := Double(0.014369729507045804812451432443580);
        WKronrod[6] := Double(0.016920889189053272627572289420322);
        WKronrod[7] := Double(0.019414141193942381173408951050128);
        WKronrod[8] := Double(0.021828035821609192297167485738339);
        WKronrod[9] := Double(0.024191162078080601365686370725232);
        WKronrod[10] := Double(0.026509954882333101610601709335075);
        WKronrod[11] := Double(0.028754048765041292843978785354334);
        WKronrod[12] := Double(0.030907257562387762472884252943092);
        WKronrod[13] := Double(0.032981447057483726031814191016854);
        WKronrod[14] := Double(0.034979338028060024137499670731468);
        WKronrod[15] := Double(0.036882364651821229223911065617136);
        WKronrod[16] := Double(0.038678945624727592950348651532281);
        WKronrod[17] := Double(0.040374538951535959111995279752468);
        WKronrod[18] := Double(0.041969810215164246147147541285970);
        WKronrod[19] := Double(0.043452539701356069316831728117073);
        WKronrod[20] := Double(0.044814800133162663192355551616723);
        WKronrod[21] := Double(0.046059238271006988116271735559374);
        WKronrod[22] := Double(0.047185546569299153945261478181099);
        WKronrod[23] := Double(0.048185861757087129140779492298305);
        WKronrod[24] := Double(0.049055434555029778887528165367238);
        WKronrod[25] := Double(0.049795683427074206357811569379942);
        WKronrod[26] := Double(0.050405921402782346840893085653585);
        WKronrod[27] := Double(0.050881795898749606492297473049805);
        WKronrod[28] := Double(0.051221547849258772170656282604944);
        WKronrod[29] := Double(0.051426128537459025933862879215781);
        WKronrod[30] := Double(0.051494729429451567558340433647099);
    end;
    
    //
    // copy nodes
    //
    I:=N-1;
    while I>=N div 2 do
    begin
        X[I] := -X[N-1-I];
        Dec(I);
    end;
    
    //
    // copy Kronrod weights
    //
    I:=N-1;
    while I>=N div 2 do
    begin
        WKronrod[I] := WKronrod[N-1-I];
        Dec(I);
    end;
    
    //
    // copy Gauss weights
    //
    I:=NG-1;
    while I>=0 do
    begin
        WGauss[N-2-2*I] := WGauss[I];
        WGauss[1+2*I] := WGauss[I];
        Dec(I);
    end;
    I:=0;
    while I<=N div 2 do
    begin
        WGauss[2*I] := 0;
        Inc(I);
    end;
    
    //
    // reorder
    //
    TagSort(X, N, P1, P2);
    I:=0;
    while I<=N-1 do
    begin
        Tmp := WKronrod[I];
        WKronrod[I] := WKronrod[P2[I]];
        WKronrod[P2[I]] := Tmp;
        Tmp := WGauss[I];
        WGauss[I] := WGauss[P2[I]];
        WGauss[P2[I]] := Tmp;
        Inc(I);
    end;
end;


end.