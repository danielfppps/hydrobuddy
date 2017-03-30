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
unit spline3;
interface
uses Math, Sysutils, Ap;

procedure BuildLinearSpline(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : TReal1DArray);
procedure BuildCubicSpline(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     BoundLType : AlglibInteger;
     BoundL : Double;
     BoundRType : AlglibInteger;
     BoundR : Double;
     var C : TReal1DArray);
procedure BuildHermiteSpline(X : TReal1DArray;
     Y : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var C : TReal1DArray);
procedure BuildAkimaSpline(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : TReal1DArray);
function SplineInterpolation(const C : TReal1DArray; X : Double):Double;
procedure SplineDifferentiation(const C : TReal1DArray;
     X : Double;
     var S : Double;
     var DS : Double;
     var D2S : Double);
procedure SplineCopy(const C : TReal1DArray; var CC : TReal1DArray);
procedure SplineUnpack(const C : TReal1DArray;
     var N : AlglibInteger;
     var Tbl : TReal2DArray);
procedure SplineLinTransX(var C : TReal1DArray; A : Double; B : Double);
procedure SplineLinTransY(var C : TReal1DArray; A : Double; B : Double);
function SplineIntegration(const C : TReal1DArray; X : Double):Double;
procedure Spline3BuildTable(N : AlglibInteger;
     const DiffN : AlglibInteger;
     x : TReal1DArray;
     y : TReal1DArray;
     const BoundL : Double;
     const BoundR : Double;
     var ctbl : TReal2DArray);
function Spline3Interpolate(N : AlglibInteger;
     const c : TReal2DArray;
     const X : Double):Double;

implementation

procedure HeapSortPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger);forward;
procedure HeapSortDPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     var D : TReal1DArray;
     N : AlglibInteger);forward;
procedure SolveTridiagonal(A : TReal1DArray;
     B : TReal1DArray;
     C : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var X : TReal1DArray);forward;
function DiffThreePoint(T : Double;
     X0 : Double;
     F0 : Double;
     X1 : Double;
     F1 : Double;
     X2 : Double;
     F2 : Double):Double;forward;


procedure BuildLinearSpline(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : TReal1DArray);
var
    I : AlglibInteger;
    TblSize : AlglibInteger;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>=2, 'BuildLinearSpline: N<2!');
    
    //
    // Sort points
    //
    HeapSortPoints(X, Y, N);
    
    //
    // Fill C:
    //  C[0]            -   length(C)
    //  C[1]            -   type(C):
    //                      3 - general cubic spline
    //  C[2]            -   N
    //  C[3]...C[3+N-1] -   x[i], i = 0...N-1
    //  C[3+N]...C[3+N+(N-1)*4-1] - coefficients table
    //
    TblSize := 3+N+(N-1)*4;
    SetLength(C, TblSize-1+1);
    C[0] := TblSize;
    C[1] := 3;
    C[2] := N;
    I:=0;
    while I<=N-1 do
    begin
        C[3+I] := X[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        C[3+N+4*I+0] := Y[I];
        C[3+N+4*I+1] := (Y[I+1]-Y[I])/(X[I+1]-X[I]);
        C[3+N+4*I+2] := 0;
        C[3+N+4*I+3] := 0;
        Inc(I);
    end;
end;


procedure BuildCubicSpline(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     BoundLType : AlglibInteger;
     BoundL : Double;
     BoundRType : AlglibInteger;
     BoundR : Double;
     var C : TReal1DArray);
var
    A1 : TReal1DArray;
    A2 : TReal1DArray;
    A3 : TReal1DArray;
    B : TReal1DArray;
    D : TReal1DArray;
    I : AlglibInteger;
    TblSize : AlglibInteger;
    Delta : Double;
    Delta2 : Double;
    Delta3 : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>=2, 'BuildCubicSpline: N<2!');
    Assert((BoundLType=0) or (BoundLType=1) or (BoundLType=2), 'BuildCubicSpline: incorrect BoundLType!');
    Assert((BoundRType=0) or (BoundRType=1) or (BoundRType=2), 'BuildCubicSpline: incorrect BoundRType!');
    SetLength(A1, N-1+1);
    SetLength(A2, N-1+1);
    SetLength(A3, N-1+1);
    SetLength(B, N-1+1);
    
    //
    // Special case:
    // * N=2
    // * parabolic terminated boundary condition on both ends
    //
    if (N=2) and (BoundLType=0) and (BoundRType=0) then
    begin
        
        //
        // Change task type
        //
        BoundLType := 2;
        BoundL := 0;
        BoundRType := 2;
        BoundR := 0;
    end;
    
    //
    //
    // Sort points
    //
    HeapSortPoints(X, Y, N);
    
    //
    // Left boundary conditions
    //
    if BoundLType=0 then
    begin
        A1[0] := 0;
        A2[0] := 1;
        A3[0] := 1;
        B[0] := 2*(Y[1]-Y[0])/(X[1]-X[0]);
    end;
    if BoundLType=1 then
    begin
        A1[0] := 0;
        A2[0] := 1;
        A3[0] := 0;
        B[0] := BoundL;
    end;
    if BoundLType=2 then
    begin
        A1[0] := 0;
        A2[0] := 2;
        A3[0] := 1;
        B[0] := 3*(Y[1]-Y[0])/(X[1]-X[0])-Double(0.5)*BoundL*(X[1]-X[0]);
    end;
    
    //
    // Central conditions
    //
    I:=1;
    while I<=N-2 do
    begin
        A1[I] := X[I+1]-X[I];
        A2[I] := 2*(X[I+1]-X[I-1]);
        A3[I] := X[I]-X[I-1];
        B[I] := 3*(Y[I]-Y[I-1])/(X[I]-X[I-1])*(X[I+1]-X[I])+3*(Y[I+1]-Y[I])/(X[I+1]-X[I])*(X[I]-X[I-1]);
        Inc(I);
    end;
    
    //
    // Right boundary conditions
    //
    if BoundRType=0 then
    begin
        A1[N-1] := 1;
        A2[N-1] := 1;
        A3[N-1] := 0;
        B[N-1] := 2*(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2]);
    end;
    if BoundRType=1 then
    begin
        A1[N-1] := 0;
        A2[N-1] := 1;
        A3[N-1] := 0;
        B[N-1] := BoundR;
    end;
    if BoundRType=2 then
    begin
        A1[N-1] := 1;
        A2[N-1] := 2;
        A3[N-1] := 0;
        B[N-1] := 3*(Y[N-1]-Y[N-2])/(X[N-1]-X[N-2])+Double(0.5)*BoundR*(X[N-1]-X[N-2]);
    end;
    
    //
    // Solve
    //
    SolveTridiagonal(A1, A2, A3, B, N, D);
    
    //
    // Now problem is reduced to the cubic Hermite spline
    //
    BuildHermiteSpline(X, Y, D, N, C);
end;


procedure BuildHermiteSpline(X : TReal1DArray;
     Y : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var C : TReal1DArray);
var
    I : AlglibInteger;
    TblSize : AlglibInteger;
    Delta : Double;
    Delta2 : Double;
    Delta3 : Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    D := DynamicArrayCopy(D);
    Assert(N>=2, 'BuildHermiteSpline: N<2!');
    
    //
    // Sort points
    //
    HeapSortDPoints(X, Y, D, N);
    
    //
    // Fill C:
    //  C[0]            -   length(C)
    //  C[1]            -   type(C):
    //                      3 - general cubic spline
    //  C[2]            -   N
    //  C[3]...C[3+N-1] -   x[i], i = 0...N-1
    //  C[3+N]...C[3+N+(N-1)*4-1] - coefficients table
    //
    TblSize := 3+N+(N-1)*4;
    SetLength(C, TblSize-1+1);
    C[0] := TblSize;
    C[1] := 3;
    C[2] := N;
    I:=0;
    while I<=N-1 do
    begin
        C[3+I] := X[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        Delta := X[I+1]-X[I];
        Delta2 := AP_Sqr(Delta);
        Delta3 := Delta*Delta2;
        C[3+N+4*I+0] := Y[I];
        C[3+N+4*I+1] := D[I];
        C[3+N+4*I+2] := (3*(Y[I+1]-Y[I])-2*D[I]*Delta-D[I+1]*Delta)/Delta2;
        C[3+N+4*I+3] := (2*(Y[I]-Y[I+1])+D[I]*Delta+D[I+1]*Delta)/Delta3;
        Inc(I);
    end;
end;


procedure BuildAkimaSpline(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger;
     var C : TReal1DArray);
var
    I : AlglibInteger;
    D : TReal1DArray;
    W : TReal1DArray;
    Diff : TReal1DArray;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    Assert(N>=5, 'BuildAkimaSpline: N<5!');
    
    //
    // Sort points
    //
    HeapSortPoints(X, Y, N);
    
    //
    // Prepare W (weights), Diff (divided differences)
    //
    SetLength(W, N-2+1);
    SetLength(Diff, N-2+1);
    I:=0;
    while I<=N-2 do
    begin
        Diff[I] := (Y[I+1]-Y[I])/(X[I+1]-X[I]);
        Inc(I);
    end;
    I:=1;
    while I<=N-2 do
    begin
        W[I] := AbsReal(Diff[I]-Diff[I-1]);
        Inc(I);
    end;
    
    //
    // Prepare Hermite interpolation scheme
    //
    SetLength(D, N-1+1);
    I:=2;
    while I<=N-3 do
    begin
        if AP_FP_Neq(AbsReal(W[I-1])+AbsReal(W[I+1]),0) then
        begin
            D[I] := (W[I+1]*Diff[I-1]+W[I-1]*Diff[I])/(W[I+1]+W[I-1]);
        end
        else
        begin
            D[I] := ((X[I+1]-X[I])*Diff[I-1]+(X[I]-X[I-1])*Diff[I])/(X[I+1]-X[I-1]);
        end;
        Inc(I);
    end;
    D[0] := DiffThreePoint(X[0], X[0], Y[0], X[1], Y[1], X[2], Y[2]);
    D[1] := DiffThreePoint(X[1], X[0], Y[0], X[1], Y[1], X[2], Y[2]);
    D[N-2] := DiffThreePoint(X[N-2], X[N-3], Y[N-3], X[N-2], Y[N-2], X[N-1], Y[N-1]);
    D[N-1] := DiffThreePoint(X[N-1], X[N-3], Y[N-3], X[N-2], Y[N-2], X[N-1], Y[N-1]);
    
    //
    // Build Akima spline using Hermite interpolation scheme
    //
    BuildHermiteSpline(X, Y, D, N, C);
end;


function SplineInterpolation(const C : TReal1DArray; X : Double):Double;
var
    N : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    M : AlglibInteger;
begin
    Assert(Round(C[1])=3, 'SplineInterpolation: incorrect C!');
    N := Round(C[2]);
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    L := 3;
    R := 3+N-2+1;
    while L<>R-1 do
    begin
        M := (L+R) div 2;
        if AP_FP_Greater_Eq(C[M],X) then
        begin
            R := M;
        end
        else
        begin
            L := M;
        end;
    end;
    
    //
    // Interpolation
    //
    X := X-C[L];
    M := 3+N+4*(L-3);
    Result := C[M]+X*(C[M+1]+X*(C[M+2]+X*C[M+3]));
end;


procedure SplineDifferentiation(const C : TReal1DArray;
     X : Double;
     var S : Double;
     var DS : Double;
     var D2S : Double);
var
    N : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    M : AlglibInteger;
begin
    Assert(Round(C[1])=3, 'SplineInterpolation: incorrect C!');
    N := Round(C[2]);
    
    //
    // Binary search
    //
    L := 3;
    R := 3+N-2+1;
    while L<>R-1 do
    begin
        M := (L+R) div 2;
        if AP_FP_Greater_Eq(C[M],X) then
        begin
            R := M;
        end
        else
        begin
            L := M;
        end;
    end;
    
    //
    // Differentiation
    //
    X := X-C[L];
    M := 3+N+4*(L-3);
    S := C[M]+X*(C[M+1]+X*(C[M+2]+X*C[M+3]));
    DS := C[M+1]+2*X*C[M+2]+3*AP_Sqr(X)*C[M+3];
    D2S := 2*C[M+2]+6*X*C[M+3];
end;


procedure SplineCopy(const C : TReal1DArray; var CC : TReal1DArray);
var
    S : AlglibInteger;
begin
    S := Round(C[0]);
    SetLength(CC, S-1+1);
    APVMove(@CC[0], 0, S-1, @C[0], 0, S-1);
end;


procedure SplineUnpack(const C : TReal1DArray;
     var N : AlglibInteger;
     var Tbl : TReal2DArray);
var
    I : AlglibInteger;
begin
    Assert(Round(C[1])=3, 'SplineUnpack: incorrect C!');
    N := Round(C[2]);
    SetLength(Tbl, N-2+1, 5+1);
    
    //
    // Fill
    //
    I:=0;
    while I<=N-2 do
    begin
        Tbl[I,0] := C[3+I];
        Tbl[I,1] := C[3+I+1];
        Tbl[I,2] := C[3+N+4*I];
        Tbl[I,3] := C[3+N+4*I+1];
        Tbl[I,4] := C[3+N+4*I+2];
        Tbl[I,5] := C[3+N+4*I+3];
        Inc(I);
    end;
end;


procedure SplineLinTransX(var C : TReal1DArray; A : Double; B : Double);
var
    I : AlglibInteger;
    N : AlglibInteger;
    V : Double;
    DV : Double;
    D2V : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    D : TReal1DArray;
begin
    Assert(Round(C[1])=3, 'SplineLinTransX: incorrect C!');
    N := Round(C[2]);
    
    //
    // Special case: A=0
    //
    if AP_FP_Eq(A,0) then
    begin
        V := SplineInterpolation(C, B);
        I:=0;
        while I<=N-2 do
        begin
            C[3+N+4*I] := V;
            C[3+N+4*I+1] := 0;
            C[3+N+4*I+2] := 0;
            C[3+N+4*I+3] := 0;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // General case: A<>0.
    // Unpack, X, Y, dY/dX.
    // Scale and pack again.
    //
    SetLength(X, N-1+1);
    SetLength(Y, N-1+1);
    SetLength(D, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := C[3+I];
        SplineDifferentiation(C, X[I], V, DV, D2V);
        X[I] := (X[I]-B)/A;
        Y[I] := V;
        D[I] := A*DV;
        Inc(I);
    end;
    BuildHermiteSpline(X, Y, D, N, C);
end;


procedure SplineLinTransY(var C : TReal1DArray; A : Double; B : Double);
var
    I : AlglibInteger;
    N : AlglibInteger;
    V : Double;
    DV : Double;
    D2V : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    D : TReal1DArray;
begin
    Assert(Round(C[1])=3, 'SplineLinTransX: incorrect C!');
    N := Round(C[2]);
    
    //
    // Special case: A=0
    //
    I:=0;
    while I<=N-2 do
    begin
        C[3+N+4*I] := A*C[3+N+4*I]+B;
        C[3+N+4*I+1] := A*C[3+N+4*I+1];
        C[3+N+4*I+2] := A*C[3+N+4*I+2];
        C[3+N+4*I+3] := A*C[3+N+4*I+3];
        Inc(I);
    end;
end;


function SplineIntegration(const C : TReal1DArray; X : Double):Double;
var
    N : AlglibInteger;
    I : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    M : AlglibInteger;
    W : Double;
begin
    Assert(Round(C[1])=3, 'SplineIntegration: incorrect C!');
    N := Round(C[2]);
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    L := 3;
    R := 3+N-2+1;
    while L<>R-1 do
    begin
        M := (L+R) div 2;
        if AP_FP_Greater_Eq(C[M],X) then
        begin
            R := M;
        end
        else
        begin
            L := M;
        end;
    end;
    
    //
    // Integration
    //
    Result := 0;
    I:=3;
    while I<=L-1 do
    begin
        W := C[I+1]-C[I];
        M := 3+N+4*(I-3);
        Result := Result+C[M]*W;
        Result := Result+C[M+1]*AP_Sqr(W)/2;
        Result := Result+C[M+2]*AP_Sqr(W)*W/3;
        Result := Result+C[M+3]*AP_Sqr(AP_Sqr(W))/4;
        Inc(I);
    end;
    W := X-C[L];
    M := 3+N+4*(L-3);
    Result := Result+C[M]*W;
    Result := Result+C[M+1]*AP_Sqr(W)/2;
    Result := Result+C[M+2]*AP_Sqr(W)*W/3;
    Result := Result+C[M+3]*AP_Sqr(AP_Sqr(W))/4;
end;


procedure Spline3BuildTable(N : AlglibInteger;
     const DiffN : AlglibInteger;
     x : TReal1DArray;
     y : TReal1DArray;
     const BoundL : Double;
     const BoundR : Double;
     var ctbl : TReal2DArray);
var
    C : Boolean;
    E : AlglibInteger;
    G : AlglibInteger;
    Tmp : Double;
    nxm1 : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    DX : Double;
    DXJ : Double;
    DYJ : Double;
    DXJP1 : Double;
    DYJP1 : Double;
    DXP : Double;
    DYP : Double;
    YPPA : Double;
    YPPB : Double;
    PJ : Double;
    b1 : Double;
    b2 : Double;
    b3 : Double;
    b4 : Double;
begin
    x := DynamicArrayCopy(x);
    y := DynamicArrayCopy(y);
    N := N-1;
    g := (n+1) div 2;
    repeat
        i := g;
        repeat
            j := i-g;
            c := True;
            repeat
                if AP_FP_Less_Eq(x[j],x[j+g]) then
                begin
                    c := False;
                end
                else
                begin
                    Tmp := x[j];
                    x[j] := x[j+g];
                    x[j+g] := Tmp;
                    Tmp := y[j];
                    y[j] := y[j+g];
                    y[j+g] := Tmp;
                end;
                j := j-1;
            until  not ((j>=0) and C);
            i := i+1;
        until  not (i<=n);
        g := g div 2;
    until  not (g>0);
    SetLength(ctbl, 4+1, N+1);
    N := N+1;
    if DiffN=1 then
    begin
        b1 := 1;
        b2 := 6/(x[1]-x[0])*((y[1]-y[0])/(x[1]-x[0])-BoundL);
        b3 := 1;
        b4 := 6/(x[n-1]-x[n-2])*(BoundR-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    end
    else
    begin
        b1 := 0;
        b2 := 2*BoundL;
        b3 := 0;
        b4 := 2*BoundR;
    end;
    nxm1 := n-1;
    if n>=2 then
    begin
        if n>2 then
        begin
            dxj := x[1]-x[0];
            dyj := y[1]-y[0];
            j := 2;
            while j<=nxm1 do
            begin
                dxjp1 := x[j]-x[j-1];
                dyjp1 := y[j]-y[j-1];
                dxp := dxj+dxjp1;
                ctbl[1,j-1] := dxjp1/dxp;
                ctbl[2,j-1] := 1-ctbl[1,j-1];
                ctbl[3,j-1] := 6*(dyjp1/dxjp1-dyj/dxj)/dxp;
                dxj := dxjp1;
                dyj := dyjp1;
                j := j+1;
            end;
        end;
        ctbl[1,0] := -b1/2;
        ctbl[2,0] := b2/2;
        if n<>2 then
        begin
            j := 2;
            while j<=nxm1 do
            begin
                pj := ctbl[2,j-1]*ctbl[1,j-2]+2;
                ctbl[1,j-1] := -ctbl[1,j-1]/pj;
                ctbl[2,j-1] := (ctbl[3,j-1]-ctbl[2,j-1]*ctbl[2,J-2])/pj;
                j := j+1;
            end;
        end;
        yppb := (b4-b3*ctbl[2,nxm1-1])/(b3*ctbl[1,nxm1-1]+2);
        i := 1;
        while i<=nxm1 do
        begin
            j := n-i;
            yppa := ctbl[1,j-1]*yppb+ctbl[2,j-1];
            dx := x[j]-x[j-1];
            ctbl[3,j-1] := (yppb-yppa)/dx/6;
            ctbl[2,j-1] := yppa/2;
            ctbl[1,j-1] := (y[j]-y[j-1])/dx-(ctbl[2,j-1]+ctbl[3,j-1]*dx)*dx;
            yppb := yppa;
            i := i+1;
        end;
        i:=1;
        while i<=n do
        begin
            ctbl[0,i-1] := y[i-1];
            ctbl[4,i-1] := x[i-1];
            Inc(i);
        end;
    end;
end;


function Spline3Interpolate(N : AlglibInteger;
     const c : TReal2DArray;
     const X : Double):Double;
var
    I : AlglibInteger;
    L : AlglibInteger;
    Half : AlglibInteger;
    First : AlglibInteger;
    Middle : AlglibInteger;
begin
    N := N-1;
    L := N;
    First := 0;
    while L>0 do
    begin
        Half := L div 2;
        Middle := First+Half;
        if AP_FP_Less(C[4,Middle],X) then
        begin
            First := Middle+1;
            L := L-Half-1;
        end
        else
        begin
            L := Half;
        end;
    end;
    I := First-1;
    if I<0 then
    begin
        I := 0;
    end;
    Result := c[0,I]+(X-c[4,i])*(C[1,I]+(X-c[4,i])*(C[2,I]+C[3,i]*(X-c[4,i])));
end;


procedure HeapSortPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    IsAscending : Boolean;
    IsDescending : Boolean;
begin
    
    //
    // Test for already sorted set
    //
    IsAscending := True;
    IsDescending := True;
    I:=1;
    while I<=N-1 do
    begin
        IsAscending := IsAscending and AP_FP_Greater(X[I],X[I-1]);
        IsDescending := IsDescending and AP_FP_Less(X[I],X[I-1]);
        Inc(I);
    end;
    if IsAscending then
    begin
        Exit;
    end;
    if IsDescending then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J := N-1-I;
            if J<=I then
            begin
                Break;
            end;
            Tmp := X[I];
            X[I] := X[J];
            X[J] := Tmp;
            Tmp := Y[I];
            Y[I] := Y[J];
            Y[J] := Tmp;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        Exit;
    end;
    
    //
    // General case
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(X[k-1],X[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := X[k-1];
                X[k-1] := X[t-1];
                X[t-1] := Tmp;
                Tmp := Y[k-1];
                Y[k-1] := Y[t-1];
                Y[t-1] := Tmp;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := X[i];
        X[i] := X[0];
        X[0] := Tmp;
        Tmp := Y[i];
        Y[i] := Y[0];
        Y[0] := Tmp;
        t := 1;
        while t<>0 do
        begin
            k := 2*t;
            if k>i then
            begin
                t := 0;
            end
            else
            begin
                if k<i then
                begin
                    if AP_FP_Greater(X[k],X[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(X[t-1],X[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := X[k-1];
                    X[k-1] := X[t-1];
                    X[t-1] := Tmp;
                    Tmp := Y[k-1];
                    Y[k-1] := Y[t-1];
                    Y[t-1] := Tmp;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


procedure HeapSortDPoints(var X : TReal1DArray;
     var Y : TReal1DArray;
     var D : TReal1DArray;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    IsAscending : Boolean;
    IsDescending : Boolean;
begin
    
    //
    // Test for already sorted set
    //
    IsAscending := True;
    IsDescending := True;
    I:=1;
    while I<=N-1 do
    begin
        IsAscending := IsAscending and AP_FP_Greater(X[I],X[I-1]);
        IsDescending := IsDescending and AP_FP_Less(X[I],X[I-1]);
        Inc(I);
    end;
    if IsAscending then
    begin
        Exit;
    end;
    if IsDescending then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J := N-1-I;
            if J<=I then
            begin
                Break;
            end;
            Tmp := X[I];
            X[I] := X[J];
            X[J] := Tmp;
            Tmp := Y[I];
            Y[I] := Y[J];
            Y[J] := Tmp;
            Tmp := D[I];
            D[I] := D[J];
            D[J] := Tmp;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Special case: N=1
    //
    if N=1 then
    begin
        Exit;
    end;
    
    //
    // General case
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(X[k-1],X[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := X[k-1];
                X[k-1] := X[t-1];
                X[t-1] := Tmp;
                Tmp := Y[k-1];
                Y[k-1] := Y[t-1];
                Y[t-1] := Tmp;
                Tmp := D[k-1];
                D[k-1] := D[t-1];
                D[t-1] := Tmp;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := X[i];
        X[i] := X[0];
        X[0] := Tmp;
        Tmp := Y[i];
        Y[i] := Y[0];
        Y[0] := Tmp;
        Tmp := D[i];
        D[i] := D[0];
        D[0] := Tmp;
        t := 1;
        while t<>0 do
        begin
            k := 2*t;
            if k>i then
            begin
                t := 0;
            end
            else
            begin
                if k<i then
                begin
                    if AP_FP_Greater(X[k],X[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(X[t-1],X[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := X[k-1];
                    X[k-1] := X[t-1];
                    X[t-1] := Tmp;
                    Tmp := Y[k-1];
                    Y[k-1] := Y[t-1];
                    Y[t-1] := Tmp;
                    Tmp := D[k-1];
                    D[k-1] := D[t-1];
                    D[t-1] := Tmp;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


procedure SolveTridiagonal(A : TReal1DArray;
     B : TReal1DArray;
     C : TReal1DArray;
     D : TReal1DArray;
     N : AlglibInteger;
     var X : TReal1DArray);
var
    K : AlglibInteger;
    T : Double;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    C := DynamicArrayCopy(C);
    D := DynamicArrayCopy(D);
    SetLength(X, N-1+1);
    A[0] := 0;
    C[N-1] := 0;
    K:=1;
    while K<=N-1 do
    begin
        T := A[K]/B[K-1];
        B[K] := B[K]-T*C[K-1];
        D[K] := D[K]-T*D[K-1];
        Inc(K);
    end;
    X[N-1] := D[N-1]/B[N-1];
    K:=N-2;
    while K>=0 do
    begin
        X[K] := (D[K]-C[K]*X[K+1])/B[K];
        Dec(K);
    end;
end;


function DiffThreePoint(T : Double;
     X0 : Double;
     F0 : Double;
     X1 : Double;
     F1 : Double;
     X2 : Double;
     F2 : Double):Double;
var
    A : Double;
    B : Double;
begin
    T := T-X0;
    X1 := X1-X0;
    X2 := X2-X0;
    A := (F2-F0-X2/X1*(F1-F0))/(AP_Sqr(X2)-X1*X2);
    B := (F1-F0-A*AP_Sqr(X1))/X1;
    Result := 2*A*T+B;
end;


end.