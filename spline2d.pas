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
unit spline2d;
interface
uses Math, Sysutils, Ap, spline3, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, apserv, spline1d;

type
(*************************************************************************
2-dimensional spline inteprolant
*************************************************************************)
Spline2DInterpolant = record
    K : AlglibInteger;
    C : TReal1DArray;
end;



procedure Spline2DBuildBilinear(X : TReal1DArray;
     Y : TReal1DArray;
     F : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var C : Spline2DInterpolant);
procedure Spline2DBuildBicubic(X : TReal1DArray;
     Y : TReal1DArray;
     F : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var C : Spline2DInterpolant);
function Spline2DCalc(const C : Spline2DInterpolant;
     X : Double;
     Y : Double):Double;
procedure Spline2DDiff(const C : Spline2DInterpolant;
     X : Double;
     Y : Double;
     var F : Double;
     var FX : Double;
     var FY : Double;
     var FXY : Double);
procedure Spline2DUnpack(const C : Spline2DInterpolant;
     var M : AlglibInteger;
     var N : AlglibInteger;
     var Tbl : TReal2DArray);
procedure Spline2DLinTransXY(var C : Spline2DInterpolant;
     AX : Double;
     BX : Double;
     AY : Double;
     BY : Double);
procedure Spline2DLinTransF(var C : Spline2DInterpolant;
     A : Double;
     B : Double);
procedure Spline2DCopy(const C : Spline2DInterpolant;
     var CC : Spline2DInterpolant);
procedure Spline2DSerialize(const C : Spline2DInterpolant;
     var RA : TReal1DArray;
     var RALen : AlglibInteger);
procedure Spline2DUnserialize(const RA : TReal1DArray;
     var C : Spline2DInterpolant);
procedure Spline2DResampleBicubic(const A : TReal2DArray;
     OldHeight : AlglibInteger;
     OldWidth : AlglibInteger;
     var B : TReal2DArray;
     NewHeight : AlglibInteger;
     NewWidth : AlglibInteger);
procedure Spline2DResampleBilinear(const A : TReal2DArray;
     OldHeight : AlglibInteger;
     OldWidth : AlglibInteger;
     var B : TReal2DArray;
     NewHeight : AlglibInteger;
     NewWidth : AlglibInteger);

implementation

const
    Spline2DVNum = 12;

procedure BicubicCalcDerivatives(const A : TReal2DArray;
     const X : TReal1DArray;
     const Y : TReal1DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var DX : TReal2DArray;
     var DY : TReal2DArray;
     var DXY : TReal2DArray);forward;


(*************************************************************************
This subroutine builds bilinear spline coefficients table.

Input parameters:
    X   -   spline abscissas, array[0..N-1]
    Y   -   spline ordinates, array[0..M-1]
    F   -   function values, array[0..M-1,0..N-1]
    M,N -   grid size, M>=2, N>=2

Output parameters:
    C   -   spline interpolant

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DBuildBilinear(X : TReal1DArray;
     Y : TReal1DArray;
     F : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var C : Spline2DInterpolant);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TblSize : AlglibInteger;
    Shift : AlglibInteger;
    T : Double;
    DX : TReal2DArray;
    DY : TReal2DArray;
    DXY : TReal2DArray;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    F := DynamicArrayCopy(F);
    Assert((N>=2) and (M>=2), 'Spline2DBuildBilinear: N<2 or M<2!');
    
    //
    // Sort points
    //
    J:=0;
    while J<=N-1 do
    begin
        K := J;
        I:=J+1;
        while I<=N-1 do
        begin
            if AP_FP_Less(X[I],X[K]) then
            begin
                K := I;
            end;
            Inc(I);
        end;
        if K<>J then
        begin
            I:=0;
            while I<=M-1 do
            begin
                T := F[I,J];
                F[I,J] := F[I,K];
                F[I,K] := T;
                Inc(I);
            end;
            T := X[J];
            X[J] := X[K];
            X[K] := T;
        end;
        Inc(J);
    end;
    I:=0;
    while I<=M-1 do
    begin
        K := I;
        J:=I+1;
        while J<=M-1 do
        begin
            if AP_FP_Less(Y[J],Y[K]) then
            begin
                K := J;
            end;
            Inc(J);
        end;
        if K<>I then
        begin
            J:=0;
            while J<=N-1 do
            begin
                T := F[I,J];
                F[I,J] := F[K,J];
                F[K,J] := T;
                Inc(J);
            end;
            T := Y[I];
            Y[I] := Y[K];
            Y[K] := T;
        end;
        Inc(I);
    end;
    
    //
    // Fill C:
    //  C[0]            -   length(C)
    //  C[1]            -   type(C):
    //                      -1 = bilinear interpolant
    //                      -3 = general cubic spline
    //                           (see BuildBicubicSpline)
    //  C[2]:
    //      N (x count)
    //  C[3]:
    //      M (y count)
    //  C[4]...C[4+N-1]:
    //      x[i], i = 0...N-1
    //  C[4+N]...C[4+N+M-1]:
    //      y[i], i = 0...M-1
    //  C[4+N+M]...C[4+N+M+(N*M-1)]:
    //      f(i,j) table. f(0,0), f(0, 1), f(0,2) and so on...
    //
    C.K := 1;
    TblSize := 4+N+M+N*M;
    SetLength(C.C, TblSize-1+1);
    C.C[0] := TblSize;
    C.C[1] := -1;
    C.C[2] := N;
    C.C[3] := M;
    I:=0;
    while I<=N-1 do
    begin
        C.C[4+I] := X[I];
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        C.C[4+N+I] := Y[I];
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Shift := I*N+J;
            C.C[4+N+M+Shift] := F[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine builds bicubic spline coefficients table.

Input parameters:
    X   -   spline abscissas, array[0..N-1]
    Y   -   spline ordinates, array[0..M-1]
    F   -   function values, array[0..M-1,0..N-1]
    M,N -   grid size, M>=2, N>=2

Output parameters:
    C   -   spline interpolant

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DBuildBicubic(X : TReal1DArray;
     Y : TReal1DArray;
     F : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var C : Spline2DInterpolant);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TblSize : AlglibInteger;
    Shift : AlglibInteger;
    T : Double;
    DX : TReal2DArray;
    DY : TReal2DArray;
    DXY : TReal2DArray;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    F := DynamicArrayCopy(F);
    Assert((N>=2) and (M>=2), 'BuildBicubicSpline: N<2 or M<2!');
    
    //
    // Sort points
    //
    J:=0;
    while J<=N-1 do
    begin
        K := J;
        I:=J+1;
        while I<=N-1 do
        begin
            if AP_FP_Less(X[I],X[K]) then
            begin
                K := I;
            end;
            Inc(I);
        end;
        if K<>J then
        begin
            I:=0;
            while I<=M-1 do
            begin
                T := F[I,J];
                F[I,J] := F[I,K];
                F[I,K] := T;
                Inc(I);
            end;
            T := X[J];
            X[J] := X[K];
            X[K] := T;
        end;
        Inc(J);
    end;
    I:=0;
    while I<=M-1 do
    begin
        K := I;
        J:=I+1;
        while J<=M-1 do
        begin
            if AP_FP_Less(Y[J],Y[K]) then
            begin
                K := J;
            end;
            Inc(J);
        end;
        if K<>I then
        begin
            J:=0;
            while J<=N-1 do
            begin
                T := F[I,J];
                F[I,J] := F[K,J];
                F[K,J] := T;
                Inc(J);
            end;
            T := Y[I];
            Y[I] := Y[K];
            Y[K] := T;
        end;
        Inc(I);
    end;
    
    //
    // Fill C:
    //  C[0]            -   length(C)
    //  C[1]            -   type(C):
    //                      -1 = bilinear interpolant
    //                           (see BuildBilinearInterpolant)
    //                      -3 = general cubic spline
    //  C[2]:
    //      N (x count)
    //  C[3]:
    //      M (y count)
    //  C[4]...C[4+N-1]:
    //      x[i], i = 0...N-1
    //  C[4+N]...C[4+N+M-1]:
    //      y[i], i = 0...M-1
    //  C[4+N+M]...C[4+N+M+(N*M-1)]:
    //      f(i,j) table. f(0,0), f(0, 1), f(0,2) and so on...
    //  C[4+N+M+N*M]...C[4+N+M+(2*N*M-1)]:
    //      df(i,j)/dx table.
    //  C[4+N+M+2*N*M]...C[4+N+M+(3*N*M-1)]:
    //      df(i,j)/dy table.
    //  C[4+N+M+3*N*M]...C[4+N+M+(4*N*M-1)]:
    //      d2f(i,j)/dxdy table.
    //
    C.K := 3;
    TblSize := 4+N+M+4*N*M;
    SetLength(C.C, TblSize-1+1);
    C.C[0] := TblSize;
    C.C[1] := -3;
    C.C[2] := N;
    C.C[3] := M;
    I:=0;
    while I<=N-1 do
    begin
        C.C[4+I] := X[I];
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        C.C[4+N+I] := Y[I];
        Inc(I);
    end;
    BicubicCalcDerivatives(F, X, Y, M, N, DX, DY, DXY);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Shift := I*N+J;
            C.C[4+N+M+Shift] := F[I,J];
            C.C[4+N+M+N*M+Shift] := DX[I,J];
            C.C[4+N+M+2*N*M+Shift] := DY[I,J];
            C.C[4+N+M+3*N*M+Shift] := DXY[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine calculates the value of the bilinear or bicubic spline  at
the given point X.

Input parameters:
    C   -   coefficients table.
            Built by BuildBilinearSpline or BuildBicubicSpline.
    X, Y-   point

Result:
    S(x,y)

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************)
function Spline2DCalc(const C : Spline2DInterpolant;
     X : Double;
     Y : Double):Double;
var
    V : Double;
    VX : Double;
    VY : Double;
    VXY : Double;
begin
    Spline2DDiff(C, X, Y, V, VX, VY, VXY);
    Result := V;
end;


(*************************************************************************
This subroutine calculates the value of the bilinear or bicubic spline  at
the given point X and its derivatives.

Input parameters:
    C   -   spline interpolant.
    X, Y-   point

Output parameters:
    F   -   S(x,y)
    FX  -   dS(x,y)/dX
    FY  -   dS(x,y)/dY
    FXY -   d2S(x,y)/dXdY

  -- ALGLIB PROJECT --
     Copyright 05.07.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DDiff(const C : Spline2DInterpolant;
     X : Double;
     Y : Double;
     var F : Double;
     var FX : Double;
     var FY : Double;
     var FXY : Double);
var
    N : AlglibInteger;
    M : AlglibInteger;
    T : Double;
    DT : Double;
    U : Double;
    DU : Double;
    IX : AlglibInteger;
    IY : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    H : AlglibInteger;
    Shift1 : AlglibInteger;
    S1 : AlglibInteger;
    S2 : AlglibInteger;
    S3 : AlglibInteger;
    S4 : AlglibInteger;
    SF : AlglibInteger;
    SFX : AlglibInteger;
    SFY : AlglibInteger;
    SFXY : AlglibInteger;
    Y1 : Double;
    Y2 : Double;
    Y3 : Double;
    Y4 : Double;
    V : Double;
    T0 : Double;
    T1 : Double;
    T2 : Double;
    T3 : Double;
    U0 : Double;
    U1 : Double;
    U2 : Double;
    U3 : Double;
begin
    Assert((Round(C.C[1])=-1) or (Round(C.C[1])=-3), 'Spline2DDiff: incorrect C!');
    N := Round(C.C[2]);
    M := Round(C.C[3]);
    
    //
    // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
    //
    L := 4;
    R := 4+N-2+1;
    while L<>R-1 do
    begin
        H := (L+R) div 2;
        if AP_FP_Greater_Eq(C.C[H],X) then
        begin
            R := H;
        end
        else
        begin
            L := H;
        end;
    end;
    T := (X-C.C[L])/(C.C[L+1]-C.C[L]);
    DT := Double(1.0)/(C.C[L+1]-C.C[L]);
    IX := L-4;
    
    //
    // Binary search in the [ y[0], ..., y[m-2] ] (y[m-1] is not included)
    //
    L := 4+N;
    R := 4+N+(M-2)+1;
    while L<>R-1 do
    begin
        H := (L+R) div 2;
        if AP_FP_Greater_Eq(C.C[H],Y) then
        begin
            R := H;
        end
        else
        begin
            L := H;
        end;
    end;
    U := (Y-C.C[L])/(C.C[L+1]-C.C[L]);
    DU := Double(1.0)/(C.C[L+1]-C.C[L]);
    IY := L-(4+N);
    
    //
    // Prepare F, dF/dX, dF/dY, d2F/dXdY
    //
    F := 0;
    FX := 0;
    FY := 0;
    FXY := 0;
    
    //
    // Bilinear interpolation
    //
    if Round(C.C[1])=-1 then
    begin
        Shift1 := 4+N+M;
        Y1 := C.C[Shift1+N*IY+IX];
        Y2 := C.C[Shift1+N*IY+(IX+1)];
        Y3 := C.C[Shift1+N*(IY+1)+(IX+1)];
        Y4 := C.C[Shift1+N*(IY+1)+IX];
        F := (1-T)*(1-U)*Y1+T*(1-U)*Y2+T*U*Y3+(1-T)*U*Y4;
        FX := (-(1-U)*Y1+(1-U)*Y2+U*Y3-U*Y4)*DT;
        FY := (-(1-T)*Y1-T*Y2+T*Y3+(1-T)*Y4)*DU;
        FXY := (Y1-Y2+Y3-Y4)*DU*DT;
        Exit;
    end;
    
    //
    // Bicubic interpolation
    //
    if Round(C.C[1])=-3 then
    begin
        
        //
        // Prepare info
        //
        T0 := 1;
        T1 := T;
        T2 := AP_Sqr(T);
        T3 := T*T2;
        U0 := 1;
        U1 := U;
        U2 := AP_Sqr(U);
        U3 := U*U2;
        SF := 4+N+M;
        SFX := 4+N+M+N*M;
        SFY := 4+N+M+2*N*M;
        SFXY := 4+N+M+3*N*M;
        S1 := N*IY+IX;
        S2 := N*IY+(IX+1);
        S3 := N*(IY+1)+(IX+1);
        S4 := N*(IY+1)+IX;
        
        //
        // Calculate
        //
        V := +1*C.C[SF+S1];
        F := F+V*T0*U0;
        V := +1*C.C[SFY+S1]/DU;
        F := F+V*T0*U1;
        FY := FY+1*V*T0*U0*DU;
        V := -3*C.C[SF+S1]+3*C.C[SF+S4]-2*C.C[SFY+S1]/DU-1*C.C[SFY+S4]/DU;
        F := F+V*T0*U2;
        FY := FY+2*V*T0*U1*DU;
        V := +2*C.C[SF+S1]-2*C.C[SF+S4]+1*C.C[SFY+S1]/DU+1*C.C[SFY+S4]/DU;
        F := F+V*T0*U3;
        FY := FY+3*V*T0*U2*DU;
        V := +1*C.C[SFX+S1]/DT;
        F := F+V*T1*U0;
        FX := FX+1*V*T0*U0*DT;
        V := +1*C.C[SFXY+S1]/(DT*DU);
        F := F+V*T1*U1;
        FX := FX+1*V*T0*U1*DT;
        FY := FY+1*V*T1*U0*DU;
        FXY := FXY+1*V*T0*U0*DT*DU;
        V := -3*C.C[SFX+S1]/DT+3*C.C[SFX+S4]/DT-2*C.C[SFXY+S1]/(DT*DU)-1*C.C[SFXY+S4]/(DT*DU);
        F := F+V*T1*U2;
        FX := FX+1*V*T0*U2*DT;
        FY := FY+2*V*T1*U1*DU;
        FXY := FXY+2*V*T0*U1*DT*DU;
        V := +2*C.C[SFX+S1]/DT-2*C.C[SFX+S4]/DT+1*C.C[SFXY+S1]/(DT*DU)+1*C.C[SFXY+S4]/(DT*DU);
        F := F+V*T1*U3;
        FX := FX+1*V*T0*U3*DT;
        FY := FY+3*V*T1*U2*DU;
        FXY := FXY+3*V*T0*U2*DT*DU;
        V := -3*C.C[SF+S1]+3*C.C[SF+S2]-2*C.C[SFX+S1]/DT-1*C.C[SFX+S2]/DT;
        F := F+V*T2*U0;
        FX := FX+2*V*T1*U0*DT;
        V := -3*C.C[SFY+S1]/DU+3*C.C[SFY+S2]/DU-2*C.C[SFXY+S1]/(DT*DU)-1*C.C[SFXY+S2]/(DT*DU);
        F := F+V*T2*U1;
        FX := FX+2*V*T1*U1*DT;
        FY := FY+1*V*T2*U0*DU;
        FXY := FXY+2*V*T1*U0*DT*DU;
        V := +9*C.C[SF+S1]-9*C.C[SF+S2]+9*C.C[SF+S3]-9*C.C[SF+S4]+6*C.C[SFX+S1]/DT+3*C.C[SFX+S2]/DT-3*C.C[SFX+S3]/DT-6*C.C[SFX+S4]/DT+6*C.C[SFY+S1]/DU-6*C.C[SFY+S2]/DU-3*C.C[SFY+S3]/DU+3*C.C[SFY+S4]/DU+4*C.C[SFXY+S1]/(DT*DU)+2*C.C[SFXY+S2]/(DT*DU)+1*C.C[SFXY+S3]/(DT*DU)+2*C.C[SFXY+S4]/(DT*DU);
        F := F+V*T2*U2;
        FX := FX+2*V*T1*U2*DT;
        FY := FY+2*V*T2*U1*DU;
        FXY := FXY+4*V*T1*U1*DT*DU;
        V := -6*C.C[SF+S1]+6*C.C[SF+S2]-6*C.C[SF+S3]+6*C.C[SF+S4]-4*C.C[SFX+S1]/DT-2*C.C[SFX+S2]/DT+2*C.C[SFX+S3]/DT+4*C.C[SFX+S4]/DT-3*C.C[SFY+S1]/DU+3*C.C[SFY+S2]/DU+3*C.C[SFY+S3]/DU-3*C.C[SFY+S4]/DU-2*C.C[SFXY+S1]/(DT*DU)-1*C.C[SFXY+S2]/(DT*DU)-1*C.C[SFXY+S3]/(DT*DU)-2*C.C[SFXY+S4]/(DT*DU);
        F := F+V*T2*U3;
        FX := FX+2*V*T1*U3*DT;
        FY := FY+3*V*T2*U2*DU;
        FXY := FXY+6*V*T1*U2*DT*DU;
        V := +2*C.C[SF+S1]-2*C.C[SF+S2]+1*C.C[SFX+S1]/DT+1*C.C[SFX+S2]/DT;
        F := F+V*T3*U0;
        FX := FX+3*V*T2*U0*DT;
        V := +2*C.C[SFY+S1]/DU-2*C.C[SFY+S2]/DU+1*C.C[SFXY+S1]/(DT*DU)+1*C.C[SFXY+S2]/(DT*DU);
        F := F+V*T3*U1;
        FX := FX+3*V*T2*U1*DT;
        FY := FY+1*V*T3*U0*DU;
        FXY := FXY+3*V*T2*U0*DT*DU;
        V := -6*C.C[SF+S1]+6*C.C[SF+S2]-6*C.C[SF+S3]+6*C.C[SF+S4]-3*C.C[SFX+S1]/DT-3*C.C[SFX+S2]/DT+3*C.C[SFX+S3]/DT+3*C.C[SFX+S4]/DT-4*C.C[SFY+S1]/DU+4*C.C[SFY+S2]/DU+2*C.C[SFY+S3]/DU-2*C.C[SFY+S4]/DU-2*C.C[SFXY+S1]/(DT*DU)-2*C.C[SFXY+S2]/(DT*DU)-1*C.C[SFXY+S3]/(DT*DU)-1*C.C[SFXY+S4]/(DT*DU);
        F := F+V*T3*U2;
        FX := FX+3*V*T2*U2*DT;
        FY := FY+2*V*T3*U1*DU;
        FXY := FXY+6*V*T2*U1*DT*DU;
        V := +4*C.C[SF+S1]-4*C.C[SF+S2]+4*C.C[SF+S3]-4*C.C[SF+S4]+2*C.C[SFX+S1]/DT+2*C.C[SFX+S2]/DT-2*C.C[SFX+S3]/DT-2*C.C[SFX+S4]/DT+2*C.C[SFY+S1]/DU-2*C.C[SFY+S2]/DU-2*C.C[SFY+S3]/DU+2*C.C[SFY+S4]/DU+1*C.C[SFXY+S1]/(DT*DU)+1*C.C[SFXY+S2]/(DT*DU)+1*C.C[SFXY+S3]/(DT*DU)+1*C.C[SFXY+S4]/(DT*DU);
        F := F+V*T3*U3;
        FX := FX+3*V*T2*U3*DT;
        FY := FY+3*V*T3*U2*DU;
        FXY := FXY+9*V*T2*U2*DT*DU;
        Exit;
    end;
end;


(*************************************************************************
This subroutine unpacks two-dimensional spline into the coefficients table

Input parameters:
    C   -   spline interpolant.

Result:
    M, N-   grid size (x-axis and y-axis)
    Tbl -   coefficients table, unpacked format,
            [0..(N-1)*(M-1)-1, 0..19].
            For I = 0...M-2, J=0..N-2:
                K =  I*(N-1)+J
                Tbl[K,0] = X[j]
                Tbl[K,1] = X[j+1]
                Tbl[K,2] = Y[i]
                Tbl[K,3] = Y[i+1]
                Tbl[K,4] = C00
                Tbl[K,5] = C01
                Tbl[K,6] = C02
                Tbl[K,7] = C03
                Tbl[K,8] = C10
                Tbl[K,9] = C11
                ...
                Tbl[K,19] = C33
            On each grid square spline is equals to:
                S(x) = SUM(c[i,j]*(x^i)*(y^j), i=0..3, j=0..3)
                t = x-x[j]
                u = y-y[i]

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DUnpack(const C : Spline2DInterpolant;
     var M : AlglibInteger;
     var N : AlglibInteger;
     var Tbl : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    CI : AlglibInteger;
    CJ : AlglibInteger;
    K : AlglibInteger;
    P : AlglibInteger;
    Shift : AlglibInteger;
    S1 : AlglibInteger;
    S2 : AlglibInteger;
    S3 : AlglibInteger;
    S4 : AlglibInteger;
    SF : AlglibInteger;
    SFX : AlglibInteger;
    SFY : AlglibInteger;
    SFXY : AlglibInteger;
    Y1 : Double;
    Y2 : Double;
    Y3 : Double;
    Y4 : Double;
    DT : Double;
    DU : Double;
begin
    Assert((Round(C.C[1])=-3) or (Round(C.C[1])=-1), 'SplineUnpack2D: incorrect C!');
    N := Round(C.C[2]);
    M := Round(C.C[3]);
    SetLength(Tbl, (N-1)*(M-1)-1+1, 19+1);
    
    //
    // Fill
    //
    I:=0;
    while I<=M-2 do
    begin
        J:=0;
        while J<=N-2 do
        begin
            P := I*(N-1)+J;
            Tbl[P,0] := C.C[4+J];
            Tbl[P,1] := C.C[4+J+1];
            Tbl[P,2] := C.C[4+N+I];
            Tbl[P,3] := C.C[4+N+I+1];
            DT := 1/(Tbl[P,1]-Tbl[P,0]);
            DU := 1/(Tbl[P,3]-Tbl[P,2]);
            
            //
            // Bilinear interpolation
            //
            if Round(C.C[1])=-1 then
            begin
                K:=4;
                while K<=19 do
                begin
                    Tbl[P,K] := 0;
                    Inc(K);
                end;
                Shift := 4+N+M;
                Y1 := C.C[Shift+N*I+J];
                Y2 := C.C[Shift+N*I+(J+1)];
                Y3 := C.C[Shift+N*(I+1)+(J+1)];
                Y4 := C.C[Shift+N*(I+1)+J];
                Tbl[P,4] := Y1;
                Tbl[P,4+1*4+0] := Y2-Y1;
                Tbl[P,4+0*4+1] := Y4-Y1;
                Tbl[P,4+1*4+1] := Y3-Y2-Y4+Y1;
            end;
            
            //
            // Bicubic interpolation
            //
            if Round(C.C[1])=-3 then
            begin
                SF := 4+N+M;
                SFX := 4+N+M+N*M;
                SFY := 4+N+M+2*N*M;
                SFXY := 4+N+M+3*N*M;
                S1 := N*I+J;
                S2 := N*I+(J+1);
                S3 := N*(I+1)+(J+1);
                S4 := N*(I+1)+J;
                Tbl[P,4+0*4+0] := +1*C.C[SF+S1];
                Tbl[P,4+0*4+1] := +1*C.C[SFY+S1]/DU;
                Tbl[P,4+0*4+2] := -3*C.C[SF+S1]+3*C.C[SF+S4]-2*C.C[SFY+S1]/DU-1*C.C[SFY+S4]/DU;
                Tbl[P,4+0*4+3] := +2*C.C[SF+S1]-2*C.C[SF+S4]+1*C.C[SFY+S1]/DU+1*C.C[SFY+S4]/DU;
                Tbl[P,4+1*4+0] := +1*C.C[SFX+S1]/DT;
                Tbl[P,4+1*4+1] := +1*C.C[SFXY+S1]/(DT*DU);
                Tbl[P,4+1*4+2] := -3*C.C[SFX+S1]/DT+3*C.C[SFX+S4]/DT-2*C.C[SFXY+S1]/(DT*DU)-1*C.C[SFXY+S4]/(DT*DU);
                Tbl[P,4+1*4+3] := +2*C.C[SFX+S1]/DT-2*C.C[SFX+S4]/DT+1*C.C[SFXY+S1]/(DT*DU)+1*C.C[SFXY+S4]/(DT*DU);
                Tbl[P,4+2*4+0] := -3*C.C[SF+S1]+3*C.C[SF+S2]-2*C.C[SFX+S1]/DT-1*C.C[SFX+S2]/DT;
                Tbl[P,4+2*4+1] := -3*C.C[SFY+S1]/DU+3*C.C[SFY+S2]/DU-2*C.C[SFXY+S1]/(DT*DU)-1*C.C[SFXY+S2]/(DT*DU);
                Tbl[P,4+2*4+2] := +9*C.C[SF+S1]-9*C.C[SF+S2]+9*C.C[SF+S3]-9*C.C[SF+S4]+6*C.C[SFX+S1]/DT+3*C.C[SFX+S2]/DT-3*C.C[SFX+S3]/DT-6*C.C[SFX+S4]/DT+6*C.C[SFY+S1]/DU-6*C.C[SFY+S2]/DU-3*C.C[SFY+S3]/DU+3*C.C[SFY+S4]/DU+4*C.C[SFXY+S1]/(DT*DU)+2*C.C[SFXY+S2]/(DT*DU)+1*C.C[SFXY+S3]/(DT*DU)+2*C.C[SFXY+S4]/(DT*DU);
                Tbl[P,4+2*4+3] := -6*C.C[SF+S1]+6*C.C[SF+S2]-6*C.C[SF+S3]+6*C.C[SF+S4]-4*C.C[SFX+S1]/DT-2*C.C[SFX+S2]/DT+2*C.C[SFX+S3]/DT+4*C.C[SFX+S4]/DT-3*C.C[SFY+S1]/DU+3*C.C[SFY+S2]/DU+3*C.C[SFY+S3]/DU-3*C.C[SFY+S4]/DU-2*C.C[SFXY+S1]/(DT*DU)-1*C.C[SFXY+S2]/(DT*DU)-1*C.C[SFXY+S3]/(DT*DU)-2*C.C[SFXY+S4]/(DT*DU);
                Tbl[P,4+3*4+0] := +2*C.C[SF+S1]-2*C.C[SF+S2]+1*C.C[SFX+S1]/DT+1*C.C[SFX+S2]/DT;
                Tbl[P,4+3*4+1] := +2*C.C[SFY+S1]/DU-2*C.C[SFY+S2]/DU+1*C.C[SFXY+S1]/(DT*DU)+1*C.C[SFXY+S2]/(DT*DU);
                Tbl[P,4+3*4+2] := -6*C.C[SF+S1]+6*C.C[SF+S2]-6*C.C[SF+S3]+6*C.C[SF+S4]-3*C.C[SFX+S1]/DT-3*C.C[SFX+S2]/DT+3*C.C[SFX+S3]/DT+3*C.C[SFX+S4]/DT-4*C.C[SFY+S1]/DU+4*C.C[SFY+S2]/DU+2*C.C[SFY+S3]/DU-2*C.C[SFY+S4]/DU-2*C.C[SFXY+S1]/(DT*DU)-2*C.C[SFXY+S2]/(DT*DU)-1*C.C[SFXY+S3]/(DT*DU)-1*C.C[SFXY+S4]/(DT*DU);
                Tbl[P,4+3*4+3] := +4*C.C[SF+S1]-4*C.C[SF+S2]+4*C.C[SF+S3]-4*C.C[SF+S4]+2*C.C[SFX+S1]/DT+2*C.C[SFX+S2]/DT-2*C.C[SFX+S3]/DT-2*C.C[SFX+S4]/DT+2*C.C[SFY+S1]/DU-2*C.C[SFY+S2]/DU-2*C.C[SFY+S3]/DU+2*C.C[SFY+S4]/DU+1*C.C[SFXY+S1]/(DT*DU)+1*C.C[SFXY+S2]/(DT*DU)+1*C.C[SFXY+S3]/(DT*DU)+1*C.C[SFXY+S4]/(DT*DU);
            end;
            
            //
            // Rescale Cij
            //
            CI:=0;
            while CI<=3 do
            begin
                CJ:=0;
                while CJ<=3 do
                begin
                    Tbl[P,4+CI*4+CJ] := Tbl[P,4+CI*4+CJ]*Power(DT, CI)*Power(DU, CJ);
                    Inc(CJ);
                end;
                Inc(CI);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
This subroutine performs linear transformation of the spline argument.

Input parameters:
    C       -   spline interpolant
    AX, BX  -   transformation coefficients: x = A*t + B
    AY, BY  -   transformation coefficients: y = A*u + B
Result:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DLinTransXY(var C : Spline2DInterpolant;
     AX : Double;
     BX : Double;
     AY : Double;
     BY : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    V : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    F : TReal2DArray;
    TypeC : AlglibInteger;
begin
    TypeC := Round(C.C[1]);
    Assert((TypeC=-3) or (TypeC=-1), 'Spline2DLinTransXY: incorrect C!');
    N := Round(C.C[2]);
    M := Round(C.C[3]);
    SetLength(X, N-1+1);
    SetLength(Y, M-1+1);
    SetLength(F, M-1+1, N-1+1);
    J:=0;
    while J<=N-1 do
    begin
        X[J] := C.C[4+J];
        Inc(J);
    end;
    I:=0;
    while I<=M-1 do
    begin
        Y[I] := C.C[4+N+I];
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            F[I,J] := C.C[4+N+M+I*N+J];
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Special case: AX=0 or AY=0
    //
    if AP_FP_Eq(AX,0) then
    begin
        I:=0;
        while I<=M-1 do
        begin
            V := Spline2DCalc(C, BX, Y[I]);
            J:=0;
            while J<=N-1 do
            begin
                F[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        if TypeC=-3 then
        begin
            Spline2DBuildBicubic(X, Y, F, M, N, C);
        end;
        if TypeC=-1 then
        begin
            Spline2DBuildBilinear(X, Y, F, M, N, C);
        end;
        AX := 1;
        BX := 0;
    end;
    if AP_FP_Eq(AY,0) then
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := Spline2DCalc(C, X[J], BY);
            I:=0;
            while I<=M-1 do
            begin
                F[I,J] := V;
                Inc(I);
            end;
            Inc(J);
        end;
        if TypeC=-3 then
        begin
            Spline2DBuildBicubic(X, Y, F, M, N, C);
        end;
        if TypeC=-1 then
        begin
            Spline2DBuildBilinear(X, Y, F, M, N, C);
        end;
        AY := 1;
        BY := 0;
    end;
    
    //
    // General case: AX<>0, AY<>0
    // Unpack, scale and pack again.
    //
    J:=0;
    while J<=N-1 do
    begin
        X[J] := (X[J]-BX)/AX;
        Inc(J);
    end;
    I:=0;
    while I<=M-1 do
    begin
        Y[I] := (Y[I]-BY)/AY;
        Inc(I);
    end;
    if TypeC=-3 then
    begin
        Spline2DBuildBicubic(X, Y, F, M, N, C);
    end;
    if TypeC=-1 then
    begin
        Spline2DBuildBilinear(X, Y, F, M, N, C);
    end;
end;


(*************************************************************************
This subroutine performs linear transformation of the spline.

Input parameters:
    C   -   spline interpolant.
    A, B-   transformation coefficients: S2(x,y) = A*S(x,y) + B
    
Output parameters:
    C   -   transformed spline

  -- ALGLIB PROJECT --
     Copyright 30.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DLinTransF(var C : Spline2DInterpolant;
     A : Double;
     B : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    X : TReal1DArray;
    Y : TReal1DArray;
    F : TReal2DArray;
    TypeC : AlglibInteger;
begin
    TypeC := Round(C.C[1]);
    Assert((TypeC=-3) or (TypeC=-1), 'Spline2DLinTransXY: incorrect C!');
    N := Round(C.C[2]);
    M := Round(C.C[3]);
    SetLength(X, N-1+1);
    SetLength(Y, M-1+1);
    SetLength(F, M-1+1, N-1+1);
    J:=0;
    while J<=N-1 do
    begin
        X[J] := C.C[4+J];
        Inc(J);
    end;
    I:=0;
    while I<=M-1 do
    begin
        Y[I] := C.C[4+N+I];
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            F[I,J] := A*C.C[4+N+M+I*N+J]+B;
            Inc(J);
        end;
        Inc(I);
    end;
    if TypeC=-3 then
    begin
        Spline2DBuildBicubic(X, Y, F, M, N, C);
    end;
    if TypeC=-1 then
    begin
        Spline2DBuildBilinear(X, Y, F, M, N, C);
    end;
end;


(*************************************************************************
This subroutine makes the copy of the spline model.

Input parameters:
    C   -   spline interpolant

Output parameters:
    CC  -   spline copy

  -- ALGLIB PROJECT --
     Copyright 29.06.2007 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DCopy(const C : Spline2DInterpolant;
     var CC : Spline2DInterpolant);
var
    N : AlglibInteger;
begin
    Assert((C.K=1) or (C.K=3), 'Spline2DCopy: incorrect C!');
    CC.K := C.K;
    N := Round(C.C[0]);
    SetLength(CC.C, N);
    APVMove(@CC.C[0], 0, N-1, @C.C[0], 0, N-1);
end;


(*************************************************************************
Serialization of the spline interpolant

INPUT PARAMETERS:
    B   -   spline interpolant

OUTPUT PARAMETERS:
    RA      -   array of real numbers which contains interpolant,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DSerialize(const C : Spline2DInterpolant;
     var RA : TReal1DArray;
     var RALen : AlglibInteger);
var
    CLen : AlglibInteger;
begin
    Assert((C.K=1) or (C.K=3), 'Spline2DSerialize: incorrect C!');
    CLen := Round(C.C[0]);
    RALen := 3+CLen;
    SetLength(RA, RALen);
    RA[0] := RALen;
    RA[1] := Spline2DVNum;
    RA[2] := C.K;
    APVMove(@RA[0], 3, 3+CLen-1, @C.C[0], 0, CLen-1);
end;


(*************************************************************************
Unserialization of the spline interpolant

INPUT PARAMETERS:
    RA  -   array of real numbers which contains interpolant,

OUTPUT PARAMETERS:
    B   -   spline interpolant

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure Spline2DUnserialize(const RA : TReal1DArray;
     var C : Spline2DInterpolant);
var
    CLen : AlglibInteger;
begin
    Assert(Round(RA[1])=Spline2DVNum, 'Spline2DUnserialize: corrupted array!');
    C.K := Round(RA[2]);
    CLen := Round(RA[3]);
    SetLength(C.C, CLen);
    APVMove(@C.C[0], 0, CLen-1, @RA[0], 3, 3+CLen-1);
end;


(*************************************************************************
Bicubic spline resampling

Input parameters:
    A           -   function values at the old grid,
                    array[0..OldHeight-1, 0..OldWidth-1]
    OldHeight   -   old grid height, OldHeight>1
    OldWidth    -   old grid width, OldWidth>1
    NewHeight   -   new grid height, NewHeight>1
    NewWidth    -   new grid width, NewWidth>1
    
Output parameters:
    B           -   function values at the new grid,
                    array[0..NewHeight-1, 0..NewWidth-1]

  -- ALGLIB routine --
     15 May, 2007
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure Spline2DResampleBicubic(const A : TReal2DArray;
     OldHeight : AlglibInteger;
     OldWidth : AlglibInteger;
     var B : TReal2DArray;
     NewHeight : AlglibInteger;
     NewWidth : AlglibInteger);
var
    Buf : TReal2DArray;
    X : TReal1DArray;
    Y : TReal1DArray;
    C : Spline1DInterpolant;
    I : AlglibInteger;
    J : AlglibInteger;
    MW : AlglibInteger;
    MH : AlglibInteger;
begin
    Assert((OldWidth>1) and (OldHeight>1), 'Spline2DResampleBicubic: width/height less than 1');
    Assert((NewWidth>1) and (NewHeight>1), 'Spline2DResampleBicubic: width/height less than 1');
    
    //
    // Prepare
    //
    MW := Max(OldWidth, NewWidth);
    MH := Max(OldHeight, NewHeight);
    SetLength(B, NewHeight-1+1, NewWidth-1+1);
    SetLength(Buf, OldHeight-1+1, NewWidth-1+1);
    SetLength(X, Max(MW, MH)-1+1);
    SetLength(Y, Max(MW, MH)-1+1);
    
    //
    // Horizontal interpolation
    //
    I:=0;
    while I<=OldHeight-1 do
    begin
        
        //
        // Fill X, Y
        //
        J:=0;
        while J<=OldWidth-1 do
        begin
            X[J] := AP_Double(J)/(OldWidth-1);
            Y[J] := A[I,J];
            Inc(J);
        end;
        
        //
        // Interpolate and place result into temporary matrix
        //
        Spline1DBuildCubic(X, Y, OldWidth, 0, Double(0.0), 0, Double(0.0), C);
        J:=0;
        while J<=NewWidth-1 do
        begin
            Buf[I,J] := Spline1DCalc(C, AP_Double(J)/(NewWidth-1));
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Vertical interpolation
    //
    J:=0;
    while J<=NewWidth-1 do
    begin
        
        //
        // Fill X, Y
        //
        I:=0;
        while I<=OldHeight-1 do
        begin
            X[I] := AP_Double(I)/(OldHeight-1);
            Y[I] := Buf[I,J];
            Inc(I);
        end;
        
        //
        // Interpolate and place result into B
        //
        Spline1DBuildCubic(X, Y, OldHeight, 0, Double(0.0), 0, Double(0.0), C);
        I:=0;
        while I<=NewHeight-1 do
        begin
            B[I,J] := Spline1DCalc(C, AP_Double(I)/(NewHeight-1));
            Inc(I);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Bilinear spline resampling

Input parameters:
    A           -   function values at the old grid,
                    array[0..OldHeight-1, 0..OldWidth-1]
    OldHeight   -   old grid height, OldHeight>1
    OldWidth    -   old grid width, OldWidth>1
    NewHeight   -   new grid height, NewHeight>1
    NewWidth    -   new grid width, NewWidth>1

Output parameters:
    B           -   function values at the new grid,
                    array[0..NewHeight-1, 0..NewWidth-1]

  -- ALGLIB routine --
     09.07.2007
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure Spline2DResampleBilinear(const A : TReal2DArray;
     OldHeight : AlglibInteger;
     OldWidth : AlglibInteger;
     var B : TReal2DArray;
     NewHeight : AlglibInteger;
     NewWidth : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    L : AlglibInteger;
    C : AlglibInteger;
    T : Double;
    U : Double;
begin
    SetLength(B, NewHeight-1+1, NewWidth-1+1);
    I:=0;
    while I<=NewHeight-1 do
    begin
        J:=0;
        while J<=NewWidth-1 do
        begin
            L := I*(OldHeight-1) div (NewHeight-1);
            if L=OldHeight-1 then
            begin
                L := OldHeight-2;
            end;
            U := AP_Double(I)/(NewHeight-1)*(OldHeight-1)-L;
            C := J*(OldWidth-1) div (NewWidth-1);
            if C=OldWidth-1 then
            begin
                C := OldWidth-2;
            end;
            T := AP_Double(J*(OldWidth-1))/(NewWidth-1)-C;
            B[I,J] := (1-T)*(1-U)*A[L,C]+T*(1-U)*A[L,C+1]+T*U*A[L+1,C+1]+(1-T)*U*A[L+1,C];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal subroutine.
Calculation of the first derivatives and the cross-derivative.
*************************************************************************)
procedure BicubicCalcDerivatives(const A : TReal2DArray;
     const X : TReal1DArray;
     const Y : TReal1DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var DX : TReal2DArray;
     var DY : TReal2DArray;
     var DXY : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    XT : TReal1DArray;
    FT : TReal1DArray;
    C : TReal1DArray;
    S : Double;
    DS : Double;
    D2S : Double;
begin
    SetLength(DX, M-1+1, N-1+1);
    SetLength(DY, M-1+1, N-1+1);
    SetLength(DXY, M-1+1, N-1+1);
    
    //
    // dF/dX
    //
    SetLength(XT, N-1+1);
    SetLength(FT, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            XT[J] := X[J];
            FT[J] := A[I,J];
            Inc(J);
        end;
        BuildCubicSpline(XT, FT, N, 0, Double(0.0), 0, Double(0.0), C);
        J:=0;
        while J<=N-1 do
        begin
            SplineDifferentiation(C, X[J], S, DS, D2S);
            DX[I,J] := DS;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // dF/dY
    //
    SetLength(XT, M-1+1);
    SetLength(FT, M-1+1);
    J:=0;
    while J<=N-1 do
    begin
        I:=0;
        while I<=M-1 do
        begin
            XT[I] := Y[I];
            FT[I] := A[I,J];
            Inc(I);
        end;
        BuildCubicSpline(XT, FT, M, 0, Double(0.0), 0, Double(0.0), C);
        I:=0;
        while I<=M-1 do
        begin
            SplineDifferentiation(C, Y[I], S, DS, D2S);
            DY[I,J] := DS;
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // d2F/dXdY
    //
    SetLength(XT, N-1+1);
    SetLength(FT, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            XT[J] := X[J];
            FT[J] := DY[I,J];
            Inc(J);
        end;
        BuildCubicSpline(XT, FT, N, 0, Double(0.0), 0, Double(0.0), C);
        J:=0;
        while J<=N-1 do
        begin
            SplineDifferentiation(C, X[J], S, DS, D2S);
            DXY[I,J] := DS;
            Inc(J);
        end;
        Inc(I);
    end;
end;


end.