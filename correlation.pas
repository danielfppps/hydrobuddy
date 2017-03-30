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
unit correlation;
interface
uses Math, Sysutils, Ap;

function PearsonCorrelation(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger):Double;
function SpearmanRankCorrelation(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger):Double;

implementation

procedure RankX(var X : TReal1DArray; N : AlglibInteger);forward;


(*************************************************************************
Pearson product-moment correlation coefficient

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   sample size.

Result:
    Pearson product-moment correlation coefficient

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************)
function PearsonCorrelation(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    XMean : Double;
    YMean : Double;
    S : Double;
    XV : Double;
    YV : Double;
    T1 : Double;
    T2 : Double;
begin
    XV := 0;
    YV := 0;
    if N<=1 then
    begin
        Result := 0;
        Exit;
    end;
    
    //
    // Mean
    //
    XMean := 0;
    YMean := 0;
    I:=0;
    while I<=N-1 do
    begin
        XMean := XMean+X[I];
        YMean := YMean+Y[I];
        Inc(I);
    end;
    XMean := XMean/N;
    YMean := YMean/N;
    
    //
    // numerator and denominator
    //
    S := 0;
    XV := 0;
    YV := 0;
    I:=0;
    while I<=N-1 do
    begin
        T1 := X[I]-XMean;
        T2 := Y[I]-YMean;
        XV := XV+AP_Sqr(T1);
        YV := YV+AP_Sqr(T2);
        S := S+T1*T2;
        Inc(I);
    end;
    if AP_FP_Eq(XV,0) or AP_FP_Eq(YV,0) then
    begin
        Result := 0;
    end
    else
    begin
        Result := S/(Sqrt(XV)*Sqrt(YV));
    end;
end;


(*************************************************************************
Spearman's rank correlation coefficient

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   sample size.

Result:
    Spearman's rank correlation coefficient

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************)
function SpearmanRankCorrelation(X : TReal1DArray;
     Y : TReal1DArray;
     N : AlglibInteger):Double;
begin
    X := DynamicArrayCopy(X);
    Y := DynamicArrayCopy(Y);
    RankX(X, N);
    RankX(Y, N);
    Result := PearsonCorrelation(X, Y, N);
end;


(*************************************************************************
Internal ranking subroutine
*************************************************************************)
procedure RankX(var X : TReal1DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    TmpI : AlglibInteger;
    R : TReal1DArray;
    C : TInteger1DArray;
begin
    
    //
    // Prepare
    //
    if N<=1 then
    begin
        X[0] := 1;
        Exit;
    end;
    SetLength(R, N-1+1);
    SetLength(C, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        R[I] := X[I];
        C[I] := I;
        Inc(I);
    end;
    
    //
    // sort {R, C}
    //
    if N<>1 then
    begin
        i := 2;
        repeat
            t := i;
            while t<>1 do
            begin
                k := t div 2;
                if AP_FP_Greater_Eq(R[k-1],R[t-1]) then
                begin
                    t := 1;
                end
                else
                begin
                    Tmp := R[k-1];
                    R[k-1] := R[t-1];
                    R[t-1] := Tmp;
                    TmpI := C[k-1];
                    C[k-1] := C[t-1];
                    C[t-1] := TmpI;
                    t := k;
                end;
            end;
            i := i+1;
        until  not (i<=N);
        i := N-1;
        repeat
            Tmp := R[i];
            R[i] := R[0];
            R[0] := Tmp;
            TmpI := C[i];
            C[i] := C[0];
            C[0] := TmpI;
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
                        if AP_FP_Greater(R[k],R[k-1]) then
                        begin
                            k := k+1;
                        end;
                    end;
                    if AP_FP_Greater_Eq(R[t-1],R[k-1]) then
                    begin
                        t := 0;
                    end
                    else
                    begin
                        Tmp := R[k-1];
                        R[k-1] := R[t-1];
                        R[t-1] := Tmp;
                        TmpI := C[k-1];
                        C[k-1] := C[t-1];
                        C[t-1] := TmpI;
                        t := k;
                    end;
                end;
            end;
            i := i-1;
        until  not (i>=1);
    end;
    
    //
    // compute tied ranks
    //
    I := 0;
    while I<=N-1 do
    begin
        J := I+1;
        while J<=N-1 do
        begin
            if AP_FP_Neq(R[J],R[I]) then
            begin
                Break;
            end;
            J := J+1;
        end;
        K:=I;
        while K<=J-1 do
        begin
            R[K] := 1+AP_Double((I+J-1))/2;
            Inc(K);
        end;
        I := J;
    end;
    
    //
    // back to x
    //
    I:=0;
    while I<=N-1 do
    begin
        X[C[I]] := R[I];
        Inc(I);
    end;
end;


end.