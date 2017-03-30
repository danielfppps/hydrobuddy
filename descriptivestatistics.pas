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
unit descriptivestatistics;
interface
uses Math, Sysutils, Ap;

procedure CalculateMoments(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var Variance : Double;
     var Skewness : Double;
     var Kurtosis : Double);
procedure CalculateADev(const X : TReal1DArray;
     N : AlglibInteger;
     var ADev : Double);
procedure CalculateMedian(X : TReal1DArray;
     N : AlglibInteger;
     var Median : Double);
procedure CalculatePercentile(X : TReal1DArray;
     N : AlglibInteger;
     P : Double;
     var V : Double);

implementation

procedure InternalStatHeapSort(var Arr : TReal1DArray;
     N : AlglibInteger);forward;


(*************************************************************************
Calculation of the distribution moments: mean, variance, slewness, kurtosis.

Input parameters:
    X       -   sample. Array with whose indexes range within [0..N-1]
    N       -   sample size.
    
Output parameters:
    Mean    -   mean.
    Variance-   variance.
    Skewness-   skewness (if variance<>0; zero otherwise).
    Kurtosis-   kurtosis (if variance<>0; zero otherwise).

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure CalculateMoments(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var Variance : Double;
     var Skewness : Double;
     var Kurtosis : Double);
var
    I : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    StdDev : Double;
begin
    Mean := 0;
    Variance := 0;
    Skewness := 0;
    Kurtosis := 0;
    StdDev := 0;
    if N<=0 then
    begin
        Exit;
    end;
    
    //
    // Mean
    //
    I:=0;
    while I<=N-1 do
    begin
        Mean := Mean+X[I];
        Inc(I);
    end;
    Mean := Mean/N;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    if N<>1 then
    begin
        V1 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V1 := V1+AP_Sqr(X[I]-Mean);
            Inc(I);
        end;
        V2 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V2 := V2+(X[I]-Mean);
            Inc(I);
        end;
        V2 := AP_Sqr(V2)/N;
        Variance := (V1-V2)/(N-1);
        if AP_FP_Less(Variance,0) then
        begin
            Variance := 0;
        end;
        StdDev := Sqrt(Variance);
    end;
    
    //
    // Skewness and kurtosis
    //
    if AP_FP_Neq(StdDev,0) then
    begin
        I:=0;
        while I<=N-1 do
        begin
            V := (X[I]-Mean)/StdDev;
            V2 := AP_Sqr(V);
            Skewness := Skewness+V2*V;
            Kurtosis := Kurtosis+AP_Sqr(V2);
            Inc(I);
        end;
        Skewness := Skewness/N;
        Kurtosis := Kurtosis/N-3;
    end;
end;


(*************************************************************************
ADev

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   sample size
    
Output parameters:
    ADev-   ADev

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure CalculateADev(const X : TReal1DArray;
     N : AlglibInteger;
     var ADev : Double);
var
    I : AlglibInteger;
    Mean : Double;
begin
    Mean := 0;
    ADev := 0;
    if N<=0 then
    begin
        Exit;
    end;
    
    //
    // Mean
    //
    I:=0;
    while I<=N-1 do
    begin
        Mean := Mean+X[I];
        Inc(I);
    end;
    Mean := Mean/N;
    
    //
    // ADev
    //
    I:=0;
    while I<=N-1 do
    begin
        ADev := ADev+AbsReal(X[I]-Mean);
        Inc(I);
    end;
    ADev := ADev/N;
end;


(*************************************************************************
Median calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   sample size

Output parameters:
    Median

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure CalculateMedian(X : TReal1DArray;
     N : AlglibInteger;
     var Median : Double);
var
    i : AlglibInteger;
    ir : AlglibInteger;
    j : AlglibInteger;
    l : AlglibInteger;
    midp : AlglibInteger;
    K : AlglibInteger;
    a : Double;
    tval : Double;
begin
    X := DynamicArrayCopy(X);
    
    //
    // Some degenerate cases
    //
    Median := 0;
    if N<=0 then
    begin
        Exit;
    end;
    if N=1 then
    begin
        Median := X[0];
        Exit;
    end;
    if N=2 then
    begin
        Median := Double(0.5)*(X[0]+X[1]);
        Exit;
    end;
    
    //
    // Common case, N>=3.
    // Choose X[(N-1)/2]
    //
    l := 0;
    ir := n-1;
    k := (N-1) div 2;
    while True do
    begin
        if ir<=l+1 then
        begin
            
            //
            // 1 or 2 elements in partition
            //
            if (ir=l+1) and AP_FP_Less(x[ir],x[l]) then
            begin
                tval := x[l];
                x[l] := x[ir];
                x[ir] := tval;
            end;
            Break;
        end
        else
        begin
            midp := (l+ir) div 2;
            tval := x[midp];
            x[midp] := x[l+1];
            x[l+1] := tval;
            if AP_FP_Greater(x[l],x[ir]) then
            begin
                tval := x[l];
                x[l] := x[ir];
                x[ir] := tval;
            end;
            if AP_FP_Greater(x[l+1],x[ir]) then
            begin
                tval := x[l+1];
                x[l+1] := x[ir];
                x[ir] := tval;
            end;
            if AP_FP_Greater(x[l],x[l+1]) then
            begin
                tval := x[l];
                x[l] := x[l+1];
                x[l+1] := tval;
            end;
            i := l+1;
            j := ir;
            a := x[l+1];
            while True do
            begin
                repeat
                    i := i+1;
                until AP_FP_Greater_Eq(x[i],a);
                repeat
                    j := j-1;
                until AP_FP_Less_Eq(x[j],a);
                if j<i then
                begin
                    Break;
                end;
                tval := x[i];
                x[i] := x[j];
                x[j] := tval;
            end;
            x[l+1] := x[j];
            x[j] := a;
            if j>=k then
            begin
                ir := j-1;
            end;
            if j<=k then
            begin
                l := i;
            end;
        end;
    end;
    
    //
    // If N is odd, return result
    //
    if N mod 2=1 then
    begin
        Median := X[k];
        Exit;
    end;
    a := x[n-1];
    i:=k+1;
    while i<=n-1 do
    begin
        if AP_FP_Less(x[i],a) then
        begin
            a := x[i];
        end;
        Inc(i);
    end;
    Median := Double(0.5)*(x[k]+a);
end;


(*************************************************************************
Percentile calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   sample size, N>1
    P   -   percentile (0<=P<=1)

Output parameters:
    V   -   percentile

  -- ALGLIB --
     Copyright 01.03.2008 by Bochkanov Sergey
*************************************************************************)
procedure CalculatePercentile(X : TReal1DArray;
     N : AlglibInteger;
     P : Double;
     var V : Double);
var
    I1 : AlglibInteger;
    T : Double;
begin
    X := DynamicArrayCopy(X);
    Assert(N>1, 'CalculatePercentile: N<=1!');
    Assert(AP_FP_Greater_Eq(P,0) and AP_FP_Less_Eq(P,1), 'CalculatePercentile: incorrect P!');
    InternalStatHeapSort(X, N);
    if AP_FP_Eq(P,0) then
    begin
        V := X[0];
        Exit;
    end;
    if AP_FP_Eq(P,1) then
    begin
        V := X[N-1];
        Exit;
    end;
    T := P*(N-1);
    I1 := Floor(T);
    T := T-Floor(T);
    V := X[I1]*(1-T)+X[I1+1]*T;
end;


procedure InternalStatHeapSort(var Arr : TReal1DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
begin
    if N=1 then
    begin
        Exit;
    end;
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(Arr[k-1],Arr[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := Arr[k-1];
                Arr[k-1] := Arr[t-1];
                Arr[t-1] := Tmp;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := Arr[i];
        Arr[i] := Arr[0];
        Arr[0] := Tmp;
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
                    if AP_FP_Greater(Arr[k],Arr[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(Arr[t-1],Arr[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := Arr[k-1];
                    Arr[k-1] := Arr[t-1];
                    Arr[t-1] := Tmp;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


end.