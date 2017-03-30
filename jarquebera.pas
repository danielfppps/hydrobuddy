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
unit jarquebera;
interface
uses Math, Sysutils, Ap;

procedure JarqueBeraTest(const X : TReal1DArray;
     N : AlglibInteger;
     var P : Double);

implementation

procedure JarqueBeraStatistic(const X : TReal1DArray;
     N : AlglibInteger;
     var S : Double);forward;
function JarqueBeraApprox(N : AlglibInteger; S : Double):Double;forward;
function JBTbl5(S : Double):Double;forward;
function JBTbl6(S : Double):Double;forward;
function JBTbl7(S : Double):Double;forward;
function JBTbl8(S : Double):Double;forward;
function JBTbl9(S : Double):Double;forward;
function JBTbl10(S : Double):Double;forward;
function JBTbl11(S : Double):Double;forward;
function JBTbl12(S : Double):Double;forward;
function JBTbl13(S : Double):Double;forward;
function JBTbl14(S : Double):Double;forward;
function JBTbl15(S : Double):Double;forward;
function JBTbl16(S : Double):Double;forward;
function JBTbl17(S : Double):Double;forward;
function JBTbl18(S : Double):Double;forward;
function JBTbl19(S : Double):Double;forward;
function JBTbl20(S : Double):Double;forward;
function JBTbl30(S : Double):Double;forward;
function JBTbl50(S : Double):Double;forward;
function JBTbl65(S : Double):Double;forward;
function JBTbl100(S : Double):Double;forward;
function JBTbl130(S : Double):Double;forward;
function JBTbl200(S : Double):Double;forward;
function JBTbl301(S : Double):Double;forward;
function JBTbl501(S : Double):Double;forward;
function JBTbl701(S : Double):Double;forward;
function JBTbl1401(S : Double):Double;forward;
procedure JBCheb(X : Double;
     C : Double;
     var TJ : Double;
     var TJ1 : Double;
     var R : Double);forward;


(*************************************************************************
Jarque-Bera test

This test checks hypotheses about the fact that a  given  sample  X  is  a
sample of normal random variable.

Requirements:
    * the number of elements in the sample is not less than 5.

Input parameters:
    X   -   sample. Array whose index goes from 0 to N-1.
    N   -   size of the sample. N>=5

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

Accuracy of the approximation used (5<=N<=1951):

p-value  	    relative error (5<=N<=1951)
[1, 0.1]            < 1%
[0.1, 0.01]         < 2%
[0.01, 0.001]       < 6%
[0.001, 0]          wasn't measured

For N>1951 accuracy wasn't measured but it shouldn't be sharply  different
from table values.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************)
procedure JarqueBeraTest(const X : TReal1DArray;
     N : AlglibInteger;
     var P : Double);
var
    S : Double;
begin
    
    //
    // N is too small
    //
    if N<5 then
    begin
        P := Double(1.0);
        Exit;
    end;
    
    //
    // N is large enough
    //
    JarqueBeraStatistic(X, N, S);
    P := JarqueBeraApprox(N, S);
end;


procedure JarqueBeraStatistic(const X : TReal1DArray;
     N : AlglibInteger;
     var S : Double);
var
    I : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    StdDev : Double;
    Mean : Double;
    Variance : Double;
    Skewness : Double;
    Kurtosis : Double;
begin
    Mean := 0;
    Variance := 0;
    Skewness := 0;
    Kurtosis := 0;
    StdDev := 0;
    Assert(N>1);
    
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
    
    //
    // Statistic
    //
    S := AP_Double(N)/6*(AP_Sqr(Skewness)+AP_Sqr(Kurtosis)/4);
end;


function JarqueBeraApprox(N : AlglibInteger; S : Double):Double;
var
    VX : TReal1DArray;
    VY : TReal1DArray;
    CTbl : TReal2DArray;
    T1 : Double;
    T2 : Double;
    T3 : Double;
    T : Double;
    F1 : Double;
    F2 : Double;
    F3 : Double;
    F12 : Double;
    F23 : Double;
    X : Double;
begin
    Result := 1;
    X := S;
    if N<5 then
    begin
        Exit;
    end;
    
    //
    // N = 5..20 are tabulated
    //
    if (N>=5) and (N<=20) then
    begin
        if N=5 then
        begin
            Result := Exp(JBTbl5(X));
        end;
        if N=6 then
        begin
            Result := Exp(JBTbl6(X));
        end;
        if N=7 then
        begin
            Result := Exp(JBTbl7(X));
        end;
        if N=8 then
        begin
            Result := Exp(JBTbl8(X));
        end;
        if N=9 then
        begin
            Result := Exp(JBTbl9(X));
        end;
        if N=10 then
        begin
            Result := Exp(JBTbl10(X));
        end;
        if N=11 then
        begin
            Result := Exp(JBTbl11(X));
        end;
        if N=12 then
        begin
            Result := Exp(JBTbl12(X));
        end;
        if N=13 then
        begin
            Result := Exp(JBTbl13(X));
        end;
        if N=14 then
        begin
            Result := Exp(JBTbl14(X));
        end;
        if N=15 then
        begin
            Result := Exp(JBTbl15(X));
        end;
        if N=16 then
        begin
            Result := Exp(JBTbl16(X));
        end;
        if N=17 then
        begin
            Result := Exp(JBTbl17(X));
        end;
        if N=18 then
        begin
            Result := Exp(JBTbl18(X));
        end;
        if N=19 then
        begin
            Result := Exp(JBTbl19(X));
        end;
        if N=20 then
        begin
            Result := Exp(JBTbl20(X));
        end;
        Exit;
    end;
    
    //
    // N = 20, 30, 50 are tabulated.
    // In-between values are interpolated
    // using interpolating polynomial of the second degree.
    //
    if (N>20) and (N<=50) then
    begin
        T1 := -Double(1.0)/Double(20.0);
        T2 := -Double(1.0)/Double(30.0);
        T3 := -Double(1.0)/Double(50.0);
        T := -Double(1.0)/N;
        F1 := JBTbl20(X);
        F2 := JBTbl30(X);
        F3 := JBTbl50(X);
        F12 := ((T-T2)*F1+(T1-T)*F2)/(T1-T2);
        F23 := ((T-T3)*F2+(T2-T)*F3)/(T2-T3);
        Result := ((T-T3)*F12+(T1-T)*F23)/(T1-T3);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Result := Exp(Result);
        Exit;
    end;
    
    //
    // N = 50, 65, 100 are tabulated.
    // In-between values are interpolated
    // using interpolating polynomial of the second degree.
    //
    if (N>50) and (N<=100) then
    begin
        T1 := -Double(1.0)/Double(50.0);
        T2 := -Double(1.0)/Double(65.0);
        T3 := -Double(1.0)/Double(100.0);
        T := -Double(1.0)/N;
        F1 := JBTbl50(X);
        F2 := JBTbl65(X);
        F3 := JBTbl100(X);
        F12 := ((T-T2)*F1+(T1-T)*F2)/(T1-T2);
        F23 := ((T-T3)*F2+(T2-T)*F3)/(T2-T3);
        Result := ((T-T3)*F12+(T1-T)*F23)/(T1-T3);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Result := Exp(Result);
        Exit;
    end;
    
    //
    // N = 100, 130, 200 are tabulated.
    // In-between values are interpolated
    // using interpolating polynomial of the second degree.
    //
    if (N>100) and (N<=200) then
    begin
        T1 := -Double(1.0)/Double(100.0);
        T2 := -Double(1.0)/Double(130.0);
        T3 := -Double(1.0)/Double(200.0);
        T := -Double(1.0)/N;
        F1 := JBTbl100(X);
        F2 := JBTbl130(X);
        F3 := JBTbl200(X);
        F12 := ((T-T2)*F1+(T1-T)*F2)/(T1-T2);
        F23 := ((T-T3)*F2+(T2-T)*F3)/(T2-T3);
        Result := ((T-T3)*F12+(T1-T)*F23)/(T1-T3);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Result := Exp(Result);
        Exit;
    end;
    
    //
    // N = 200, 301, 501 are tabulated.
    // In-between values are interpolated
    // using interpolating polynomial of the second degree.
    //
    if (N>200) and (N<=501) then
    begin
        T1 := -Double(1.0)/Double(200.0);
        T2 := -Double(1.0)/Double(301.0);
        T3 := -Double(1.0)/Double(501.0);
        T := -Double(1.0)/N;
        F1 := JBTbl200(X);
        F2 := JBTbl301(X);
        F3 := JBTbl501(X);
        F12 := ((T-T2)*F1+(T1-T)*F2)/(T1-T2);
        F23 := ((T-T3)*F2+(T2-T)*F3)/(T2-T3);
        Result := ((T-T3)*F12+(T1-T)*F23)/(T1-T3);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Result := Exp(Result);
        Exit;
    end;
    
    //
    // N = 501, 701, 1401 are tabulated.
    // In-between values are interpolated
    // using interpolating polynomial of the second degree.
    //
    if (N>501) and (N<=1401) then
    begin
        T1 := -Double(1.0)/Double(501.0);
        T2 := -Double(1.0)/Double(701.0);
        T3 := -Double(1.0)/Double(1401.0);
        T := -Double(1.0)/N;
        F1 := JBTbl501(X);
        F2 := JBTbl701(X);
        F3 := JBTbl1401(X);
        F12 := ((T-T2)*F1+(T1-T)*F2)/(T1-T2);
        F23 := ((T-T3)*F2+(T2-T)*F3)/(T2-T3);
        Result := ((T-T3)*F12+(T1-T)*F23)/(T1-T3);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Result := Exp(Result);
        Exit;
    end;
    
    //
    // Asymptotic expansion
    //
    if N>1401 then
    begin
        Result := -Double(0.5)*X+(JBTbl1401(X)+Double(0.5)*X)*Sqrt(AP_Double(1401)/N);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Result := Exp(Result);
        Exit;
    end;
end;


function JBTbl5(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(0.4000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(0.400000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.097885e-20), TJ, TJ1, Result);
        JBCheb(X, -Double(2.854501e-20), TJ, TJ1, Result);
        JBCheb(X, -Double(1.756616e-20), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(1.1000)) then
    begin
        X := 2*(S-Double(0.400000))/Double(0.700000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.324545e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.075941e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(9.772272e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.175686e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.576162e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.126861e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.434425e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.790359e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.809178e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.479704e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.717040e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.294170e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.880632e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(3.023344e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.601531e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(7.920403e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(5.188419e+02)*(S-Double(1.100000e+00))-Double(4.767297e+00);
end;


function JBTbl6(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(0.2500)) then
    begin
        X := 2*(S-Double(0.000000))/Double(0.250000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.274707e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(5.700471e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(3.425764e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(1.3000)) then
    begin
        X := 2*(S-Double(0.250000))/Double(1.050000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.339000e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.011104e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.168177e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.085666e-01), TJ, TJ1, Result);
        JBCheb(X, Double(7.738606e-02), TJ, TJ1, Result);
        JBCheb(X, Double(7.022876e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.462402e-02), TJ, TJ1, Result);
        JBCheb(X, Double(6.908270e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(8.230772e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.006996e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.410222e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.893768e-03), TJ, TJ1, Result);
        JBCheb(X, Double(8.114564e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(1.8500)) then
    begin
        X := 2*(S-Double(1.300000))/Double(0.550000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.794311e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(3.578700e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.394664e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.928290e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.813273e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.076063e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.835380e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.013013e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(5.058903e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.856915e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(6.710887e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.770029e+02)*(S-Double(1.850000e+00))-Double(1.371015e+01);
end;


function JBTbl7(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.4000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.400000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.093681e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.695911e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.473192e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.203236e-01), TJ, TJ1, Result);
        JBCheb(X, Double(6.590379e-02), TJ, TJ1, Result);
        JBCheb(X, Double(6.291876e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.132007e-02), TJ, TJ1, Result);
        JBCheb(X, Double(9.411147e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.180067e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(3.487610e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.436561e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(1.400000))/Double(1.600000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.947854e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.772675e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(4.707912e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.691171e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.132795e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.481310e-02), TJ, TJ1, Result);
        JBCheb(X, Double(2.867536e-03), TJ, TJ1, Result);
        JBCheb(X, Double(8.772327e-04), TJ, TJ1, Result);
        JBCheb(X, Double(5.033387e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.378277e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.497964e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(3.636814e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(9.581640e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(3.2000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(0.200000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(7.511008e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.140472e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.682053e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.568561e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.933930e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.140472e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.895025e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.140472e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.933930e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.568561e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.682053e+00), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.824116e+03)*(S-Double(3.200000e+00))-Double(1.440330e+01);
end;


function JBTbl8(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.3000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.300000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(7.199015e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.095921e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(4.736828e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.047438e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.484320e-03), TJ, TJ1, Result);
        JBCheb(X, Double(7.937923e-03), TJ, TJ1, Result);
        JBCheb(X, Double(4.810470e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.139780e-03), TJ, TJ1, Result);
        JBCheb(X, Double(6.708443e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(2.0000)) then
    begin
        X := 2*(S-Double(1.300000))/Double(0.700000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(3.378966e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.802461e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.547593e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(6.241042e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.203274e-02), TJ, TJ1, Result);
        JBCheb(X, Double(5.201990e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(5.125597e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.584426e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.546069e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(5.0000)) then
    begin
        X := 2*(S-Double(2.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.828366e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(3.137533e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(5.016671e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.745637e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(5.189801e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.621610e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(6.741122e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.516368e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.552085e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.787029e-03), TJ, TJ1, Result);
        JBCheb(X, Double(5.359774e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(5.087028e+00)*(S-Double(5.000000e+00))-Double(1.071300e+01);
end;


function JBTbl9(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.3000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.300000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.279320e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(9.277151e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.669339e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(7.086149e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.333816e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.871249e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.007048e-03), TJ, TJ1, Result);
        JBCheb(X, Double(7.482245e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.355615e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(2.0000)) then
    begin
        X := 2*(S-Double(1.300000))/Double(0.700000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.981430e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.972248e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.747737e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.808530e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(7.888305e-03), TJ, TJ1, Result);
        JBCheb(X, Double(9.001302e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.378767e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.108510e-03), TJ, TJ1, Result);
        JBCheb(X, Double(5.915372e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(7.0000)) then
    begin
        X := 2*(S-Double(2.000000))/Double(5.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.387463e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.845231e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.809956e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(7.543461e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.880397e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.160074e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(7.356527e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.394428e-03), TJ, TJ1, Result);
        JBCheb(X, Double(9.619892e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(2.758763e-04), TJ, TJ1, Result);
        JBCheb(X, Double(4.790977e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.020952e+00)*(S-Double(7.000000e+00))-Double(9.516623e+00);
end;


function JBTbl10(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.2000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.200000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.590993e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(6.562730e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.353934e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.069933e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.849151e-03), TJ, TJ1, Result);
        JBCheb(X, Double(8.931406e-04), TJ, TJ1, Result);
        JBCheb(X, Double(3.636295e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.178340e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(8.917749e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(2.0000)) then
    begin
        X := 2*(S-Double(1.200000))/Double(0.800000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.537658e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(9.962401e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.838715e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.055792e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.580316e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.781701e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.770362e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.838983e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(6.999052e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(7.0000)) then
    begin
        X := 2*(S-Double(2.000000))/Double(5.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.337524e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.877029e+00), TJ, TJ1, Result);
        JBCheb(X, Double(4.734650e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.249254e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.320250e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(6.432266e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(8.711035e-01)*(S-Double(7.000000e+00))-Double(7.212811e+00);
end;


function JBTbl11(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.2000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.200000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.339517e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(6.051558e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.000992e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.022547e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(9.808401e-04), TJ, TJ1, Result);
        JBCheb(X, Double(5.592870e-04), TJ, TJ1, Result);
        JBCheb(X, Double(3.575081e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.086173e-04), TJ, TJ1, Result);
        JBCheb(X, Double(6.089011e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(2.2500)) then
    begin
        X := 2*(S-Double(1.200000))/Double(1.050000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.523221e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.068388e+00), TJ, TJ1, Result);
        JBCheb(X, Double(2.179661e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.555524e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(3.238964e-02), TJ, TJ1, Result);
        JBCheb(X, Double(7.364320e-03), TJ, TJ1, Result);
        JBCheb(X, Double(4.895771e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.762774e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(8.201340e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(8.0000)) then
    begin
        X := 2*(S-Double(2.250000))/Double(5.750000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.212179e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.684579e+00), TJ, TJ1, Result);
        JBCheb(X, Double(8.299519e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(3.606261e-02), TJ, TJ1, Result);
        JBCheb(X, Double(7.310869e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(3.320115e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(5.715445e-01)*(S-Double(8.000000e+00))-Double(6.845834e+00);
end;


function JBTbl12(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.736742e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.657836e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.047209e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.319599e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.545631e-04), TJ, TJ1, Result);
        JBCheb(X, Double(9.280445e-05), TJ, TJ1, Result);
        JBCheb(X, Double(2.815679e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(2.213519e-05), TJ, TJ1, Result);
        JBCheb(X, Double(1.256838e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(1.000000))/Double(2.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.573947e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.515287e+00), TJ, TJ1, Result);
        JBCheb(X, Double(3.611880e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.271311e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(6.495815e-02), TJ, TJ1, Result);
        JBCheb(X, Double(4.141186e-02), TJ, TJ1, Result);
        JBCheb(X, Double(7.180886e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.388211e-02), TJ, TJ1, Result);
        JBCheb(X, Double(4.890761e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.233175e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.946156e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(12.0000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(9.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.947819e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.034157e+00), TJ, TJ1, Result);
        JBCheb(X, Double(6.878986e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.078603e-02), TJ, TJ1, Result);
        JBCheb(X, Double(6.990977e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.866215e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.897866e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.512252e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.073743e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.022621e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.501343e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.877243e-01)*(S-Double(1.200000e+01))-Double(7.936839e+00);
end;


function JBTbl13(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.713276e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.557541e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(9.459092e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.044145e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.546132e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.002374e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.349456e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(7.025669e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(1.590242e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(1.000000))/Double(2.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.454383e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.467539e+00), TJ, TJ1, Result);
        JBCheb(X, Double(3.270774e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(8.075763e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(6.611647e-02), TJ, TJ1, Result);
        JBCheb(X, Double(2.990785e-02), TJ, TJ1, Result);
        JBCheb(X, Double(8.109212e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.135031e-02), TJ, TJ1, Result);
        JBCheb(X, Double(5.915919e-04), TJ, TJ1, Result);
        JBCheb(X, Double(3.522390e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.144701e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(13.0000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.736127e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.920809e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.175858e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.002049e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.158966e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(3.157781e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.762172e-03), TJ, TJ1, Result);
        JBCheb(X, Double(5.780347e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.193310e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.442421e-05), TJ, TJ1, Result);
        JBCheb(X, Double(2.547756e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.799944e-01)*(S-Double(1.300000e+01))-Double(7.566269e+00);
end;


function JBTbl14(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(1.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(1.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.698527e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.479081e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(8.640733e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(8.466899e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.469485e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.150009e-05), TJ, TJ1, Result);
        JBCheb(X, Double(1.965975e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(4.710210e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(1.327808e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(1.000000))/Double(2.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(2.350359e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.421365e+00), TJ, TJ1, Result);
        JBCheb(X, Double(2.960468e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.149167e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(6.361109e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.976022e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.082700e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(8.563328e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.453123e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.917559e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.151067e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(12.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.746892e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.010441e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.566146e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(5.129690e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.929724e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.524227e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.192933e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.254730e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.620685e-03), TJ, TJ1, Result);
        JBCheb(X, Double(7.289618e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(2.112350e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.590621e-01)*(S-Double(1.500000e+01))-Double(7.632238e+00);
end;


function JBTbl15(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(2.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(2.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.043660e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.361653e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(3.009497e-01), TJ, TJ1, Result);
        JBCheb(X, Double(4.951784e-02), TJ, TJ1, Result);
        JBCheb(X, Double(4.377903e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.003253e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.271309e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(5.0000)) then
    begin
        X := 2*(S-Double(2.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(3.582778e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.349578e-01), TJ, TJ1, Result);
        JBCheb(X, Double(9.476514e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.717385e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.222591e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(6.635124e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.815993e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(17.0000)) then
    begin
        X := 2*(S-Double(5.000000))/Double(12.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.115476e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.655936e+00), TJ, TJ1, Result);
        JBCheb(X, Double(8.404310e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.663794e-02), TJ, TJ1, Result);
        JBCheb(X, Double(8.868618e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.381447e-03), TJ, TJ1, Result);
        JBCheb(X, Double(9.444801e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.581503e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(9.468696e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.728509e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.206470e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.927937e-01)*(S-Double(1.700000e+01))-Double(7.700983e+00);
end;


function JBTbl16(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(2.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(2.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.002570e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.298141e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.832803e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.877026e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.539436e-02), TJ, TJ1, Result);
        JBCheb(X, Double(8.439658e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.756911e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(5.0000)) then
    begin
        X := 2*(S-Double(2.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(3.486198e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.242944e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.020002e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.130531e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.512373e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(8.054876e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.556839e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(20.0000)) then
    begin
        X := 2*(S-Double(5.000000))/Double(15.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.241608e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.832655e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.340545e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.361143e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.283219e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.484549e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.805968e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.057243e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.454439e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.177513e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.819209e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.391580e-01)*(S-Double(2.000000e+01))-Double(7.963205e+00);
end;


function JBTbl17(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.566973e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.810330e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(4.840039e-02), TJ, TJ1, Result);
        JBCheb(X, Double(2.337294e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(5.383549e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(5.556515e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(8.656965e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.404569e-02), TJ, TJ1, Result);
        JBCheb(X, Double(6.447867e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(6.0000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(3.905684e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.222920e-01), TJ, TJ1, Result);
        JBCheb(X, Double(4.146667e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.809176e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.057028e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.211838e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(4.099683e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.161105e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.225465e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(24.0000)) then
    begin
        X := 2*(S-Double(6.000000))/Double(18.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.594282e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.917838e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.455980e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.999589e-02), TJ, TJ1, Result);
        JBCheb(X, Double(5.604263e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(3.484445e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.819937e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.930390e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.771761e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(6.232581e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(7.029083e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.127771e-01)*(S-Double(2.400000e+01))-Double(8.400197e+00);
end;


function JBTbl18(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.526802e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.762373e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(5.598890e-02), TJ, TJ1, Result);
        JBCheb(X, Double(2.189437e-01), TJ, TJ1, Result);
        JBCheb(X, Double(5.971721e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.823067e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.064501e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.014932e-02), TJ, TJ1, Result);
        JBCheb(X, Double(5.953513e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(6.0000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(3.818669e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.070918e-01), TJ, TJ1, Result);
        JBCheb(X, Double(4.277196e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.879817e-03), TJ, TJ1, Result);
        JBCheb(X, Double(6.887357e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.638451e-05), TJ, TJ1, Result);
        JBCheb(X, Double(1.502800e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(3.165796e-05), TJ, TJ1, Result);
        JBCheb(X, Double(5.034960e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(20.0000)) then
    begin
        X := 2*(S-Double(6.000000))/Double(14.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.010656e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.496296e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.002227e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.338250e-02), TJ, TJ1, Result);
        JBCheb(X, Double(4.137036e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.586202e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(9.736384e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.332251e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.877982e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.160963e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(2.547247e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.684623e-01)*(S-Double(2.000000e+01))-Double(7.428883e+00);
end;


function JBTbl19(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(3.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.490213e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.719633e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.459123e-02), TJ, TJ1, Result);
        JBCheb(X, Double(2.034878e-01), TJ, TJ1, Result);
        JBCheb(X, Double(1.113868e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.030922e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.054022e-02), TJ, TJ1, Result);
        JBCheb(X, Double(7.525623e-03), TJ, TJ1, Result);
        JBCheb(X, Double(5.277360e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(6.0000)) then
    begin
        X := 2*(S-Double(3.000000))/Double(3.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(3.744750e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(5.977749e-01), TJ, TJ1, Result);
        JBCheb(X, Double(4.223716e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.363889e-03), TJ, TJ1, Result);
        JBCheb(X, Double(5.711774e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(5.557257e-04), TJ, TJ1, Result);
        JBCheb(X, Double(4.254794e-04), TJ, TJ1, Result);
        JBCheb(X, Double(9.034207e-05), TJ, TJ1, Result);
        JBCheb(X, Double(5.498107e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(20.0000)) then
    begin
        X := 2*(S-Double(6.000000))/Double(14.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.872768e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.430689e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.136575e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.726627e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.421110e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.581510e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(5.559520e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(6.838208e-04), TJ, TJ1, Result);
        JBCheb(X, Double(8.428839e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(7.170682e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(6.006647e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.539373e-01)*(S-Double(2.000000e+01))-Double(7.206941e+00);
end;


function JBTbl20(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.854794e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.948947e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.632184e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.139397e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(1.006237e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.810031e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.573620e-02), TJ, TJ1, Result);
        JBCheb(X, Double(9.951242e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.274092e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(3.464196e-03), TJ, TJ1, Result);
        JBCheb(X, Double(4.882139e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.575144e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.822804e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(7.061348e-04), TJ, TJ1, Result);
        JBCheb(X, Double(5.908404e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.978353e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.030989e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.327151e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.346404e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.840051e-02), TJ, TJ1, Result);
        JBCheb(X, Double(7.578551e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(9.813886e-04), TJ, TJ1, Result);
        JBCheb(X, Double(5.905973e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(5.358489e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(3.450795e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(6.941157e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(7.432418e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(2.070537e-04), TJ, TJ1, Result);
        JBCheb(X, Double(9.375654e-04), TJ, TJ1, Result);
        JBCheb(X, Double(5.367378e-04), TJ, TJ1, Result);
        JBCheb(X, Double(9.890859e-04), TJ, TJ1, Result);
        JBCheb(X, Double(6.679782e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(7.015854e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.487737e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.244254e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.318007e-01)*(S-Double(2.500000e+01))-Double(7.742185e+00);
end;


function JBTbl30(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.630822e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.724298e+00), TJ, TJ1, Result);
        JBCheb(X, Double(7.872756e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.658268e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.573597e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.994157e-02), TJ, TJ1, Result);
        JBCheb(X, Double(5.994825e-03), TJ, TJ1, Result);
        JBCheb(X, Double(7.394303e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(5.785029e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.990264e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.037838e-04), TJ, TJ1, Result);
        JBCheb(X, Double(6.755546e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.774473e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(2.821395e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.392603e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.353313e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.539322e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.197018e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.396848e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.804293e-02), TJ, TJ1, Result);
        JBCheb(X, Double(6.867928e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(2.768758e-03), TJ, TJ1, Result);
        JBCheb(X, Double(5.211792e-04), TJ, TJ1, Result);
        JBCheb(X, Double(4.925799e-04), TJ, TJ1, Result);
        JBCheb(X, Double(5.046235e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(9.536469e-05), TJ, TJ1, Result);
        JBCheb(X, -Double(6.489642e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.263462e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.177316e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.590637e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.028212e-01)*(S-Double(2.500000e+01))-Double(6.855288e+00);
end;


function JBTbl50(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.436279e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.519711e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.148699e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.001204e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.207620e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.034778e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(1.220322e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.033260e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.588280e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.851653e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.287733e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.234645e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.189127e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.429738e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.058822e-02), TJ, TJ1, Result);
        JBCheb(X, Double(9.086776e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.445783e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.311671e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(7.261298e-04), TJ, TJ1, Result);
        JBCheb(X, Double(6.496987e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.605249e-04), TJ, TJ1, Result);
        JBCheb(X, Double(8.162282e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.921095e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(5.888603e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.080113e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(9.313116e-02)*(S-Double(2.500000e+01))-Double(6.479154e+00);
end;


function JBTbl65(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.360024e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.434631e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.514580e-03), TJ, TJ1, Result);
        JBCheb(X, Double(7.332038e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.158197e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(5.121233e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.051056e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.148601e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.214233e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.487977e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.424720e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.116715e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(4.043152e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.718149e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.313701e-03), TJ, TJ1, Result);
        JBCheb(X, Double(3.097305e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.181031e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.256975e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.858951e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(5.895179e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.933237e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(9.443768e-02)*(S-Double(2.500000e+01))-Double(6.419137e+00);
end;


function JBTbl100(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.257021e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.313418e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.628931e-02), TJ, TJ1, Result);
        JBCheb(X, Double(4.264287e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.518487e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.499826e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.836044e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.056508e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.279690e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.665746e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.290012e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.487632e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.704465e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.211669e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.866099e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.399767e-01), TJ, TJ1, Result);
        JBCheb(X, Double(2.498208e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.080097e-01)*(S-Double(2.500000e+01))-Double(6.481094e+00);
end;


function JBTbl130(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.207999e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.253864e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.618032e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.112729e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.210546e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(4.732602e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(2.410527e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.026324e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.331990e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.779129e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.674749e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.669077e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(5.679136e-03), TJ, TJ1, Result);
        JBCheb(X, Double(8.833221e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(5.893951e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(6.475304e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.116734e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.045722e-01)*(S-Double(2.500000e+01))-Double(6.510314e+00);
end;


function JBTbl200(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.146155e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.177398e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.297970e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.869745e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.717288e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(1.982108e-04), TJ, TJ1, Result);
        JBCheb(X, Double(6.427636e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.034235e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.455006e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.942996e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.973795e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.418812e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(3.156778e-03), TJ, TJ1, Result);
        JBCheb(X, Double(4.896705e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.086071e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.152176e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.725393e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.132404e-01)*(S-Double(2.500000e+01))-Double(6.764034e+00);
end;


function JBTbl301(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.104290e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.125800e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(9.595847e-03), TJ, TJ1, Result);
        JBCheb(X, Double(1.219666e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.502210e-04), TJ, TJ1, Result);
        JBCheb(X, -Double(6.414543e-05), TJ, TJ1, Result);
        JBCheb(X, Double(6.754115e-05), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.065955e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.582060e+00), TJ, TJ1, Result);
        JBCheb(X, Double(2.004472e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(4.709092e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.105779e-02), TJ, TJ1, Result);
        JBCheb(X, Double(1.197391e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(8.386780e-04), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.311384e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(7.918763e-01), TJ, TJ1, Result);
        JBCheb(X, Double(3.626584e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.293626e-01)*(S-Double(2.500000e+01))-Double(7.066995e+00);
end;


function JBTbl501(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.067426e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.079765e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(5.463005e-03), TJ, TJ1, Result);
        JBCheb(X, Double(6.875659e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.127574e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.740694e+00), TJ, TJ1, Result);
        JBCheb(X, Double(2.044502e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(3.746714e-02), TJ, TJ1, Result);
        JBCheb(X, Double(3.810594e-04), TJ, TJ1, Result);
        JBCheb(X, Double(1.197111e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.628194e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(8.846221e-01), TJ, TJ1, Result);
        JBCheb(X, Double(4.386405e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.418332e-01)*(S-Double(2.500000e+01))-Double(7.468952e+00);
end;


function JBTbl701(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.050999e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.059769e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(3.922680e-03), TJ, TJ1, Result);
        JBCheb(X, Double(4.847054e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.192182e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.860007e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.963942e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(2.838711e-02), TJ, TJ1, Result);
        JBCheb(X, -Double(2.893112e-04), TJ, TJ1, Result);
        JBCheb(X, Double(2.159788e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(6.917851e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(9.817020e-01), TJ, TJ1, Result);
        JBCheb(X, Double(5.383727e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(1.532706e-01)*(S-Double(2.500000e+01))-Double(7.845715e+00);
end;


function JBTbl1401(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    if AP_FP_Less_Eq(S,Double(4.0000)) then
    begin
        X := 2*(S-Double(0.000000))/Double(4.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(1.026266e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.030061e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.259222e-03), TJ, TJ1, Result);
        JBCheb(X, Double(2.536254e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(15.0000)) then
    begin
        X := 2*(S-Double(4.000000))/Double(11.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(4.329849e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(2.095443e+00), TJ, TJ1, Result);
        JBCheb(X, Double(1.759363e-01), TJ, TJ1, Result);
        JBCheb(X, -Double(7.751359e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(6.124368e-03), TJ, TJ1, Result);
        JBCheb(X, -Double(1.793114e-03), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(S,Double(25.0000)) then
    begin
        X := 2*(S-Double(15.000000))/Double(10.000000)-1;
        TJ := 1;
        TJ1 := X;
        JBCheb(X, -Double(7.544330e+00), TJ, TJ1, Result);
        JBCheb(X, -Double(1.225382e+00), TJ, TJ1, Result);
        JBCheb(X, Double(5.392349e-02), TJ, TJ1, Result);
        if AP_FP_Greater(Result,0) then
        begin
            Result := 0;
        end;
        Exit;
    end;
    Result := -Double(2.019375e-01)*(S-Double(2.500000e+01))-Double(8.715788e+00);
end;


procedure JBCheb(X : Double;
     C : Double;
     var TJ : Double;
     var TJ1 : Double;
     var R : Double);
var
    T : Double;
begin
    R := R+C*TJ;
    T := 2*X*TJ1-TJ;
    TJ := TJ1;
    TJ1 := T;
end;


end.