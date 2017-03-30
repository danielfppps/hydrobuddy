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
unit wsr;
interface
uses Math, Sysutils, Ap;

procedure WilcoxonSignedRankTest(X : TReal1DArray;
     N : AlglibInteger;
     E : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);

implementation

procedure WCheb(X : Double;
     C : Double;
     var TJ : Double;
     var TJ1 : Double;
     var R : Double);forward;
function W5(S : Double):Double;forward;
function W6(S : Double):Double;forward;
function W7(S : Double):Double;forward;
function W8(S : Double):Double;forward;
function W9(S : Double):Double;forward;
function W10(S : Double):Double;forward;
function W11(S : Double):Double;forward;
function W12(S : Double):Double;forward;
function W13(S : Double):Double;forward;
function W14(S : Double):Double;forward;
function W15(S : Double):Double;forward;
function W16(S : Double):Double;forward;
function W17(S : Double):Double;forward;
function W18(S : Double):Double;forward;
function W19(S : Double):Double;forward;
function W20(S : Double):Double;forward;
function W21(S : Double):Double;forward;
function W22(S : Double):Double;forward;
function W23(S : Double):Double;forward;
function W24(S : Double):Double;forward;
function W25(S : Double):Double;forward;
function W26(S : Double):Double;forward;
function W27(S : Double):Double;forward;
function W28(S : Double):Double;forward;
function W29(S : Double):Double;forward;
function W30(S : Double):Double;forward;
function W40(S : Double):Double;forward;
function W60(S : Double):Double;forward;
function W120(S : Double):Double;forward;
function W200(S : Double):Double;forward;
function WSigma(S : Double; N : AlglibInteger):Double;forward;


(*************************************************************************
Wilcoxon signed-rank test

This test checks three hypotheses about the median  of  the  given sample.
The following tests are performed:
    * two-tailed test (null hypothesis - the median is equal to the  given
      value)
    * left-tailed test (null hypothesis - the median is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis  -  the  median  is  less than or
      equal to the given value)

Requirements:
    * the scale of measurement should be ordinal, interval or  ratio (i.e.
      the test could not be applied to nominal variables).
    * the distribution should be continuous and symmetric relative to  its
      median.
    * number of distinct values in the X array should be greater than 4

The test is non-parametric and doesn't require distribution X to be normal

Input parameters:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of the sample.
    Median  -   assumed median value.

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

To calculate p-values, special approximation is used. This method lets  us
calculate p-values with two decimal places in interval [0.0001, 1].

"Two decimal places" does not sound very impressive, but in  practice  the
relative error of less than 1% is enough to make a decision.

There is no approximation outside the [0.0001, 1] interval. Therefore,  if
the significance level outlies this interval, the test returns 0.0001.

  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure WilcoxonSignedRankTest(X : TReal1DArray;
     N : AlglibInteger;
     E : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    TmpI : AlglibInteger;
    NS : AlglibInteger;
    R : TReal1DArray;
    C : TInteger1DArray;
    W : Double;
    P : Double;
    MP : Double;
    S : Double;
    Sigma : Double;
    Mu : Double;
begin
    X := DynamicArrayCopy(X);
    
    //
    // Prepare
    //
    if N<5 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    NS := 0;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Eq(X[I],E) then
        begin
            Inc(I);
            Continue;
        end;
        X[NS] := X[I];
        NS := NS+1;
        Inc(I);
    end;
    if NS<5 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    SetLength(R, NS-1+1);
    SetLength(C, NS-1+1);
    I:=0;
    while I<=NS-1 do
    begin
        R[I] := AbsReal(X[I]-E);
        C[I] := I;
        Inc(I);
    end;
    
    //
    // sort {R, C}
    //
    if NS<>1 then
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
        until  not (i<=NS);
        i := NS-1;
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
    while I<=NS-1 do
    begin
        J := I+1;
        while J<=NS-1 do
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
    // Compute W+
    //
    W := 0;
    I:=0;
    while I<=NS-1 do
    begin
        if AP_FP_Greater(X[C[I]],E) then
        begin
            W := W+R[I];
        end;
        Inc(I);
    end;
    
    //
    // Result
    //
    Mu := AP_Double(NS*(NS+1))/4;
    Sigma := Sqrt(AP_Double(NS*(NS+1)*(2*NS+1))/24);
    S := (W-Mu)/Sigma;
    if AP_FP_Less_Eq(S,0) then
    begin
        P := Exp(WSigma(-(W-Mu)/Sigma, NS));
        MP := 1-Exp(WSigma(-(W-1-Mu)/Sigma, NS));
    end
    else
    begin
        MP := Exp(WSigma((W-Mu)/Sigma, NS));
        P := 1-Exp(WSigma((W+1-Mu)/Sigma, NS));
    end;
    BothTails := Max(2*Min(P, MP), Double(1.0E-4));
    LeftTail := Max(P, Double(1.0E-4));
    RightTail := Max(MP, Double(1.0E-4));
end;


(*************************************************************************
Sequential Chebyshev interpolation.
*************************************************************************)
procedure WCheb(X : Double;
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


(*************************************************************************
Tail(S, 5)
*************************************************************************)
function W5(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(3.708099e+00)*S+Double(7.500000e+00));
    if W>=7 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=6 then
    begin
        R := -Double(9.008e-01);
    end;
    if W=5 then
    begin
        R := -Double(1.163e+00);
    end;
    if W=4 then
    begin
        R := -Double(1.520e+00);
    end;
    if W=3 then
    begin
        R := -Double(1.856e+00);
    end;
    if W=2 then
    begin
        R := -Double(2.367e+00);
    end;
    if W=1 then
    begin
        R := -Double(2.773e+00);
    end;
    if W<=0 then
    begin
        R := -Double(3.466e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 6)
*************************************************************************)
function W6(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(4.769696e+00)*S+Double(1.050000e+01));
    if W>=10 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=9 then
    begin
        R := -Double(8.630e-01);
    end;
    if W=8 then
    begin
        R := -Double(1.068e+00);
    end;
    if W=7 then
    begin
        R := -Double(1.269e+00);
    end;
    if W=6 then
    begin
        R := -Double(1.520e+00);
    end;
    if W=5 then
    begin
        R := -Double(1.856e+00);
    end;
    if W=4 then
    begin
        R := -Double(2.213e+00);
    end;
    if W=3 then
    begin
        R := -Double(2.549e+00);
    end;
    if W=2 then
    begin
        R := -Double(3.060e+00);
    end;
    if W=1 then
    begin
        R := -Double(3.466e+00);
    end;
    if W<=0 then
    begin
        R := -Double(4.159e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 7)
*************************************************************************)
function W7(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(5.916080e+00)*S+Double(1.400000e+01));
    if W>=14 then
    begin
        R := -Double(6.325e-01);
    end;
    if W=13 then
    begin
        R := -Double(7.577e-01);
    end;
    if W=12 then
    begin
        R := -Double(9.008e-01);
    end;
    if W=11 then
    begin
        R := -Double(1.068e+00);
    end;
    if W=10 then
    begin
        R := -Double(1.241e+00);
    end;
    if W=9 then
    begin
        R := -Double(1.451e+00);
    end;
    if W=8 then
    begin
        R := -Double(1.674e+00);
    end;
    if W=7 then
    begin
        R := -Double(1.908e+00);
    end;
    if W=6 then
    begin
        R := -Double(2.213e+00);
    end;
    if W=5 then
    begin
        R := -Double(2.549e+00);
    end;
    if W=4 then
    begin
        R := -Double(2.906e+00);
    end;
    if W=3 then
    begin
        R := -Double(3.243e+00);
    end;
    if W=2 then
    begin
        R := -Double(3.753e+00);
    end;
    if W=1 then
    begin
        R := -Double(4.159e+00);
    end;
    if W<=0 then
    begin
        R := -Double(4.852e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 8)
*************************************************************************)
function W8(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(7.141428e+00)*S+Double(1.800000e+01));
    if W>=18 then
    begin
        R := -Double(6.399e-01);
    end;
    if W=17 then
    begin
        R := -Double(7.494e-01);
    end;
    if W=16 then
    begin
        R := -Double(8.630e-01);
    end;
    if W=15 then
    begin
        R := -Double(9.913e-01);
    end;
    if W=14 then
    begin
        R := -Double(1.138e+00);
    end;
    if W=13 then
    begin
        R := -Double(1.297e+00);
    end;
    if W=12 then
    begin
        R := -Double(1.468e+00);
    end;
    if W=11 then
    begin
        R := -Double(1.653e+00);
    end;
    if W=10 then
    begin
        R := -Double(1.856e+00);
    end;
    if W=9 then
    begin
        R := -Double(2.079e+00);
    end;
    if W=8 then
    begin
        R := -Double(2.326e+00);
    end;
    if W=7 then
    begin
        R := -Double(2.601e+00);
    end;
    if W=6 then
    begin
        R := -Double(2.906e+00);
    end;
    if W=5 then
    begin
        R := -Double(3.243e+00);
    end;
    if W=4 then
    begin
        R := -Double(3.599e+00);
    end;
    if W=3 then
    begin
        R := -Double(3.936e+00);
    end;
    if W=2 then
    begin
        R := -Double(4.447e+00);
    end;
    if W=1 then
    begin
        R := -Double(4.852e+00);
    end;
    if W<=0 then
    begin
        R := -Double(5.545e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 9)
*************************************************************************)
function W9(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(8.440972e+00)*S+Double(2.250000e+01));
    if W>=22 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=21 then
    begin
        R := -Double(7.873e-01);
    end;
    if W=20 then
    begin
        R := -Double(8.912e-01);
    end;
    if W=19 then
    begin
        R := -Double(1.002e+00);
    end;
    if W=18 then
    begin
        R := -Double(1.120e+00);
    end;
    if W=17 then
    begin
        R := -Double(1.255e+00);
    end;
    if W=16 then
    begin
        R := -Double(1.394e+00);
    end;
    if W=15 then
    begin
        R := -Double(1.547e+00);
    end;
    if W=14 then
    begin
        R := -Double(1.717e+00);
    end;
    if W=13 then
    begin
        R := -Double(1.895e+00);
    end;
    if W=12 then
    begin
        R := -Double(2.079e+00);
    end;
    if W=11 then
    begin
        R := -Double(2.287e+00);
    end;
    if W=10 then
    begin
        R := -Double(2.501e+00);
    end;
    if W=9 then
    begin
        R := -Double(2.742e+00);
    end;
    if W=8 then
    begin
        R := -Double(3.019e+00);
    end;
    if W=7 then
    begin
        R := -Double(3.294e+00);
    end;
    if W=6 then
    begin
        R := -Double(3.599e+00);
    end;
    if W=5 then
    begin
        R := -Double(3.936e+00);
    end;
    if W=4 then
    begin
        R := -Double(4.292e+00);
    end;
    if W=3 then
    begin
        R := -Double(4.629e+00);
    end;
    if W=2 then
    begin
        R := -Double(5.140e+00);
    end;
    if W=1 then
    begin
        R := -Double(5.545e+00);
    end;
    if W<=0 then
    begin
        R := -Double(6.238e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 10)
*************************************************************************)
function W10(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(9.810708e+00)*S+Double(2.750000e+01));
    if W>=27 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=26 then
    begin
        R := -Double(7.745e-01);
    end;
    if W=25 then
    begin
        R := -Double(8.607e-01);
    end;
    if W=24 then
    begin
        R := -Double(9.551e-01);
    end;
    if W=23 then
    begin
        R := -Double(1.057e+00);
    end;
    if W=22 then
    begin
        R := -Double(1.163e+00);
    end;
    if W=21 then
    begin
        R := -Double(1.279e+00);
    end;
    if W=20 then
    begin
        R := -Double(1.402e+00);
    end;
    if W=19 then
    begin
        R := -Double(1.533e+00);
    end;
    if W=18 then
    begin
        R := -Double(1.674e+00);
    end;
    if W=17 then
    begin
        R := -Double(1.826e+00);
    end;
    if W=16 then
    begin
        R := -Double(1.983e+00);
    end;
    if W=15 then
    begin
        R := -Double(2.152e+00);
    end;
    if W=14 then
    begin
        R := -Double(2.336e+00);
    end;
    if W=13 then
    begin
        R := -Double(2.525e+00);
    end;
    if W=12 then
    begin
        R := -Double(2.727e+00);
    end;
    if W=11 then
    begin
        R := -Double(2.942e+00);
    end;
    if W=10 then
    begin
        R := -Double(3.170e+00);
    end;
    if W=9 then
    begin
        R := -Double(3.435e+00);
    end;
    if W=8 then
    begin
        R := -Double(3.713e+00);
    end;
    if W=7 then
    begin
        R := -Double(3.987e+00);
    end;
    if W=6 then
    begin
        R := -Double(4.292e+00);
    end;
    if W=5 then
    begin
        R := -Double(4.629e+00);
    end;
    if W=4 then
    begin
        R := -Double(4.986e+00);
    end;
    if W=3 then
    begin
        R := -Double(5.322e+00);
    end;
    if W=2 then
    begin
        R := -Double(5.833e+00);
    end;
    if W=1 then
    begin
        R := -Double(6.238e+00);
    end;
    if W<=0 then
    begin
        R := -Double(6.931e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 11)
*************************************************************************)
function W11(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(1.124722e+01)*S+Double(3.300000e+01));
    if W>=33 then
    begin
        R := -Double(6.595e-01);
    end;
    if W=32 then
    begin
        R := -Double(7.279e-01);
    end;
    if W=31 then
    begin
        R := -Double(8.002e-01);
    end;
    if W=30 then
    begin
        R := -Double(8.782e-01);
    end;
    if W=29 then
    begin
        R := -Double(9.615e-01);
    end;
    if W=28 then
    begin
        R := -Double(1.050e+00);
    end;
    if W=27 then
    begin
        R := -Double(1.143e+00);
    end;
    if W=26 then
    begin
        R := -Double(1.243e+00);
    end;
    if W=25 then
    begin
        R := -Double(1.348e+00);
    end;
    if W=24 then
    begin
        R := -Double(1.459e+00);
    end;
    if W=23 then
    begin
        R := -Double(1.577e+00);
    end;
    if W=22 then
    begin
        R := -Double(1.700e+00);
    end;
    if W=21 then
    begin
        R := -Double(1.832e+00);
    end;
    if W=20 then
    begin
        R := -Double(1.972e+00);
    end;
    if W=19 then
    begin
        R := -Double(2.119e+00);
    end;
    if W=18 then
    begin
        R := -Double(2.273e+00);
    end;
    if W=17 then
    begin
        R := -Double(2.437e+00);
    end;
    if W=16 then
    begin
        R := -Double(2.607e+00);
    end;
    if W=15 then
    begin
        R := -Double(2.788e+00);
    end;
    if W=14 then
    begin
        R := -Double(2.980e+00);
    end;
    if W=13 then
    begin
        R := -Double(3.182e+00);
    end;
    if W=12 then
    begin
        R := -Double(3.391e+00);
    end;
    if W=11 then
    begin
        R := -Double(3.617e+00);
    end;
    if W=10 then
    begin
        R := -Double(3.863e+00);
    end;
    if W=9 then
    begin
        R := -Double(4.128e+00);
    end;
    if W=8 then
    begin
        R := -Double(4.406e+00);
    end;
    if W=7 then
    begin
        R := -Double(4.680e+00);
    end;
    if W=6 then
    begin
        R := -Double(4.986e+00);
    end;
    if W=5 then
    begin
        R := -Double(5.322e+00);
    end;
    if W=4 then
    begin
        R := -Double(5.679e+00);
    end;
    if W=3 then
    begin
        R := -Double(6.015e+00);
    end;
    if W=2 then
    begin
        R := -Double(6.526e+00);
    end;
    if W=1 then
    begin
        R := -Double(6.931e+00);
    end;
    if W<=0 then
    begin
        R := -Double(7.625e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 12)
*************************************************************************)
function W12(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(1.274755e+01)*S+Double(3.900000e+01));
    if W>=39 then
    begin
        R := -Double(6.633e-01);
    end;
    if W=38 then
    begin
        R := -Double(7.239e-01);
    end;
    if W=37 then
    begin
        R := -Double(7.878e-01);
    end;
    if W=36 then
    begin
        R := -Double(8.556e-01);
    end;
    if W=35 then
    begin
        R := -Double(9.276e-01);
    end;
    if W=34 then
    begin
        R := -Double(1.003e+00);
    end;
    if W=33 then
    begin
        R := -Double(1.083e+00);
    end;
    if W=32 then
    begin
        R := -Double(1.168e+00);
    end;
    if W=31 then
    begin
        R := -Double(1.256e+00);
    end;
    if W=30 then
    begin
        R := -Double(1.350e+00);
    end;
    if W=29 then
    begin
        R := -Double(1.449e+00);
    end;
    if W=28 then
    begin
        R := -Double(1.552e+00);
    end;
    if W=27 then
    begin
        R := -Double(1.660e+00);
    end;
    if W=26 then
    begin
        R := -Double(1.774e+00);
    end;
    if W=25 then
    begin
        R := -Double(1.893e+00);
    end;
    if W=24 then
    begin
        R := -Double(2.017e+00);
    end;
    if W=23 then
    begin
        R := -Double(2.148e+00);
    end;
    if W=22 then
    begin
        R := -Double(2.285e+00);
    end;
    if W=21 then
    begin
        R := -Double(2.429e+00);
    end;
    if W=20 then
    begin
        R := -Double(2.581e+00);
    end;
    if W=19 then
    begin
        R := -Double(2.738e+00);
    end;
    if W=18 then
    begin
        R := -Double(2.902e+00);
    end;
    if W=17 then
    begin
        R := -Double(3.076e+00);
    end;
    if W=16 then
    begin
        R := -Double(3.255e+00);
    end;
    if W=15 then
    begin
        R := -Double(3.443e+00);
    end;
    if W=14 then
    begin
        R := -Double(3.645e+00);
    end;
    if W=13 then
    begin
        R := -Double(3.852e+00);
    end;
    if W=12 then
    begin
        R := -Double(4.069e+00);
    end;
    if W=11 then
    begin
        R := -Double(4.310e+00);
    end;
    if W=10 then
    begin
        R := -Double(4.557e+00);
    end;
    if W=9 then
    begin
        R := -Double(4.821e+00);
    end;
    if W=8 then
    begin
        R := -Double(5.099e+00);
    end;
    if W=7 then
    begin
        R := -Double(5.373e+00);
    end;
    if W=6 then
    begin
        R := -Double(5.679e+00);
    end;
    if W=5 then
    begin
        R := -Double(6.015e+00);
    end;
    if W=4 then
    begin
        R := -Double(6.372e+00);
    end;
    if W=3 then
    begin
        R := -Double(6.708e+00);
    end;
    if W=2 then
    begin
        R := -Double(7.219e+00);
    end;
    if W=1 then
    begin
        R := -Double(7.625e+00);
    end;
    if W<=0 then
    begin
        R := -Double(8.318e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 13)
*************************************************************************)
function W13(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(1.430909e+01)*S+Double(4.550000e+01));
    if W>=45 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=44 then
    begin
        R := -Double(7.486e-01);
    end;
    if W=43 then
    begin
        R := -Double(8.068e-01);
    end;
    if W=42 then
    begin
        R := -Double(8.683e-01);
    end;
    if W=41 then
    begin
        R := -Double(9.328e-01);
    end;
    if W=40 then
    begin
        R := -Double(1.001e+00);
    end;
    if W=39 then
    begin
        R := -Double(1.072e+00);
    end;
    if W=38 then
    begin
        R := -Double(1.146e+00);
    end;
    if W=37 then
    begin
        R := -Double(1.224e+00);
    end;
    if W=36 then
    begin
        R := -Double(1.306e+00);
    end;
    if W=35 then
    begin
        R := -Double(1.392e+00);
    end;
    if W=34 then
    begin
        R := -Double(1.481e+00);
    end;
    if W=33 then
    begin
        R := -Double(1.574e+00);
    end;
    if W=32 then
    begin
        R := -Double(1.672e+00);
    end;
    if W=31 then
    begin
        R := -Double(1.773e+00);
    end;
    if W=30 then
    begin
        R := -Double(1.879e+00);
    end;
    if W=29 then
    begin
        R := -Double(1.990e+00);
    end;
    if W=28 then
    begin
        R := -Double(2.104e+00);
    end;
    if W=27 then
    begin
        R := -Double(2.224e+00);
    end;
    if W=26 then
    begin
        R := -Double(2.349e+00);
    end;
    if W=25 then
    begin
        R := -Double(2.479e+00);
    end;
    if W=24 then
    begin
        R := -Double(2.614e+00);
    end;
    if W=23 then
    begin
        R := -Double(2.755e+00);
    end;
    if W=22 then
    begin
        R := -Double(2.902e+00);
    end;
    if W=21 then
    begin
        R := -Double(3.055e+00);
    end;
    if W=20 then
    begin
        R := -Double(3.215e+00);
    end;
    if W=19 then
    begin
        R := -Double(3.380e+00);
    end;
    if W=18 then
    begin
        R := -Double(3.551e+00);
    end;
    if W=17 then
    begin
        R := -Double(3.733e+00);
    end;
    if W=16 then
    begin
        R := -Double(3.917e+00);
    end;
    if W=15 then
    begin
        R := -Double(4.113e+00);
    end;
    if W=14 then
    begin
        R := -Double(4.320e+00);
    end;
    if W=13 then
    begin
        R := -Double(4.534e+00);
    end;
    if W=12 then
    begin
        R := -Double(4.762e+00);
    end;
    if W=11 then
    begin
        R := -Double(5.004e+00);
    end;
    if W=10 then
    begin
        R := -Double(5.250e+00);
    end;
    if W=9 then
    begin
        R := -Double(5.514e+00);
    end;
    if W=8 then
    begin
        R := -Double(5.792e+00);
    end;
    if W=7 then
    begin
        R := -Double(6.066e+00);
    end;
    if W=6 then
    begin
        R := -Double(6.372e+00);
    end;
    if W=5 then
    begin
        R := -Double(6.708e+00);
    end;
    if W=4 then
    begin
        R := -Double(7.065e+00);
    end;
    if W=3 then
    begin
        R := -Double(7.401e+00);
    end;
    if W=2 then
    begin
        R := -Double(7.912e+00);
    end;
    if W=1 then
    begin
        R := -Double(8.318e+00);
    end;
    if W<=0 then
    begin
        R := -Double(9.011e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 14)
*************************************************************************)
function W14(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(1.592953e+01)*S+Double(5.250000e+01));
    if W>=52 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=51 then
    begin
        R := -Double(7.428e-01);
    end;
    if W=50 then
    begin
        R := -Double(7.950e-01);
    end;
    if W=49 then
    begin
        R := -Double(8.495e-01);
    end;
    if W=48 then
    begin
        R := -Double(9.067e-01);
    end;
    if W=47 then
    begin
        R := -Double(9.664e-01);
    end;
    if W=46 then
    begin
        R := -Double(1.029e+00);
    end;
    if W=45 then
    begin
        R := -Double(1.094e+00);
    end;
    if W=44 then
    begin
        R := -Double(1.162e+00);
    end;
    if W=43 then
    begin
        R := -Double(1.233e+00);
    end;
    if W=42 then
    begin
        R := -Double(1.306e+00);
    end;
    if W=41 then
    begin
        R := -Double(1.383e+00);
    end;
    if W=40 then
    begin
        R := -Double(1.463e+00);
    end;
    if W=39 then
    begin
        R := -Double(1.546e+00);
    end;
    if W=38 then
    begin
        R := -Double(1.632e+00);
    end;
    if W=37 then
    begin
        R := -Double(1.722e+00);
    end;
    if W=36 then
    begin
        R := -Double(1.815e+00);
    end;
    if W=35 then
    begin
        R := -Double(1.911e+00);
    end;
    if W=34 then
    begin
        R := -Double(2.011e+00);
    end;
    if W=33 then
    begin
        R := -Double(2.115e+00);
    end;
    if W=32 then
    begin
        R := -Double(2.223e+00);
    end;
    if W=31 then
    begin
        R := -Double(2.334e+00);
    end;
    if W=30 then
    begin
        R := -Double(2.450e+00);
    end;
    if W=29 then
    begin
        R := -Double(2.570e+00);
    end;
    if W=28 then
    begin
        R := -Double(2.694e+00);
    end;
    if W=27 then
    begin
        R := -Double(2.823e+00);
    end;
    if W=26 then
    begin
        R := -Double(2.956e+00);
    end;
    if W=25 then
    begin
        R := -Double(3.095e+00);
    end;
    if W=24 then
    begin
        R := -Double(3.238e+00);
    end;
    if W=23 then
    begin
        R := -Double(3.387e+00);
    end;
    if W=22 then
    begin
        R := -Double(3.541e+00);
    end;
    if W=21 then
    begin
        R := -Double(3.700e+00);
    end;
    if W=20 then
    begin
        R := -Double(3.866e+00);
    end;
    if W=19 then
    begin
        R := -Double(4.038e+00);
    end;
    if W=18 then
    begin
        R := -Double(4.215e+00);
    end;
    if W=17 then
    begin
        R := -Double(4.401e+00);
    end;
    if W=16 then
    begin
        R := -Double(4.592e+00);
    end;
    if W=15 then
    begin
        R := -Double(4.791e+00);
    end;
    if W=14 then
    begin
        R := -Double(5.004e+00);
    end;
    if W=13 then
    begin
        R := -Double(5.227e+00);
    end;
    if W=12 then
    begin
        R := -Double(5.456e+00);
    end;
    if W=11 then
    begin
        R := -Double(5.697e+00);
    end;
    if W=10 then
    begin
        R := -Double(5.943e+00);
    end;
    if W=9 then
    begin
        R := -Double(6.208e+00);
    end;
    if W=8 then
    begin
        R := -Double(6.485e+00);
    end;
    if W=7 then
    begin
        R := -Double(6.760e+00);
    end;
    if W=6 then
    begin
        R := -Double(7.065e+00);
    end;
    if W=5 then
    begin
        R := -Double(7.401e+00);
    end;
    if W=4 then
    begin
        R := -Double(7.758e+00);
    end;
    if W=3 then
    begin
        R := -Double(8.095e+00);
    end;
    if W=2 then
    begin
        R := -Double(8.605e+00);
    end;
    if W=1 then
    begin
        R := -Double(9.011e+00);
    end;
    if W<=0 then
    begin
        R := -Double(9.704e+00);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 15)
*************************************************************************)
function W15(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(1.760682e+01)*S+Double(6.000000e+01));
    if W>=60 then
    begin
        R := -Double(6.714e-01);
    end;
    if W=59 then
    begin
        R := -Double(7.154e-01);
    end;
    if W=58 then
    begin
        R := -Double(7.613e-01);
    end;
    if W=57 then
    begin
        R := -Double(8.093e-01);
    end;
    if W=56 then
    begin
        R := -Double(8.593e-01);
    end;
    if W=55 then
    begin
        R := -Double(9.114e-01);
    end;
    if W=54 then
    begin
        R := -Double(9.656e-01);
    end;
    if W=53 then
    begin
        R := -Double(1.022e+00);
    end;
    if W=52 then
    begin
        R := -Double(1.081e+00);
    end;
    if W=51 then
    begin
        R := -Double(1.142e+00);
    end;
    if W=50 then
    begin
        R := -Double(1.205e+00);
    end;
    if W=49 then
    begin
        R := -Double(1.270e+00);
    end;
    if W=48 then
    begin
        R := -Double(1.339e+00);
    end;
    if W=47 then
    begin
        R := -Double(1.409e+00);
    end;
    if W=46 then
    begin
        R := -Double(1.482e+00);
    end;
    if W=45 then
    begin
        R := -Double(1.558e+00);
    end;
    if W=44 then
    begin
        R := -Double(1.636e+00);
    end;
    if W=43 then
    begin
        R := -Double(1.717e+00);
    end;
    if W=42 then
    begin
        R := -Double(1.801e+00);
    end;
    if W=41 then
    begin
        R := -Double(1.888e+00);
    end;
    if W=40 then
    begin
        R := -Double(1.977e+00);
    end;
    if W=39 then
    begin
        R := -Double(2.070e+00);
    end;
    if W=38 then
    begin
        R := -Double(2.166e+00);
    end;
    if W=37 then
    begin
        R := -Double(2.265e+00);
    end;
    if W=36 then
    begin
        R := -Double(2.366e+00);
    end;
    if W=35 then
    begin
        R := -Double(2.472e+00);
    end;
    if W=34 then
    begin
        R := -Double(2.581e+00);
    end;
    if W=33 then
    begin
        R := -Double(2.693e+00);
    end;
    if W=32 then
    begin
        R := -Double(2.809e+00);
    end;
    if W=31 then
    begin
        R := -Double(2.928e+00);
    end;
    if W=30 then
    begin
        R := -Double(3.051e+00);
    end;
    if W=29 then
    begin
        R := -Double(3.179e+00);
    end;
    if W=28 then
    begin
        R := -Double(3.310e+00);
    end;
    if W=27 then
    begin
        R := -Double(3.446e+00);
    end;
    if W=26 then
    begin
        R := -Double(3.587e+00);
    end;
    if W=25 then
    begin
        R := -Double(3.732e+00);
    end;
    if W=24 then
    begin
        R := -Double(3.881e+00);
    end;
    if W=23 then
    begin
        R := -Double(4.036e+00);
    end;
    if W=22 then
    begin
        R := -Double(4.195e+00);
    end;
    if W=21 then
    begin
        R := -Double(4.359e+00);
    end;
    if W=20 then
    begin
        R := -Double(4.531e+00);
    end;
    if W=19 then
    begin
        R := -Double(4.707e+00);
    end;
    if W=18 then
    begin
        R := -Double(4.888e+00);
    end;
    if W=17 then
    begin
        R := -Double(5.079e+00);
    end;
    if W=16 then
    begin
        R := -Double(5.273e+00);
    end;
    if W=15 then
    begin
        R := -Double(5.477e+00);
    end;
    if W=14 then
    begin
        R := -Double(5.697e+00);
    end;
    if W=13 then
    begin
        R := -Double(5.920e+00);
    end;
    if W=12 then
    begin
        R := -Double(6.149e+00);
    end;
    if W=11 then
    begin
        R := -Double(6.390e+00);
    end;
    if W=10 then
    begin
        R := -Double(6.636e+00);
    end;
    if W=9 then
    begin
        R := -Double(6.901e+00);
    end;
    if W=8 then
    begin
        R := -Double(7.178e+00);
    end;
    if W=7 then
    begin
        R := -Double(7.453e+00);
    end;
    if W=6 then
    begin
        R := -Double(7.758e+00);
    end;
    if W=5 then
    begin
        R := -Double(8.095e+00);
    end;
    if W=4 then
    begin
        R := -Double(8.451e+00);
    end;
    if W=3 then
    begin
        R := -Double(8.788e+00);
    end;
    if W=2 then
    begin
        R := -Double(9.299e+00);
    end;
    if W=1 then
    begin
        R := -Double(9.704e+00);
    end;
    if W<=0 then
    begin
        R := -Double(1.040e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 16)
*************************************************************************)
function W16(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(1.933908e+01)*S+Double(6.800000e+01));
    if W>=68 then
    begin
        R := -Double(6.733e-01);
    end;
    if W=67 then
    begin
        R := -Double(7.134e-01);
    end;
    if W=66 then
    begin
        R := -Double(7.551e-01);
    end;
    if W=65 then
    begin
        R := -Double(7.986e-01);
    end;
    if W=64 then
    begin
        R := -Double(8.437e-01);
    end;
    if W=63 then
    begin
        R := -Double(8.905e-01);
    end;
    if W=62 then
    begin
        R := -Double(9.391e-01);
    end;
    if W=61 then
    begin
        R := -Double(9.895e-01);
    end;
    if W=60 then
    begin
        R := -Double(1.042e+00);
    end;
    if W=59 then
    begin
        R := -Double(1.096e+00);
    end;
    if W=58 then
    begin
        R := -Double(1.152e+00);
    end;
    if W=57 then
    begin
        R := -Double(1.210e+00);
    end;
    if W=56 then
    begin
        R := -Double(1.270e+00);
    end;
    if W=55 then
    begin
        R := -Double(1.331e+00);
    end;
    if W=54 then
    begin
        R := -Double(1.395e+00);
    end;
    if W=53 then
    begin
        R := -Double(1.462e+00);
    end;
    if W=52 then
    begin
        R := -Double(1.530e+00);
    end;
    if W=51 then
    begin
        R := -Double(1.600e+00);
    end;
    if W=50 then
    begin
        R := -Double(1.673e+00);
    end;
    if W=49 then
    begin
        R := -Double(1.748e+00);
    end;
    if W=48 then
    begin
        R := -Double(1.825e+00);
    end;
    if W=47 then
    begin
        R := -Double(1.904e+00);
    end;
    if W=46 then
    begin
        R := -Double(1.986e+00);
    end;
    if W=45 then
    begin
        R := -Double(2.071e+00);
    end;
    if W=44 then
    begin
        R := -Double(2.158e+00);
    end;
    if W=43 then
    begin
        R := -Double(2.247e+00);
    end;
    if W=42 then
    begin
        R := -Double(2.339e+00);
    end;
    if W=41 then
    begin
        R := -Double(2.434e+00);
    end;
    if W=40 then
    begin
        R := -Double(2.532e+00);
    end;
    if W=39 then
    begin
        R := -Double(2.632e+00);
    end;
    if W=38 then
    begin
        R := -Double(2.735e+00);
    end;
    if W=37 then
    begin
        R := -Double(2.842e+00);
    end;
    if W=36 then
    begin
        R := -Double(2.951e+00);
    end;
    if W=35 then
    begin
        R := -Double(3.064e+00);
    end;
    if W=34 then
    begin
        R := -Double(3.179e+00);
    end;
    if W=33 then
    begin
        R := -Double(3.298e+00);
    end;
    if W=32 then
    begin
        R := -Double(3.420e+00);
    end;
    if W=31 then
    begin
        R := -Double(3.546e+00);
    end;
    if W=30 then
    begin
        R := -Double(3.676e+00);
    end;
    if W=29 then
    begin
        R := -Double(3.810e+00);
    end;
    if W=28 then
    begin
        R := -Double(3.947e+00);
    end;
    if W=27 then
    begin
        R := -Double(4.088e+00);
    end;
    if W=26 then
    begin
        R := -Double(4.234e+00);
    end;
    if W=25 then
    begin
        R := -Double(4.383e+00);
    end;
    if W=24 then
    begin
        R := -Double(4.538e+00);
    end;
    if W=23 then
    begin
        R := -Double(4.697e+00);
    end;
    if W=22 then
    begin
        R := -Double(4.860e+00);
    end;
    if W=21 then
    begin
        R := -Double(5.029e+00);
    end;
    if W=20 then
    begin
        R := -Double(5.204e+00);
    end;
    if W=19 then
    begin
        R := -Double(5.383e+00);
    end;
    if W=18 then
    begin
        R := -Double(5.569e+00);
    end;
    if W=17 then
    begin
        R := -Double(5.762e+00);
    end;
    if W=16 then
    begin
        R := -Double(5.960e+00);
    end;
    if W=15 then
    begin
        R := -Double(6.170e+00);
    end;
    if W=14 then
    begin
        R := -Double(6.390e+00);
    end;
    if W=13 then
    begin
        R := -Double(6.613e+00);
    end;
    if W=12 then
    begin
        R := -Double(6.842e+00);
    end;
    if W=11 then
    begin
        R := -Double(7.083e+00);
    end;
    if W=10 then
    begin
        R := -Double(7.329e+00);
    end;
    if W=9 then
    begin
        R := -Double(7.594e+00);
    end;
    if W=8 then
    begin
        R := -Double(7.871e+00);
    end;
    if W=7 then
    begin
        R := -Double(8.146e+00);
    end;
    if W=6 then
    begin
        R := -Double(8.451e+00);
    end;
    if W=5 then
    begin
        R := -Double(8.788e+00);
    end;
    if W=4 then
    begin
        R := -Double(9.144e+00);
    end;
    if W=3 then
    begin
        R := -Double(9.481e+00);
    end;
    if W=2 then
    begin
        R := -Double(9.992e+00);
    end;
    if W=1 then
    begin
        R := -Double(1.040e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.109e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 17)
*************************************************************************)
function W17(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(2.112463e+01)*S+Double(7.650000e+01));
    if W>=76 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=75 then
    begin
        R := -Double(7.306e-01);
    end;
    if W=74 then
    begin
        R := -Double(7.695e-01);
    end;
    if W=73 then
    begin
        R := -Double(8.097e-01);
    end;
    if W=72 then
    begin
        R := -Double(8.514e-01);
    end;
    if W=71 then
    begin
        R := -Double(8.946e-01);
    end;
    if W=70 then
    begin
        R := -Double(9.392e-01);
    end;
    if W=69 then
    begin
        R := -Double(9.853e-01);
    end;
    if W=68 then
    begin
        R := -Double(1.033e+00);
    end;
    if W=67 then
    begin
        R := -Double(1.082e+00);
    end;
    if W=66 then
    begin
        R := -Double(1.133e+00);
    end;
    if W=65 then
    begin
        R := -Double(1.185e+00);
    end;
    if W=64 then
    begin
        R := -Double(1.240e+00);
    end;
    if W=63 then
    begin
        R := -Double(1.295e+00);
    end;
    if W=62 then
    begin
        R := -Double(1.353e+00);
    end;
    if W=61 then
    begin
        R := -Double(1.412e+00);
    end;
    if W=60 then
    begin
        R := -Double(1.473e+00);
    end;
    if W=59 then
    begin
        R := -Double(1.536e+00);
    end;
    if W=58 then
    begin
        R := -Double(1.600e+00);
    end;
    if W=57 then
    begin
        R := -Double(1.666e+00);
    end;
    if W=56 then
    begin
        R := -Double(1.735e+00);
    end;
    if W=55 then
    begin
        R := -Double(1.805e+00);
    end;
    if W=54 then
    begin
        R := -Double(1.877e+00);
    end;
    if W=53 then
    begin
        R := -Double(1.951e+00);
    end;
    if W=52 then
    begin
        R := -Double(2.028e+00);
    end;
    if W=51 then
    begin
        R := -Double(2.106e+00);
    end;
    if W=50 then
    begin
        R := -Double(2.186e+00);
    end;
    if W=49 then
    begin
        R := -Double(2.269e+00);
    end;
    if W=48 then
    begin
        R := -Double(2.353e+00);
    end;
    if W=47 then
    begin
        R := -Double(2.440e+00);
    end;
    if W=46 then
    begin
        R := -Double(2.530e+00);
    end;
    if W=45 then
    begin
        R := -Double(2.621e+00);
    end;
    if W=44 then
    begin
        R := -Double(2.715e+00);
    end;
    if W=43 then
    begin
        R := -Double(2.812e+00);
    end;
    if W=42 then
    begin
        R := -Double(2.911e+00);
    end;
    if W=41 then
    begin
        R := -Double(3.012e+00);
    end;
    if W=40 then
    begin
        R := -Double(3.116e+00);
    end;
    if W=39 then
    begin
        R := -Double(3.223e+00);
    end;
    if W=38 then
    begin
        R := -Double(3.332e+00);
    end;
    if W=37 then
    begin
        R := -Double(3.445e+00);
    end;
    if W=36 then
    begin
        R := -Double(3.560e+00);
    end;
    if W=35 then
    begin
        R := -Double(3.678e+00);
    end;
    if W=34 then
    begin
        R := -Double(3.799e+00);
    end;
    if W=33 then
    begin
        R := -Double(3.924e+00);
    end;
    if W=32 then
    begin
        R := -Double(4.052e+00);
    end;
    if W=31 then
    begin
        R := -Double(4.183e+00);
    end;
    if W=30 then
    begin
        R := -Double(4.317e+00);
    end;
    if W=29 then
    begin
        R := -Double(4.456e+00);
    end;
    if W=28 then
    begin
        R := -Double(4.597e+00);
    end;
    if W=27 then
    begin
        R := -Double(4.743e+00);
    end;
    if W=26 then
    begin
        R := -Double(4.893e+00);
    end;
    if W=25 then
    begin
        R := -Double(5.047e+00);
    end;
    if W=24 then
    begin
        R := -Double(5.204e+00);
    end;
    if W=23 then
    begin
        R := -Double(5.367e+00);
    end;
    if W=22 then
    begin
        R := -Double(5.534e+00);
    end;
    if W=21 then
    begin
        R := -Double(5.706e+00);
    end;
    if W=20 then
    begin
        R := -Double(5.884e+00);
    end;
    if W=19 then
    begin
        R := -Double(6.066e+00);
    end;
    if W=18 then
    begin
        R := -Double(6.254e+00);
    end;
    if W=17 then
    begin
        R := -Double(6.451e+00);
    end;
    if W=16 then
    begin
        R := -Double(6.654e+00);
    end;
    if W=15 then
    begin
        R := -Double(6.864e+00);
    end;
    if W=14 then
    begin
        R := -Double(7.083e+00);
    end;
    if W=13 then
    begin
        R := -Double(7.306e+00);
    end;
    if W=12 then
    begin
        R := -Double(7.535e+00);
    end;
    if W=11 then
    begin
        R := -Double(7.776e+00);
    end;
    if W=10 then
    begin
        R := -Double(8.022e+00);
    end;
    if W=9 then
    begin
        R := -Double(8.287e+00);
    end;
    if W=8 then
    begin
        R := -Double(8.565e+00);
    end;
    if W=7 then
    begin
        R := -Double(8.839e+00);
    end;
    if W=6 then
    begin
        R := -Double(9.144e+00);
    end;
    if W=5 then
    begin
        R := -Double(9.481e+00);
    end;
    if W=4 then
    begin
        R := -Double(9.838e+00);
    end;
    if W=3 then
    begin
        R := -Double(1.017e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.068e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.109e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.178e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 18)
*************************************************************************)
function W18(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(2.296193e+01)*S+Double(8.550000e+01));
    if W>=85 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=84 then
    begin
        R := -Double(7.276e-01);
    end;
    if W=83 then
    begin
        R := -Double(7.633e-01);
    end;
    if W=82 then
    begin
        R := -Double(8.001e-01);
    end;
    if W=81 then
    begin
        R := -Double(8.381e-01);
    end;
    if W=80 then
    begin
        R := -Double(8.774e-01);
    end;
    if W=79 then
    begin
        R := -Double(9.179e-01);
    end;
    if W=78 then
    begin
        R := -Double(9.597e-01);
    end;
    if W=77 then
    begin
        R := -Double(1.003e+00);
    end;
    if W=76 then
    begin
        R := -Double(1.047e+00);
    end;
    if W=75 then
    begin
        R := -Double(1.093e+00);
    end;
    if W=74 then
    begin
        R := -Double(1.140e+00);
    end;
    if W=73 then
    begin
        R := -Double(1.188e+00);
    end;
    if W=72 then
    begin
        R := -Double(1.238e+00);
    end;
    if W=71 then
    begin
        R := -Double(1.289e+00);
    end;
    if W=70 then
    begin
        R := -Double(1.342e+00);
    end;
    if W=69 then
    begin
        R := -Double(1.396e+00);
    end;
    if W=68 then
    begin
        R := -Double(1.452e+00);
    end;
    if W=67 then
    begin
        R := -Double(1.509e+00);
    end;
    if W=66 then
    begin
        R := -Double(1.568e+00);
    end;
    if W=65 then
    begin
        R := -Double(1.628e+00);
    end;
    if W=64 then
    begin
        R := -Double(1.690e+00);
    end;
    if W=63 then
    begin
        R := -Double(1.753e+00);
    end;
    if W=62 then
    begin
        R := -Double(1.818e+00);
    end;
    if W=61 then
    begin
        R := -Double(1.885e+00);
    end;
    if W=60 then
    begin
        R := -Double(1.953e+00);
    end;
    if W=59 then
    begin
        R := -Double(2.023e+00);
    end;
    if W=58 then
    begin
        R := -Double(2.095e+00);
    end;
    if W=57 then
    begin
        R := -Double(2.168e+00);
    end;
    if W=56 then
    begin
        R := -Double(2.244e+00);
    end;
    if W=55 then
    begin
        R := -Double(2.321e+00);
    end;
    if W=54 then
    begin
        R := -Double(2.400e+00);
    end;
    if W=53 then
    begin
        R := -Double(2.481e+00);
    end;
    if W=52 then
    begin
        R := -Double(2.564e+00);
    end;
    if W=51 then
    begin
        R := -Double(2.648e+00);
    end;
    if W=50 then
    begin
        R := -Double(2.735e+00);
    end;
    if W=49 then
    begin
        R := -Double(2.824e+00);
    end;
    if W=48 then
    begin
        R := -Double(2.915e+00);
    end;
    if W=47 then
    begin
        R := -Double(3.008e+00);
    end;
    if W=46 then
    begin
        R := -Double(3.104e+00);
    end;
    if W=45 then
    begin
        R := -Double(3.201e+00);
    end;
    if W=44 then
    begin
        R := -Double(3.301e+00);
    end;
    if W=43 then
    begin
        R := -Double(3.403e+00);
    end;
    if W=42 then
    begin
        R := -Double(3.508e+00);
    end;
    if W=41 then
    begin
        R := -Double(3.615e+00);
    end;
    if W=40 then
    begin
        R := -Double(3.724e+00);
    end;
    if W=39 then
    begin
        R := -Double(3.836e+00);
    end;
    if W=38 then
    begin
        R := -Double(3.950e+00);
    end;
    if W=37 then
    begin
        R := -Double(4.068e+00);
    end;
    if W=36 then
    begin
        R := -Double(4.188e+00);
    end;
    if W=35 then
    begin
        R := -Double(4.311e+00);
    end;
    if W=34 then
    begin
        R := -Double(4.437e+00);
    end;
    if W=33 then
    begin
        R := -Double(4.565e+00);
    end;
    if W=32 then
    begin
        R := -Double(4.698e+00);
    end;
    if W=31 then
    begin
        R := -Double(4.833e+00);
    end;
    if W=30 then
    begin
        R := -Double(4.971e+00);
    end;
    if W=29 then
    begin
        R := -Double(5.113e+00);
    end;
    if W=28 then
    begin
        R := -Double(5.258e+00);
    end;
    if W=27 then
    begin
        R := -Double(5.408e+00);
    end;
    if W=26 then
    begin
        R := -Double(5.561e+00);
    end;
    if W=25 then
    begin
        R := -Double(5.717e+00);
    end;
    if W=24 then
    begin
        R := -Double(5.878e+00);
    end;
    if W=23 then
    begin
        R := -Double(6.044e+00);
    end;
    if W=22 then
    begin
        R := -Double(6.213e+00);
    end;
    if W=21 then
    begin
        R := -Double(6.388e+00);
    end;
    if W=20 then
    begin
        R := -Double(6.569e+00);
    end;
    if W=19 then
    begin
        R := -Double(6.753e+00);
    end;
    if W=18 then
    begin
        R := -Double(6.943e+00);
    end;
    if W=17 then
    begin
        R := -Double(7.144e+00);
    end;
    if W=16 then
    begin
        R := -Double(7.347e+00);
    end;
    if W=15 then
    begin
        R := -Double(7.557e+00);
    end;
    if W=14 then
    begin
        R := -Double(7.776e+00);
    end;
    if W=13 then
    begin
        R := -Double(7.999e+00);
    end;
    if W=12 then
    begin
        R := -Double(8.228e+00);
    end;
    if W=11 then
    begin
        R := -Double(8.469e+00);
    end;
    if W=10 then
    begin
        R := -Double(8.715e+00);
    end;
    if W=9 then
    begin
        R := -Double(8.980e+00);
    end;
    if W=8 then
    begin
        R := -Double(9.258e+00);
    end;
    if W=7 then
    begin
        R := -Double(9.532e+00);
    end;
    if W=6 then
    begin
        R := -Double(9.838e+00);
    end;
    if W=5 then
    begin
        R := -Double(1.017e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.053e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.087e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.138e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.178e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.248e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 19)
*************************************************************************)
function W19(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(2.484955e+01)*S+Double(9.500000e+01));
    if W>=95 then
    begin
        R := -Double(6.776e-01);
    end;
    if W=94 then
    begin
        R := -Double(7.089e-01);
    end;
    if W=93 then
    begin
        R := -Double(7.413e-01);
    end;
    if W=92 then
    begin
        R := -Double(7.747e-01);
    end;
    if W=91 then
    begin
        R := -Double(8.090e-01);
    end;
    if W=90 then
    begin
        R := -Double(8.445e-01);
    end;
    if W=89 then
    begin
        R := -Double(8.809e-01);
    end;
    if W=88 then
    begin
        R := -Double(9.185e-01);
    end;
    if W=87 then
    begin
        R := -Double(9.571e-01);
    end;
    if W=86 then
    begin
        R := -Double(9.968e-01);
    end;
    if W=85 then
    begin
        R := -Double(1.038e+00);
    end;
    if W=84 then
    begin
        R := -Double(1.080e+00);
    end;
    if W=83 then
    begin
        R := -Double(1.123e+00);
    end;
    if W=82 then
    begin
        R := -Double(1.167e+00);
    end;
    if W=81 then
    begin
        R := -Double(1.213e+00);
    end;
    if W=80 then
    begin
        R := -Double(1.259e+00);
    end;
    if W=79 then
    begin
        R := -Double(1.307e+00);
    end;
    if W=78 then
    begin
        R := -Double(1.356e+00);
    end;
    if W=77 then
    begin
        R := -Double(1.407e+00);
    end;
    if W=76 then
    begin
        R := -Double(1.458e+00);
    end;
    if W=75 then
    begin
        R := -Double(1.511e+00);
    end;
    if W=74 then
    begin
        R := -Double(1.565e+00);
    end;
    if W=73 then
    begin
        R := -Double(1.621e+00);
    end;
    if W=72 then
    begin
        R := -Double(1.678e+00);
    end;
    if W=71 then
    begin
        R := -Double(1.736e+00);
    end;
    if W=70 then
    begin
        R := -Double(1.796e+00);
    end;
    if W=69 then
    begin
        R := -Double(1.857e+00);
    end;
    if W=68 then
    begin
        R := -Double(1.919e+00);
    end;
    if W=67 then
    begin
        R := -Double(1.983e+00);
    end;
    if W=66 then
    begin
        R := -Double(2.048e+00);
    end;
    if W=65 then
    begin
        R := -Double(2.115e+00);
    end;
    if W=64 then
    begin
        R := -Double(2.183e+00);
    end;
    if W=63 then
    begin
        R := -Double(2.253e+00);
    end;
    if W=62 then
    begin
        R := -Double(2.325e+00);
    end;
    if W=61 then
    begin
        R := -Double(2.398e+00);
    end;
    if W=60 then
    begin
        R := -Double(2.472e+00);
    end;
    if W=59 then
    begin
        R := -Double(2.548e+00);
    end;
    if W=58 then
    begin
        R := -Double(2.626e+00);
    end;
    if W=57 then
    begin
        R := -Double(2.706e+00);
    end;
    if W=56 then
    begin
        R := -Double(2.787e+00);
    end;
    if W=55 then
    begin
        R := -Double(2.870e+00);
    end;
    if W=54 then
    begin
        R := -Double(2.955e+00);
    end;
    if W=53 then
    begin
        R := -Double(3.042e+00);
    end;
    if W=52 then
    begin
        R := -Double(3.130e+00);
    end;
    if W=51 then
    begin
        R := -Double(3.220e+00);
    end;
    if W=50 then
    begin
        R := -Double(3.313e+00);
    end;
    if W=49 then
    begin
        R := -Double(3.407e+00);
    end;
    if W=48 then
    begin
        R := -Double(3.503e+00);
    end;
    if W=47 then
    begin
        R := -Double(3.601e+00);
    end;
    if W=46 then
    begin
        R := -Double(3.702e+00);
    end;
    if W=45 then
    begin
        R := -Double(3.804e+00);
    end;
    if W=44 then
    begin
        R := -Double(3.909e+00);
    end;
    if W=43 then
    begin
        R := -Double(4.015e+00);
    end;
    if W=42 then
    begin
        R := -Double(4.125e+00);
    end;
    if W=41 then
    begin
        R := -Double(4.236e+00);
    end;
    if W=40 then
    begin
        R := -Double(4.350e+00);
    end;
    if W=39 then
    begin
        R := -Double(4.466e+00);
    end;
    if W=38 then
    begin
        R := -Double(4.585e+00);
    end;
    if W=37 then
    begin
        R := -Double(4.706e+00);
    end;
    if W=36 then
    begin
        R := -Double(4.830e+00);
    end;
    if W=35 then
    begin
        R := -Double(4.957e+00);
    end;
    if W=34 then
    begin
        R := -Double(5.086e+00);
    end;
    if W=33 then
    begin
        R := -Double(5.219e+00);
    end;
    if W=32 then
    begin
        R := -Double(5.355e+00);
    end;
    if W=31 then
    begin
        R := -Double(5.493e+00);
    end;
    if W=30 then
    begin
        R := -Double(5.634e+00);
    end;
    if W=29 then
    begin
        R := -Double(5.780e+00);
    end;
    if W=28 then
    begin
        R := -Double(5.928e+00);
    end;
    if W=27 then
    begin
        R := -Double(6.080e+00);
    end;
    if W=26 then
    begin
        R := -Double(6.235e+00);
    end;
    if W=25 then
    begin
        R := -Double(6.394e+00);
    end;
    if W=24 then
    begin
        R := -Double(6.558e+00);
    end;
    if W=23 then
    begin
        R := -Double(6.726e+00);
    end;
    if W=22 then
    begin
        R := -Double(6.897e+00);
    end;
    if W=21 then
    begin
        R := -Double(7.074e+00);
    end;
    if W=20 then
    begin
        R := -Double(7.256e+00);
    end;
    if W=19 then
    begin
        R := -Double(7.443e+00);
    end;
    if W=18 then
    begin
        R := -Double(7.636e+00);
    end;
    if W=17 then
    begin
        R := -Double(7.837e+00);
    end;
    if W=16 then
    begin
        R := -Double(8.040e+00);
    end;
    if W=15 then
    begin
        R := -Double(8.250e+00);
    end;
    if W=14 then
    begin
        R := -Double(8.469e+00);
    end;
    if W=13 then
    begin
        R := -Double(8.692e+00);
    end;
    if W=12 then
    begin
        R := -Double(8.921e+00);
    end;
    if W=11 then
    begin
        R := -Double(9.162e+00);
    end;
    if W=10 then
    begin
        R := -Double(9.409e+00);
    end;
    if W=9 then
    begin
        R := -Double(9.673e+00);
    end;
    if W=8 then
    begin
        R := -Double(9.951e+00);
    end;
    if W=7 then
    begin
        R := -Double(1.023e+01);
    end;
    if W=6 then
    begin
        R := -Double(1.053e+01);
    end;
    if W=5 then
    begin
        R := -Double(1.087e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.122e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.156e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.207e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.248e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.317e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 20)
*************************************************************************)
function W20(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(2.678619e+01)*S+Double(1.050000e+02));
    if W>=105 then
    begin
        R := -Double(6.787e-01);
    end;
    if W=104 then
    begin
        R := -Double(7.078e-01);
    end;
    if W=103 then
    begin
        R := -Double(7.378e-01);
    end;
    if W=102 then
    begin
        R := -Double(7.686e-01);
    end;
    if W=101 then
    begin
        R := -Double(8.004e-01);
    end;
    if W=100 then
    begin
        R := -Double(8.330e-01);
    end;
    if W=99 then
    begin
        R := -Double(8.665e-01);
    end;
    if W=98 then
    begin
        R := -Double(9.010e-01);
    end;
    if W=97 then
    begin
        R := -Double(9.363e-01);
    end;
    if W=96 then
    begin
        R := -Double(9.726e-01);
    end;
    if W=95 then
    begin
        R := -Double(1.010e+00);
    end;
    if W=94 then
    begin
        R := -Double(1.048e+00);
    end;
    if W=93 then
    begin
        R := -Double(1.087e+00);
    end;
    if W=92 then
    begin
        R := -Double(1.128e+00);
    end;
    if W=91 then
    begin
        R := -Double(1.169e+00);
    end;
    if W=90 then
    begin
        R := -Double(1.211e+00);
    end;
    if W=89 then
    begin
        R := -Double(1.254e+00);
    end;
    if W=88 then
    begin
        R := -Double(1.299e+00);
    end;
    if W=87 then
    begin
        R := -Double(1.344e+00);
    end;
    if W=86 then
    begin
        R := -Double(1.390e+00);
    end;
    if W=85 then
    begin
        R := -Double(1.438e+00);
    end;
    if W=84 then
    begin
        R := -Double(1.486e+00);
    end;
    if W=83 then
    begin
        R := -Double(1.536e+00);
    end;
    if W=82 then
    begin
        R := -Double(1.587e+00);
    end;
    if W=81 then
    begin
        R := -Double(1.639e+00);
    end;
    if W=80 then
    begin
        R := -Double(1.692e+00);
    end;
    if W=79 then
    begin
        R := -Double(1.746e+00);
    end;
    if W=78 then
    begin
        R := -Double(1.802e+00);
    end;
    if W=77 then
    begin
        R := -Double(1.859e+00);
    end;
    if W=76 then
    begin
        R := -Double(1.916e+00);
    end;
    if W=75 then
    begin
        R := -Double(1.976e+00);
    end;
    if W=74 then
    begin
        R := -Double(2.036e+00);
    end;
    if W=73 then
    begin
        R := -Double(2.098e+00);
    end;
    if W=72 then
    begin
        R := -Double(2.161e+00);
    end;
    if W=71 then
    begin
        R := -Double(2.225e+00);
    end;
    if W=70 then
    begin
        R := -Double(2.290e+00);
    end;
    if W=69 then
    begin
        R := -Double(2.357e+00);
    end;
    if W=68 then
    begin
        R := -Double(2.426e+00);
    end;
    if W=67 then
    begin
        R := -Double(2.495e+00);
    end;
    if W=66 then
    begin
        R := -Double(2.566e+00);
    end;
    if W=65 then
    begin
        R := -Double(2.639e+00);
    end;
    if W=64 then
    begin
        R := -Double(2.713e+00);
    end;
    if W=63 then
    begin
        R := -Double(2.788e+00);
    end;
    if W=62 then
    begin
        R := -Double(2.865e+00);
    end;
    if W=61 then
    begin
        R := -Double(2.943e+00);
    end;
    if W=60 then
    begin
        R := -Double(3.023e+00);
    end;
    if W=59 then
    begin
        R := -Double(3.104e+00);
    end;
    if W=58 then
    begin
        R := -Double(3.187e+00);
    end;
    if W=57 then
    begin
        R := -Double(3.272e+00);
    end;
    if W=56 then
    begin
        R := -Double(3.358e+00);
    end;
    if W=55 then
    begin
        R := -Double(3.446e+00);
    end;
    if W=54 then
    begin
        R := -Double(3.536e+00);
    end;
    if W=53 then
    begin
        R := -Double(3.627e+00);
    end;
    if W=52 then
    begin
        R := -Double(3.721e+00);
    end;
    if W=51 then
    begin
        R := -Double(3.815e+00);
    end;
    if W=50 then
    begin
        R := -Double(3.912e+00);
    end;
    if W=49 then
    begin
        R := -Double(4.011e+00);
    end;
    if W=48 then
    begin
        R := -Double(4.111e+00);
    end;
    if W=47 then
    begin
        R := -Double(4.214e+00);
    end;
    if W=46 then
    begin
        R := -Double(4.318e+00);
    end;
    if W=45 then
    begin
        R := -Double(4.425e+00);
    end;
    if W=44 then
    begin
        R := -Double(4.534e+00);
    end;
    if W=43 then
    begin
        R := -Double(4.644e+00);
    end;
    if W=42 then
    begin
        R := -Double(4.757e+00);
    end;
    if W=41 then
    begin
        R := -Double(4.872e+00);
    end;
    if W=40 then
    begin
        R := -Double(4.990e+00);
    end;
    if W=39 then
    begin
        R := -Double(5.109e+00);
    end;
    if W=38 then
    begin
        R := -Double(5.232e+00);
    end;
    if W=37 then
    begin
        R := -Double(5.356e+00);
    end;
    if W=36 then
    begin
        R := -Double(5.484e+00);
    end;
    if W=35 then
    begin
        R := -Double(5.614e+00);
    end;
    if W=34 then
    begin
        R := -Double(5.746e+00);
    end;
    if W=33 then
    begin
        R := -Double(5.882e+00);
    end;
    if W=32 then
    begin
        R := -Double(6.020e+00);
    end;
    if W=31 then
    begin
        R := -Double(6.161e+00);
    end;
    if W=30 then
    begin
        R := -Double(6.305e+00);
    end;
    if W=29 then
    begin
        R := -Double(6.453e+00);
    end;
    if W=28 then
    begin
        R := -Double(6.603e+00);
    end;
    if W=27 then
    begin
        R := -Double(6.757e+00);
    end;
    if W=26 then
    begin
        R := -Double(6.915e+00);
    end;
    if W=25 then
    begin
        R := -Double(7.076e+00);
    end;
    if W=24 then
    begin
        R := -Double(7.242e+00);
    end;
    if W=23 then
    begin
        R := -Double(7.411e+00);
    end;
    if W=22 then
    begin
        R := -Double(7.584e+00);
    end;
    if W=21 then
    begin
        R := -Double(7.763e+00);
    end;
    if W=20 then
    begin
        R := -Double(7.947e+00);
    end;
    if W=19 then
    begin
        R := -Double(8.136e+00);
    end;
    if W=18 then
    begin
        R := -Double(8.330e+00);
    end;
    if W=17 then
    begin
        R := -Double(8.530e+00);
    end;
    if W=16 then
    begin
        R := -Double(8.733e+00);
    end;
    if W=15 then
    begin
        R := -Double(8.943e+00);
    end;
    if W=14 then
    begin
        R := -Double(9.162e+00);
    end;
    if W=13 then
    begin
        R := -Double(9.386e+00);
    end;
    if W=12 then
    begin
        R := -Double(9.614e+00);
    end;
    if W=11 then
    begin
        R := -Double(9.856e+00);
    end;
    if W=10 then
    begin
        R := -Double(1.010e+01);
    end;
    if W=9 then
    begin
        R := -Double(1.037e+01);
    end;
    if W=8 then
    begin
        R := -Double(1.064e+01);
    end;
    if W=7 then
    begin
        R := -Double(1.092e+01);
    end;
    if W=6 then
    begin
        R := -Double(1.122e+01);
    end;
    if W=5 then
    begin
        R := -Double(1.156e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.192e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.225e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.276e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.317e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.386e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 21)
*************************************************************************)
function W21(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(2.877064e+01)*S+Double(1.155000e+02));
    if W>=115 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=114 then
    begin
        R := -Double(7.207e-01);
    end;
    if W=113 then
    begin
        R := -Double(7.489e-01);
    end;
    if W=112 then
    begin
        R := -Double(7.779e-01);
    end;
    if W=111 then
    begin
        R := -Double(8.077e-01);
    end;
    if W=110 then
    begin
        R := -Double(8.383e-01);
    end;
    if W=109 then
    begin
        R := -Double(8.697e-01);
    end;
    if W=108 then
    begin
        R := -Double(9.018e-01);
    end;
    if W=107 then
    begin
        R := -Double(9.348e-01);
    end;
    if W=106 then
    begin
        R := -Double(9.685e-01);
    end;
    if W=105 then
    begin
        R := -Double(1.003e+00);
    end;
    if W=104 then
    begin
        R := -Double(1.039e+00);
    end;
    if W=103 then
    begin
        R := -Double(1.075e+00);
    end;
    if W=102 then
    begin
        R := -Double(1.112e+00);
    end;
    if W=101 then
    begin
        R := -Double(1.150e+00);
    end;
    if W=100 then
    begin
        R := -Double(1.189e+00);
    end;
    if W=99 then
    begin
        R := -Double(1.229e+00);
    end;
    if W=98 then
    begin
        R := -Double(1.269e+00);
    end;
    if W=97 then
    begin
        R := -Double(1.311e+00);
    end;
    if W=96 then
    begin
        R := -Double(1.353e+00);
    end;
    if W=95 then
    begin
        R := -Double(1.397e+00);
    end;
    if W=94 then
    begin
        R := -Double(1.441e+00);
    end;
    if W=93 then
    begin
        R := -Double(1.486e+00);
    end;
    if W=92 then
    begin
        R := -Double(1.533e+00);
    end;
    if W=91 then
    begin
        R := -Double(1.580e+00);
    end;
    if W=90 then
    begin
        R := -Double(1.628e+00);
    end;
    if W=89 then
    begin
        R := -Double(1.677e+00);
    end;
    if W=88 then
    begin
        R := -Double(1.728e+00);
    end;
    if W=87 then
    begin
        R := -Double(1.779e+00);
    end;
    if W=86 then
    begin
        R := -Double(1.831e+00);
    end;
    if W=85 then
    begin
        R := -Double(1.884e+00);
    end;
    if W=84 then
    begin
        R := -Double(1.939e+00);
    end;
    if W=83 then
    begin
        R := -Double(1.994e+00);
    end;
    if W=82 then
    begin
        R := -Double(2.051e+00);
    end;
    if W=81 then
    begin
        R := -Double(2.108e+00);
    end;
    if W=80 then
    begin
        R := -Double(2.167e+00);
    end;
    if W=79 then
    begin
        R := -Double(2.227e+00);
    end;
    if W=78 then
    begin
        R := -Double(2.288e+00);
    end;
    if W=77 then
    begin
        R := -Double(2.350e+00);
    end;
    if W=76 then
    begin
        R := -Double(2.414e+00);
    end;
    if W=75 then
    begin
        R := -Double(2.478e+00);
    end;
    if W=74 then
    begin
        R := -Double(2.544e+00);
    end;
    if W=73 then
    begin
        R := -Double(2.611e+00);
    end;
    if W=72 then
    begin
        R := -Double(2.679e+00);
    end;
    if W=71 then
    begin
        R := -Double(2.748e+00);
    end;
    if W=70 then
    begin
        R := -Double(2.819e+00);
    end;
    if W=69 then
    begin
        R := -Double(2.891e+00);
    end;
    if W=68 then
    begin
        R := -Double(2.964e+00);
    end;
    if W=67 then
    begin
        R := -Double(3.039e+00);
    end;
    if W=66 then
    begin
        R := -Double(3.115e+00);
    end;
    if W=65 then
    begin
        R := -Double(3.192e+00);
    end;
    if W=64 then
    begin
        R := -Double(3.270e+00);
    end;
    if W=63 then
    begin
        R := -Double(3.350e+00);
    end;
    if W=62 then
    begin
        R := -Double(3.432e+00);
    end;
    if W=61 then
    begin
        R := -Double(3.515e+00);
    end;
    if W=60 then
    begin
        R := -Double(3.599e+00);
    end;
    if W=59 then
    begin
        R := -Double(3.685e+00);
    end;
    if W=58 then
    begin
        R := -Double(3.772e+00);
    end;
    if W=57 then
    begin
        R := -Double(3.861e+00);
    end;
    if W=56 then
    begin
        R := -Double(3.952e+00);
    end;
    if W=55 then
    begin
        R := -Double(4.044e+00);
    end;
    if W=54 then
    begin
        R := -Double(4.138e+00);
    end;
    if W=53 then
    begin
        R := -Double(4.233e+00);
    end;
    if W=52 then
    begin
        R := -Double(4.330e+00);
    end;
    if W=51 then
    begin
        R := -Double(4.429e+00);
    end;
    if W=50 then
    begin
        R := -Double(4.530e+00);
    end;
    if W=49 then
    begin
        R := -Double(4.632e+00);
    end;
    if W=48 then
    begin
        R := -Double(4.736e+00);
    end;
    if W=47 then
    begin
        R := -Double(4.842e+00);
    end;
    if W=46 then
    begin
        R := -Double(4.950e+00);
    end;
    if W=45 then
    begin
        R := -Double(5.060e+00);
    end;
    if W=44 then
    begin
        R := -Double(5.172e+00);
    end;
    if W=43 then
    begin
        R := -Double(5.286e+00);
    end;
    if W=42 then
    begin
        R := -Double(5.402e+00);
    end;
    if W=41 then
    begin
        R := -Double(5.520e+00);
    end;
    if W=40 then
    begin
        R := -Double(5.641e+00);
    end;
    if W=39 then
    begin
        R := -Double(5.763e+00);
    end;
    if W=38 then
    begin
        R := -Double(5.889e+00);
    end;
    if W=37 then
    begin
        R := -Double(6.016e+00);
    end;
    if W=36 then
    begin
        R := -Double(6.146e+00);
    end;
    if W=35 then
    begin
        R := -Double(6.278e+00);
    end;
    if W=34 then
    begin
        R := -Double(6.413e+00);
    end;
    if W=33 then
    begin
        R := -Double(6.551e+00);
    end;
    if W=32 then
    begin
        R := -Double(6.692e+00);
    end;
    if W=31 then
    begin
        R := -Double(6.835e+00);
    end;
    if W=30 then
    begin
        R := -Double(6.981e+00);
    end;
    if W=29 then
    begin
        R := -Double(7.131e+00);
    end;
    if W=28 then
    begin
        R := -Double(7.283e+00);
    end;
    if W=27 then
    begin
        R := -Double(7.439e+00);
    end;
    if W=26 then
    begin
        R := -Double(7.599e+00);
    end;
    if W=25 then
    begin
        R := -Double(7.762e+00);
    end;
    if W=24 then
    begin
        R := -Double(7.928e+00);
    end;
    if W=23 then
    begin
        R := -Double(8.099e+00);
    end;
    if W=22 then
    begin
        R := -Double(8.274e+00);
    end;
    if W=21 then
    begin
        R := -Double(8.454e+00);
    end;
    if W=20 then
    begin
        R := -Double(8.640e+00);
    end;
    if W=19 then
    begin
        R := -Double(8.829e+00);
    end;
    if W=18 then
    begin
        R := -Double(9.023e+00);
    end;
    if W=17 then
    begin
        R := -Double(9.223e+00);
    end;
    if W=16 then
    begin
        R := -Double(9.426e+00);
    end;
    if W=15 then
    begin
        R := -Double(9.636e+00);
    end;
    if W=14 then
    begin
        R := -Double(9.856e+00);
    end;
    if W=13 then
    begin
        R := -Double(1.008e+01);
    end;
    if W=12 then
    begin
        R := -Double(1.031e+01);
    end;
    if W=11 then
    begin
        R := -Double(1.055e+01);
    end;
    if W=10 then
    begin
        R := -Double(1.079e+01);
    end;
    if W=9 then
    begin
        R := -Double(1.106e+01);
    end;
    if W=8 then
    begin
        R := -Double(1.134e+01);
    end;
    if W=7 then
    begin
        R := -Double(1.161e+01);
    end;
    if W=6 then
    begin
        R := -Double(1.192e+01);
    end;
    if W=5 then
    begin
        R := -Double(1.225e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.261e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.295e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.346e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.386e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.456e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 22)
*************************************************************************)
function W22(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(3.080179e+01)*S+Double(1.265000e+02));
    if W>=126 then
    begin
        R := -Double(6.931e-01);
    end;
    if W=125 then
    begin
        R := -Double(7.189e-01);
    end;
    if W=124 then
    begin
        R := -Double(7.452e-01);
    end;
    if W=123 then
    begin
        R := -Double(7.722e-01);
    end;
    if W=122 then
    begin
        R := -Double(7.999e-01);
    end;
    if W=121 then
    begin
        R := -Double(8.283e-01);
    end;
    if W=120 then
    begin
        R := -Double(8.573e-01);
    end;
    if W=119 then
    begin
        R := -Double(8.871e-01);
    end;
    if W=118 then
    begin
        R := -Double(9.175e-01);
    end;
    if W=117 then
    begin
        R := -Double(9.486e-01);
    end;
    if W=116 then
    begin
        R := -Double(9.805e-01);
    end;
    if W=115 then
    begin
        R := -Double(1.013e+00);
    end;
    if W=114 then
    begin
        R := -Double(1.046e+00);
    end;
    if W=113 then
    begin
        R := -Double(1.080e+00);
    end;
    if W=112 then
    begin
        R := -Double(1.115e+00);
    end;
    if W=111 then
    begin
        R := -Double(1.151e+00);
    end;
    if W=110 then
    begin
        R := -Double(1.187e+00);
    end;
    if W=109 then
    begin
        R := -Double(1.224e+00);
    end;
    if W=108 then
    begin
        R := -Double(1.262e+00);
    end;
    if W=107 then
    begin
        R := -Double(1.301e+00);
    end;
    if W=106 then
    begin
        R := -Double(1.340e+00);
    end;
    if W=105 then
    begin
        R := -Double(1.381e+00);
    end;
    if W=104 then
    begin
        R := -Double(1.422e+00);
    end;
    if W=103 then
    begin
        R := -Double(1.464e+00);
    end;
    if W=102 then
    begin
        R := -Double(1.506e+00);
    end;
    if W=101 then
    begin
        R := -Double(1.550e+00);
    end;
    if W=100 then
    begin
        R := -Double(1.594e+00);
    end;
    if W=99 then
    begin
        R := -Double(1.640e+00);
    end;
    if W=98 then
    begin
        R := -Double(1.686e+00);
    end;
    if W=97 then
    begin
        R := -Double(1.733e+00);
    end;
    if W=96 then
    begin
        R := -Double(1.781e+00);
    end;
    if W=95 then
    begin
        R := -Double(1.830e+00);
    end;
    if W=94 then
    begin
        R := -Double(1.880e+00);
    end;
    if W=93 then
    begin
        R := -Double(1.930e+00);
    end;
    if W=92 then
    begin
        R := -Double(1.982e+00);
    end;
    if W=91 then
    begin
        R := -Double(2.034e+00);
    end;
    if W=90 then
    begin
        R := -Double(2.088e+00);
    end;
    if W=89 then
    begin
        R := -Double(2.142e+00);
    end;
    if W=88 then
    begin
        R := -Double(2.198e+00);
    end;
    if W=87 then
    begin
        R := -Double(2.254e+00);
    end;
    if W=86 then
    begin
        R := -Double(2.312e+00);
    end;
    if W=85 then
    begin
        R := -Double(2.370e+00);
    end;
    if W=84 then
    begin
        R := -Double(2.429e+00);
    end;
    if W=83 then
    begin
        R := -Double(2.490e+00);
    end;
    if W=82 then
    begin
        R := -Double(2.551e+00);
    end;
    if W=81 then
    begin
        R := -Double(2.614e+00);
    end;
    if W=80 then
    begin
        R := -Double(2.677e+00);
    end;
    if W=79 then
    begin
        R := -Double(2.742e+00);
    end;
    if W=78 then
    begin
        R := -Double(2.808e+00);
    end;
    if W=77 then
    begin
        R := -Double(2.875e+00);
    end;
    if W=76 then
    begin
        R := -Double(2.943e+00);
    end;
    if W=75 then
    begin
        R := -Double(3.012e+00);
    end;
    if W=74 then
    begin
        R := -Double(3.082e+00);
    end;
    if W=73 then
    begin
        R := -Double(3.153e+00);
    end;
    if W=72 then
    begin
        R := -Double(3.226e+00);
    end;
    if W=71 then
    begin
        R := -Double(3.300e+00);
    end;
    if W=70 then
    begin
        R := -Double(3.375e+00);
    end;
    if W=69 then
    begin
        R := -Double(3.451e+00);
    end;
    if W=68 then
    begin
        R := -Double(3.529e+00);
    end;
    if W=67 then
    begin
        R := -Double(3.607e+00);
    end;
    if W=66 then
    begin
        R := -Double(3.687e+00);
    end;
    if W=65 then
    begin
        R := -Double(3.769e+00);
    end;
    if W=64 then
    begin
        R := -Double(3.851e+00);
    end;
    if W=63 then
    begin
        R := -Double(3.935e+00);
    end;
    if W=62 then
    begin
        R := -Double(4.021e+00);
    end;
    if W=61 then
    begin
        R := -Double(4.108e+00);
    end;
    if W=60 then
    begin
        R := -Double(4.196e+00);
    end;
    if W=59 then
    begin
        R := -Double(4.285e+00);
    end;
    if W=58 then
    begin
        R := -Double(4.376e+00);
    end;
    if W=57 then
    begin
        R := -Double(4.469e+00);
    end;
    if W=56 then
    begin
        R := -Double(4.563e+00);
    end;
    if W=55 then
    begin
        R := -Double(4.659e+00);
    end;
    if W=54 then
    begin
        R := -Double(4.756e+00);
    end;
    if W=53 then
    begin
        R := -Double(4.855e+00);
    end;
    if W=52 then
    begin
        R := -Double(4.955e+00);
    end;
    if W=51 then
    begin
        R := -Double(5.057e+00);
    end;
    if W=50 then
    begin
        R := -Double(5.161e+00);
    end;
    if W=49 then
    begin
        R := -Double(5.266e+00);
    end;
    if W=48 then
    begin
        R := -Double(5.374e+00);
    end;
    if W=47 then
    begin
        R := -Double(5.483e+00);
    end;
    if W=46 then
    begin
        R := -Double(5.594e+00);
    end;
    if W=45 then
    begin
        R := -Double(5.706e+00);
    end;
    if W=44 then
    begin
        R := -Double(5.821e+00);
    end;
    if W=43 then
    begin
        R := -Double(5.938e+00);
    end;
    if W=42 then
    begin
        R := -Double(6.057e+00);
    end;
    if W=41 then
    begin
        R := -Double(6.177e+00);
    end;
    if W=40 then
    begin
        R := -Double(6.300e+00);
    end;
    if W=39 then
    begin
        R := -Double(6.426e+00);
    end;
    if W=38 then
    begin
        R := -Double(6.553e+00);
    end;
    if W=37 then
    begin
        R := -Double(6.683e+00);
    end;
    if W=36 then
    begin
        R := -Double(6.815e+00);
    end;
    if W=35 then
    begin
        R := -Double(6.949e+00);
    end;
    if W=34 then
    begin
        R := -Double(7.086e+00);
    end;
    if W=33 then
    begin
        R := -Double(7.226e+00);
    end;
    if W=32 then
    begin
        R := -Double(7.368e+00);
    end;
    if W=31 then
    begin
        R := -Double(7.513e+00);
    end;
    if W=30 then
    begin
        R := -Double(7.661e+00);
    end;
    if W=29 then
    begin
        R := -Double(7.813e+00);
    end;
    if W=28 then
    begin
        R := -Double(7.966e+00);
    end;
    if W=27 then
    begin
        R := -Double(8.124e+00);
    end;
    if W=26 then
    begin
        R := -Double(8.285e+00);
    end;
    if W=25 then
    begin
        R := -Double(8.449e+00);
    end;
    if W=24 then
    begin
        R := -Double(8.617e+00);
    end;
    if W=23 then
    begin
        R := -Double(8.789e+00);
    end;
    if W=22 then
    begin
        R := -Double(8.965e+00);
    end;
    if W=21 then
    begin
        R := -Double(9.147e+00);
    end;
    if W=20 then
    begin
        R := -Double(9.333e+00);
    end;
    if W=19 then
    begin
        R := -Double(9.522e+00);
    end;
    if W=18 then
    begin
        R := -Double(9.716e+00);
    end;
    if W=17 then
    begin
        R := -Double(9.917e+00);
    end;
    if W=16 then
    begin
        R := -Double(1.012e+01);
    end;
    if W=15 then
    begin
        R := -Double(1.033e+01);
    end;
    if W=14 then
    begin
        R := -Double(1.055e+01);
    end;
    if W=13 then
    begin
        R := -Double(1.077e+01);
    end;
    if W=12 then
    begin
        R := -Double(1.100e+01);
    end;
    if W=11 then
    begin
        R := -Double(1.124e+01);
    end;
    if W=10 then
    begin
        R := -Double(1.149e+01);
    end;
    if W=9 then
    begin
        R := -Double(1.175e+01);
    end;
    if W=8 then
    begin
        R := -Double(1.203e+01);
    end;
    if W=7 then
    begin
        R := -Double(1.230e+01);
    end;
    if W=6 then
    begin
        R := -Double(1.261e+01);
    end;
    if W=5 then
    begin
        R := -Double(1.295e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.330e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.364e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.415e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.456e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.525e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 23)
*************************************************************************)
function W23(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(3.287856e+01)*S+Double(1.380000e+02));
    if W>=138 then
    begin
        R := -Double(6.813e-01);
    end;
    if W=137 then
    begin
        R := -Double(7.051e-01);
    end;
    if W=136 then
    begin
        R := -Double(7.295e-01);
    end;
    if W=135 then
    begin
        R := -Double(7.544e-01);
    end;
    if W=134 then
    begin
        R := -Double(7.800e-01);
    end;
    if W=133 then
    begin
        R := -Double(8.061e-01);
    end;
    if W=132 then
    begin
        R := -Double(8.328e-01);
    end;
    if W=131 then
    begin
        R := -Double(8.601e-01);
    end;
    if W=130 then
    begin
        R := -Double(8.880e-01);
    end;
    if W=129 then
    begin
        R := -Double(9.166e-01);
    end;
    if W=128 then
    begin
        R := -Double(9.457e-01);
    end;
    if W=127 then
    begin
        R := -Double(9.755e-01);
    end;
    if W=126 then
    begin
        R := -Double(1.006e+00);
    end;
    if W=125 then
    begin
        R := -Double(1.037e+00);
    end;
    if W=124 then
    begin
        R := -Double(1.069e+00);
    end;
    if W=123 then
    begin
        R := -Double(1.101e+00);
    end;
    if W=122 then
    begin
        R := -Double(1.134e+00);
    end;
    if W=121 then
    begin
        R := -Double(1.168e+00);
    end;
    if W=120 then
    begin
        R := -Double(1.202e+00);
    end;
    if W=119 then
    begin
        R := -Double(1.237e+00);
    end;
    if W=118 then
    begin
        R := -Double(1.273e+00);
    end;
    if W=117 then
    begin
        R := -Double(1.309e+00);
    end;
    if W=116 then
    begin
        R := -Double(1.347e+00);
    end;
    if W=115 then
    begin
        R := -Double(1.384e+00);
    end;
    if W=114 then
    begin
        R := -Double(1.423e+00);
    end;
    if W=113 then
    begin
        R := -Double(1.462e+00);
    end;
    if W=112 then
    begin
        R := -Double(1.502e+00);
    end;
    if W=111 then
    begin
        R := -Double(1.543e+00);
    end;
    if W=110 then
    begin
        R := -Double(1.585e+00);
    end;
    if W=109 then
    begin
        R := -Double(1.627e+00);
    end;
    if W=108 then
    begin
        R := -Double(1.670e+00);
    end;
    if W=107 then
    begin
        R := -Double(1.714e+00);
    end;
    if W=106 then
    begin
        R := -Double(1.758e+00);
    end;
    if W=105 then
    begin
        R := -Double(1.804e+00);
    end;
    if W=104 then
    begin
        R := -Double(1.850e+00);
    end;
    if W=103 then
    begin
        R := -Double(1.897e+00);
    end;
    if W=102 then
    begin
        R := -Double(1.944e+00);
    end;
    if W=101 then
    begin
        R := -Double(1.993e+00);
    end;
    if W=100 then
    begin
        R := -Double(2.042e+00);
    end;
    if W=99 then
    begin
        R := -Double(2.093e+00);
    end;
    if W=98 then
    begin
        R := -Double(2.144e+00);
    end;
    if W=97 then
    begin
        R := -Double(2.195e+00);
    end;
    if W=96 then
    begin
        R := -Double(2.248e+00);
    end;
    if W=95 then
    begin
        R := -Double(2.302e+00);
    end;
    if W=94 then
    begin
        R := -Double(2.356e+00);
    end;
    if W=93 then
    begin
        R := -Double(2.412e+00);
    end;
    if W=92 then
    begin
        R := -Double(2.468e+00);
    end;
    if W=91 then
    begin
        R := -Double(2.525e+00);
    end;
    if W=90 then
    begin
        R := -Double(2.583e+00);
    end;
    if W=89 then
    begin
        R := -Double(2.642e+00);
    end;
    if W=88 then
    begin
        R := -Double(2.702e+00);
    end;
    if W=87 then
    begin
        R := -Double(2.763e+00);
    end;
    if W=86 then
    begin
        R := -Double(2.825e+00);
    end;
    if W=85 then
    begin
        R := -Double(2.888e+00);
    end;
    if W=84 then
    begin
        R := -Double(2.951e+00);
    end;
    if W=83 then
    begin
        R := -Double(3.016e+00);
    end;
    if W=82 then
    begin
        R := -Double(3.082e+00);
    end;
    if W=81 then
    begin
        R := -Double(3.149e+00);
    end;
    if W=80 then
    begin
        R := -Double(3.216e+00);
    end;
    if W=79 then
    begin
        R := -Double(3.285e+00);
    end;
    if W=78 then
    begin
        R := -Double(3.355e+00);
    end;
    if W=77 then
    begin
        R := -Double(3.426e+00);
    end;
    if W=76 then
    begin
        R := -Double(3.498e+00);
    end;
    if W=75 then
    begin
        R := -Double(3.571e+00);
    end;
    if W=74 then
    begin
        R := -Double(3.645e+00);
    end;
    if W=73 then
    begin
        R := -Double(3.721e+00);
    end;
    if W=72 then
    begin
        R := -Double(3.797e+00);
    end;
    if W=71 then
    begin
        R := -Double(3.875e+00);
    end;
    if W=70 then
    begin
        R := -Double(3.953e+00);
    end;
    if W=69 then
    begin
        R := -Double(4.033e+00);
    end;
    if W=68 then
    begin
        R := -Double(4.114e+00);
    end;
    if W=67 then
    begin
        R := -Double(4.197e+00);
    end;
    if W=66 then
    begin
        R := -Double(4.280e+00);
    end;
    if W=65 then
    begin
        R := -Double(4.365e+00);
    end;
    if W=64 then
    begin
        R := -Double(4.451e+00);
    end;
    if W=63 then
    begin
        R := -Double(4.539e+00);
    end;
    if W=62 then
    begin
        R := -Double(4.628e+00);
    end;
    if W=61 then
    begin
        R := -Double(4.718e+00);
    end;
    if W=60 then
    begin
        R := -Double(4.809e+00);
    end;
    if W=59 then
    begin
        R := -Double(4.902e+00);
    end;
    if W=58 then
    begin
        R := -Double(4.996e+00);
    end;
    if W=57 then
    begin
        R := -Double(5.092e+00);
    end;
    if W=56 then
    begin
        R := -Double(5.189e+00);
    end;
    if W=55 then
    begin
        R := -Double(5.287e+00);
    end;
    if W=54 then
    begin
        R := -Double(5.388e+00);
    end;
    if W=53 then
    begin
        R := -Double(5.489e+00);
    end;
    if W=52 then
    begin
        R := -Double(5.592e+00);
    end;
    if W=51 then
    begin
        R := -Double(5.697e+00);
    end;
    if W=50 then
    begin
        R := -Double(5.804e+00);
    end;
    if W=49 then
    begin
        R := -Double(5.912e+00);
    end;
    if W=48 then
    begin
        R := -Double(6.022e+00);
    end;
    if W=47 then
    begin
        R := -Double(6.133e+00);
    end;
    if W=46 then
    begin
        R := -Double(6.247e+00);
    end;
    if W=45 then
    begin
        R := -Double(6.362e+00);
    end;
    if W=44 then
    begin
        R := -Double(6.479e+00);
    end;
    if W=43 then
    begin
        R := -Double(6.598e+00);
    end;
    if W=42 then
    begin
        R := -Double(6.719e+00);
    end;
    if W=41 then
    begin
        R := -Double(6.842e+00);
    end;
    if W=40 then
    begin
        R := -Double(6.967e+00);
    end;
    if W=39 then
    begin
        R := -Double(7.094e+00);
    end;
    if W=38 then
    begin
        R := -Double(7.224e+00);
    end;
    if W=37 then
    begin
        R := -Double(7.355e+00);
    end;
    if W=36 then
    begin
        R := -Double(7.489e+00);
    end;
    if W=35 then
    begin
        R := -Double(7.625e+00);
    end;
    if W=34 then
    begin
        R := -Double(7.764e+00);
    end;
    if W=33 then
    begin
        R := -Double(7.905e+00);
    end;
    if W=32 then
    begin
        R := -Double(8.049e+00);
    end;
    if W=31 then
    begin
        R := -Double(8.196e+00);
    end;
    if W=30 then
    begin
        R := -Double(8.345e+00);
    end;
    if W=29 then
    begin
        R := -Double(8.498e+00);
    end;
    if W=28 then
    begin
        R := -Double(8.653e+00);
    end;
    if W=27 then
    begin
        R := -Double(8.811e+00);
    end;
    if W=26 then
    begin
        R := -Double(8.974e+00);
    end;
    if W=25 then
    begin
        R := -Double(9.139e+00);
    end;
    if W=24 then
    begin
        R := -Double(9.308e+00);
    end;
    if W=23 then
    begin
        R := -Double(9.481e+00);
    end;
    if W=22 then
    begin
        R := -Double(9.658e+00);
    end;
    if W=21 then
    begin
        R := -Double(9.840e+00);
    end;
    if W=20 then
    begin
        R := -Double(1.003e+01);
    end;
    if W=19 then
    begin
        R := -Double(1.022e+01);
    end;
    if W=18 then
    begin
        R := -Double(1.041e+01);
    end;
    if W=17 then
    begin
        R := -Double(1.061e+01);
    end;
    if W=16 then
    begin
        R := -Double(1.081e+01);
    end;
    if W=15 then
    begin
        R := -Double(1.102e+01);
    end;
    if W=14 then
    begin
        R := -Double(1.124e+01);
    end;
    if W=13 then
    begin
        R := -Double(1.147e+01);
    end;
    if W=12 then
    begin
        R := -Double(1.169e+01);
    end;
    if W=11 then
    begin
        R := -Double(1.194e+01);
    end;
    if W=10 then
    begin
        R := -Double(1.218e+01);
    end;
    if W=9 then
    begin
        R := -Double(1.245e+01);
    end;
    if W=8 then
    begin
        R := -Double(1.272e+01);
    end;
    if W=7 then
    begin
        R := -Double(1.300e+01);
    end;
    if W=6 then
    begin
        R := -Double(1.330e+01);
    end;
    if W=5 then
    begin
        R := -Double(1.364e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.400e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.433e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.484e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.525e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.594e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 24)
*************************************************************************)
function W24(S : Double):Double;
var
    W : AlglibInteger;
    R : Double;
begin
    W := Round(-Double(3.500000e+01)*S+Double(1.500000e+02));
    if W>=150 then
    begin
        R := -Double(6.820e-01);
    end;
    if W=149 then
    begin
        R := -Double(7.044e-01);
    end;
    if W=148 then
    begin
        R := -Double(7.273e-01);
    end;
    if W=147 then
    begin
        R := -Double(7.507e-01);
    end;
    if W=146 then
    begin
        R := -Double(7.746e-01);
    end;
    if W=145 then
    begin
        R := -Double(7.990e-01);
    end;
    if W=144 then
    begin
        R := -Double(8.239e-01);
    end;
    if W=143 then
    begin
        R := -Double(8.494e-01);
    end;
    if W=142 then
    begin
        R := -Double(8.754e-01);
    end;
    if W=141 then
    begin
        R := -Double(9.020e-01);
    end;
    if W=140 then
    begin
        R := -Double(9.291e-01);
    end;
    if W=139 then
    begin
        R := -Double(9.567e-01);
    end;
    if W=138 then
    begin
        R := -Double(9.849e-01);
    end;
    if W=137 then
    begin
        R := -Double(1.014e+00);
    end;
    if W=136 then
    begin
        R := -Double(1.043e+00);
    end;
    if W=135 then
    begin
        R := -Double(1.073e+00);
    end;
    if W=134 then
    begin
        R := -Double(1.103e+00);
    end;
    if W=133 then
    begin
        R := -Double(1.135e+00);
    end;
    if W=132 then
    begin
        R := -Double(1.166e+00);
    end;
    if W=131 then
    begin
        R := -Double(1.198e+00);
    end;
    if W=130 then
    begin
        R := -Double(1.231e+00);
    end;
    if W=129 then
    begin
        R := -Double(1.265e+00);
    end;
    if W=128 then
    begin
        R := -Double(1.299e+00);
    end;
    if W=127 then
    begin
        R := -Double(1.334e+00);
    end;
    if W=126 then
    begin
        R := -Double(1.369e+00);
    end;
    if W=125 then
    begin
        R := -Double(1.405e+00);
    end;
    if W=124 then
    begin
        R := -Double(1.441e+00);
    end;
    if W=123 then
    begin
        R := -Double(1.479e+00);
    end;
    if W=122 then
    begin
        R := -Double(1.517e+00);
    end;
    if W=121 then
    begin
        R := -Double(1.555e+00);
    end;
    if W=120 then
    begin
        R := -Double(1.594e+00);
    end;
    if W=119 then
    begin
        R := -Double(1.634e+00);
    end;
    if W=118 then
    begin
        R := -Double(1.675e+00);
    end;
    if W=117 then
    begin
        R := -Double(1.716e+00);
    end;
    if W=116 then
    begin
        R := -Double(1.758e+00);
    end;
    if W=115 then
    begin
        R := -Double(1.800e+00);
    end;
    if W=114 then
    begin
        R := -Double(1.844e+00);
    end;
    if W=113 then
    begin
        R := -Double(1.888e+00);
    end;
    if W=112 then
    begin
        R := -Double(1.932e+00);
    end;
    if W=111 then
    begin
        R := -Double(1.978e+00);
    end;
    if W=110 then
    begin
        R := -Double(2.024e+00);
    end;
    if W=109 then
    begin
        R := -Double(2.070e+00);
    end;
    if W=108 then
    begin
        R := -Double(2.118e+00);
    end;
    if W=107 then
    begin
        R := -Double(2.166e+00);
    end;
    if W=106 then
    begin
        R := -Double(2.215e+00);
    end;
    if W=105 then
    begin
        R := -Double(2.265e+00);
    end;
    if W=104 then
    begin
        R := -Double(2.316e+00);
    end;
    if W=103 then
    begin
        R := -Double(2.367e+00);
    end;
    if W=102 then
    begin
        R := -Double(2.419e+00);
    end;
    if W=101 then
    begin
        R := -Double(2.472e+00);
    end;
    if W=100 then
    begin
        R := -Double(2.526e+00);
    end;
    if W=99 then
    begin
        R := -Double(2.580e+00);
    end;
    if W=98 then
    begin
        R := -Double(2.636e+00);
    end;
    if W=97 then
    begin
        R := -Double(2.692e+00);
    end;
    if W=96 then
    begin
        R := -Double(2.749e+00);
    end;
    if W=95 then
    begin
        R := -Double(2.806e+00);
    end;
    if W=94 then
    begin
        R := -Double(2.865e+00);
    end;
    if W=93 then
    begin
        R := -Double(2.925e+00);
    end;
    if W=92 then
    begin
        R := -Double(2.985e+00);
    end;
    if W=91 then
    begin
        R := -Double(3.046e+00);
    end;
    if W=90 then
    begin
        R := -Double(3.108e+00);
    end;
    if W=89 then
    begin
        R := -Double(3.171e+00);
    end;
    if W=88 then
    begin
        R := -Double(3.235e+00);
    end;
    if W=87 then
    begin
        R := -Double(3.300e+00);
    end;
    if W=86 then
    begin
        R := -Double(3.365e+00);
    end;
    if W=85 then
    begin
        R := -Double(3.432e+00);
    end;
    if W=84 then
    begin
        R := -Double(3.499e+00);
    end;
    if W=83 then
    begin
        R := -Double(3.568e+00);
    end;
    if W=82 then
    begin
        R := -Double(3.637e+00);
    end;
    if W=81 then
    begin
        R := -Double(3.708e+00);
    end;
    if W=80 then
    begin
        R := -Double(3.779e+00);
    end;
    if W=79 then
    begin
        R := -Double(3.852e+00);
    end;
    if W=78 then
    begin
        R := -Double(3.925e+00);
    end;
    if W=77 then
    begin
        R := -Double(4.000e+00);
    end;
    if W=76 then
    begin
        R := -Double(4.075e+00);
    end;
    if W=75 then
    begin
        R := -Double(4.151e+00);
    end;
    if W=74 then
    begin
        R := -Double(4.229e+00);
    end;
    if W=73 then
    begin
        R := -Double(4.308e+00);
    end;
    if W=72 then
    begin
        R := -Double(4.387e+00);
    end;
    if W=71 then
    begin
        R := -Double(4.468e+00);
    end;
    if W=70 then
    begin
        R := -Double(4.550e+00);
    end;
    if W=69 then
    begin
        R := -Double(4.633e+00);
    end;
    if W=68 then
    begin
        R := -Double(4.718e+00);
    end;
    if W=67 then
    begin
        R := -Double(4.803e+00);
    end;
    if W=66 then
    begin
        R := -Double(4.890e+00);
    end;
    if W=65 then
    begin
        R := -Double(4.978e+00);
    end;
    if W=64 then
    begin
        R := -Double(5.067e+00);
    end;
    if W=63 then
    begin
        R := -Double(5.157e+00);
    end;
    if W=62 then
    begin
        R := -Double(5.249e+00);
    end;
    if W=61 then
    begin
        R := -Double(5.342e+00);
    end;
    if W=60 then
    begin
        R := -Double(5.436e+00);
    end;
    if W=59 then
    begin
        R := -Double(5.531e+00);
    end;
    if W=58 then
    begin
        R := -Double(5.628e+00);
    end;
    if W=57 then
    begin
        R := -Double(5.727e+00);
    end;
    if W=56 then
    begin
        R := -Double(5.826e+00);
    end;
    if W=55 then
    begin
        R := -Double(5.927e+00);
    end;
    if W=54 then
    begin
        R := -Double(6.030e+00);
    end;
    if W=53 then
    begin
        R := -Double(6.134e+00);
    end;
    if W=52 then
    begin
        R := -Double(6.240e+00);
    end;
    if W=51 then
    begin
        R := -Double(6.347e+00);
    end;
    if W=50 then
    begin
        R := -Double(6.456e+00);
    end;
    if W=49 then
    begin
        R := -Double(6.566e+00);
    end;
    if W=48 then
    begin
        R := -Double(6.678e+00);
    end;
    if W=47 then
    begin
        R := -Double(6.792e+00);
    end;
    if W=46 then
    begin
        R := -Double(6.907e+00);
    end;
    if W=45 then
    begin
        R := -Double(7.025e+00);
    end;
    if W=44 then
    begin
        R := -Double(7.144e+00);
    end;
    if W=43 then
    begin
        R := -Double(7.265e+00);
    end;
    if W=42 then
    begin
        R := -Double(7.387e+00);
    end;
    if W=41 then
    begin
        R := -Double(7.512e+00);
    end;
    if W=40 then
    begin
        R := -Double(7.639e+00);
    end;
    if W=39 then
    begin
        R := -Double(7.768e+00);
    end;
    if W=38 then
    begin
        R := -Double(7.899e+00);
    end;
    if W=37 then
    begin
        R := -Double(8.032e+00);
    end;
    if W=36 then
    begin
        R := -Double(8.167e+00);
    end;
    if W=35 then
    begin
        R := -Double(8.305e+00);
    end;
    if W=34 then
    begin
        R := -Double(8.445e+00);
    end;
    if W=33 then
    begin
        R := -Double(8.588e+00);
    end;
    if W=32 then
    begin
        R := -Double(8.733e+00);
    end;
    if W=31 then
    begin
        R := -Double(8.881e+00);
    end;
    if W=30 then
    begin
        R := -Double(9.031e+00);
    end;
    if W=29 then
    begin
        R := -Double(9.185e+00);
    end;
    if W=28 then
    begin
        R := -Double(9.341e+00);
    end;
    if W=27 then
    begin
        R := -Double(9.501e+00);
    end;
    if W=26 then
    begin
        R := -Double(9.664e+00);
    end;
    if W=25 then
    begin
        R := -Double(9.830e+00);
    end;
    if W=24 then
    begin
        R := -Double(1.000e+01);
    end;
    if W=23 then
    begin
        R := -Double(1.017e+01);
    end;
    if W=22 then
    begin
        R := -Double(1.035e+01);
    end;
    if W=21 then
    begin
        R := -Double(1.053e+01);
    end;
    if W=20 then
    begin
        R := -Double(1.072e+01);
    end;
    if W=19 then
    begin
        R := -Double(1.091e+01);
    end;
    if W=18 then
    begin
        R := -Double(1.110e+01);
    end;
    if W=17 then
    begin
        R := -Double(1.130e+01);
    end;
    if W=16 then
    begin
        R := -Double(1.151e+01);
    end;
    if W=15 then
    begin
        R := -Double(1.172e+01);
    end;
    if W=14 then
    begin
        R := -Double(1.194e+01);
    end;
    if W=13 then
    begin
        R := -Double(1.216e+01);
    end;
    if W=12 then
    begin
        R := -Double(1.239e+01);
    end;
    if W=11 then
    begin
        R := -Double(1.263e+01);
    end;
    if W=10 then
    begin
        R := -Double(1.287e+01);
    end;
    if W=9 then
    begin
        R := -Double(1.314e+01);
    end;
    if W=8 then
    begin
        R := -Double(1.342e+01);
    end;
    if W=7 then
    begin
        R := -Double(1.369e+01);
    end;
    if W=6 then
    begin
        R := -Double(1.400e+01);
    end;
    if W=5 then
    begin
        R := -Double(1.433e+01);
    end;
    if W=4 then
    begin
        R := -Double(1.469e+01);
    end;
    if W=3 then
    begin
        R := -Double(1.503e+01);
    end;
    if W=2 then
    begin
        R := -Double(1.554e+01);
    end;
    if W=1 then
    begin
        R := -Double(1.594e+01);
    end;
    if W<=0 then
    begin
        R := -Double(1.664e+01);
    end;
    Result := R;
end;


(*************************************************************************
Tail(S, 25)
*************************************************************************)
function W25(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(5.150509e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.695528e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.437637e+00), TJ, TJ1, Result);
    WCheb(X, -Double(2.611906e-01), TJ, TJ1, Result);
    WCheb(X, -Double(7.625722e-02), TJ, TJ1, Result);
    WCheb(X, -Double(2.579892e-02), TJ, TJ1, Result);
    WCheb(X, -Double(1.086876e-02), TJ, TJ1, Result);
    WCheb(X, -Double(2.906543e-03), TJ, TJ1, Result);
    WCheb(X, -Double(2.354881e-03), TJ, TJ1, Result);
    WCheb(X, Double(1.007195e-04), TJ, TJ1, Result);
    WCheb(X, -Double(8.437327e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 26)
*************************************************************************)
function W26(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(5.117622e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.635159e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.395167e+00), TJ, TJ1, Result);
    WCheb(X, -Double(2.382823e-01), TJ, TJ1, Result);
    WCheb(X, -Double(6.531987e-02), TJ, TJ1, Result);
    WCheb(X, -Double(2.060112e-02), TJ, TJ1, Result);
    WCheb(X, -Double(8.203697e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.516523e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.431364e-03), TJ, TJ1, Result);
    WCheb(X, Double(6.384553e-04), TJ, TJ1, Result);
    WCheb(X, -Double(3.238369e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 27)
*************************************************************************)
function W27(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(5.089731e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.584248e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.359966e+00), TJ, TJ1, Result);
    WCheb(X, -Double(2.203696e-01), TJ, TJ1, Result);
    WCheb(X, -Double(5.753344e-02), TJ, TJ1, Result);
    WCheb(X, -Double(1.761891e-02), TJ, TJ1, Result);
    WCheb(X, -Double(7.096897e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.419108e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.581214e-03), TJ, TJ1, Result);
    WCheb(X, Double(3.033766e-04), TJ, TJ1, Result);
    WCheb(X, -Double(5.901441e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 28)
*************************************************************************)
function W28(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(5.065046e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.539163e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.328939e+00), TJ, TJ1, Result);
    WCheb(X, -Double(2.046376e-01), TJ, TJ1, Result);
    WCheb(X, -Double(5.061515e-02), TJ, TJ1, Result);
    WCheb(X, -Double(1.469271e-02), TJ, TJ1, Result);
    WCheb(X, -Double(5.711578e-03), TJ, TJ1, Result);
    WCheb(X, -Double(8.389153e-04), TJ, TJ1, Result);
    WCheb(X, -Double(1.250575e-03), TJ, TJ1, Result);
    WCheb(X, Double(4.047245e-04), TJ, TJ1, Result);
    WCheb(X, -Double(5.128555e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 29)
*************************************************************************)
function W29(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(5.043413e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.499756e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.302137e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.915129e-01), TJ, TJ1, Result);
    WCheb(X, -Double(4.516329e-02), TJ, TJ1, Result);
    WCheb(X, -Double(1.260064e-02), TJ, TJ1, Result);
    WCheb(X, -Double(4.817269e-03), TJ, TJ1, Result);
    WCheb(X, -Double(5.478130e-04), TJ, TJ1, Result);
    WCheb(X, -Double(1.111668e-03), TJ, TJ1, Result);
    WCheb(X, Double(4.093451e-04), TJ, TJ1, Result);
    WCheb(X, -Double(5.135860e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 30)
*************************************************************************)
function W30(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(5.024071e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.464515e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.278342e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.800030e-01), TJ, TJ1, Result);
    WCheb(X, -Double(4.046294e-02), TJ, TJ1, Result);
    WCheb(X, -Double(1.076162e-02), TJ, TJ1, Result);
    WCheb(X, -Double(3.968677e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.911679e-04), TJ, TJ1, Result);
    WCheb(X, -Double(8.619185e-04), TJ, TJ1, Result);
    WCheb(X, Double(5.125362e-04), TJ, TJ1, Result);
    WCheb(X, -Double(3.984370e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 40)
*************************************************************************)
function W40(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(4.904809e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.248327e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.136698e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.170982e-01), TJ, TJ1, Result);
    WCheb(X, -Double(1.824427e-02), TJ, TJ1, Result);
    WCheb(X, -Double(3.888648e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.344929e-03), TJ, TJ1, Result);
    WCheb(X, Double(2.790407e-04), TJ, TJ1, Result);
    WCheb(X, -Double(4.619858e-04), TJ, TJ1, Result);
    WCheb(X, Double(3.359121e-04), TJ, TJ1, Result);
    WCheb(X, -Double(2.883026e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 60)
*************************************************************************)
function W60(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(4.809656e+00), TJ, TJ1, Result);
    WCheb(X, -Double(5.077191e+00), TJ, TJ1, Result);
    WCheb(X, -Double(1.029402e+00), TJ, TJ1, Result);
    WCheb(X, -Double(7.507931e-02), TJ, TJ1, Result);
    WCheb(X, -Double(6.506226e-03), TJ, TJ1, Result);
    WCheb(X, -Double(1.391278e-03), TJ, TJ1, Result);
    WCheb(X, -Double(4.263635e-04), TJ, TJ1, Result);
    WCheb(X, Double(2.302271e-04), TJ, TJ1, Result);
    WCheb(X, -Double(2.384348e-04), TJ, TJ1, Result);
    WCheb(X, Double(1.865587e-04), TJ, TJ1, Result);
    WCheb(X, -Double(1.622355e-04), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 120)
*************************************************************************)
function W120(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(4.729426e+00), TJ, TJ1, Result);
    WCheb(X, -Double(4.934426e+00), TJ, TJ1, Result);
    WCheb(X, -Double(9.433231e-01), TJ, TJ1, Result);
    WCheb(X, -Double(4.492504e-02), TJ, TJ1, Result);
    WCheb(X, Double(1.673948e-05), TJ, TJ1, Result);
    WCheb(X, -Double(6.077014e-04), TJ, TJ1, Result);
    WCheb(X, -Double(7.215768e-05), TJ, TJ1, Result);
    WCheb(X, Double(9.086734e-05), TJ, TJ1, Result);
    WCheb(X, -Double(8.447980e-05), TJ, TJ1, Result);
    WCheb(X, Double(6.705028e-05), TJ, TJ1, Result);
    WCheb(X, -Double(5.828507e-05), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S, 200)
*************************************************************************)
function W200(S : Double):Double;
var
    X : Double;
    TJ : Double;
    TJ1 : Double;
begin
    Result := 0;
    X := Min(2*(S-Double(0.000000e+00))/Double(4.000000e+00)-1, Double(1.0));
    TJ := 1;
    TJ1 := X;
    WCheb(X, -Double(4.700240e+00), TJ, TJ1, Result);
    WCheb(X, -Double(4.883080e+00), TJ, TJ1, Result);
    WCheb(X, -Double(9.132168e-01), TJ, TJ1, Result);
    WCheb(X, -Double(3.512684e-02), TJ, TJ1, Result);
    WCheb(X, Double(1.726342e-03), TJ, TJ1, Result);
    WCheb(X, -Double(5.189796e-04), TJ, TJ1, Result);
    WCheb(X, -Double(1.628659e-06), TJ, TJ1, Result);
    WCheb(X, Double(4.261786e-05), TJ, TJ1, Result);
    WCheb(X, -Double(4.002498e-05), TJ, TJ1, Result);
    WCheb(X, Double(3.146287e-05), TJ, TJ1, Result);
    WCheb(X, -Double(2.727576e-05), TJ, TJ1, Result);
end;


(*************************************************************************
Tail(S,N), S>=0
*************************************************************************)
function WSigma(S : Double; N : AlglibInteger):Double;
var
    F0 : Double;
    F1 : Double;
    F2 : Double;
    F3 : Double;
    F4 : Double;
    X0 : Double;
    X1 : Double;
    X2 : Double;
    X3 : Double;
    X4 : Double;
    X : Double;
begin
    if N=5 then
    begin
        Result := W5(S);
    end;
    if N=6 then
    begin
        Result := W6(S);
    end;
    if N=7 then
    begin
        Result := W7(S);
    end;
    if N=8 then
    begin
        Result := W8(S);
    end;
    if N=9 then
    begin
        Result := W9(S);
    end;
    if N=10 then
    begin
        Result := W10(S);
    end;
    if N=11 then
    begin
        Result := W11(S);
    end;
    if N=12 then
    begin
        Result := W12(S);
    end;
    if N=13 then
    begin
        Result := W13(S);
    end;
    if N=14 then
    begin
        Result := W14(S);
    end;
    if N=15 then
    begin
        Result := W15(S);
    end;
    if N=16 then
    begin
        Result := W16(S);
    end;
    if N=17 then
    begin
        Result := W17(S);
    end;
    if N=18 then
    begin
        Result := W18(S);
    end;
    if N=19 then
    begin
        Result := W19(S);
    end;
    if N=20 then
    begin
        Result := W20(S);
    end;
    if N=21 then
    begin
        Result := W21(S);
    end;
    if N=22 then
    begin
        Result := W22(S);
    end;
    if N=23 then
    begin
        Result := W23(S);
    end;
    if N=24 then
    begin
        Result := W24(S);
    end;
    if N=25 then
    begin
        Result := W25(S);
    end;
    if N=26 then
    begin
        Result := W26(S);
    end;
    if N=27 then
    begin
        Result := W27(S);
    end;
    if N=28 then
    begin
        Result := W28(S);
    end;
    if N=29 then
    begin
        Result := W29(S);
    end;
    if N=30 then
    begin
        Result := W30(S);
    end;
    if N>30 then
    begin
        X := Double(1.0)/N;
        X0 := Double(1.0)/30;
        F0 := W30(S);
        X1 := Double(1.0)/40;
        F1 := W40(S);
        X2 := Double(1.0)/60;
        F2 := W60(S);
        X3 := Double(1.0)/120;
        F3 := W120(S);
        X4 := Double(1.0)/200;
        F4 := W200(S);
        F1 := ((X-X0)*F1-(X-X1)*F0)/(X1-X0);
        F2 := ((X-X0)*F2-(X-X2)*F0)/(X2-X0);
        F3 := ((X-X0)*F3-(X-X3)*F0)/(X3-X0);
        F4 := ((X-X0)*F4-(X-X4)*F0)/(X4-X0);
        F2 := ((X-X1)*F2-(X-X2)*F1)/(X2-X1);
        F3 := ((X-X1)*F3-(X-X3)*F1)/(X3-X1);
        F4 := ((X-X1)*F4-(X-X4)*F1)/(X4-X1);
        F3 := ((X-X2)*F3-(X-X3)*F2)/(X3-X2);
        F4 := ((X-X2)*F4-(X-X4)*F2)/(X4-X2);
        F4 := ((X-X3)*F4-(X-X4)*F3)/(X4-X3);
        Result := F4;
    end;
end;


end.