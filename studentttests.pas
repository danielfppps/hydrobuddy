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
unit studentttests;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, ibetaf, studenttdistr;

procedure StudentTTest1(const X : TReal1DArray;
     N : AlglibInteger;
     Mean : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
procedure StudentTTest2(const X : TReal1DArray;
     N : AlglibInteger;
     const Y : TReal1DArray;
     M : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
procedure UnequalVarianceTTest(const X : TReal1DArray;
     N : AlglibInteger;
     const Y : TReal1DArray;
     M : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);

implementation

(*************************************************************************
One-sample t-test

This test checks three hypotheses about the mean of the given sample.  The
following tests are performed:
    * two-tailed test (null hypothesis - the mean is equal  to  the  given
      value)
    * left-tailed test (null hypothesis - the  mean  is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis - the mean is less than or  equal
      to the given value).

The test is based on the assumption that  a  given  sample  has  a  normal
distribution and  an  unknown  dispersion.  If  the  distribution  sharply
differs from normal, the test will work incorrectly.

Input parameters:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of sample.
    Mean    -   assumed value of the mean.

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

  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure StudentTTest1(const X : TReal1DArray;
     N : AlglibInteger;
     Mean : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    XMean : Double;
    XVariance : Double;
    XStdDev : Double;
    V1 : Double;
    V2 : Double;
    Stat : Double;
    S : Double;
begin
    if N<=1 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Mean
    //
    XMean := 0;
    I:=0;
    while I<=N-1 do
    begin
        XMean := XMean+X[I];
        Inc(I);
    end;
    XMean := XMean/N;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    XVariance := 0;
    XStdDev := 0;
    if N<>1 then
    begin
        V1 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V1 := V1+AP_Sqr(X[I]-XMean);
            Inc(I);
        end;
        V2 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V2 := V2+(X[I]-XMean);
            Inc(I);
        end;
        V2 := AP_Sqr(V2)/N;
        XVariance := (V1-V2)/(N-1);
        if AP_FP_Less(XVariance,0) then
        begin
            XVariance := 0;
        end;
        XStdDev := Sqrt(XVariance);
    end;
    if AP_FP_Eq(XStdDev,0) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Statistic
    //
    Stat := (XMean-Mean)/(XStdDev/Sqrt(N));
    S := StudentTDistribution(N-1, Stat);
    BothTails := 2*Min(S, 1-S);
    LeftTail := S;
    RightTail := 1-S;
end;


(*************************************************************************
Two-sample pooled test

This test checks three hypotheses about the mean of the given samples. The
following tests are performed:
    * two-tailed test (null hypothesis - the means are equal)
    * left-tailed test (null hypothesis - the mean of the first sample  is
      greater than or equal to the mean of the second sample)
    * right-tailed test (null hypothesis - the mean of the first sample is
      less than or equal to the mean of the second sample).

Test is based on the following assumptions:
    * given samples have normal distributions
    * dispersions are equal
    * samples are independent.

Input parameters:
    X       -   sample 1. Array whose index goes from 0 to N-1.
    N       -   size of sample.
    Y       -   sample 2. Array whose index goes from 0 to M-1.
    M       -   size of sample.

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

  -- ALGLIB --
     Copyright 18.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure StudentTTest2(const X : TReal1DArray;
     N : AlglibInteger;
     const Y : TReal1DArray;
     M : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    XMean : Double;
    YMean : Double;
    Stat : Double;
    S : Double;
    P : Double;
begin
    if (N<=1) or (M<=1) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Mean
    //
    XMean := 0;
    I:=0;
    while I<=N-1 do
    begin
        XMean := XMean+X[I];
        Inc(I);
    end;
    XMean := XMean/N;
    YMean := 0;
    I:=0;
    while I<=M-1 do
    begin
        YMean := YMean+Y[I];
        Inc(I);
    end;
    YMean := YMean/M;
    
    //
    // S
    //
    S := 0;
    I:=0;
    while I<=N-1 do
    begin
        S := S+AP_Sqr(X[I]-XMean);
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        S := S+AP_Sqr(Y[I]-YMean);
        Inc(I);
    end;
    S := Sqrt(S*(AP_Double(1)/N+AP_Double(1)/M)/(N+M-2));
    if AP_FP_Eq(S,0) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Statistic
    //
    Stat := (XMean-YMean)/S;
    P := StudentTDistribution(N+M-2, Stat);
    BothTails := 2*Min(P, 1-P);
    LeftTail := P;
    RightTail := 1-P;
end;


(*************************************************************************
Two-sample unpooled test

This test checks three hypotheses about the mean of the given samples. The
following tests are performed:
    * two-tailed test (null hypothesis - the means are equal)
    * left-tailed test (null hypothesis - the mean of the first sample  is
      greater than or equal to the mean of the second sample)
    * right-tailed test (null hypothesis - the mean of the first sample is
      less than or equal to the mean of the second sample).

Test is based on the following assumptions:
    * given samples have normal distributions
    * samples are independent.
Dispersion equality is not required

Input parameters:
    X - sample 1. Array whose index goes from 0 to N-1.
    N - size of the sample.
    Y - sample 2. Array whose index goes from 0 to M-1.
    M - size of the sample.

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

  -- ALGLIB --
     Copyright 18.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure UnequalVarianceTTest(const X : TReal1DArray;
     N : AlglibInteger;
     const Y : TReal1DArray;
     M : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    XMean : Double;
    YMean : Double;
    XVar : Double;
    YVar : Double;
    DF : Double;
    P : Double;
    Stat : Double;
    C : Double;
begin
    if (N<=1) or (M<=1) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Mean
    //
    XMean := 0;
    I:=0;
    while I<=N-1 do
    begin
        XMean := XMean+X[I];
        Inc(I);
    end;
    XMean := XMean/N;
    YMean := 0;
    I:=0;
    while I<=M-1 do
    begin
        YMean := YMean+Y[I];
        Inc(I);
    end;
    YMean := YMean/M;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    XVar := 0;
    I:=0;
    while I<=N-1 do
    begin
        XVar := XVar+AP_Sqr(X[I]-XMean);
        Inc(I);
    end;
    XVar := XVar/(N-1);
    YVar := 0;
    I:=0;
    while I<=M-1 do
    begin
        YVar := YVar+AP_Sqr(Y[I]-YMean);
        Inc(I);
    end;
    YVar := YVar/(M-1);
    if AP_FP_Eq(XVar,0) or AP_FP_Eq(YVar,0) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Statistic
    //
    Stat := (XMean-YMean)/Sqrt(XVar/N+YVar/M);
    C := XVar/N/(XVar/N+YVar/M);
    DF := (N-1)*(M-1)/((M-1)*AP_Sqr(C)+(N-1)*(1-AP_Sqr(C)));
    if AP_FP_Greater(Stat,0) then
    begin
        P := 1-Double(0.5)*IncompleteBeta(DF/2, Double(0.5), DF/(DF+AP_Sqr(Stat)));
    end
    else
    begin
        P := Double(0.5)*IncompleteBeta(DF/2, Double(0.5), DF/(DF+AP_Sqr(Stat)));
    end;
    BothTails := 2*Min(P, 1-P);
    LeftTail := P;
    RightTail := 1-P;
end;


end.