{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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
unit linreg;
interface
uses Math, Sysutils, Ap, descriptivestatistics, gammafunc, normaldistr, igammaf, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd;

type
LinearModel = record
    W : TReal1DArray;
end;


(*************************************************************************
LRReport structure contains additional information about linear model:
* C             -   covariation matrix,  array[0..NVars,0..NVars].
                    C[i,j] = Cov(A[i],A[j])
* RMSError      -   root mean square error on a training set
* AvgError      -   average error on a training set
* AvgRelError   -   average relative error on a training set (excluding
                    observations with zero function value).
* CVRMSError    -   leave-one-out cross-validation estimate of
                    generalization error. Calculated using fast algorithm
                    with O(NVars*NPoints) complexity.
* CVAvgError    -   cross-validation estimate of average error
* CVAvgRelError -   cross-validation estimate of average relative error

All other fields of the structure are intended for internal use and should
not be used outside ALGLIB.
*************************************************************************)
LRReport = record
    C : TReal2DArray;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    CVRMSError : Double;
    CVAvgError : Double;
    CVAvgRelError : Double;
    NCVDefects : AlglibInteger;
    CVDefects : TInteger1DArray;
end;



procedure LRBuild(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
procedure LRBuildS(const XY : TReal2DArray;
     const S : TReal1DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
procedure LRBuildZS(const XY : TReal2DArray;
     const S : TReal1DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
procedure LRBuildZ(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
procedure LRUnpack(const LM : LinearModel;
     var V : TReal1DArray;
     var NVars : AlglibInteger);
procedure LRPack(const V : TReal1DArray;
     NVars : AlglibInteger;
     var LM : LinearModel);
function LRProcess(const LM : LinearModel; const X : TReal1DArray):Double;
function LRRMSError(const LM : LinearModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function LRAvgError(const LM : LinearModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function LRAvgRelError(const LM : LinearModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
procedure LRCopy(const LM1 : LinearModel; var LM2 : LinearModel);
procedure LRSerialize(const LM : LinearModel;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
procedure LRUnserialize(const RA : TReal1DArray; var LM : LinearModel);
procedure LRLineS(const XY : TReal2DArray;
     const S : TReal1DArray;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var A : Double;
     var B : Double;
     var VarA : Double;
     var VarB : Double;
     var CovAB : Double;
     var CorrAB : Double;
     var P : Double);
procedure LRLine(const XY : TReal2DArray;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var A : Double;
     var B : Double);

implementation

const
    LRVNum = 5;

procedure LRInternal(const XY : TReal2DArray;
     const S : TReal1DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);forward;


(*************************************************************************
Linear regression

Subroutine builds model:

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)

and model found in ALGLIB format, covariation matrix, training set  errors
(rms,  average,  average  relative)   and  leave-one-out  cross-validation
estimate of the generalization error. CV  estimate calculated  using  fast
algorithm with O(NPoints*NVars) complexity.

When  covariation  matrix  is  calculated  standard deviations of function
values are assumed to be equal to RMS error on the training set.

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************)
procedure LRBuild(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
var
    S : TReal1DArray;
    I : AlglibInteger;
    Sigma2 : Double;
begin
    if (NPoints<=NVars+1) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(S, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        S[I] := 1;
        Inc(I);
    end;
    LRBuildS(XY, S, NPoints, NVars, Info, LM, AR);
    if Info<0 then
    begin
        Exit;
    end;
    Sigma2 := AP_Sqr(AR.RMSError)*NPoints/(NPoints-NVars-1);
    I:=0;
    while I<=NVars do
    begin
        APVMul(@AR.C[I][0], 0, NVars, Sigma2);
        Inc(I);
    end;
end;


(*************************************************************************
Linear regression

Variant of LRBuild which uses vector of standatd deviations (errors in
function values).

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    S           -   standard deviations (errors in function values)
                    array[0..NPoints-1], S[i]>0.
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    * -2, if S[I]<=0
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************)
procedure LRBuildS(const XY : TReal2DArray;
     const S : TReal1DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
var
    XYI : TReal2DArray;
    X : TReal1DArray;
    Means : TReal1DArray;
    Sigmas : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Offs : AlglibInteger;
    Mean : Double;
    Variance : Double;
    Skewness : Double;
    Kurtosis : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Test parameters
    //
    if (NPoints<=NVars+1) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // Copy data, add one more column (constant term)
    //
    SetLength(XYI, NPoints-1+1, NVars+1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@XYI[I][0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        XYI[I,NVars] := 1;
        XYI[I,NVars+1] := XY[I,NVars];
        Inc(I);
    end;
    
    //
    // Standartization
    //
    SetLength(X, NPoints-1+1);
    SetLength(Means, NVars-1+1);
    SetLength(Sigmas, NVars-1+1);
    J:=0;
    while J<=NVars-1 do
    begin
        for i_ := 0 to NPoints-1 do
        begin
            X[i_] := XY[i_,J];
        end;
        CalculateMoments(X, NPoints, Mean, Variance, Skewness, Kurtosis);
        Means[J] := Mean;
        Sigmas[J] := Sqrt(Variance);
        if AP_FP_Eq(Sigmas[J],0) then
        begin
            Sigmas[J] := 1;
        end;
        I:=0;
        while I<=NPoints-1 do
        begin
            XYI[I,J] := (XYI[I,J]-Means[J])/Sigmas[J];
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // Internal processing
    //
    LRInternal(XYI, S, NPoints, NVars+1, Info, LM, AR);
    if Info<0 then
    begin
        Exit;
    end;
    
    //
    // Un-standartization
    //
    Offs := Round(LM.W[3]);
    J:=0;
    while J<=NVars-1 do
    begin
        
        //
        // Constant term is updated (and its covariance too,
        // since it gets some variance from J-th component)
        //
        LM.W[Offs+NVars] := LM.W[Offs+NVars]-LM.W[Offs+J]*Means[J]/Sigmas[J];
        V := Means[J]/Sigmas[J];
        APVSub(@AR.C[NVars][0], 0, NVars, @AR.C[J][0], 0, NVars, V);
        for i_ := 0 to NVars do
        begin
            AR.C[i_,NVars] := AR.C[i_,NVars] - V*AR.C[i_,J];
        end;
        
        //
        // J-th term is updated
        //
        LM.W[Offs+J] := LM.W[Offs+J]/Sigmas[J];
        V := 1/Sigmas[J];
        APVMul(@AR.C[J][0], 0, NVars, V);
        for i_ := 0 to NVars do
        begin
            AR.C[i_,J] := V*AR.C[i_,J];
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************)
procedure LRBuildZS(const XY : TReal2DArray;
     const S : TReal1DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
var
    XYI : TReal2DArray;
    X : TReal1DArray;
    C : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Offs : AlglibInteger;
    Mean : Double;
    Variance : Double;
    Skewness : Double;
    Kurtosis : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Test parameters
    //
    if (NPoints<=NVars+1) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // Copy data, add one more column (constant term)
    //
    SetLength(XYI, NPoints-1+1, NVars+1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@XYI[I][0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        XYI[I,NVars] := 0;
        XYI[I,NVars+1] := XY[I,NVars];
        Inc(I);
    end;
    
    //
    // Standartization: unusual scaling
    //
    SetLength(X, NPoints-1+1);
    SetLength(C, NVars-1+1);
    J:=0;
    while J<=NVars-1 do
    begin
        for i_ := 0 to NPoints-1 do
        begin
            X[i_] := XY[i_,J];
        end;
        CalculateMoments(X, NPoints, Mean, Variance, Skewness, Kurtosis);
        if AP_FP_Greater(AbsReal(Mean),Sqrt(Variance)) then
        begin
            
            //
            // variation is relatively small, it is better to
            // bring mean value to 1
            //
            C[J] := Mean;
        end
        else
        begin
            
            //
            // variation is large, it is better to bring variance to 1
            //
            if AP_FP_Eq(Variance,0) then
            begin
                Variance := 1;
            end;
            C[J] := Sqrt(Variance);
        end;
        I:=0;
        while I<=NPoints-1 do
        begin
            XYI[I,J] := XYI[I,J]/C[J];
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // Internal processing
    //
    LRInternal(XYI, S, NPoints, NVars+1, Info, LM, AR);
    if Info<0 then
    begin
        Exit;
    end;
    
    //
    // Un-standartization
    //
    Offs := Round(LM.W[3]);
    J:=0;
    while J<=NVars-1 do
    begin
        
        //
        // J-th term is updated
        //
        LM.W[Offs+J] := LM.W[Offs+J]/C[J];
        V := 1/C[J];
        APVMul(@AR.C[J][0], 0, NVars, V);
        for i_ := 0 to NVars do
        begin
            AR.C[i_,J] := V*AR.C[i_,J];
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************)
procedure LRBuildZ(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
var
    S : TReal1DArray;
    I : AlglibInteger;
    Sigma2 : Double;
begin
    if (NPoints<=NVars+1) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(S, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        S[I] := 1;
        Inc(I);
    end;
    LRBuildZS(XY, S, NPoints, NVars, Info, LM, AR);
    if Info<0 then
    begin
        Exit;
    end;
    Sigma2 := AP_Sqr(AR.RMSError)*NPoints/(NPoints-NVars-1);
    I:=0;
    while I<=NVars do
    begin
        APVMul(@AR.C[I][0], 0, NVars, Sigma2);
        Inc(I);
    end;
end;


(*************************************************************************
Unpacks coefficients of linear model.

INPUT PARAMETERS:
    LM          -   linear model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables (one less than number
                    of coefficients)

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
procedure LRUnpack(const LM : LinearModel;
     var V : TReal1DArray;
     var NVars : AlglibInteger);
var
    Offs : AlglibInteger;
begin
    Assert(Round(LM.W[1])=LRVNum, 'LINREG: Incorrect LINREG version!');
    NVars := Round(LM.W[2]);
    Offs := Round(LM.W[3]);
    SetLength(V, NVars+1);
    APVMove(@V[0], 0, NVars, @LM.W[0], Offs, Offs+NVars);
end;


(*************************************************************************
"Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
reversed).

INPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables

OUTPUT PAREMETERS:
    LM          -   linear model.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
procedure LRPack(const V : TReal1DArray;
     NVars : AlglibInteger;
     var LM : LinearModel);
var
    Offs : AlglibInteger;
begin
    SetLength(LM.W, 4+NVars+1);
    Offs := 4;
    LM.W[0] := 4+NVars+1;
    LM.W[1] := LRVNum;
    LM.W[2] := NVars;
    LM.W[3] := Offs;
    APVMove(@LM.W[0], Offs, Offs+NVars, @V[0], 0, NVars);
end;


(*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   linear model
    X       -   input vector,  array[0..NVars-1].

Result:
    value of linear model regression estimate

  -- ALGLIB --
     Copyright 03.09.2008 by Bochkanov Sergey
*************************************************************************)
function LRProcess(const LM : LinearModel; const X : TReal1DArray):Double;
var
    V : Double;
    Offs : AlglibInteger;
    NVars : AlglibInteger;
begin
    Assert(Round(LM.W[1])=LRVNum, 'LINREG: Incorrect LINREG version!');
    NVars := Round(LM.W[2]);
    Offs := Round(LM.W[3]);
    V := APVDotProduct(@X[0], 0, NVars-1, @LM.W[0], Offs, Offs+NVars-1);
    Result := V+LM.W[Offs+NVars];
end;


(*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
function LRRMSError(const LM : LinearModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    I : AlglibInteger;
    V : Double;
    Offs : AlglibInteger;
    NVars : AlglibInteger;
begin
    Assert(Round(LM.W[1])=LRVNum, 'LINREG: Incorrect LINREG version!');
    NVars := Round(LM.W[2]);
    Offs := Round(LM.W[3]);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        V := APVDotProduct(@XY[I][0], 0, NVars-1, @LM.W[0], Offs, Offs+NVars-1);
        V := V+LM.W[Offs+NVars];
        Result := Result+AP_Sqr(V-XY[I,NVars]);
        Inc(I);
    end;
    Result := Sqrt(Result/NPoints);
end;


(*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
function LRAvgError(const LM : LinearModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    I : AlglibInteger;
    V : Double;
    Offs : AlglibInteger;
    NVars : AlglibInteger;
begin
    Assert(Round(LM.W[1])=LRVNum, 'LINREG: Incorrect LINREG version!');
    NVars := Round(LM.W[2]);
    Offs := Round(LM.W[3]);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        V := APVDotProduct(@XY[I][0], 0, NVars-1, @LM.W[0], Offs, Offs+NVars-1);
        V := V+LM.W[Offs+NVars];
        Result := Result+AbsReal(V-XY[I,NVars]);
        Inc(I);
    end;
    Result := Result/NPoints;
end;


(*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
function LRAvgRelError(const LM : LinearModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    Offs : AlglibInteger;
    NVars : AlglibInteger;
begin
    Assert(Round(LM.W[1])=LRVNum, 'LINREG: Incorrect LINREG version!');
    NVars := Round(LM.W[2]);
    Offs := Round(LM.W[3]);
    Result := 0;
    K := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        if AP_FP_Neq(XY[I,NVars],0) then
        begin
            V := APVDotProduct(@XY[I][0], 0, NVars-1, @LM.W[0], Offs, Offs+NVars-1);
            V := V+LM.W[Offs+NVars];
            Result := Result+AbsReal((V-XY[I,NVars])/XY[I,NVars]);
            K := K+1;
        end;
        Inc(I);
    end;
    if K<>0 then
    begin
        Result := Result/K;
    end;
end;


(*************************************************************************
Copying of LinearModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure LRCopy(const LM1 : LinearModel; var LM2 : LinearModel);
var
    K : AlglibInteger;
begin
    K := Round(LM1.W[0]);
    SetLength(LM2.W, K-1+1);
    APVMove(@LM2.W[0], 0, K-1, @LM1.W[0], 0, K-1);
end;


(*************************************************************************
Serialization of LinearModel strucure

INPUT PARAMETERS:
    LM      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores model,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure LRSerialize(const LM : LinearModel;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
begin
    RLen := Round(LM.W[0])+1;
    SetLength(RA, RLen-1+1);
    RA[0] := LRVNum;
    APVMove(@RA[0], 1, RLen-1, @LM.W[0], 0, RLen-2);
end;


(*************************************************************************
Unserialization of DecisionForest strucure

INPUT PARAMETERS:
    RA      -   real array which stores decision forest

OUTPUT PARAMETERS:
    LM      -   unserialized structure

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure LRUnserialize(const RA : TReal1DArray; var LM : LinearModel);
begin
    Assert(Round(RA[0])=LRVNum, 'LRUnserialize: incorrect array!');
    SetLength(LM.W, Round(RA[1])-1+1);
    APVMove(@LM.W[0], 0, Round(RA[1])-1, @RA[0], 1, Round(RA[1]));
end;


procedure LRLineS(const XY : TReal2DArray;
     const S : TReal1DArray;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var A : Double;
     var B : Double;
     var VarA : Double;
     var VarB : Double;
     var CovAB : Double;
     var CorrAB : Double;
     var P : Double);
var
    I : AlglibInteger;
    SS : Double;
    SX : Double;
    SXX : Double;
    SY : Double;
    STT : Double;
    E1 : Double;
    E2 : Double;
    T : Double;
    Chi2 : Double;
begin
    if N<2 then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Less_Eq(S[I],0) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Calculate S, SX, SY, SXX
    //
    SS := 0;
    SX := 0;
    SY := 0;
    SXX := 0;
    I:=0;
    while I<=N-1 do
    begin
        T := AP_Sqr(S[I]);
        SS := SS+1/T;
        SX := SX+XY[I,0]/T;
        SY := SY+XY[I,1]/T;
        SXX := SXX+AP_Sqr(XY[I,0])/T;
        Inc(I);
    end;
    
    //
    // Test for condition number
    //
    T := Sqrt(4*AP_Sqr(SX)+AP_Sqr(SS-SXX));
    E1 := Double(0.5)*(SS+SXX+T);
    E2 := Double(0.5)*(SS+SXX-T);
    if AP_FP_Less_Eq(Min(E1, E2),1000*MachineEpsilon*Max(E1, E2)) then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // Calculate A, B
    //
    A := 0;
    B := 0;
    STT := 0;
    I:=0;
    while I<=N-1 do
    begin
        T := (XY[I,0]-SX/SS)/S[I];
        B := B+T*XY[I,1]/S[I];
        STT := STT+AP_Sqr(T);
        Inc(I);
    end;
    B := B/STT;
    A := (SY-SX*B)/SS;
    
    //
    // Calculate goodness-of-fit
    //
    if N>2 then
    begin
        Chi2 := 0;
        I:=0;
        while I<=N-1 do
        begin
            Chi2 := Chi2+AP_Sqr((XY[I,1]-A-B*XY[I,0])/S[I]);
            Inc(I);
        end;
        P := IncompleteGammaC(AP_Double((N-2))/2, Chi2/2);
    end
    else
    begin
        P := 1;
    end;
    
    //
    // Calculate other parameters
    //
    VarA := (1+AP_Sqr(SX)/(SS*STT))/SS;
    VarB := 1/STT;
    CovAB := -SX/(SS*STT);
    CorrAB := CovAB/Sqrt(VarA*VarB);
end;


procedure LRLine(const XY : TReal2DArray;
     N : AlglibInteger;
     var Info : AlglibInteger;
     var A : Double;
     var B : Double);
var
    S : TReal1DArray;
    I : AlglibInteger;
    VarA : Double;
    VarB : Double;
    CovAB : Double;
    CorrAB : Double;
    P : Double;
begin
    if N<2 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(S, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        S[I] := 1;
        Inc(I);
    end;
    LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
end;


(*************************************************************************
Internal linear regression subroutine
*************************************************************************)
procedure LRInternal(const XY : TReal2DArray;
     const S : TReal1DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LinearModel;
     var AR : LRReport);
var
    A : TReal2DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    VM : TReal2DArray;
    XYM : TReal2DArray;
    B : TReal1DArray;
    SV : TReal1DArray;
    T : TReal1DArray;
    SVI : TReal1DArray;
    WORK : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    NCV : AlglibInteger;
    NA : AlglibInteger;
    NACV : AlglibInteger;
    R : Double;
    P : Double;
    EpsTol : Double;
    AR2 : LRReport;
    Offs : AlglibInteger;
    TLM : LinearModel;
    i_ : AlglibInteger;
begin
    EpsTol := 1000;
    
    //
    // Check for errors in data
    //
    if (NPoints<NVars) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        if AP_FP_Less_Eq(S[I],0) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Create design matrix
    //
    SetLength(A, NPoints-1+1, NVars-1+1);
    SetLength(B, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        R := 1/S[I];
        APVMove(@A[I][0], 0, NVars-1, @XY[I][0], 0, NVars-1, R);
        B[I] := XY[I,NVars]/S[I];
        Inc(I);
    end;
    
    //
    // Allocate W:
    // W[0]     array size
    // W[1]     version number, 0
    // W[2]     NVars (minus 1, to be compatible with external representation)
    // W[3]     coefficients offset
    //
    SetLength(LM.W, 4+NVars-1+1);
    Offs := 4;
    LM.W[0] := 4+NVars;
    LM.W[1] := LRVNum;
    LM.W[2] := NVars-1;
    LM.W[3] := Offs;
    
    //
    // Solve problem using SVD:
    //
    // 0. check for degeneracy (different types)
    // 1. A = U*diag(sv)*V'
    // 2. T = b'*U
    // 3. w = SUM((T[i]/sv[i])*V[..,i])
    // 4. cov(wi,wj) = SUM(Vji*Vjk/sv[i]^2,K=1..M)
    //
    // see $15.4 of "Numerical Recipes in C" for more information
    //
    SetLength(T, NVars-1+1);
    SetLength(SVI, NVars-1+1);
    SetLength(AR.C, NVars-1+1, NVars-1+1);
    SetLength(VM, NVars-1+1, NVars-1+1);
    if  not RMatrixSVD(A, NPoints, NVars, 1, 1, 2, SV, U, VT) then
    begin
        Info := -4;
        Exit;
    end;
    if AP_FP_Less_Eq(SV[0],0) then
    begin
        
        //
        // Degenerate case: zero design matrix.
        //
        I:=Offs;
        while I<=Offs+NVars-1 do
        begin
            LM.W[I] := 0;
            Inc(I);
        end;
        AR.RMSError := LRRMSError(LM, XY, NPoints);
        AR.AvgError := LRAvgError(LM, XY, NPoints);
        AR.AvgRelError := LRAvgRelError(LM, XY, NPoints);
        AR.CVRMSError := AR.RMSError;
        AR.CVAvgError := AR.AvgError;
        AR.CVAvgRelError := AR.AvgRelError;
        AR.NCVDefects := 0;
        SetLength(AR.CVDefects, NVars-1+1);
        SetLength(AR.C, NVars-1+1, NVars-1+1);
        I:=0;
        while I<=NVars-1 do
        begin
            J:=0;
            while J<=NVars-1 do
            begin
                AR.C[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        Exit;
    end;
    if AP_FP_Less_Eq(SV[NVars-1],EpsTol*MachineEpsilon*SV[0]) then
    begin
        
        //
        // Degenerate case, non-zero design matrix.
        //
        // We can leave it and solve task in SVD least squares fashion.
        // Solution and covariance matrix will be obtained correctly,
        // but CV error estimates - will not. It is better to reduce
        // it to non-degenerate task and to obtain correct CV estimates.
        //
        K:=NVars;
        while K>=1 do
        begin
            if AP_FP_Greater(SV[K-1],EpsTol*MachineEpsilon*SV[0]) then
            begin
                
                //
                // Reduce
                //
                SetLength(XYM, NPoints-1+1, K+1);
                I:=0;
                while I<=NPoints-1 do
                begin
                    J:=0;
                    while J<=K-1 do
                    begin
                        R := APVDotProduct(@XY[I][0], 0, NVars-1, @VT[J][0], 0, NVars-1);
                        XYM[I,J] := R;
                        Inc(J);
                    end;
                    XYM[I,K] := XY[I,NVars];
                    Inc(I);
                end;
                
                //
                // Solve
                //
                LRInternal(XYM, S, NPoints, K, Info, TLM, AR2);
                if Info<>1 then
                begin
                    Exit;
                end;
                
                //
                // Convert back to un-reduced format
                //
                J:=0;
                while J<=NVars-1 do
                begin
                    LM.W[Offs+J] := 0;
                    Inc(J);
                end;
                J:=0;
                while J<=K-1 do
                begin
                    R := TLM.W[Offs+J];
                    APVAdd(@LM.W[0], Offs, Offs+NVars-1, @VT[J][0], 0, NVars-1, R);
                    Inc(J);
                end;
                AR.RMSError := AR2.RMSError;
                AR.AvgError := AR2.AvgError;
                AR.AvgRelError := AR2.AvgRelError;
                AR.CVRMSError := AR2.CVRMSError;
                AR.CVAvgError := AR2.CVAvgError;
                AR.CVAvgRelError := AR2.CVAvgRelError;
                AR.NCVDefects := AR2.NCVDefects;
                SetLength(AR.CVDefects, NVars-1+1);
                J:=0;
                while J<=AR.NCVDefects-1 do
                begin
                    AR.CVDefects[J] := AR2.CVDefects[J];
                    Inc(J);
                end;
                SetLength(AR.C, NVars-1+1, NVars-1+1);
                SetLength(WORK, NVars+1);
                MatrixMatrixMultiply(AR2.C, 0, K-1, 0, K-1, False, VT, 0, K-1, 0, NVars-1, False, Double(1.0), VM, 0, K-1, 0, NVars-1, Double(0.0), WORK);
                MatrixMatrixMultiply(VT, 0, K-1, 0, NVars-1, True, VM, 0, K-1, 0, NVars-1, False, Double(1.0), AR.C, 0, NVars-1, 0, NVars-1, Double(0.0), WORK);
                Exit;
            end;
            Dec(K);
        end;
        Info := -255;
        Exit;
    end;
    I:=0;
    while I<=NVars-1 do
    begin
        if AP_FP_Greater(SV[I],EpsTol*MachineEpsilon*SV[0]) then
        begin
            SVI[I] := 1/SV[I];
        end
        else
        begin
            SVI[I] := 0;
        end;
        Inc(I);
    end;
    I:=0;
    while I<=NVars-1 do
    begin
        T[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        R := B[I];
        APVAdd(@T[0], 0, NVars-1, @U[I][0], 0, NVars-1, R);
        Inc(I);
    end;
    I:=0;
    while I<=NVars-1 do
    begin
        LM.W[Offs+I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=NVars-1 do
    begin
        R := T[I]*SVI[I];
        APVAdd(@LM.W[0], Offs, Offs+NVars-1, @VT[I][0], 0, NVars-1, R);
        Inc(I);
    end;
    J:=0;
    while J<=NVars-1 do
    begin
        R := SVI[J];
        for i_ := 0 to NVars-1 do
        begin
            VM[i_,J] := R*VT[J,i_];
        end;
        Inc(J);
    end;
    I:=0;
    while I<=NVars-1 do
    begin
        J:=I;
        while J<=NVars-1 do
        begin
            R := APVDotProduct(@VM[I][0], 0, NVars-1, @VM[J][0], 0, NVars-1);
            AR.C[I,J] := R;
            AR.C[J,I] := R;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Leave-1-out cross-validation error.
    //
    // NOTATIONS:
    // A            design matrix
    // A*x = b      original linear least squares task
    // U*S*V'       SVD of A
    // ai           i-th row of the A
    // bi           i-th element of the b
    // xf           solution of the original LLS task
    //
    // Cross-validation error of i-th element from a sample is
    // calculated using following formula:
    //
    //     ERRi = ai*xf - (ai*xf-bi*(ui*ui'))/(1-ui*ui')     (1)
    //
    // This formula can be derived from normal equations of the
    // original task
    //
    //     (A'*A)x = A'*b                                    (2)
    //
    // by applying modification (zeroing out i-th row of A) to (2):
    //
    //     (A-ai)'*(A-ai) = (A-ai)'*b
    //
    // and using Sherman-Morrison formula for updating matrix inverse
    //
    // NOTE 1: b is not zeroed out since it is much simpler and
    // does not influence final result.
    //
    // NOTE 2: some design matrices A have such ui that 1-ui*ui'=0.
    // Formula (1) can't be applied for such cases and they are skipped
    // from CV calculation (which distorts resulting CV estimate).
    // But from the properties of U we can conclude that there can
    // be no more than NVars such vectors. Usually
    // NVars << NPoints, so in a normal case it only slightly
    // influences result.
    //
    NCV := 0;
    NA := 0;
    NACV := 0;
    AR.RMSError := 0;
    AR.AvgError := 0;
    AR.AvgRelError := 0;
    AR.CVRMSError := 0;
    AR.CVAvgError := 0;
    AR.CVAvgRelError := 0;
    AR.NCVDefects := 0;
    SetLength(AR.CVDefects, NVars-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        
        //
        // Error on a training set
        //
        R := APVDotProduct(@XY[I][0], 0, NVars-1, @LM.W[0], Offs, Offs+NVars-1);
        AR.RMSError := AR.RMSError+AP_Sqr(R-XY[I,NVars]);
        AR.AvgError := AR.AvgError+AbsReal(R-XY[I,NVars]);
        if AP_FP_Neq(XY[I,NVars],0) then
        begin
            AR.AvgRelError := AR.AvgRelError+AbsReal((R-XY[I,NVars])/XY[I,NVars]);
            NA := NA+1;
        end;
        
        //
        // Error using fast leave-one-out cross-validation
        //
        P := APVDotProduct(@U[I][0], 0, NVars-1, @U[I][0], 0, NVars-1);
        if AP_FP_Greater(P,1-EpsTol*MachineEpsilon) then
        begin
            AR.CVDefects[AR.NCVDefects] := I;
            AR.NCVDefects := AR.NCVDefects+1;
            Inc(I);
            Continue;
        end;
        R := S[I]*(R/S[I]-B[I]*P)/(1-P);
        AR.CVRMSError := AR.CVRMSError+AP_Sqr(R-XY[I,NVars]);
        AR.CVAvgError := AR.CVAvgError+AbsReal(R-XY[I,NVars]);
        if AP_FP_Neq(XY[I,NVars],0) then
        begin
            AR.CVAvgRelError := AR.CVAvgRelError+AbsReal((R-XY[I,NVars])/XY[I,NVars]);
            NACV := NACV+1;
        end;
        NCV := NCV+1;
        Inc(I);
    end;
    if NCV=0 then
    begin
        
        //
        // Something strange: ALL ui are degenerate.
        // Unexpected...
        //
        Info := -255;
        Exit;
    end;
    AR.RMSError := Sqrt(AR.RMSError/NPoints);
    AR.AvgError := AR.AvgError/NPoints;
    if NA<>0 then
    begin
        AR.AvgRelError := AR.AvgRelError/NA;
    end;
    AR.CVRMSError := Sqrt(AR.CVRMSError/NCV);
    AR.CVAvgError := AR.CVAvgError/NCV;
    if NACV<>0 then
    begin
        AR.CVAvgRelError := AR.CVAvgRelError/NACV;
    end;
end;


end.