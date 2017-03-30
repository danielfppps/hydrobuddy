{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2008, Sergey Bochkanov (ALGLIB project).

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
unit logit;
interface
uses Math, Sysutils, Ap, descriptivestatistics, mlpbase, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd, hqrnd, matgen, trfac, trlinsolve, safesolve, rcond, xblas, densesolver, tsort, bdss;

type
LogitModel = record
    W : TReal1DArray;
end;


LOGITMCState = record
    BRACKT : Boolean;
    STAGE1 : Boolean;
    INFOC : AlglibInteger;
    DG : Double;
    DGM : Double;
    DGINIT : Double;
    DGTEST : Double;
    DGX : Double;
    DGXM : Double;
    DGY : Double;
    DGYM : Double;
    FINIT : Double;
    FTEST1 : Double;
    FM : Double;
    FX : Double;
    FXM : Double;
    FY : Double;
    FYM : Double;
    STX : Double;
    STY : Double;
    STMIN : Double;
    STMAX : Double;
    WIDTH : Double;
    WIDTH1 : Double;
    XTRAPF : Double;
end;


(*************************************************************************
MNLReport structure contains information about training process:
* NGrad     -   number of gradient calculations
* NHess     -   number of Hessian calculations
*************************************************************************)
MNLReport = record
    NGrad : AlglibInteger;
    NHess : AlglibInteger;
end;



procedure MNLTrainH(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LogitModel;
     var Rep : MNLReport);
procedure MNLProcess(var LM : LogitModel;
     const X : TReal1DArray;
     var Y : TReal1DArray);
procedure MNLUnpack(const LM : LogitModel;
     var A : TReal2DArray;
     var NVars : AlglibInteger;
     var NClasses : AlglibInteger);
procedure MNLPack(const A : TReal2DArray;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var LM : LogitModel);
procedure MNLCopy(const LM1 : LogitModel; var LM2 : LogitModel);
procedure MNLSerialize(const LM : LogitModel;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
procedure MNLUnserialize(const RA : TReal1DArray; var LM : LogitModel);
function MNLAvgCE(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MNLRelClsError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MNLRMSError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MNLAvgError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
function MNLAvgRelError(var LM : LogitModel;
     const XY : TReal2DArray;
     SSize : AlglibInteger):Double;
function MNLClsError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):AlglibInteger;

implementation

const
    XTOL = 100*MachineEpsilon;
    FTOL = Double(0.0001);
    GTOL = Double(0.3);
    MAXFEV = 20;
    STPMIN = Double(1.0E-2);
    STPMAX = Double(1.0E5);
    LogitVNum = 6;

procedure MNLIExp(var W : TReal1DArray; const X : TReal1DArray);forward;
procedure MNLAllErrors(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     var RelCls : Double;
     var AvgCE : Double;
     var RMS : Double;
     var Avg : Double;
     var AvgRel : Double);forward;
procedure MNLMCSRCH(const N : AlglibInteger;
     var X : TReal1DArray;
     var F : Double;
     var G : TReal1DArray;
     const S : TReal1DArray;
     var STP : Double;
     var INFO : AlglibInteger;
     var NFEV : AlglibInteger;
     var WA : TReal1DArray;
     var State : LOGITMCState;
     var Stage : AlglibInteger);forward;
procedure MNLMCSTEP(var STX : Double;
     var FX : Double;
     var DX : Double;
     var STY : Double;
     var FY : Double;
     var DY : Double;
     var STP : Double;
     const FP : Double;
     const DP : Double;
     var BRACKT : Boolean;
     const STMIN : Double;
     const STMAX : Double;
     var INFO : AlglibInteger);forward;


(*************************************************************************
This subroutine trains logit model.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars]
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<NVars+2, NVars<1, NClasses<2).
                    *  1, if task has been solved
    LM          -   model built
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
procedure MNLTrainH(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var LM : LogitModel;
     var Rep : MNLReport);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    SSize : AlglibInteger;
    AllSame : Boolean;
    Offs : AlglibInteger;
    Threshold : Double;
    WMinStep : Double;
    Decay : Double;
    WDim : AlglibInteger;
    ExpOffs : AlglibInteger;
    V : Double;
    S : Double;
    Network : MultiLayerPerceptron;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    WCount : AlglibInteger;
    E : Double;
    G : TReal1DArray;
    H : TReal2DArray;
    SPD : Boolean;
    X : TReal1DArray;
    Y : TReal1DArray;
    WBase : TReal1DArray;
    WStep : Double;
    WDir : TReal1DArray;
    WORK : TReal1DArray;
    MCStage : AlglibInteger;
    MCState : LOGITMCState;
    MCInfo : AlglibInteger;
    MCNFEV : AlglibInteger;
    SolverInfo : AlglibInteger;
    SolverRep : DenseSolverReport;
begin
    Threshold := 1000*MachineEpsilon;
    WMinStep := Double(0.001);
    Decay := Double(0.001);
    
    //
    // Test for inputs
    //
    if (NPoints<NVars+2) or (NVars<1) or (NClasses<2) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        if (Round(XY[I,NVars])<0) or (Round(XY[I,NVars])>=NClasses) then
        begin
            Info := -2;
            Exit;
        end;
        Inc(I);
    end;
    Info := 1;
    
    //
    // Initialize data
    //
    Rep.NGrad := 0;
    Rep.NHess := 0;
    
    //
    // Allocate array
    //
    WDim := (NVars+1)*(NClasses-1);
    Offs := 5;
    ExpOffs := Offs+WDim;
    SSize := 5+(NVars+1)*(NClasses-1)+NClasses;
    SetLength(LM.W, SSize-1+1);
    LM.W[0] := SSize;
    LM.W[1] := LogitVNum;
    LM.W[2] := NVars;
    LM.W[3] := NClasses;
    LM.W[4] := Offs;
    
    //
    // Degenerate case: all outputs are equal
    //
    AllSame := True;
    I:=1;
    while I<=NPoints-1 do
    begin
        if Round(XY[I,NVars])<>Round(XY[I-1,NVars]) then
        begin
            AllSame := False;
        end;
        Inc(I);
    end;
    if AllSame then
    begin
        I:=0;
        while I<=(NVars+1)*(NClasses-1)-1 do
        begin
            LM.W[Offs+I] := 0;
            Inc(I);
        end;
        V := -2*Ln(MinRealNumber);
        K := Round(XY[0,NVars]);
        if K=NClasses-1 then
        begin
            I:=0;
            while I<=NClasses-2 do
            begin
                LM.W[Offs+I*(NVars+1)+NVars] := -V;
                Inc(I);
            end;
        end
        else
        begin
            I:=0;
            while I<=NClasses-2 do
            begin
                if I=K then
                begin
                    LM.W[Offs+I*(NVars+1)+NVars] := +V;
                end
                else
                begin
                    LM.W[Offs+I*(NVars+1)+NVars] := 0;
                end;
                Inc(I);
            end;
        end;
        Exit;
    end;
    
    //
    // General case.
    // Prepare task and network. Allocate space.
    //
    MLPCreateC0(NVars, NClasses, Network);
    MLPInitPreprocessor(Network, XY, NPoints);
    MLPProperties(Network, NIn, NOut, WCount);
    I:=0;
    while I<=WCount-1 do
    begin
        Network.Weights[I] := (2*RandomReal-1)/NVars;
        Inc(I);
    end;
    SetLength(G, WCount-1+1);
    SetLength(H, WCount-1+1, WCount-1+1);
    SetLength(WBase, WCount-1+1);
    SetLength(WDir, WCount-1+1);
    SetLength(WORK, WCount-1+1);
    
    //
    // First stage: optimize in gradient direction.
    //
    K:=0;
    while K<=WCount div 3+10 do
    begin
        
        //
        // Calculate gradient in starting point
        //
        MLPGradNBatch(Network, XY, NPoints, E, G);
        V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        E := E+Double(0.5)*Decay*V;
        APVAdd(@G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
        Rep.NGrad := Rep.NGrad+1;
        
        //
        // Setup optimization scheme
        //
        APVMoveNeg(@WDir[0], 0, WCount-1, @G[0], 0, WCount-1);
        V := APVDotProduct(@WDir[0], 0, WCount-1, @WDir[0], 0, WCount-1);
        WStep := Sqrt(V);
        V := 1/Sqrt(V);
        APVMul(@WDir[0], 0, WCount-1, V);
        MCStage := 0;
        MNLMCSRCH(WCount, Network.Weights, E, G, WDir, WStep, MCInfo, MCNFEV, WORK, MCState, MCStage);
        while MCStage<>0 do
        begin
            MLPGradNBatch(Network, XY, NPoints, E, G);
            V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            E := E+Double(0.5)*Decay*V;
            APVAdd(@G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
            Rep.NGrad := Rep.NGrad+1;
            MNLMCSRCH(WCount, Network.Weights, E, G, WDir, WStep, MCInfo, MCNFEV, WORK, MCState, MCStage);
        end;
        Inc(K);
    end;
    
    //
    // Second stage: use Hessian when we are close to the minimum
    //
    while True do
    begin
        
        //
        // Calculate and update E/G/H
        //
        MLPHessianNBatch(Network, XY, NPoints, E, G, H);
        V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        E := E+Double(0.5)*Decay*V;
        APVAdd(@G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
        K:=0;
        while K<=WCount-1 do
        begin
            H[K,K] := H[K,K]+Decay;
            Inc(K);
        end;
        Rep.NHess := Rep.NHess+1;
        
        //
        // Select step direction
        // NOTE: it is important to use lower-triangle Cholesky
        // factorization since it is much faster than higher-triangle version.
        //
        SPD := SPDMatrixCholesky(H, WCount, False);
        SPDMatrixCholeskySolve(H, WCount, False, G, SolverInfo, SolverRep, WDir);
        SPD := SolverInfo>0;
        if SPD then
        begin
            
            //
            // H is positive definite.
            // Step in Newton direction.
            //
            APVMul(@WDir[0], 0, WCount-1, -1);
            SPD := True;
        end
        else
        begin
            
            //
            // H is indefinite.
            // Step in gradient direction.
            //
            APVMoveNeg(@WDir[0], 0, WCount-1, @G[0], 0, WCount-1);
            SPD := False;
        end;
        
        //
        // Optimize in WDir direction
        //
        V := APVDotProduct(@WDir[0], 0, WCount-1, @WDir[0], 0, WCount-1);
        WStep := Sqrt(V);
        V := 1/Sqrt(V);
        APVMul(@WDir[0], 0, WCount-1, V);
        MCStage := 0;
        MNLMCSRCH(WCount, Network.Weights, E, G, WDir, WStep, MCInfo, MCNFEV, WORK, MCState, MCStage);
        while MCStage<>0 do
        begin
            MLPGradNBatch(Network, XY, NPoints, E, G);
            V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            E := E+Double(0.5)*Decay*V;
            APVAdd(@G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
            Rep.NGrad := Rep.NGrad+1;
            MNLMCSRCH(WCount, Network.Weights, E, G, WDir, WStep, MCInfo, MCNFEV, WORK, MCState, MCStage);
        end;
        if SPD and ((MCInfo=2) or (MCInfo=4) or (MCInfo=6)) then
        begin
            Break;
        end;
    end;
    
    //
    // Convert from NN format to MNL format
    //
    APVMove(@LM.W[0], Offs, Offs+WCount-1, @Network.Weights[0], 0, WCount-1);
    K:=0;
    while K<=NVars-1 do
    begin
        I:=0;
        while I<=NClasses-2 do
        begin
            S := Network.ColumnSigmas[K];
            if AP_FP_Eq(S,0) then
            begin
                S := 1;
            end;
            J := Offs+(NVars+1)*I;
            V := LM.W[J+K];
            LM.W[J+K] := V/S;
            LM.W[J+NVars] := LM.W[J+NVars]+V*Network.ColumnMeans[K]/S;
            Inc(I);
        end;
        Inc(K);
    end;
    K:=0;
    while K<=NClasses-2 do
    begin
        LM.W[Offs+(NVars+1)*K+NVars] := -LM.W[Offs+(NVars+1)*K+NVars];
        Inc(K);
    end;
end;


(*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   logit model, passed by non-constant reference
                (some fields of structure are used as temporaries
                when calculating model output).
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result, array[0..NClasses-1]
                Vector of posterior probabilities for classification task.
                Subroutine does not allocate memory for this vector, it is
                responsibility of a caller to allocate it. Array  must  be
                at least [0..NClasses-1].

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
procedure MNLProcess(var LM : LogitModel;
     const X : TReal1DArray;
     var Y : TReal1DArray);
var
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    Offs : AlglibInteger;
    I : AlglibInteger;
    I1 : AlglibInteger;
    S : Double;
begin
    Assert(AP_FP_Eq(LM.W[1],LogitVNum), 'MNLProcess: unexpected model version');
    NVars := Round(LM.W[2]);
    NClasses := Round(LM.W[3]);
    Offs := Round(LM.W[4]);
    MNLIExp(LM.W, X);
    S := 0;
    I1 := Offs+(NVars+1)*(NClasses-1);
    I:=I1;
    while I<=I1+NClasses-1 do
    begin
        S := S+LM.W[I];
        Inc(I);
    end;
    I:=0;
    while I<=NClasses-1 do
    begin
        Y[I] := LM.W[I1+I]/S;
        Inc(I);
    end;
end;


(*************************************************************************
Unpacks coefficients of logit model. Logit model have form:

    P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))
          S(i) = Exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]), when i<M-1
        S(M-1) = 1

INPUT PARAMETERS:
    LM          -   logit model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NClasses-2,0..NVars]
    NVars       -   number of independent variables
    NClasses    -   number of classes

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
procedure MNLUnpack(const LM : LogitModel;
     var A : TReal2DArray;
     var NVars : AlglibInteger;
     var NClasses : AlglibInteger);
var
    Offs : AlglibInteger;
    I : AlglibInteger;
begin
    Assert(AP_FP_Eq(LM.W[1],LogitVNum), 'MNLUnpack: unexpected model version');
    NVars := Round(LM.W[2]);
    NClasses := Round(LM.W[3]);
    Offs := Round(LM.W[4]);
    SetLength(A, NClasses-2+1, NVars+1);
    I:=0;
    while I<=NClasses-2 do
    begin
        APVMove(@A[I][0], 0, NVars, @LM.W[0], Offs+I*(NVars+1), Offs+I*(NVars+1)+NVars);
        Inc(I);
    end;
end;


(*************************************************************************
"Packs" coefficients and creates logit model in ALGLIB format (MNLUnpack
reversed).

INPUT PARAMETERS:
    A           -   model (see MNLUnpack)
    NVars       -   number of independent variables
    NClasses    -   number of classes

OUTPUT PARAMETERS:
    LM          -   logit model.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
procedure MNLPack(const A : TReal2DArray;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var LM : LogitModel);
var
    Offs : AlglibInteger;
    I : AlglibInteger;
    WDim : AlglibInteger;
    SSize : AlglibInteger;
begin
    WDim := (NVars+1)*(NClasses-1);
    Offs := 5;
    SSize := 5+(NVars+1)*(NClasses-1)+NClasses;
    SetLength(LM.W, SSize-1+1);
    LM.W[0] := SSize;
    LM.W[1] := LogitVNum;
    LM.W[2] := NVars;
    LM.W[3] := NClasses;
    LM.W[4] := Offs;
    I:=0;
    while I<=NClasses-2 do
    begin
        APVMove(@LM.W[0], Offs+I*(NVars+1), Offs+I*(NVars+1)+NVars, @A[I][0], 0, NVars);
        Inc(I);
    end;
end;


(*************************************************************************
Copying of LogitModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MNLCopy(const LM1 : LogitModel; var LM2 : LogitModel);
var
    K : AlglibInteger;
begin
    K := Round(LM1.W[0]);
    SetLength(LM2.W, K-1+1);
    APVMove(@LM2.W[0], 0, K-1, @LM1.W[0], 0, K-1);
end;


(*************************************************************************
Serialization of LogitModel strucure

INPUT PARAMETERS:
    LM      -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores model,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MNLSerialize(const LM : LogitModel;
     var RA : TReal1DArray;
     var RLen : AlglibInteger);
begin
    RLen := Round(LM.W[0])+1;
    SetLength(RA, RLen-1+1);
    RA[0] := LogitVNum;
    APVMove(@RA[0], 1, RLen-1, @LM.W[0], 0, RLen-2);
end;


(*************************************************************************
Unserialization of LogitModel strucure

INPUT PARAMETERS:
    RA      -   real array which stores model

OUTPUT PARAMETERS:
    LM      -   restored model

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MNLUnserialize(const RA : TReal1DArray; var LM : LogitModel);
begin
    Assert(Round(RA[0])=LogitVNum, 'MNLUnserialize: incorrect array!');
    SetLength(LM.W, Round(RA[1])-1+1);
    APVMove(@LM.W[0], 0, Round(RA[1])-1, @RA[0], 1, Round(RA[1]));
end;


(*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*ln(2)).

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
function MNLAvgCE(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    I : AlglibInteger;
    WorkX : TReal1DArray;
    WorkY : TReal1DArray;
begin
    Assert(AP_FP_Eq(LM.W[1],LogitVNum), 'MNLClsError: unexpected model version');
    NVars := Round(LM.W[2]);
    NClasses := Round(LM.W[3]);
    SetLength(WorkX, NVars-1+1);
    SetLength(WorkY, NClasses-1+1);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        Assert((Round(XY[I,NVars])>=0) and (Round(XY[I,NVars])<NClasses), 'MNLAvgCE: incorrect class number!');
        
        //
        // Process
        //
        APVMove(@WorkX[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        MNLProcess(LM, WorkX, WorkY);
        if AP_FP_Greater(WorkY[Round(XY[I,NVars])],0) then
        begin
            Result := Result-Ln(WorkY[Round(XY[I,NVars])]);
        end
        else
        begin
            Result := Result-Ln(MinRealNumber);
        end;
        Inc(I);
    end;
    Result := Result/(NPoints*Ln(2));
end;


(*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
function MNLRelClsError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
begin
    Result := AP_Double(MNLClsError(LM, XY, NPoints))/NPoints;
end;


(*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
function MNLRMSError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    Assert(Round(LM.W[1])=LogitVNum, 'MNLRMSError: Incorrect MNL version!');
    MNLAllErrors(LM, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := RMS;
end;


(*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
function MNLAvgError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    Assert(Round(LM.W[1])=LogitVNum, 'MNLRMSError: Incorrect MNL version!');
    MNLAllErrors(LM, XY, NPoints, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := Avg;
end;


(*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
function MNLAvgRelError(var LM : LogitModel;
     const XY : TReal2DArray;
     SSize : AlglibInteger):Double;
var
    RelCls : Double;
    AvgCE : Double;
    RMS : Double;
    Avg : Double;
    AvgRel : Double;
begin
    Assert(Round(LM.W[1])=LogitVNum, 'MNLRMSError: Incorrect MNL version!');
    MNLAllErrors(LM, XY, SSize, RelCls, AvgCE, RMS, Avg, AvgRel);
    Result := AvgRel;
end;


(*************************************************************************
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************)
function MNLClsError(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger):AlglibInteger;
var
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    WorkX : TReal1DArray;
    WorkY : TReal1DArray;
    NMAX : AlglibInteger;
begin
    Assert(AP_FP_Eq(LM.W[1],LogitVNum), 'MNLClsError: unexpected model version');
    NVars := Round(LM.W[2]);
    NClasses := Round(LM.W[3]);
    SetLength(WorkX, NVars-1+1);
    SetLength(WorkY, NClasses-1+1);
    Result := 0;
    I:=0;
    while I<=NPoints-1 do
    begin
        
        //
        // Process
        //
        APVMove(@WorkX[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        MNLProcess(LM, WorkX, WorkY);
        
        //
        // Logit version of the answer
        //
        NMAX := 0;
        J:=0;
        while J<=NClasses-1 do
        begin
            if AP_FP_Greater(WorkY[J],WorkY[NMAX]) then
            begin
                NMAX := J;
            end;
            Inc(J);
        end;
        
        //
        // compare
        //
        if NMAX<>Round(XY[I,NVars]) then
        begin
            Result := Result+1;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal subroutine. Places exponents of the anti-overflow shifted
internal linear outputs into the service part of the W array.
*************************************************************************)
procedure MNLIExp(var W : TReal1DArray; const X : TReal1DArray);
var
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    Offs : AlglibInteger;
    I : AlglibInteger;
    I1 : AlglibInteger;
    V : Double;
    MX : Double;
begin
    Assert(AP_FP_Eq(W[1],LogitVNum), 'LOGIT: unexpected model version');
    NVars := Round(W[2]);
    NClasses := Round(W[3]);
    Offs := Round(W[4]);
    I1 := Offs+(NVars+1)*(NClasses-1);
    I:=0;
    while I<=NClasses-2 do
    begin
        V := APVDotProduct(@W[0], Offs+I*(NVars+1), Offs+I*(NVars+1)+NVars-1, @X[0], 0, NVars-1);
        W[I1+I] := V+W[Offs+I*(NVars+1)+NVars];
        Inc(I);
    end;
    W[I1+NClasses-1] := 0;
    MX := 0;
    I:=I1;
    while I<=I1+NClasses-1 do
    begin
        MX := Max(MX, W[I]);
        Inc(I);
    end;
    I:=I1;
    while I<=I1+NClasses-1 do
    begin
        W[I] := Exp(W[I]-MX);
        Inc(I);
    end;
end;


(*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************)
procedure MNLAllErrors(var LM : LogitModel;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     var RelCls : Double;
     var AvgCE : Double;
     var RMS : Double;
     var Avg : Double;
     var AvgRel : Double);
var
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    I : AlglibInteger;
    Buf : TReal1DArray;
    WorkX : TReal1DArray;
    Y : TReal1DArray;
    DY : TReal1DArray;
begin
    Assert(Round(LM.W[1])=LogitVNum, 'MNL unit: Incorrect MNL version!');
    NVars := Round(LM.W[2]);
    NClasses := Round(LM.W[3]);
    SetLength(WorkX, NVars-1+1);
    SetLength(Y, NClasses-1+1);
    SetLength(DY, 0+1);
    DSErrAllocate(NClasses, Buf);
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@WorkX[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        MNLProcess(LM, WorkX, Y);
        DY[0] := XY[I,NVars];
        DSErrAccumulate(Buf, Y, DY);
        Inc(I);
    end;
    DSErrFinish(Buf);
    RelCls := Buf[0];
    AvgCE := Buf[1];
    Rms := Buf[2];
    Avg := Buf[3];
    AvgRel := Buf[4];
end;


(*************************************************************************
THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
DECREASE CONDITION AND A CURVATURE CONDITION.

AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION

    F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).

THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
DECREASE CONDITION

    F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

AND THE CURVATURE CONDITION

    ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.

PARAMETERS DESCRIPRION

N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.

X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.

F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
IT CONTAINS THE VALUE OF F AT X + STP*S.

G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.

S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.

STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.

FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
SATISFIED.

XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.

STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
UPPER BOUNDS FOR THE STEP.

MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.

INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    INFO = 0  IMPROPER INPUT PARAMETERS.

    INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
              DIRECTIONAL DERIVATIVE CONDITION HOLD.

    INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
              IS AT MOST XTOL.

    INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

    INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

    INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

    INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
              THERE MAY NOT BE A STEP WHICH SATISFIES THE
              SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
              TOLERANCES MAY BE TOO SMALL.

NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

WA IS A WORK ARRAY OF LENGTH N.

ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE
*************************************************************************)
procedure MNLMCSRCH(const N : AlglibInteger;
     var X : TReal1DArray;
     var F : Double;
     var G : TReal1DArray;
     const S : TReal1DArray;
     var STP : Double;
     var INFO : AlglibInteger;
     var NFEV : AlglibInteger;
     var WA : TReal1DArray;
     var State : LOGITMCState;
     var Stage : AlglibInteger);
var
    V : Double;
    P5 : Double;
    P66 : Double;
    ZERO : Double;
begin
    
    //
    // init
    //
    P5 := Double(0.5);
    P66 := Double(0.66);
    State.XTRAPF := Double(4.0);
    ZERO := 0;
    
    //
    // Main cycle
    //
    while True do
    begin
        if Stage=0 then
        begin
            
            //
            // NEXT
            //
            Stage := 2;
            Continue;
        end;
        if Stage=2 then
        begin
            State.INFOC := 1;
            INFO := 0;
            
            //
            //     CHECK THE INPUT PARAMETERS FOR ERRORS.
            //
            if (N<=0) or AP_FP_Less_Eq(STP,0) or AP_FP_Less(FTOL,0) or AP_FP_Less(GTOL,ZERO) or AP_FP_Less(XTOL,ZERO) or AP_FP_Less(STPMIN,ZERO) or AP_FP_Less(STPMAX,STPMIN) or (MAXFEV<=0) then
            begin
                Stage := 0;
                Exit;
            end;
            
            //
            //     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
            //     AND CHECK THAT S IS A DESCENT DIRECTION.
            //
            V := APVDotProduct(@G[0], 0, N-1, @S[0], 0, N-1);
            State.DGINIT := V;
            if AP_FP_Greater_Eq(State.DGINIT,0) then
            begin
                Stage := 0;
                Exit;
            end;
            
            //
            //     INITIALIZE LOCAL VARIABLES.
            //
            State.BRACKT := False;
            State.STAGE1 := True;
            NFEV := 0;
            State.FINIT := F;
            State.DGTEST := FTOL*State.DGINIT;
            State.WIDTH := STPMAX-STPMIN;
            State.WIDTH1 := State.WIDTH/P5;
            APVMove(@WA[0], 0, N-1, @X[0], 0, N-1);
            
            //
            //     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
            //     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
            //     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
            //     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
            //     THE INTERVAL OF UNCERTAINTY.
            //     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
            //     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
            //
            State.STX := 0;
            State.FX := State.FINIT;
            State.DGX := State.DGINIT;
            State.STY := 0;
            State.FY := State.FINIT;
            State.DGY := State.DGINIT;
            
            //
            // NEXT
            //
            Stage := 3;
            Continue;
        end;
        if Stage=3 then
        begin
            
            //
            //     START OF ITERATION.
            //
            //     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
            //     TO THE PRESENT INTERVAL OF UNCERTAINTY.
            //
            if State.BRACKT then
            begin
                if AP_FP_Less(State.STX,State.STY) then
                begin
                    State.STMIN := State.STX;
                    State.STMAX := State.STY;
                end
                else
                begin
                    State.STMIN := State.STY;
                    State.STMAX := State.STX;
                end;
            end
            else
            begin
                State.STMIN := State.STX;
                State.STMAX := STP+State.XTRAPF*(STP-State.STX);
            end;
            
            //
            //        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
            //
            if AP_FP_Greater(STP,STPMAX) then
            begin
                STP := STPMAX;
            end;
            if AP_FP_Less(STP,STPMIN) then
            begin
                STP := STPMIN;
            end;
            
            //
            //        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
            //        STP BE THE LOWEST POINT OBTAINED SO FAR.
            //
            if State.BRACKT and (AP_FP_Less_Eq(STP,State.STMIN) or AP_FP_Greater_Eq(STP,State.STMAX)) or (NFEV>=MAXFEV-1) or (State.INFOC=0) or State.BRACKT and AP_FP_Less_Eq(State.STMAX-State.STMIN,XTOL*State.STMAX) then
            begin
                STP := State.STX;
            end;
            
            //
            //        EVALUATE THE FUNCTION AND GRADIENT AT STP
            //        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
            //
            APVMove(@X[0], 0, N-1, @WA[0], 0, N-1);
            APVAdd(@X[0], 0, N-1, @S[0], 0, N-1, STP);
            
            //
            // NEXT
            //
            Stage := 4;
            Exit;
        end;
        if Stage=4 then
        begin
            INFO := 0;
            NFEV := NFEV+1;
            V := APVDotProduct(@G[0], 0, N-1, @S[0], 0, N-1);
            State.DG := V;
            State.FTEST1 := State.FINIT+STP*State.DGTEST;
            
            //
            //        TEST FOR CONVERGENCE.
            //
            if State.BRACKT and (AP_FP_Less_Eq(STP,State.STMIN) or AP_FP_Greater_Eq(STP,State.STMAX)) or (State.INFOC=0) then
            begin
                INFO := 6;
            end;
            if AP_FP_Eq(STP,STPMAX) and AP_FP_Less_Eq(F,State.FTEST1) and AP_FP_Less_Eq(State.DG,State.DGTEST) then
            begin
                INFO := 5;
            end;
            if AP_FP_Eq(STP,STPMIN) and (AP_FP_Greater(F,State.FTEST1) or AP_FP_Greater_Eq(State.DG,State.DGTEST)) then
            begin
                INFO := 4;
            end;
            if NFEV>=MAXFEV then
            begin
                INFO := 3;
            end;
            if State.BRACKT and AP_FP_Less_Eq(State.STMAX-State.STMIN,XTOL*State.STMAX) then
            begin
                INFO := 2;
            end;
            if AP_FP_Less_Eq(F,State.FTEST1) and AP_FP_Less_Eq(AbsReal(State.DG),-GTOL*State.DGINIT) then
            begin
                INFO := 1;
            end;
            
            //
            //        CHECK FOR TERMINATION.
            //
            if INFO<>0 then
            begin
                Stage := 0;
                Exit;
            end;
            
            //
            //        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
            //        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
            //
            if State.STAGE1 and AP_FP_Less_Eq(F,State.FTEST1) and AP_FP_Greater_Eq(State.DG,Min(FTOL, GTOL)*State.DGINIT) then
            begin
                State.STAGE1 := False;
            end;
            
            //
            //        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
            //        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
            //        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
            //        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
            //        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
            //
            if State.STAGE1 and AP_FP_Less_Eq(F,State.FX) and AP_FP_Greater(F,State.FTEST1) then
            begin
                
                //
                //           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
                //
                State.FM := F-STP*State.DGTEST;
                State.FXM := State.FX-State.STX*State.DGTEST;
                State.FYM := State.FY-State.STY*State.DGTEST;
                State.DGM := State.DG-State.DGTEST;
                State.DGXM := State.DGX-State.DGTEST;
                State.DGYM := State.DGY-State.DGTEST;
                
                //
                //           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                //           AND TO COMPUTE THE NEW STEP.
                //
                MNLMCSTEP(State.STX, State.FXM, State.DGXM, State.STY, State.FYM, State.DGYM, STP, State.FM, State.DGM, State.BRACKT, State.STMIN, State.STMAX, State.INFOC);
                
                //
                //           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
                //
                State.FX := State.FXM+State.STX*State.DGTEST;
                State.FY := State.FYM+State.STY*State.DGTEST;
                State.DGX := State.DGXM+State.DGTEST;
                State.DGY := State.DGYM+State.DGTEST;
            end
            else
            begin
                
                //
                //           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                //           AND TO COMPUTE THE NEW STEP.
                //
                MNLMCSTEP(State.STX, State.FX, State.DGX, State.STY, State.FY, State.DGY, STP, F, State.DG, State.BRACKT, State.STMIN, State.STMAX, State.INFOC);
            end;
            
            //
            //        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
            //        INTERVAL OF UNCERTAINTY.
            //
            if State.BRACKT then
            begin
                if AP_FP_Greater_Eq(AbsReal(State.STY-State.STX),P66*State.WIDTH1) then
                begin
                    STP := State.STX+P5*(State.STY-State.STX);
                end;
                State.WIDTH1 := State.WIDTH;
                State.WIDTH := AbsReal(State.STY-State.STX);
            end;
            
            //
            //  NEXT.
            //
            Stage := 3;
            Continue;
        end;
    end;
end;


procedure MNLMCSTEP(var STX : Double;
     var FX : Double;
     var DX : Double;
     var STY : Double;
     var FY : Double;
     var DY : Double;
     var STP : Double;
     const FP : Double;
     const DP : Double;
     var BRACKT : Boolean;
     const STMIN : Double;
     const STMAX : Double;
     var INFO : AlglibInteger);
var
    BOUND : Boolean;
    GAMMA : Double;
    P : Double;
    Q : Double;
    R : Double;
    S : Double;
    SGND : Double;
    STPC : Double;
    STPF : Double;
    STPQ : Double;
    THETA : Double;
begin
    INFO := 0;
    
    //
    //     CHECK THE INPUT PARAMETERS FOR ERRORS.
    //
    if BRACKT and (AP_FP_Less_Eq(STP,Min(STX, STY)) or AP_FP_Greater_Eq(STP,Max(STX, STY))) or AP_FP_Greater_Eq(DX*(STP-STX),0) or AP_FP_Less(STMAX,STMIN) then
    begin
        Exit;
    end;
    
    //
    //     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
    //
    SGND := DP*(DX/AbsReal(DX));
    
    //
    //     FIRST CASE. A HIGHER FUNCTION VALUE.
    //     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
    //     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
    //     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
    //
    if AP_FP_Greater(FP,FX) then
    begin
        INFO := 1;
        BOUND := True;
        THETA := 3*(FX-FP)/(STP-STX)+DX+DP;
        S := Max(AbsReal(THETA), Max(AbsReal(DX), AbsReal(DP)));
        GAMMA := S*Sqrt(AP_Sqr(THETA/S)-DX/S*(DP/S));
        if AP_FP_Less(STP,STX) then
        begin
            GAMMA := -GAMMA;
        end;
        P := GAMMA-DX+THETA;
        Q := GAMMA-DX+GAMMA+DP;
        R := P/Q;
        STPC := STX+R*(STP-STX);
        STPQ := STX+DX/((FX-FP)/(STP-STX)+DX)/2*(STP-STX);
        if AP_FP_Less(AbsReal(STPC-STX),AbsReal(STPQ-STX)) then
        begin
            STPF := STPC;
        end
        else
        begin
            STPF := STPC+(STPQ-STPC)/2;
        end;
        BRACKT := True;
    end
    else
    begin
        if AP_FP_Less(SGND,0) then
        begin
            
            //
            //     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
            //     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
            //     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
            //     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
            //
            INFO := 2;
            BOUND := False;
            THETA := 3*(FX-FP)/(STP-STX)+DX+DP;
            S := Max(AbsReal(THETA), Max(AbsReal(DX), AbsReal(DP)));
            GAMMA := S*SQRT(AP_Sqr(THETA/S)-DX/S*(DP/S));
            if AP_FP_Greater(STP,STX) then
            begin
                GAMMA := -GAMMA;
            end;
            P := GAMMA-DP+THETA;
            Q := GAMMA-DP+GAMMA+DX;
            R := P/Q;
            STPC := STP+R*(STX-STP);
            STPQ := STP+DP/(DP-DX)*(STX-STP);
            if AP_FP_Greater(AbsReal(STPC-STP),AbsReal(STPQ-STP)) then
            begin
                STPF := STPC;
            end
            else
            begin
                STPF := STPQ;
            end;
            BRACKT := True;
        end
        else
        begin
            if AP_FP_Less(AbsReal(DP),AbsReal(DX)) then
            begin
                
                //
                //     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
                //     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
                //     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
                //     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
                //     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
                //     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
                //     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
                //
                INFO := 3;
                BOUND := True;
                THETA := 3*(FX-FP)/(STP-STX)+DX+DP;
                S := Max(AbsReal(THETA), Max(AbsReal(DX), AbsReal(DP)));
                
                //
                //        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
                //        TO INFINITY IN THE DIRECTION OF THE STEP.
                //
                GAMMA := S*SQRT(Max(0, AP_Sqr(THETA/S)-DX/S*(DP/S)));
                if AP_FP_Greater(STP,STX) then
                begin
                    GAMMA := -GAMMA;
                end;
                P := GAMMA-DP+THETA;
                Q := GAMMA+(DX-DP)+GAMMA;
                R := P/Q;
                if AP_FP_Less(R,0) and AP_FP_Neq(GAMMA,0) then
                begin
                    STPC := STP+R*(STX-STP);
                end
                else
                begin
                    if AP_FP_Greater(STP,STX) then
                    begin
                        STPC := STMAX;
                    end
                    else
                    begin
                        STPC := STMIN;
                    end;
                end;
                STPQ := STP+DP/(DP-DX)*(STX-STP);
                if BRACKT then
                begin
                    if AP_FP_Less(AbsReal(STP-STPC),AbsReal(STP-STPQ)) then
                    begin
                        STPF := STPC;
                    end
                    else
                    begin
                        STPF := STPQ;
                    end;
                end
                else
                begin
                    if AP_FP_Greater(AbsReal(STP-STPC),AbsReal(STP-STPQ)) then
                    begin
                        STPF := STPC;
                    end
                    else
                    begin
                        STPF := STPQ;
                    end;
                end;
            end
            else
            begin
                
                //
                //     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
                //     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
                //     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
                //
                INFO := 4;
                BOUND := False;
                if BRACKT then
                begin
                    THETA := 3*(FP-FY)/(STY-STP)+DY+DP;
                    S := Max(ABSReal(THETA), Max(ABSReal(DY), ABSReal(DP)));
                    GAMMA := S*SQRT(AP_Sqr(THETA/S)-DY/S*(DP/S));
                    if AP_FP_Greater(STP,STY) then
                    begin
                        GAMMA := -GAMMA;
                    end;
                    P := GAMMA-DP+THETA;
                    Q := GAMMA-DP+GAMMA+DY;
                    R := P/Q;
                    STPC := STP+R*(STY-STP);
                    STPF := STPC;
                end
                else
                begin
                    if AP_FP_Greater(STP,STX) then
                    begin
                        STPF := STMAX;
                    end
                    else
                    begin
                        STPF := STMIN;
                    end;
                end;
            end;
        end;
    end;
    
    //
    //     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
    //     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
    //
    if AP_FP_Greater(FP,FX) then
    begin
        STY := STP;
        FY := FP;
        DY := DP;
    end
    else
    begin
        if AP_FP_Less(SGND,Double(0.0)) then
        begin
            STY := STX;
            FY := FX;
            DY := DX;
        end;
        STX := STP;
        FX := FP;
        DX := DP;
    end;
    
    //
    //     COMPUTE THE NEW STEP AND SAFEGUARD IT.
    //
    STPF := Min(STMAX, STPF);
    STPF := Max(STMIN, STPF);
    STP := STPF;
    if BRACKT and BOUND then
    begin
        if AP_FP_Greater(STY,STX) then
        begin
            STP := Min(STX+Double(0.66)*(STY-STX), STP);
        end
        else
        begin
            STP := Max(STX+Double(0.66)*(STY-STX), STP);
        end;
    end;
end;


end.