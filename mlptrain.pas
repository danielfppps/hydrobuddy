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
unit mlptrain;
interface
uses Math, Sysutils, Ap, mlpbase, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, linmin, minlbfgs, hblas, sblas, ortfac, blas, rotations, bdsvd, svd, xblas, densesolver;

type
(*************************************************************************
Training report:
    * NGrad     - number of gradient calculations
    * NHess     - number of Hessian calculations
    * NCholesky - number of Cholesky decompositions
*************************************************************************)
MLPReport = record
    NGrad : AlglibInteger;
    NHess : AlglibInteger;
    NCholesky : AlglibInteger;
end;


(*************************************************************************
Cross-validation estimates of generalization error
*************************************************************************)
MLPCVReport = record
    RelCLSError : Double;
    AvgCE : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
end;



procedure MLPTrainLM(var Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
procedure MLPTrainLBFGS(var Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
procedure MLPTrainES(var Network : MultiLayerPerceptron;
     const TrnXY : TReal2DArray;
     TrnSize : AlglibInteger;
     const ValXY : TReal2DArray;
     ValSize : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
procedure MLPKFoldCVLBFGS(const Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     FoldsCount : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var CVRep : MLPCVReport);
procedure MLPKFoldCVLM(const Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     FoldsCount : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var CVRep : MLPCVReport);

implementation

const
    MinDecay = Double(0.001);

procedure MLPKFoldCVGeneral(const N : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     FoldsCount : AlglibInteger;
     LMAlgorithm : Boolean;
     WStep : Double;
     MaxIts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var CVRep : MLPCVReport);forward;
procedure MLPKFoldSplit(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NClasses : AlglibInteger;
     FoldsCount : AlglibInteger;
     StratifiedSplits : Boolean;
     var Folds : TInteger1DArray);forward;


(*************************************************************************
Neural network training  using  modified  Levenberg-Marquardt  with  exact
Hessian calculation and regularization. Subroutine trains  neural  network
with restarts from random positions. Algorithm is well  suited  for  small
and medium scale problems (hundreds of weights).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -9, if internal matrix inverse subroutine failed
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPTrainLM(var Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
var
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    WCount : AlglibInteger;
    LMFTol : Double;
    LMStepTol : Double;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    E : Double;
    ENew : Double;
    XNorm2 : Double;
    StepNorm : Double;
    G : TReal1DArray;
    D : TReal1DArray;
    H : TReal2DArray;
    HMod : TReal2DArray;
    Z : TReal2DArray;
    SPD : Boolean;
    Nu : Double;
    Lambda : Double;
    LambdaUp : Double;
    LambdaDown : Double;
    InternalRep : MinLBFGSReport;
    State : MinLBFGSState;
    X : TReal1DArray;
    Y : TReal1DArray;
    WBase : TReal1DArray;
    WDir : TReal1DArray;
    WT : TReal1DArray;
    WX : TReal1DArray;
    Pass : AlglibInteger;
    WBest : TReal1DArray;
    EBest : Double;
    InvInfo : AlglibInteger;
    InvRep : MatInvReport;
    SolverInfo : AlglibInteger;
    SolverRep : DenseSolverReport;
begin
    MLPProperties(Network, NIn, NOut, WCount);
    LambdaUp := 10;
    LambdaDown := Double(0.3);
    LMFTol := Double(0.001);
    LMStepTol := Double(0.001);
    
    //
    // Test for inputs
    //
    if (NPoints<=0) or (Restarts<1) then
    begin
        Info := -1;
        Exit;
    end;
    if MLPIsSoftmax(Network) then
    begin
        I:=0;
        while I<=NPoints-1 do
        begin
            if (Round(XY[I,NIn])<0) or (Round(XY[I,NIn])>=NOut) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
    end;
    Decay := Max(Decay, MinDecay);
    Info := 2;
    
    //
    // Initialize data
    //
    Rep.NGrad := 0;
    Rep.NHess := 0;
    Rep.NCholesky := 0;
    
    //
    // General case.
    // Prepare task and network. Allocate space.
    //
    MLPInitPreprocessor(Network, XY, NPoints);
    SetLength(G, WCount-1+1);
    SetLength(H, WCount-1+1, WCount-1+1);
    SetLength(HMod, WCount-1+1, WCount-1+1);
    SetLength(WBase, WCount-1+1);
    SetLength(WDir, WCount-1+1);
    SetLength(WBest, WCount-1+1);
    SetLength(WT, WCount-1+1);
    SetLength(WX, WCount-1+1);
    EBest := MaxRealNumber;
    
    //
    // Multiple passes
    //
    Pass:=1;
    while Pass<=Restarts do
    begin
        
        //
        // Initialize weights
        //
        MLPRandomize(Network);
        
        //
        // First stage of the hybrid algorithm: LBFGS
        //
        APVMove(@WBase[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        MinLBFGSCreate(WCount, Min(WCount, 5), WBase, State);
        MinLBFGSSetCond(State, 0, 0, 0, Max(25, WCount));
        while MinLBFGSIteration(State) do
        begin
            
            //
            // gradient
            //
            APVMove(@Network.Weights[0], 0, WCount-1, @State.X[0], 0, WCount-1);
            MLPGradBatch(Network, XY, NPoints, State.F, State.G);
            
            //
            // weight decay
            //
            V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            State.F := State.F+Double(0.5)*Decay*V;
            APVAdd(@State.G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
            
            //
            // next iteration
            //
            Rep.NGrad := Rep.NGrad+1;
        end;
        MinLBFGSResults(State, WBase, InternalRep);
        APVMove(@Network.Weights[0], 0, WCount-1, @WBase[0], 0, WCount-1);
        
        //
        // Second stage of the hybrid algorithm: LM
        //
        // Initialize H with identity matrix,
        // G with gradient,
        // E with regularized error.
        //
        MLPHessianBatch(Network, XY, NPoints, E, G, H);
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
        Lambda := Double(0.001);
        Nu := 2;
        while True do
        begin
            
            //
            // 1. HMod = H+lambda*I
            // 2. Try to solve (H+Lambda*I)*dx = -g.
            //    Increase lambda if left part is not positive definite.
            //
            I:=0;
            while I<=WCount-1 do
            begin
                APVMove(@HMod[I][0], 0, WCount-1, @H[I][0], 0, WCount-1);
                HMod[I,I] := HMod[I,I]+Lambda;
                Inc(I);
            end;
            SPD := SPDMatrixCholesky(HMod, WCount, True);
            Rep.NCholesky := Rep.NCholesky+1;
            if  not SPD then
            begin
                Lambda := Lambda*LambdaUp*Nu;
                Nu := Nu*2;
                Continue;
            end;
            SPDMatrixCholeskySolve(HMod, WCount, True, G, SolverInfo, SolverRep, WDir);
            if SolverInfo<0 then
            begin
                Lambda := Lambda*LambdaUp*Nu;
                Nu := Nu*2;
                Continue;
            end;
            APVMul(@WDir[0], 0, WCount-1, -1);
            
            //
            // Lambda found.
            // 1. Save old w in WBase
            // 1. Test some stopping criterions
            // 2. If error(w+wdir)>error(w), increase lambda
            //
            APVAdd(@Network.Weights[0], 0, WCount-1, @WDir[0], 0, WCount-1);
            XNorm2 := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            StepNorm := APVDotProduct(@WDir[0], 0, WCount-1, @WDir[0], 0, WCount-1);
            StepNorm := Sqrt(StepNorm);
            ENew := MLPError(Network, XY, NPoints)+Double(0.5)*Decay*XNorm2;
            if AP_FP_Less(StepNorm,LMStepTol*(1+Sqrt(XNorm2))) then
            begin
                Break;
            end;
            if AP_FP_Greater(ENew,E) then
            begin
                Lambda := Lambda*LambdaUp*Nu;
                Nu := Nu*2;
                Continue;
            end;
            
            //
            // Optimize using inv(cholesky(H)) as preconditioner
            //
            RMatrixTRInverse(HMod, WCount, True, False, InvInfo, InvRep);
            if InvInfo<=0 then
            begin
                
                //
                // if matrix can't be inverted then exit with errors
                // TODO: make WCount steps in direction suggested by HMod
                //
                Info := -9;
                Exit;
            end;
            APVMove(@WBase[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            I:=0;
            while I<=WCount-1 do
            begin
                WT[I] := 0;
                Inc(I);
            end;
            MinLBFGSCreateX(WCount, WCount, WT, 1, State);
            MinLBFGSSetCond(State, 0, 0, 0, 5);
            while MinLBFGSIteration(State) do
            begin
                
                //
                // gradient
                //
                I:=0;
                while I<=WCount-1 do
                begin
                    V := APVDotProduct(@State.X[0], I, WCount-1, @HMod[I][0], I, WCount-1);
                    Network.Weights[I] := WBase[I]+V;
                    Inc(I);
                end;
                MLPGradBatch(Network, XY, NPoints, State.F, G);
                I:=0;
                while I<=WCount-1 do
                begin
                    State.G[I] := 0;
                    Inc(I);
                end;
                I:=0;
                while I<=WCount-1 do
                begin
                    V := G[I];
                    APVAdd(@State.G[0], I, WCount-1, @HMod[I][0], I, WCount-1, V);
                    Inc(I);
                end;
                
                //
                // weight decay
                // grad(x'*x) = A'*(x0+A*t)
                //
                V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
                State.F := State.F+Double(0.5)*Decay*V;
                I:=0;
                while I<=WCount-1 do
                begin
                    V := Decay*Network.Weights[I];
                    APVAdd(@State.G[0], I, WCount-1, @HMod[I][0], I, WCount-1, V);
                    Inc(I);
                end;
                
                //
                // next iteration
                //
                Rep.NGrad := Rep.NGrad+1;
            end;
            MinLBFGSResults(State, WT, InternalRep);
            
            //
            // Accept new position.
            // Calculate Hessian
            //
            I:=0;
            while I<=WCount-1 do
            begin
                V := APVDotProduct(@WT[0], I, WCount-1, @HMod[I][0], I, WCount-1);
                Network.Weights[I] := WBase[I]+V;
                Inc(I);
            end;
            MLPHessianBatch(Network, XY, NPoints, E, G, H);
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
            // Update lambda
            //
            Lambda := Lambda*LambdaDown;
            Nu := 2;
        end;
        
        //
        // update WBest
        //
        V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        E := Double(0.5)*Decay*V+MLPError(Network, XY, NPoints);
        if AP_FP_Less(E,EBest) then
        begin
            EBest := E;
            APVMove(@WBest[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        end;
        Inc(Pass);
    end;
    
    //
    // copy WBest to output
    //
    APVMove(@Network.Weights[0], 0, WCount-1, @WBest[0], 0, WCount-1);
end;


(*************************************************************************
Neural  network  training  using  L-BFGS  algorithm  with  regularization.
Subroutine  trains  neural  network  with  restarts from random positions.
Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
requirements and step complexity are linear by weights number).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    iterations (NOT gradient  calculations).  Zero  MaxIts
                    means stopping when step is sufficiently small.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************)
procedure MLPTrainLBFGS(var Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
var
    I : AlglibInteger;
    Pass : AlglibInteger;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    WCount : AlglibInteger;
    W : TReal1DArray;
    WBest : TReal1DArray;
    E : Double;
    V : Double;
    EBest : Double;
    InternalRep : MinLBFGSReport;
    State : MinLBFGSState;
begin
    
    //
    // Test inputs, parse flags, read network geometry
    //
    if AP_FP_Eq(WStep,0) and (MaxIts=0) then
    begin
        Info := -8;
        Exit;
    end;
    if (NPoints<=0) or (Restarts<1) or AP_FP_Less(WStep,0) or (MaxIts<0) then
    begin
        Info := -1;
        Exit;
    end;
    MLPProperties(Network, NIn, NOut, WCount);
    if MLPIsSoftmax(Network) then
    begin
        I:=0;
        while I<=NPoints-1 do
        begin
            if (Round(XY[I,NIn])<0) or (Round(XY[I,NIn])>=NOut) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
    end;
    Decay := Max(Decay, MinDecay);
    Info := 2;
    
    //
    // Prepare
    //
    MLPInitPreprocessor(Network, XY, NPoints);
    SetLength(W, WCount-1+1);
    SetLength(WBest, WCount-1+1);
    EBest := MaxRealNumber;
    
    //
    // Multiple starts
    //
    Rep.NCholesky := 0;
    Rep.NHess := 0;
    Rep.NGrad := 0;
    Pass:=1;
    while Pass<=Restarts do
    begin
        
        //
        // Process
        //
        MLPRandomize(Network);
        APVMove(@W[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        MinLBFGSCreate(WCount, Min(WCount, 10), W, State);
        MinLBFGSSetCond(State, Double(0.0), Double(0.0), WStep, MaxIts);
        while MinLBFGSIteration(State) do
        begin
            APVMove(@Network.Weights[0], 0, WCount-1, @State.X[0], 0, WCount-1);
            MLPGradNBatch(Network, XY, NPoints, State.F, State.G);
            V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            State.F := State.F+Double(0.5)*Decay*V;
            APVAdd(@State.G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
            Rep.NGrad := Rep.NGrad+1;
        end;
        MinLBFGSResults(State, W, InternalRep);
        APVMove(@Network.Weights[0], 0, WCount-1, @W[0], 0, WCount-1);
        
        //
        // Compare with best
        //
        V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        E := MLPErrorN(Network, XY, NPoints)+Double(0.5)*Decay*V;
        if AP_FP_Less(E,EBest) then
        begin
            APVMove(@WBest[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            EBest := E;
        end;
        Inc(Pass);
    end;
    
    //
    // The best network
    //
    APVMove(@Network.Weights[0], 0, WCount-1, @WBest[0], 0, WCount-1);
end;


(*************************************************************************
Neural network training using early stopping (base algorithm - L-BFGS with
regularization).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    TrnXY       -   training set
    TrnSize     -   training set size
    ValXY       -   validation set
    ValSize     -   validation set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1, ...).
                    *  2, task has been solved, stopping  criterion  met -
                          sufficiently small step size.  Not expected  (we
                          use  EARLY  stopping)  but  possible  and not an
                          error.
                    *  6, task has been solved, stopping  criterion  met -
                          increasing of validation set error.
    Rep         -   training report

NOTE:

Algorithm stops if validation set error increases for  a  long  enough  or
step size is small enought  (there  are  task  where  validation  set  may
decrease for eternity). In any case solution returned corresponds  to  the
minimum of validation set error.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MLPTrainES(var Network : MultiLayerPerceptron;
     const TrnXY : TReal2DArray;
     TrnSize : AlglibInteger;
     const ValXY : TReal2DArray;
     ValSize : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport);
var
    I : AlglibInteger;
    Pass : AlglibInteger;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    WCount : AlglibInteger;
    W : TReal1DArray;
    WBest : TReal1DArray;
    E : Double;
    V : Double;
    EBest : Double;
    WFinal : TReal1DArray;
    EFinal : Double;
    ItBest : AlglibInteger;
    InternalRep : MinLBFGSReport;
    State : MinLBFGSState;
    WStep : Double;
begin
    WStep := Double(0.001);
    
    //
    // Test inputs, parse flags, read network geometry
    //
    if (TrnSize<=0) or (ValSize<=0) or (Restarts<1) or AP_FP_Less(Decay,0) then
    begin
        Info := -1;
        Exit;
    end;
    MLPProperties(Network, NIn, NOut, WCount);
    if MLPIsSoftmax(Network) then
    begin
        I:=0;
        while I<=TrnSize-1 do
        begin
            if (Round(TrnXY[I,NIn])<0) or (Round(TrnXY[I,NIn])>=NOut) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
        I:=0;
        while I<=ValSize-1 do
        begin
            if (Round(ValXY[I,NIn])<0) or (Round(ValXY[I,NIn])>=NOut) then
            begin
                Info := -2;
                Exit;
            end;
            Inc(I);
        end;
    end;
    Info := 2;
    
    //
    // Prepare
    //
    MLPInitPreprocessor(Network, TrnXY, TrnSize);
    SetLength(W, WCount-1+1);
    SetLength(WBest, WCount-1+1);
    SetLength(WFinal, WCount-1+1);
    EFinal := MaxRealNumber;
    I:=0;
    while I<=WCount-1 do
    begin
        WFinal[I] := 0;
        Inc(I);
    end;
    
    //
    // Multiple starts
    //
    Rep.NCholesky := 0;
    Rep.NHess := 0;
    Rep.NGrad := 0;
    Pass:=1;
    while Pass<=Restarts do
    begin
        
        //
        // Process
        //
        MLPRandomize(Network);
        EBest := MLPError(Network, ValXY, ValSize);
        APVMove(@WBest[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        ItBest := 0;
        APVMove(@W[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
        MinLBFGSCreate(WCount, Min(WCount, 10), W, State);
        MinLBFGSSetCond(State, Double(0.0), Double(0.0), WStep, 0);
        MinLBFGSSetXRep(State, True);
        while MinLBFGSIteration(State) do
        begin
            
            //
            // Calculate gradient
            //
            APVMove(@Network.Weights[0], 0, WCount-1, @State.X[0], 0, WCount-1);
            MLPGradNBatch(Network, TrnXY, TrnSize, State.F, State.G);
            V := APVDotProduct(@Network.Weights[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
            State.F := State.F+Double(0.5)*Decay*V;
            APVAdd(@State.G[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1, Decay);
            Rep.NGrad := Rep.NGrad+1;
            
            //
            // Validation set
            //
            if State.XUpdated then
            begin
                APVMove(@Network.Weights[0], 0, WCount-1, @W[0], 0, WCount-1);
                E := MLPError(Network, ValXY, ValSize);
                if AP_FP_Less(E,EBest) then
                begin
                    EBest := E;
                    APVMove(@WBest[0], 0, WCount-1, @Network.Weights[0], 0, WCount-1);
                    ItBest := InternalRep.IterationsCount;
                end;
                if (InternalRep.IterationsCount>30) and AP_FP_Greater(InternalRep.IterationsCount,Double(1.5)*ItBest) then
                begin
                    Info := 6;
                    Break;
                end;
            end;
        end;
        MinLBFGSResults(State, W, InternalRep);
        
        //
        // Compare with final answer
        //
        if AP_FP_Less(EBest,EFinal) then
        begin
            APVMove(@WFinal[0], 0, WCount-1, @WBest[0], 0, WCount-1);
            EFinal := EBest;
        end;
        Inc(Pass);
    end;
    
    //
    // The best network
    //
    APVMove(@Network.Weights[0], 0, WCount-1, @WFinal[0], 0, WCount-1);
end;


(*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - L-BFGS.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************)
procedure MLPKFoldCVLBFGS(const Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     WStep : Double;
     MaxIts : AlglibInteger;
     FoldsCount : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var CVRep : MLPCVReport);
begin
    MLPKFoldCVGeneral(Network, XY, NPoints, Decay, Restarts, FoldsCount, False, WStep, MaxIts, Info, Rep, CVRep);
end;


(*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - Levenberg-Marquardt.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************)
procedure MLPKFoldCVLM(const Network : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     FoldsCount : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var CVRep : MLPCVReport);
begin
    MLPKFoldCVGeneral(Network, XY, NPoints, Decay, Restarts, FoldsCount, True, Double(0.0), 0, Info, Rep, CVRep);
end;


(*************************************************************************
Internal cross-validation subroutine
*************************************************************************)
procedure MLPKFoldCVGeneral(const N : MultiLayerPerceptron;
     const XY : TReal2DArray;
     NPoints : AlglibInteger;
     Decay : Double;
     Restarts : AlglibInteger;
     FoldsCount : AlglibInteger;
     LMAlgorithm : Boolean;
     WStep : Double;
     MaxIts : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : MLPReport;
     var CVRep : MLPCVReport);
var
    I : AlglibInteger;
    Fold : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    Network : MultiLayerPerceptron;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    RowLen : AlglibInteger;
    WCount : AlglibInteger;
    NClasses : AlglibInteger;
    TSSize : AlglibInteger;
    CVSSize : AlglibInteger;
    CVSet : TReal2DArray;
    TestSet : TReal2DArray;
    Folds : TInteger1DArray;
    RelCnt : AlglibInteger;
    InternalRep : MLPReport;
    X : TReal1DArray;
    Y : TReal1DArray;
begin
    
    //
    // Read network geometry, test parameters
    //
    MLPProperties(N, NIn, NOut, WCount);
    if MLPIsSoftmax(N) then
    begin
        NClasses := NOut;
        RowLen := NIn+1;
    end
    else
    begin
        NClasses := -NOut;
        RowLen := NIn+NOut;
    end;
    if (NPoints<=0) or (FoldsCount<2) or (FoldsCount>NPoints) then
    begin
        Info := -1;
        Exit;
    end;
    MLPCopy(N, Network);
    
    //
    // K-fold out cross-validation.
    // First, estimate generalization error
    //
    SetLength(TestSet, NPoints-1+1, RowLen-1+1);
    SetLength(CVSet, NPoints-1+1, RowLen-1+1);
    SetLength(X, NIn-1+1);
    SetLength(Y, NOut-1+1);
    MLPKFoldSplit(XY, NPoints, NClasses, FoldsCount, False, Folds);
    CVRep.RelCLSError := 0;
    CVRep.AvgCE := 0;
    CVRep.RMSError := 0;
    CVRep.AvgError := 0;
    CVRep.AvgRelError := 0;
    Rep.NGrad := 0;
    Rep.NHess := 0;
    Rep.NCholesky := 0;
    RelCnt := 0;
    Fold:=0;
    while Fold<=FoldsCount-1 do
    begin
        
        //
        // Separate set
        //
        TSSize := 0;
        CVSSize := 0;
        I:=0;
        while I<=NPoints-1 do
        begin
            if Folds[I]=Fold then
            begin
                APVMove(@TestSet[TSSize][0], 0, RowLen-1, @XY[I][0], 0, RowLen-1);
                TSSize := TSSize+1;
            end
            else
            begin
                APVMove(@CVSet[CVSSize][0], 0, RowLen-1, @XY[I][0], 0, RowLen-1);
                CVSSize := CVSSize+1;
            end;
            Inc(I);
        end;
        
        //
        // Train on CV training set
        //
        if LMAlgorithm then
        begin
            MLPTrainLM(Network, CVSet, CVSSize, Decay, Restarts, Info, InternalRep);
        end
        else
        begin
            MLPTrainLBFGS(Network, CVSet, CVSSize, Decay, Restarts, WStep, MaxIts, Info, InternalRep);
        end;
        if Info<0 then
        begin
            CVRep.RelCLSError := 0;
            CVRep.AvgCE := 0;
            CVRep.RMSError := 0;
            CVRep.AvgError := 0;
            CVRep.AvgRelError := 0;
            Exit;
        end;
        Rep.NGrad := Rep.NGrad+InternalRep.NGrad;
        Rep.NHess := Rep.NHess+InternalRep.NHess;
        Rep.NCholesky := Rep.NCholesky+InternalRep.NCholesky;
        
        //
        // Estimate error using CV test set
        //
        if MLPIsSoftmax(Network) then
        begin
            
            //
            // classification-only code
            //
            CVRep.RelCLSError := CVRep.RelCLSError+MLPClsError(Network, TestSet, TSSize);
            CVRep.AvgCE := CVRep.AvgCE+MLPErrorN(Network, TestSet, TSSize);
        end;
        I:=0;
        while I<=TSSize-1 do
        begin
            APVMove(@X[0], 0, NIn-1, @TestSet[I][0], 0, NIn-1);
            MLPProcess(Network, X, Y);
            if MLPIsSoftmax(Network) then
            begin
                
                //
                // Classification-specific code
                //
                K := Round(TestSet[I,NIn]);
                J:=0;
                while J<=NOut-1 do
                begin
                    if J=K then
                    begin
                        CVRep.RMSError := CVRep.RMSError+AP_Sqr(Y[J]-1);
                        CVRep.AvgError := CVRep.AvgError+AbsReal(Y[J]-1);
                        CVRep.AvgRelError := CVRep.AvgRelError+AbsReal(Y[J]-1);
                        RelCnt := RelCnt+1;
                    end
                    else
                    begin
                        CVRep.RMSError := CVRep.RMSError+AP_Sqr(Y[J]);
                        CVRep.AvgError := CVRep.AvgError+AbsReal(Y[J]);
                    end;
                    Inc(J);
                end;
            end
            else
            begin
                
                //
                // Regression-specific code
                //
                J:=0;
                while J<=NOut-1 do
                begin
                    CVRep.RMSError := CVRep.RMSError+AP_Sqr(Y[J]-TestSet[I,NIn+J]);
                    CVRep.AvgError := CVRep.AvgError+AbsReal(Y[J]-TestSet[I,NIn+J]);
                    if AP_FP_Neq(TestSet[I,NIn+J],0) then
                    begin
                        CVRep.AvgRelError := CVRep.AvgRelError+AbsReal((Y[J]-TestSet[I,NIn+J])/TestSet[I,NIn+J]);
                        RelCnt := RelCnt+1;
                    end;
                    Inc(J);
                end;
            end;
            Inc(I);
        end;
        Inc(Fold);
    end;
    if MLPIsSoftMax(Network) then
    begin
        CVRep.RelCLSError := CVRep.RelCLSError/NPoints;
        CVRep.AvgCE := CVRep.AvgCE/(Ln(2)*NPoints);
    end;
    CVRep.RMSError := Sqrt(CVRep.RMSError/(NPoints*NOut));
    CVRep.AvgError := CVRep.AvgError/(NPoints*NOut);
    CVRep.AvgRelError := CVRep.AvgRelError/RelCnt;
    Info := 1;
end;


(*************************************************************************
Subroutine prepares K-fold split of the training set.

NOTES:
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.
*************************************************************************)
procedure MLPKFoldSplit(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NClasses : AlglibInteger;
     FoldsCount : AlglibInteger;
     StratifiedSplits : Boolean;
     var Folds : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
begin
    
    //
    // test parameters
    //
    Assert(NPoints>0, 'MLPKFoldSplit: wrong NPoints!');
    Assert((NClasses>1) or (NClasses<0), 'MLPKFoldSplit: wrong NClasses!');
    Assert((FoldsCount>=2) and (FoldsCount<=NPoints), 'MLPKFoldSplit: wrong FoldsCount!');
    Assert( not StratifiedSplits, 'MLPKFoldSplit: stratified splits are not supported!');
    
    //
    // Folds
    //
    SetLength(Folds, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        Folds[I] := I*FoldsCount div NPoints;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-2 do
    begin
        J := I+RandomInteger(NPoints-I);
        if J<>I then
        begin
            K := Folds[I];
            Folds[I] := Folds[J];
            Folds[J] := K;
        end;
        Inc(I);
    end;
end;


end.