{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
unit minlm;
interface
uses Math, Sysutils, Ap, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs;

type
MinLMState = record
    WrongParams : Boolean;
    N : AlglibInteger;
    M : AlglibInteger;
    EpsG : Double;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    XRep : Boolean;
    StpMax : Double;
    Flags : AlglibInteger;
    UserMode : AlglibInteger;
    X : TReal1DArray;
    F : Double;
    FI : TReal1DArray;
    J : TReal2DArray;
    H : TReal2DArray;
    G : TReal1DArray;
    NeedF : Boolean;
    NeedFG : Boolean;
    NeedFGH : Boolean;
    NeedFiJ : Boolean;
    XUpdated : Boolean;
    InternalState : MinLBFGSState;
    InternalRep : MinLBFGSReport;
    XPrec : TReal1DArray;
    XBase : TReal1DArray;
    XDir : TReal1DArray;
    GBase : TReal1DArray;
    XPrev : TReal1DArray;
    FPrev : Double;
    RawModel : TReal2DArray;
    Model : TReal2DArray;
    WORK : TReal1DArray;
    RState : RCommState;
    RepIterationsCount : AlglibInteger;
    RepTerminationType : AlglibInteger;
    RepNFunc : AlglibInteger;
    RepNJac : AlglibInteger;
    RepNGrad : AlglibInteger;
    RepNHess : AlglibInteger;
    RepNCholesky : AlglibInteger;
    SolverInfo : AlglibInteger;
    SolverRep : DenseSolverReport;
    InvInfo : AlglibInteger;
    InvRep : MatInvReport;
end;


MinLMReport = record
    IterationsCount : AlglibInteger;
    TerminationType : AlglibInteger;
    NFunc : AlglibInteger;
    NJac : AlglibInteger;
    NGrad : AlglibInteger;
    NHess : AlglibInteger;
    NCholesky : AlglibInteger;
end;



procedure MinLMCreateFGH(const N : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLMState);
procedure MinLMCreateFGJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLMState);
procedure MinLMCreateFJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLMState);
procedure MinLMSetCond(var State : MinLMState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
procedure MinLMSetXRep(var State : MinLMState; NeedXRep : Boolean);
procedure MinLMSetStpMax(var State : MinLMState; StpMax : Double);
function MinLMIteration(var State : MinLMState):Boolean;
procedure MinLMResults(const State : MinLMState;
     var X : TReal1DArray;
     var Rep : MinLMReport);

implementation

const
    LMModeFJ = 0;
    LMModeFGJ = 1;
    LMModeFGH = 2;
    LMFlagNoPreLBFGS = 1;
    LMFlagNoIntLBFGS = 2;
    LMPreLBFGSM = 5;
    LMIntLBFGSIts = 5;
    LBFGSNoRealloc = 1;

procedure LMPrepare(N : AlglibInteger;
     M : AlglibInteger;
     HaveGrad : Boolean;
     var State : MinLMState);forward;
procedure LMClearRequestFields(var State : MinLMState);forward;
function IncreaseLambda(var Lambda : Double;
     var Nu : Double;
     LambdaUp : Double):Boolean;forward;
procedure DecreaseLambda(var Lambda : Double;
     var Nu : Double;
     LambdaDown : Double);forward;


(*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Hessian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F has general form (not "sum-of-squares"):

    F = F(x[0], ..., x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTES:

1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMCreateFGH(const N : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLMState);
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 7+1);
    State.RState.Stage := -1;
    
    //
    // prepare internal structures
    //
    LMPrepare(N, 0, True, State);
    
    //
    // initialize, check parameters
    //
    MinLMSetCond(State, 0, 0, 0, 0);
    MinLMSetXRep(State, False);
    MinLMSetStpMax(State, 0);
    State.N := N;
    State.M := 0;
    State.Flags := 0;
    State.UserMode := LMModeFGH;
    State.WrongParams := False;
    if N<1 then
    begin
        State.WrongParams := True;
        Exit;
    end;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
end;


(*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Jacobian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTES:

1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMCreateFGJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLMState);
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 7+1);
    State.RState.Stage := -1;
    
    //
    // prepare internal structures
    //
    LMPrepare(N, M, True, State);
    
    //
    // initialize, check parameters
    //
    MinLMSetCond(State, 0, 0, 0, 0);
    MinLMSetXRep(State, False);
    MinLMSetStpMax(State, 0);
    State.N := N;
    State.M := M;
    State.Flags := 0;
    State.UserMode := LMModeFGJ;
    State.WrongParams := False;
    if N<1 then
    begin
        State.WrongParams := True;
        Exit;
    end;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
end;


(*************************************************************************
    CLASSIC LEVENBERG-MARQUARDT METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using Jacobi matrix. Algorithm  -  classic Levenberg-Marquardt
method.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTES:

1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMCreateFJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLMState);
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 7+1);
    State.RState.Stage := -1;
    
    //
    // prepare internal structures
    //
    LMPrepare(N, M, True, State);
    
    //
    // initialize, check parameters
    //
    MinLMSetCond(State, 0, 0, 0, 0);
    MinLMSetXRep(State, False);
    MinLMSetStpMax(State, 0);
    State.N := N;
    State.M := M;
    State.Flags := 0;
    State.UserMode := LMModeFJ;
    State.WrongParams := False;
    if N<1 then
    begin
        State.WrongParams := True;
        Exit;
    end;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
end;


(*************************************************************************
This function sets stopping conditions for Levenberg-Marquardt optimization
algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLMCreate???()
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                ||G||<EpsG is satisfied, where ||.|| means Euclidian norm,
                G - gradient.
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations   is    unlimited.   Only   Levenberg-Marquardt
                iterations  are  counted  (L-BFGS/CG  iterations  are  NOT
                counted  because their cost is very low copared to that of
                LM).

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLMSetCond(var State : MinLMState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
begin
    Assert(AP_FP_Greater_Eq(EpsG,0), 'MinLMSetCond: negative EpsG!');
    Assert(AP_FP_Greater_Eq(EpsF,0), 'MinLMSetCond: negative EpsF!');
    Assert(AP_FP_Greater_Eq(EpsX,0), 'MinLMSetCond: negative EpsX!');
    Assert(MaxIts>=0, 'MinLMSetCond: negative MaxIts!');
    if AP_FP_Eq(EpsG,0) and AP_FP_Eq(EpsF,0) and AP_FP_Eq(EpsX,0) and (MaxIts=0) then
    begin
        EpsX := Double(1.0E-6);
    end;
    State.EpsG := EpsG;
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
end;


(*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinLMCreate???()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns  from  MinLMIteration()  only  when  it  needs
function/gradient/Hessian. However, with this function we can let it  stop
after  each  iteration  (one iteration may include  more than one function
evaluation), which is indicated by XUpdated field.

Both Levenberg-Marquardt and L-BFGS iterations are reported.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLMSetXRep(var State : MinLMState; NeedXRep : Boolean);
begin
    State.XRep := NeedXRep;
end;


(*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate???()
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLMSetStpMax(var State : MinLMState; StpMax : Double);
begin
    Assert(AP_FP_Greater_Eq(StpMax,0), 'MinLMSetStpMax: StpMax<0!');
    State.StpMax := StpMax;
end;


(*************************************************************************
One Levenberg-Marquardt iteration.

Called after inialization of State structure with MinLMXXX subroutine.
See HTML docs for examples.

Input parameters:
    State   -   structure which stores algorithm state between subsequent
                calls and which is used for reverse communication. Must be
                initialized with MinLMXXX call first.

If subroutine returned False, iterative algorithm has converged.

If subroutine returned True, then:
* if State.NeedF=True,      -   function value F at State.X[0..N-1]
                                is required
* if State.NeedFG=True      -   function value F and gradient G
                                are required
* if State.NeedFiJ=True     -   function vector f[i] and Jacobi matrix J
                                are required
* if State.NeedFGH=True     -   function value F, gradient G and Hesian H
                                are required
* if State.XUpdated=True    -   algorithm reports about new iteration,
                                State.X contains current point,
                                State.F contains function value.

One and only one of this fields can be set at time.

Results are stored:
* function value            -   in MinLMState.F
* gradient                  -   in MinLMState.G[0..N-1]
* Jacobi matrix             -   in MinLMState.J[0..M-1,0..N-1]
* Hessian                   -   in MinLMState.H[0..N-1,0..N-1]

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
function MinLMIteration(var State : MinLMState):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    I : AlglibInteger;
    StepNorm : Double;
    SPD : Boolean;
    FBase : Double;
    FNew : Double;
    Lambda : Double;
    Nu : Double;
    LambdaUp : Double;
    LambdaDown : Double;
    LBFGSFlags : AlglibInteger;
    V : Double;
label
lbl_18, lbl_0, lbl_20, lbl_1, lbl_22, lbl_19, lbl_16, lbl_2, lbl_3, lbl_24, lbl_17, lbl_4, lbl_26, lbl_5, lbl_28, lbl_30, lbl_34, lbl_6, lbl_35, lbl_32, lbl_38, lbl_7, lbl_39, lbl_36, lbl_42, lbl_8, lbl_43, lbl_40, lbl_9, lbl_46, lbl_10, lbl_47, lbl_44, lbl_52, lbl_11, lbl_53, lbl_50, lbl_48, lbl_12, lbl_54, lbl_13, lbl_56, lbl_14, lbl_58, lbl_31, lbl_15, lbl_60, lbl_rcomm;
begin
    
    //
    // Reverse communication preparations
    // I know it looks ugly, but it works the same way
    // anywhere from C++ to Python.
    //
    // This code initializes locals by:
    // * random values determined during code
    //   generation - on first subroutine call
    // * values from previous call - on subsequent calls
    //
    if State.RState.Stage>=0 then
    begin
        N := State.RState.IA[0];
        M := State.RState.IA[1];
        I := State.RState.IA[2];
        LBFGSFlags := State.RState.IA[3];
        SPD := State.RState.BA[0];
        StepNorm := State.RState.RA[0];
        FBase := State.RState.RA[1];
        FNew := State.RState.RA[2];
        Lambda := State.RState.RA[3];
        Nu := State.RState.RA[4];
        LambdaUp := State.RState.RA[5];
        LambdaDown := State.RState.RA[6];
        V := State.RState.RA[7];
    end
    else
    begin
        N := -983;
        M := -989;
        I := -834;
        LBFGSFlags := 900;
        SPD := True;
        StepNorm := 364;
        FBase := 214;
        FNew := -338;
        Lambda := -686;
        Nu := 912;
        LambdaUp := 585;
        LambdaDown := 497;
        V := -271;
    end;
    if State.RState.Stage=0 then
    begin
        goto lbl_0;
    end;
    if State.RState.Stage=1 then
    begin
        goto lbl_1;
    end;
    if State.RState.Stage=2 then
    begin
        goto lbl_2;
    end;
    if State.RState.Stage=3 then
    begin
        goto lbl_3;
    end;
    if State.RState.Stage=4 then
    begin
        goto lbl_4;
    end;
    if State.RState.Stage=5 then
    begin
        goto lbl_5;
    end;
    if State.RState.Stage=6 then
    begin
        goto lbl_6;
    end;
    if State.RState.Stage=7 then
    begin
        goto lbl_7;
    end;
    if State.RState.Stage=8 then
    begin
        goto lbl_8;
    end;
    if State.RState.Stage=9 then
    begin
        goto lbl_9;
    end;
    if State.RState.Stage=10 then
    begin
        goto lbl_10;
    end;
    if State.RState.Stage=11 then
    begin
        goto lbl_11;
    end;
    if State.RState.Stage=12 then
    begin
        goto lbl_12;
    end;
    if State.RState.Stage=13 then
    begin
        goto lbl_13;
    end;
    if State.RState.Stage=14 then
    begin
        goto lbl_14;
    end;
    if State.RState.Stage=15 then
    begin
        goto lbl_15;
    end;
    
    //
    // Routine body
    //
    Assert((State.UserMode=LMModeFJ) or (State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFGH), 'LM: internal error');
    if State.WrongParams then
    begin
        State.RepTerminationType := -1;
        Result := False;
        Exit;
    end;
    
    //
    // prepare params
    //
    N := State.N;
    M := State.M;
    LambdaUp := 20;
    LambdaDown := Double(0.5);
    Nu := 1;
    LBFGSFlags := 0;
    
    //
    // if we have F and G
    //
    if not(((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFGH)) and (State.Flags div LMFlagNoPreLBFGS mod 2=0)) then
    begin
        goto lbl_16;
    end;
    
    //
    // First stage of the hybrid algorithm: LBFGS
    //
    MinLBFGSCreate(N, Min(N, LMPreLBFGSM), State.X, State.InternalState);
    MinLBFGSSetCond(State.InternalState, 0, 0, 0, Max(5, N));
    MinLBFGSSetXRep(State.InternalState, State.XRep);
    MinLBFGSSetStpMax(State.InternalState, State.StpMax);
lbl_18:
    if not MinLBFGSIteration(State.InternalState) then
    begin
        goto lbl_19;
    end;
    if not State.InternalState.NeedFG then
    begin
        goto lbl_20;
    end;
    
    //
    // RComm
    //
    APVMove(@State.X[0], 0, N-1, @State.InternalState.X[0], 0, N-1);
    LMClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    
    //
    // Call LBFGS
    //
    State.InternalState.F := State.F;
    APVMove(@State.InternalState.G[0], 0, N-1, @State.G[0], 0, N-1);
lbl_20:
    if not(State.InternalState.XUpdated and State.XRep) then
    begin
        goto lbl_22;
    end;
    LMClearRequestFields(State);
    State.F := State.InternalState.F;
    APVMove(@State.X[0], 0, N-1, @State.InternalState.X[0], 0, N-1);
    State.XUpdated := True;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
lbl_22:
    goto lbl_18;
lbl_19:
    MinLBFGSResults(State.InternalState, State.X, State.InternalRep);
    goto lbl_17;
lbl_16:
    
    //
    // No first stage.
    // However, we may need to report initial point
    //
    if not State.XRep then
    begin
        goto lbl_24;
    end;
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    LMClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 3;
    goto lbl_rcomm;
lbl_3:
lbl_24:
lbl_17:
    
    //
    // Second stage of the hybrid algorithm: LM
    // Initialize quadratic model.
    //
    if State.UserMode<>LMModeFGH then
    begin
        goto lbl_26;
    end;
    
    //
    // RComm
    //
    LMClearRequestFields(State);
    State.NeedFGH := True;
    State.RState.Stage := 4;
    goto lbl_rcomm;
lbl_4:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    State.RepNHess := State.RepNHess+1;
    
    //
    // generate raw quadratic model
    //
    RMatrixCopy(N, N, State.H, 0, 0, State.RawModel, 0, 0);
    APVMove(@State.GBase[0], 0, N-1, @State.G[0], 0, N-1);
    FBase := State.F;
lbl_26:
    if not((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFJ)) then
    begin
        goto lbl_28;
    end;
    
    //
    // RComm
    //
    LMClearRequestFields(State);
    State.NeedFiJ := True;
    State.RState.Stage := 5;
    goto lbl_rcomm;
lbl_5:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNJac := State.RepNJac+1;
    
    //
    // generate raw quadratic model
    //
    RMatrixGEMM(N, N, M, Double(2.0), State.J, 0, 0, 1, State.J, 0, 0, 0, Double(0.0), State.RawModel, 0, 0);
    RMatrixMV(N, M, State.J, 0, 0, 1, State.FI, 0, State.GBase, 0);
    APVMul(@State.GBase[0], 0, N-1, 2);
    FBase := APVDotProduct(@State.FI[0], 0, M-1, @State.FI[0], 0, M-1);
lbl_28:
    Lambda := Double(0.001);
lbl_30:
    if False then
    begin
        goto lbl_31;
    end;
    
    //
    // 1. Model = RawModel+lambda*I
    // 2. Try to solve (RawModel+Lambda*I)*dx = -g.
    //    Increase lambda if left part is not positive definite.
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@State.Model[I][0], 0, N-1, @State.RawModel[I][0], 0, N-1);
        State.Model[I,I] := State.Model[I,I]+Lambda;
        Inc(I);
    end;
    SPD := SPDMatrixCholesky(State.Model, N, True);
    State.RepNCholesky := State.RepNCholesky+1;
    if SPD then
    begin
        goto lbl_32;
    end;
    if not IncreaseLambda(Lambda, Nu, LambdaUp) then
    begin
        goto lbl_34;
    end;
    goto lbl_30;
    goto lbl_35;
lbl_34:
    State.RepTerminationType := 7;
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 6;
    goto lbl_rcomm;
lbl_6:
    goto lbl_31;
lbl_35:
lbl_32:
    SPDMatrixCholeskySolve(State.Model, N, True, State.GBase, State.SolverInfo, State.SolverRep, State.XDir);
    if State.SolverInfo>=0 then
    begin
        goto lbl_36;
    end;
    if not IncreaseLambda(Lambda, Nu, LambdaUp) then
    begin
        goto lbl_38;
    end;
    goto lbl_30;
    goto lbl_39;
lbl_38:
    State.RepTerminationType := 7;
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 7;
    goto lbl_rcomm;
lbl_7:
    goto lbl_31;
lbl_39:
lbl_36:
    APVMul(@State.XDir[0], 0, N-1, -1);
    
    //
    // Candidate lambda is found.
    // 1. Save old w in WBase
    // 1. Test some stopping criterions
    // 2. If error(w+wdir)>error(w), increase lambda
    //
    APVMove(@State.XPrev[0], 0, N-1, @State.X[0], 0, N-1);
    State.FPrev := State.F;
    APVMove(@State.XBase[0], 0, N-1, @State.X[0], 0, N-1);
    APVAdd(@State.X[0], 0, N-1, @State.XDir[0], 0, N-1);
    StepNorm := APVDotProduct(@State.XDir[0], 0, N-1, @State.XDir[0], 0, N-1);
    StepNorm := Sqrt(StepNorm);
    if not(AP_FP_Greater(State.StpMax,0) and AP_FP_Greater(StepNorm,State.StpMax)) then
    begin
        goto lbl_40;
    end;
    
    //
    // Step is larger than the limit,
    // larger lambda is needed
    //
    APVMove(@State.X[0], 0, N-1, @State.XBase[0], 0, N-1);
    if not IncreaseLambda(Lambda, Nu, LambdaUp) then
    begin
        goto lbl_42;
    end;
    goto lbl_30;
    goto lbl_43;
lbl_42:
    State.RepTerminationType := 7;
    APVMove(@State.X[0], 0, N-1, @State.XPrev[0], 0, N-1);
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 8;
    goto lbl_rcomm;
lbl_8:
    goto lbl_31;
lbl_43:
lbl_40:
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 9;
    goto lbl_rcomm;
lbl_9:
    State.RepNFunc := State.RepNFunc+1;
    FNew := State.F;
    if AP_FP_Less_Eq(FNew,FBase) then
    begin
        goto lbl_44;
    end;
    
    //
    // restore state and continue search for lambda
    //
    APVMove(@State.X[0], 0, N-1, @State.XBase[0], 0, N-1);
    if not IncreaseLambda(Lambda, Nu, LambdaUp) then
    begin
        goto lbl_46;
    end;
    goto lbl_30;
    goto lbl_47;
lbl_46:
    State.RepTerminationType := 7;
    APVMove(@State.X[0], 0, N-1, @State.XPrev[0], 0, N-1);
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 10;
    goto lbl_rcomm;
lbl_10:
    goto lbl_31;
lbl_47:
lbl_44:
    if not(AP_FP_Eq(State.StpMax,0) and ((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFGH)) and (State.Flags div LMFlagNoIntLBFGS mod 2=0)) then
    begin
        goto lbl_48;
    end;
    
    //
    // Optimize using LBFGS, with inv(cholesky(H)) as preconditioner.
    //
    // It is possible only when StpMax=0, because we can't guarantee
    // that step remains bounded when preconditioner is used (we need
    // SVD decomposition to do that, which is too slow).
    //
    RMatrixTRInverse(State.Model, N, True, False, State.InvInfo, State.InvRep);
    if State.InvInfo<=0 then
    begin
        goto lbl_50;
    end;
    
    //
    // if matrix can be inverted, use it.
    // just silently move to next iteration otherwise.
    // (will be very, very rare, mostly for specially
    // designed near-degenerate tasks)
    //
    APVMove(@State.XBase[0], 0, N-1, @State.X[0], 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        State.XPrec[I] := 0;
        Inc(I);
    end;
    MinLBFGSCreateX(N, Min(N, LMIntLBFGSIts), State.XPrec, LBFGSFlags, State.InternalState);
    MinLBFGSSetCond(State.InternalState, 0, 0, 0, LMIntLBFGSIts);
lbl_52:
    if not MinLBFGSIteration(State.InternalState) then
    begin
        goto lbl_53;
    end;
    
    //
    // convert XPrec to unpreconditioned form, then call RComm.
    //
    I:=0;
    while I<=N-1 do
    begin
        V := APVDotProduct(@State.InternalState.X[0], I, N-1, @State.Model[I][0], I, N-1);
        State.X[I] := State.XBase[I]+V;
        Inc(I);
    end;
    LMClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 11;
    goto lbl_rcomm;
lbl_11:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    
    //
    // 1. pass State.F to State.InternalState.F
    // 2. convert gradient back to preconditioned form
    //
    State.InternalState.F := State.F;
    I:=0;
    while I<=N-1 do
    begin
        State.InternalState.G[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        V := State.G[I];
        APVAdd(@State.InternalState.G[0], I, N-1, @State.Model[I][0], I, N-1, V);
        Inc(I);
    end;
    
    //
    // next iteration
    //
    goto lbl_52;
lbl_53:
    
    //
    // change LBFGS flags to NoRealloc.
    // L-BFGS subroutine will use memory allocated from previous run.
    // it is possible since all subsequent calls will be with same N/M.
    //
    LBFGSFlags := LBFGSNoRealloc;
    
    //
    // back to unpreconditioned X
    //
    MinLBFGSResults(State.InternalState, State.XPrec, State.InternalRep);
    I:=0;
    while I<=N-1 do
    begin
        V := APVDotProduct(@State.XPrec[0], I, N-1, @State.Model[I][0], I, N-1);
        State.X[I] := State.XBase[I]+V;
        Inc(I);
    end;
lbl_50:
lbl_48:
    
    //
    // Composite iteration is almost over:
    // * accept new position.
    // * rebuild quadratic model
    //
    State.RepIterationsCount := State.RepIterationsCount+1;
    if State.UserMode<>LMModeFGH then
    begin
        goto lbl_54;
    end;
    LMClearRequestFields(State);
    State.NeedFGH := True;
    State.RState.Stage := 12;
    goto lbl_rcomm;
lbl_12:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    State.RepNHess := State.RepNHess+1;
    RMatrixCopy(N, N, State.H, 0, 0, State.RawModel, 0, 0);
    APVMove(@State.GBase[0], 0, N-1, @State.G[0], 0, N-1);
    FNew := State.F;
lbl_54:
    if not((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFJ)) then
    begin
        goto lbl_56;
    end;
    LMClearRequestFields(State);
    State.NeedFiJ := True;
    State.RState.Stage := 13;
    goto lbl_rcomm;
lbl_13:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNJac := State.RepNJac+1;
    RMatrixGEMM(N, N, M, Double(2.0), State.J, 0, 0, 1, State.J, 0, 0, 0, Double(0.0), State.RawModel, 0, 0);
    RMatrixMV(N, M, State.J, 0, 0, 1, State.FI, 0, State.GBase, 0);
    APVMul(@State.GBase[0], 0, N-1, 2);
    FNew := APVDotProduct(@State.FI[0], 0, M-1, @State.FI[0], 0, M-1);
lbl_56:
    
    //
    // Stopping conditions
    //
    APVMove(@State.WORK[0], 0, N-1, @State.XPrev[0], 0, N-1);
    APVSub(@State.WORK[0], 0, N-1, @State.X[0], 0, N-1);
    StepNorm := APVDotProduct(@State.WORK[0], 0, N-1, @State.WORK[0], 0, N-1);
    StepNorm := Sqrt(StepNorm);
    if AP_FP_Less_Eq(StepNorm,State.EpsX) then
    begin
        State.RepTerminationType := 2;
        goto lbl_31;
    end;
    if (State.RepIterationsCount>=State.MaxIts) and (State.MaxIts>0) then
    begin
        State.RepTerminationType := 5;
        goto lbl_31;
    end;
    V := APVDotProduct(@State.GBase[0], 0, N-1, @State.GBase[0], 0, N-1);
    V := Sqrt(V);
    if AP_FP_Less_Eq(V,State.EpsG) then
    begin
        State.RepTerminationType := 4;
        goto lbl_31;
    end;
    if AP_FP_Less_Eq(AbsReal(FNew-FBase),State.EpsF*Max(1, Max(AbsReal(FNew), AbsReal(FBase)))) then
    begin
        State.RepTerminationType := 1;
        goto lbl_31;
    end;
    
    //
    // Now, iteration is finally over:
    // * update FBase
    // * decrease lambda
    // * report new iteration
    //
    if not State.XRep then
    begin
        goto lbl_58;
    end;
    LMClearRequestFields(State);
    State.XUpdated := True;
    State.F := FNew;
    State.RState.Stage := 14;
    goto lbl_rcomm;
lbl_14:
lbl_58:
    FBase := FNew;
    DecreaseLambda(Lambda, Nu, LambdaDown);
    goto lbl_30;
lbl_31:
    
    //
    // final point is reported
    //
    if not State.XRep then
    begin
        goto lbl_60;
    end;
    LMClearRequestFields(State);
    State.XUpdated := True;
    State.F := FNew;
    State.RState.Stage := 15;
    goto lbl_rcomm;
lbl_15:
lbl_60:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := M;
    State.RState.IA[2] := I;
    State.RState.IA[3] := LBFGSFlags;
    State.RState.BA[0] := SPD;
    State.RState.RA[0] := StepNorm;
    State.RState.RA[1] := FBase;
    State.RState.RA[2] := FNew;
    State.RState.RA[3] := Lambda;
    State.RState.RA[4] := Nu;
    State.RState.RA[5] := LambdaUp;
    State.RState.RA[6] := LambdaDown;
    State.RState.RA[7] := V;
end;


(*************************************************************************
Levenberg-Marquardt algorithm results

Called after MinLMIteration returned False.

Input parameters:
    State   -   algorithm state (used by MinLMIteration).

Output parameters:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient is no more than EpsG.
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * Rep.NFunc     - number of function calculations
                * Rep.NJac      - number of Jacobi matrix calculations
                * Rep.NGrad     - number of gradient calculations
                * Rep.NHess     - number of Hessian calculations
                * Rep.NCholesky - number of Cholesky decomposition calculations

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMResults(const State : MinLMState;
     var X : TReal1DArray;
     var Rep : MinLMReport);
begin
    SetLength(X, State.N-1+1);
    APVMove(@X[0], 0, State.N-1, @State.X[0], 0, State.N-1);
    Rep.IterationsCount := State.RepIterationsCount;
    Rep.TerminationType := State.RepTerminationType;
    Rep.NFunc := State.RepNFunc;
    Rep.NJac := State.RepNJac;
    Rep.NGrad := State.RepNGrad;
    Rep.NHess := State.RepNHess;
    Rep.NCholesky := State.RepNCholesky;
end;


(*************************************************************************
Prepare internal structures (except for RComm).

Note: M must be zero for FGH mode, non-zero for FJ/FGJ mode.
*************************************************************************)
procedure LMPrepare(N : AlglibInteger;
     M : AlglibInteger;
     HaveGrad : Boolean;
     var State : MinLMState);
begin
    State.RepIterationsCount := 0;
    State.RepTerminationType := 0;
    State.RepNFunc := 0;
    State.RepNJac := 0;
    State.RepNGrad := 0;
    State.RepNHess := 0;
    State.RepNCholesky := 0;
    if (N<=0) or (M<0) then
    begin
        Exit;
    end;
    if HaveGrad then
    begin
        SetLength(State.G, N-1+1);
    end;
    if M<>0 then
    begin
        SetLength(State.J, M-1+1, N-1+1);
        SetLength(State.FI, M-1+1);
        SetLength(State.H, 0+1, 0+1);
    end
    else
    begin
        SetLength(State.J, 0+1, 0+1);
        SetLength(State.FI, 0+1);
        SetLength(State.H, N-1+1, N-1+1);
    end;
    SetLength(State.X, N-1+1);
    SetLength(State.RawModel, N-1+1, N-1+1);
    SetLength(State.Model, N-1+1, N-1+1);
    SetLength(State.XBase, N-1+1);
    SetLength(State.XPrec, N-1+1);
    SetLength(State.GBase, N-1+1);
    SetLength(State.XDir, N-1+1);
    SetLength(State.XPrev, N-1+1);
    SetLength(State.Work, Max(N, M)+1);
end;


(*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************)
procedure LMClearRequestFields(var State : MinLMState);
begin
    State.NeedF := False;
    State.NeedFG := False;
    State.NeedFGH := False;
    State.NeedFiJ := False;
    State.XUpdated := False;
end;


(*************************************************************************
Increases lambda, returns False when there is a danger of overflow
*************************************************************************)
function IncreaseLambda(var Lambda : Double;
     var Nu : Double;
     LambdaUp : Double):Boolean;
var
    LnLambda : Double;
    LnNu : Double;
    LnLambdaUp : Double;
    LnMax : Double;
begin
    Result := False;
    LnLambda := Ln(Lambda);
    LnLambdaUp := Ln(LambdaUp);
    LnNu := Ln(Nu);
    LnMax := Ln(MaxRealNumber);
    if AP_FP_Greater(LnLambda+LnLambdaUp+LnNu,LnMax) then
    begin
        Exit;
    end;
    if AP_FP_Greater(LnNu+Ln(2),LnMax) then
    begin
        Exit;
    end;
    Lambda := Lambda*LambdaUp*Nu;
    Nu := Nu*2;
    Result := True;
end;


(*************************************************************************
Decreases lambda, but leaves it unchanged when there is danger of underflow.
*************************************************************************)
procedure DecreaseLambda(var Lambda : Double;
     var Nu : Double;
     LambdaDown : Double);
begin
    Nu := 1;
    if AP_FP_Less(Ln(Lambda)+Ln(LambdaDown),Ln(MinRealNumber)) then
    begin
        Lambda := MinRealNumber;
    end
    else
    begin
        Lambda := Lambda*LambdaDown;
    end;
end;


end.