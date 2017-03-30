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
unit minlbfgs;
interface
uses Math, Sysutils, Ap, linmin;

type
MinLBFGSState = record
    N : AlglibInteger;
    M : AlglibInteger;
    EpsG : Double;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    Flags : AlglibInteger;
    XRep : Boolean;
    StpMax : Double;
    NFEV : AlglibInteger;
    MCStage : AlglibInteger;
    K : AlglibInteger;
    Q : AlglibInteger;
    P : AlglibInteger;
    Rho : TReal1DArray;
    Y : TReal2DArray;
    S : TReal2DArray;
    Theta : TReal1DArray;
    D : TReal1DArray;
    Stp : Double;
    WORK : TReal1DArray;
    FOld : Double;
    GammaK : Double;
    X : TReal1DArray;
    F : Double;
    G : TReal1DArray;
    NeedFG : Boolean;
    XUpdated : Boolean;
    RState : RCommState;
    RepIterationsCount : AlglibInteger;
    RepNFEV : AlglibInteger;
    RepTerminationType : AlglibInteger;
    LState : LINMINState;
end;


MinLBFGSReport = record
    IterationsCount : AlglibInteger;
    NFEV : AlglibInteger;
    TerminationType : AlglibInteger;
end;



procedure MinLBFGSCreate(N : AlglibInteger;
     M : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLBFGSState);
procedure MinLBFGSSetCond(var State : MinLBFGSState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
procedure MinLBFGSSetXRep(var State : MinLBFGSState; NeedXRep : Boolean);
procedure MinLBFGSSetStpMax(var State : MinLBFGSState; StpMax : Double);
procedure MinLBFGSCreateX(N : AlglibInteger;
     M : AlglibInteger;
     const X : TReal1DArray;
     Flags : AlglibInteger;
     var State : MinLBFGSState);
function MinLBFGSIteration(var State : MinLBFGSState):Boolean;
procedure MinLBFGSResults(const State : MinLBFGSState;
     var X : TReal1DArray;
     var Rep : MinLBFGSReport);

implementation

procedure ClearRequestFields(var State : MinLBFGSState);forward;


(*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.

The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.

INPUT PARAMETERS:
    N       -   problem dimension. N>0
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   initial solution approximation, array[0..N-1].

OUTPUT PARAMETERS:
    State   -   structure used for reverse communication.
    
This function  initializes  State   structure  with  default  optimization
parameters (stopping conditions, step size, etc.). Use MinLBFGSSet??????()
functions to tune optimization parameters.

After   all   optimization   parameters   are   tuned,   you   should  use
MinLBFGSIteration() function to advance algorithm iterations.

NOTES:

1. you may tune stopping conditions with MinLBFGSSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLBFGSSetStpMax() function to bound algorithm's  steps.  However,
   L-BFGS rarely needs such a tuning.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSCreate(N : AlglibInteger;
     M : AlglibInteger;
     const X : TReal1DArray;
     var State : MinLBFGSState);
begin
    MinLBFGSCreateX(N, M, X, 0, State);
end;


(*************************************************************************
This function sets stopping conditions for L-BFGS optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLBFGSCreate()
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
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSSetCond(var State : MinLBFGSState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
begin
    Assert(AP_FP_Greater_Eq(EpsG,0), 'MinLBFGSSetCond: negative EpsG!');
    Assert(AP_FP_Greater_Eq(EpsF,0), 'MinLBFGSSetCond: negative EpsF!');
    Assert(AP_FP_Greater_Eq(EpsX,0), 'MinLBFGSSetCond: negative EpsX!');
    Assert(MaxIts>=0, 'MinLBFGSSetCond: negative MaxIts!');
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
                initialized with MinLBFGSCreate()
    NeedXRep-   whether iteration reports are needed or not

Usually algorithm returns  from  MinLBFGSIteration()  only when  it  needs
function/gradient/ (which is indicated by NeedFG field. However, with this
function we can let it  stop  after  each  iteration  (one  iteration  may
include more than one function evaluation), which is indicated by XUpdated
field.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSSetXRep(var State : MinLBFGSState; NeedXRep : Boolean);
begin
    State.XRep := NeedXRep;
end;


(*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinLBFGSCreate()
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSSetStpMax(var State : MinLBFGSState; StpMax : Double);
begin
    Assert(AP_FP_Greater_Eq(StpMax,0), 'MinLBFGSSetStpMax: StpMax<0!');
    State.StpMax := StpMax;
end;


(*************************************************************************
Extended subroutine for internal use only.

Accepts additional parameters:

    Flags - additional settings:
            * Flags = 0     means no additional settings
            * Flags = 1     "do not allocate memory". used when solving
                            a many subsequent tasks with  same N/M  values.
                            First  call MUST  be without this flag bit set,
                            subsequent  calls   of   MinLBFGS   with   same
                            MinLBFGSState structure can set Flags to 1.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSCreateX(N : AlglibInteger;
     M : AlglibInteger;
     const X : TReal1DArray;
     Flags : AlglibInteger;
     var State : MinLBFGSState);
var
    AllocateMem : Boolean;
begin
    Assert(N>=1, 'MinLBFGS: N too small!');
    Assert(M>=1, 'MinLBFGS: M too small!');
    Assert(M<=N, 'MinLBFGS: M too large!');
    
    //
    // Initialize
    //
    State.N := N;
    State.M := M;
    State.Flags := Flags;
    AllocateMem := Flags mod 2=0;
    Flags := Flags div 2;
    if AllocateMem then
    begin
        SetLength(State.Rho, M-1+1);
        SetLength(State.Theta, M-1+1);
        SetLength(State.Y, M-1+1, N-1+1);
        SetLength(State.S, M-1+1, N-1+1);
        SetLength(State.D, N-1+1);
        SetLength(State.X, N-1+1);
        SetLength(State.G, N-1+1);
        SetLength(State.WORK, N-1+1);
    end;
    MinLBFGSSetCond(State, 0, 0, 0, 0);
    MinLBFGSSetXRep(State, False);
    MinLBFGSSetStpMax(State, 0);
    
    //
    // Prepare first run
    //
    State.K := 0;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
    SetLength(State.RState.IA, 6+1);
    SetLength(State.RState.RA, 4+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
L-BFGS iterations

Called after initialization with MinLBFGSCreate() function.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLBFGSCreate()

RESULT:
* if function returned False, iterative proces has converged.
  Use MinLBFGSResults() to obtain optimization results.
* if subroutine returned True, then, depending on structure fields, we
  have one of the following situations


=== FUNC/GRAD REQUEST ===
State.NeedFG is True => function value/gradient are needed.
Caller should calculate function value State.F and gradient
State.G[0..N-1] at State.X[0..N-1] and call MinLBFGSIteration() again.

=== NEW INTERATION IS REPORTED ===
State.XUpdated is True => one more iteration was made.
State.X contains current position, State.F contains function value at X.
You can read info from these fields, but never modify  them  because  they
contain the only copy of optimization algorithm state.


One and only one of these fields (NeedFG, XUpdated) is true on return. New
iterations are reported only when reports  are  explicitly  turned  on  by
MinLBFGSSetXRep() function, so if you never called it, you can expect that
NeedFG is always True.


  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function MinLBFGSIteration(var State : MinLBFGSState):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    MaxIts : AlglibInteger;
    EpsF : Double;
    EpsG : Double;
    EpsX : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    IC : AlglibInteger;
    MCINFO : AlglibInteger;
    V : Double;
    VV : Double;
label
lbl_0, lbl_1, lbl_4, lbl_6, lbl_8, lbl_2, lbl_9, lbl_3, lbl_10, lbl_7, lbl_rcomm;
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
        MaxIts := State.RState.IA[2];
        I := State.RState.IA[3];
        J := State.RState.IA[4];
        IC := State.RState.IA[5];
        MCINFO := State.RState.IA[6];
        EpsF := State.RState.RA[0];
        EpsG := State.RState.RA[1];
        EpsX := State.RState.RA[2];
        V := State.RState.RA[3];
        VV := State.RState.RA[4];
    end
    else
    begin
        N := -983;
        M := -989;
        MaxIts := -834;
        I := 900;
        J := -287;
        IC := 364;
        MCINFO := 214;
        EpsF := -338;
        EpsG := -686;
        EpsX := 912;
        V := 585;
        VV := 497;
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
    
    //
    // Routine body
    //
    
    //
    // Unload frequently used variables from State structure
    // (just for typing convinience)
    //
    N := State.N;
    M := State.M;
    EpsG := State.EpsG;
    EpsF := State.EpsF;
    EpsX := State.EpsX;
    MaxIts := State.MaxIts;
    State.RepTerminationType := 0;
    State.RepIterationsCount := 0;
    State.RepNFEV := 0;
    
    //
    // Calculate F/G at the initial point
    //
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    if not State.XRep then
    begin
        goto lbl_4;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
lbl_4:
    State.RepNFEV := 1;
    State.FOld := State.F;
    V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
    V := Sqrt(V);
    if AP_FP_Less_Eq(V,EpsG) then
    begin
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    
    //
    // Choose initial step
    //
    if AP_FP_Eq(State.StpMax,0) then
    begin
        State.Stp := Min(Double(1.0)/V, 1);
    end
    else
    begin
        State.Stp := Min(Double(1.0)/V, State.StpMax);
    end;
    APVMoveNeg(@State.D[0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Main cycle
    //
lbl_6:
    if False then
    begin
        goto lbl_7;
    end;
    
    //
    // Main cycle: prepare to 1-D line search
    //
    State.P := State.K mod M;
    State.Q := Min(State.K, M-1);
    
    //
    // Store X[k], G[k]
    //
    APVMoveNeg(@State.S[State.P][0], 0, N-1, @State.X[0], 0, N-1);
    APVMoveNeg(@State.Y[State.P][0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Minimize F(x+alpha*d)
    // Calculate S[k], Y[k]
    //
    State.MCStage := 0;
    if State.K<>0 then
    begin
        State.Stp := Double(1.0);
    end;
    LinMinNormalizeD(State.D, State.Stp, N);
    MCSRCH(N, State.X, State.F, State.G, State.D, State.Stp, State.StpMax, MCINFO, State.NFEV, State.WORK, State.LState, State.MCStage);
lbl_8:
    if State.MCStage=0 then
    begin
        goto lbl_9;
    end;
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    MCSRCH(N, State.X, State.F, State.G, State.D, State.Stp, State.StpMax, MCINFO, State.NFEV, State.WORK, State.LState, State.MCStage);
    goto lbl_8;
lbl_9:
    if not State.XRep then
    begin
        goto lbl_10;
    end;
    
    //
    // report
    //
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 3;
    goto lbl_rcomm;
lbl_3:
lbl_10:
    State.RepNFEV := State.RepNFEV+State.NFEV;
    State.RepIterationsCount := State.RepIterationsCount+1;
    APVAdd(@State.S[State.P][0], 0, N-1, @State.X[0], 0, N-1);
    APVAdd(@State.Y[State.P][0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Stopping conditions
    //
    if (State.RepIterationsCount>=MaxIts) and (MaxIts>0) then
    begin
        
        //
        // Too many iterations
        //
        State.RepTerminationType := 5;
        Result := False;
        Exit;
    end;
    V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
    if AP_FP_Less_Eq(Sqrt(V),EpsG) then
    begin
        
        //
        // Gradient is small enough
        //
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    if AP_FP_Less_Eq(State.FOld-State.F,EpsF*Max(AbsReal(State.FOld), Max(AbsReal(State.F), Double(1.0)))) then
    begin
        
        //
        // F(k+1)-F(k) is small enough
        //
        State.RepTerminationType := 1;
        Result := False;
        Exit;
    end;
    V := APVDotProduct(@State.S[State.P][0], 0, N-1, @State.S[State.P][0], 0, N-1);
    if AP_FP_Less_Eq(Sqrt(V),EpsX) then
    begin
        
        //
        // X(k+1)-X(k) is small enough
        //
        State.RepTerminationType := 2;
        Result := False;
        Exit;
    end;
    
    //
    // If Wolfe conditions are satisfied, we can update
    // limited memory model.
    //
    // However, if conditions are not satisfied (NFEV limit is met,
    // function is too wild, ...), we'll skip L-BFGS update
    //
    if MCINFO<>1 then
    begin
        
        //
        // Skip update.
        //
        // In such cases we'll initialize search direction by
        // antigradient vector, because it  leads to more
        // transparent code with less number of special cases
        //
        State.FOld := State.F;
        APVMoveNeg(@State.D[0], 0, N-1, @State.G[0], 0, N-1);
    end
    else
    begin
        
        //
        // Calculate Rho[k], GammaK
        //
        V := APVDotProduct(@State.Y[State.P][0], 0, N-1, @State.S[State.P][0], 0, N-1);
        VV := APVDotProduct(@State.Y[State.P][0], 0, N-1, @State.Y[State.P][0], 0, N-1);
        if AP_FP_Eq(V,0) or AP_FP_Eq(VV,0) then
        begin
            
            //
            // Rounding errors make further iterations impossible.
            //
            State.RepTerminationType := -2;
            Result := False;
            Exit;
        end;
        State.Rho[State.P] := 1/V;
        State.GammaK := V/VV;
        
        //
        //  Calculate d(k+1) = -H(k+1)*g(k+1)
        //
        //  for I:=K downto K-Q do
        //      V = s(i)^T * work(iteration:I)
        //      theta(i) = V
        //      work(iteration:I+1) = work(iteration:I) - V*Rho(i)*y(i)
        //  work(last iteration) = H0*work(last iteration)
        //  for I:=K-Q to K do
        //      V = y(i)^T*work(iteration:I)
        //      work(iteration:I+1) = work(iteration:I) +(-V+theta(i))*Rho(i)*s(i)
        //
        //  NOW WORK CONTAINS d(k+1)
        //
        APVMove(@State.WORK[0], 0, N-1, @State.G[0], 0, N-1);
        I:=State.K;
        while I>=State.K-State.Q do
        begin
            IC := I mod M;
            V := APVDotProduct(@State.S[IC][0], 0, N-1, @State.Work[0], 0, N-1);
            State.Theta[IC] := V;
            VV := V*State.Rho[IC];
            APVSub(@State.Work[0], 0, N-1, @State.Y[IC][0], 0, N-1, VV);
            Dec(I);
        end;
        V := State.GammaK;
        APVMul(@State.Work[0], 0, N-1, V);
        I:=State.K-State.Q;
        while I<=State.K do
        begin
            IC := I mod M;
            V := APVDotProduct(@State.Y[IC][0], 0, N-1, @State.Work[0], 0, N-1);
            VV := State.Rho[IC]*(-V+State.Theta[IC]);
            APVAdd(@State.Work[0], 0, N-1, @State.S[IC][0], 0, N-1, VV);
            Inc(I);
        end;
        APVMoveNeg(@State.D[0], 0, N-1, @State.Work[0], 0, N-1);
        
        //
        // Next step
        //
        State.FOld := State.F;
        State.K := State.K+1;
    end;
    goto lbl_6;
lbl_7:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := M;
    State.RState.IA[2] := MaxIts;
    State.RState.IA[3] := I;
    State.RState.IA[4] := J;
    State.RState.IA[5] := IC;
    State.RState.IA[6] := MCINFO;
    State.RState.RA[0] := EpsF;
    State.RState.RA[1] := EpsG;
    State.RState.RA[2] := EpsX;
    State.RState.RA[3] := V;
    State.RState.RA[4] := VV;
end;


(*************************************************************************
L-BFGS algorithm results

Called after MinLBFGSIteration() returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by MinLBFGSIteration).

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSResults(const State : MinLBFGSState;
     var X : TReal1DArray;
     var Rep : MinLBFGSReport);
begin
    SetLength(X, State.N-1+1);
    APVMove(@X[0], 0, State.N-1, @State.X[0], 0, State.N-1);
    Rep.IterationsCount := State.RepIterationsCount;
    Rep.NFEV := State.RepNFEV;
    Rep.TerminationType := State.RepTerminationType;
end;


(*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************)
procedure ClearRequestFields(var State : MinLBFGSState);
begin
    State.NeedFG := False;
    State.XUpdated := False;
end;


end.