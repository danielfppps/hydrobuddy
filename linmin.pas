{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE

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
unit linmin;
interface
uses Math, Sysutils, Ap;

type
LINMINState = record
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



procedure LinMinNormalizeD(var D : TReal1DArray;
     var Stp : Double;
     N : AlglibInteger);
procedure MCSRCH(const N : AlglibInteger;
     var X : TReal1DArray;
     var F : Double;
     var G : TReal1DArray;
     const S : TReal1DArray;
     var STP : Double;
     STPMAX : Double;
     var INFO : AlglibInteger;
     var NFEV : AlglibInteger;
     var WA : TReal1DArray;
     var State : LINMINState;
     var Stage : AlglibInteger);

implementation

const
    FTOL = Double(0.001);
    XTOL = 100*MachineEpsilon;
    GTOL = Double(0.3);
    MAXFEV = 20;
    STPMIN = Double(1.0E-50);
    DefSTPMAX = Double(1.0E+50);

procedure MCSTEP(var STX : Double;
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
Normalizes direction/step pair: makes |D|=1, scales Stp.
If |D|=0, it returns, leavind D/Stp unchanged.

  -- ALGLIB --
     Copyright 01.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure LinMinNormalizeD(var D : TReal1DArray;
     var Stp : Double;
     N : AlglibInteger);
var
    MX : Double;
    S : Double;
    I : AlglibInteger;
begin
    
    //
    // first, scale D to avoid underflow/overflow durng squaring
    //
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        MX := Max(MX, AbsReal(D[I]));
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        Exit;
    end;
    S := 1/MX;
    APVMul(@D[0], 0, N-1, S);
    Stp := Stp/S;
    
    //
    // normalize D
    //
    S := APVDotProduct(@D[0], 0, N-1, @D[0], 0, N-1);
    S := 1/Sqrt(S);
    APVMul(@D[0], 0, N-1, S);
    Stp := Stp/S;
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
procedure MCSRCH(const N : AlglibInteger;
     var X : TReal1DArray;
     var F : Double;
     var G : TReal1DArray;
     const S : TReal1DArray;
     var STP : Double;
     STPMAX : Double;
     var INFO : AlglibInteger;
     var NFEV : AlglibInteger;
     var WA : TReal1DArray;
     var State : LINMINState;
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
    if AP_FP_Eq(STPMAX,0) then
    begin
        STPMAX := DefSTPMAX;
    end;
    if AP_FP_Less(STP,STPMIN) then
    begin
        STP := STPMIN;
    end;
    if AP_FP_Greater(STP,STPMAX) then
    begin
        STP := STPMAX;
    end;
    
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
                MCSTEP(State.STX, State.FXM, State.DGXM, State.STY, State.FYM, State.DGYM, STP, State.FM, State.DGM, State.BRACKT, State.STMIN, State.STMAX, State.INFOC);
                
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
                MCSTEP(State.STX, State.FX, State.DGX, State.STY, State.FY, State.DGY, STP, F, State.DG, State.BRACKT, State.STMIN, State.STMAX, State.INFOC);
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


procedure MCSTEP(var STX : Double;
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