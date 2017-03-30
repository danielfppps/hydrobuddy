{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
This file is a part of ALGLIB project.

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
unit safesolve;
interface
uses Math, Sysutils, Ap;

function RMatrixScaledTRSafeSolve(const A : TReal2DArray;
     SA : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     IsUpper : Boolean;
     Trans : AlglibInteger;
     IsUnit : Boolean;
     MaxGrowth : Double):Boolean;
function CMatrixScaledTRSafeSolve(const A : TComplex2DArray;
     SA : Double;
     N : AlglibInteger;
     var X : TComplex1DArray;
     IsUpper : Boolean;
     Trans : AlglibInteger;
     IsUnit : Boolean;
     MaxGrowth : Double):Boolean;

implementation

function CBasicSolveAndUpdate(Alpha : Complex;
     Beta : Complex;
     LnMax : Double;
     BNorm : Double;
     MaxGrowth : Double;
     var XNorm : Double;
     var X : Complex):Boolean;forward;


(*************************************************************************
Real implementation of CMatrixScaledTRSafeSolve

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************)
function RMatrixScaledTRSafeSolve(const A : TReal2DArray;
     SA : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     IsUpper : Boolean;
     Trans : AlglibInteger;
     IsUnit : Boolean;
     MaxGrowth : Double):Boolean;
var
    LnMax : Double;
    NrmB : Double;
    NrmX : Double;
    I : AlglibInteger;
    Alpha : Complex;
    Beta : Complex;
    VR : Double;
    CX : Complex;
    Tmp : TReal1DArray;
begin
    Assert(N>0, 'RMatrixTRSafeSolve: incorrect N!');
    Assert((Trans=0) or (Trans=1), 'RMatrixTRSafeSolve: incorrect Trans!');
    Result := True;
    LnMax := Ln(MaxRealNumber);
    
    //
    // Quick return if possible
    //
    if N<=0 then
    begin
        Exit;
    end;
    
    //
    // Load norms: right part and X
    //
    NrmB := 0;
    I:=0;
    while I<=N-1 do
    begin
        NrmB := Max(NrmB, AbsReal(X[I]));
        Inc(I);
    end;
    NrmX := 0;
    
    //
    // Solve
    //
    SetLength(Tmp, N);
    Result := True;
    if IsUpper and (Trans=0) then
    begin
        
        //
        // U*x = b
        //
        I:=N-1;
        while I>=0 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_Complex(A[I,I]*SA);
            end;
            if I<N-1 then
            begin
                APVMove(@Tmp[0], I+1, N-1, @A[I][0], I+1, N-1, SA);
                VR := APVDotProduct(@Tmp[0], I+1, N-1, @X[0], I+1, N-1);
                Beta := C_Complex(X[I]-VR);
            end
            else
            begin
                Beta := C_Complex(X[I]);
            end;
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, CX);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := CX.X;
            Dec(I);
        end;
        Exit;
    end;
    if  not IsUpper and (Trans=0) then
    begin
        
        //
        // L*x = b
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_Complex(A[I,I]*SA);
            end;
            if I>0 then
            begin
                APVMove(@Tmp[0], 0, I-1, @A[I][0], 0, I-1, SA);
                VR := APVDotProduct(@Tmp[0], 0, I-1, @X[0], 0, I-1);
                Beta := C_Complex(X[I]-VR);
            end
            else
            begin
                Beta := C_Complex(X[I]);
            end;
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, CX);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := CX.X;
            Inc(I);
        end;
        Exit;
    end;
    if IsUpper and (Trans=1) then
    begin
        
        //
        // U^T*x = b
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_Complex(A[I,I]*SA);
            end;
            Beta := C_Complex(X[I]);
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, CX);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := CX.X;
            
            //
            // update the rest of right part
            //
            if I<N-1 then
            begin
                VR := CX.X;
                APVMove(@Tmp[0], I+1, N-1, @A[I][0], I+1, N-1, SA);
                APVSub(@X[0], I+1, N-1, @Tmp[0], I+1, N-1, VR);
            end;
            Inc(I);
        end;
        Exit;
    end;
    if  not IsUpper and (Trans=1) then
    begin
        
        //
        // L^T*x = b
        //
        I:=N-1;
        while I>=0 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_Complex(A[I,I]*SA);
            end;
            Beta := C_Complex(X[I]);
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, CX);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := CX.X;
            
            //
            // update the rest of right part
            //
            if I>0 then
            begin
                VR := CX.X;
                APVMove(@Tmp[0], 0, I-1, @A[I][0], 0, I-1, SA);
                APVSub(@X[0], 0, I-1, @Tmp[0], 0, I-1, VR);
            end;
            Dec(I);
        end;
        Exit;
    end;
    Result := False;
end;


(*************************************************************************
Internal subroutine for safe solution of

    SA*op(A)=b
    
where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
either identity transform, transposition or Hermitian transposition, SA is
a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.

This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
near-degenerate matrices are handled correctly (False is returned) as long
as MaxGrowth is significantly less than MaxRealNumber/norm(b).

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************)
function CMatrixScaledTRSafeSolve(const A : TComplex2DArray;
     SA : Double;
     N : AlglibInteger;
     var X : TComplex1DArray;
     IsUpper : Boolean;
     Trans : AlglibInteger;
     IsUnit : Boolean;
     MaxGrowth : Double):Boolean;
var
    LnMax : Double;
    NrmB : Double;
    NrmX : Double;
    I : AlglibInteger;
    Alpha : Complex;
    Beta : Complex;
    VC : Complex;
    Tmp : TComplex1DArray;
    i_ : AlglibInteger;
begin
    Assert(N>0, 'CMatrixTRSafeSolve: incorrect N!');
    Assert((Trans=0) or (Trans=1) or (Trans=2), 'CMatrixTRSafeSolve: incorrect Trans!');
    Result := True;
    LnMax := Ln(MaxRealNumber);
    
    //
    // Quick return if possible
    //
    if N<=0 then
    begin
        Exit;
    end;
    
    //
    // Load norms: right part and X
    //
    NrmB := 0;
    I:=0;
    while I<=N-1 do
    begin
        NrmB := Max(NrmB, AbsComplex(X[I]));
        Inc(I);
    end;
    NrmX := 0;
    
    //
    // Solve
    //
    SetLength(Tmp, N);
    Result := True;
    if IsUpper and (Trans=0) then
    begin
        
        //
        // U*x = b
        //
        I:=N-1;
        while I>=0 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_MulR(A[I,I],SA);
            end;
            if I<N-1 then
            begin
                for i_ := I+1 to N-1 do
                begin
                    Tmp[i_] := C_MulR(A[I,i_],SA);
                end;
                VC := C_Complex(0.0);
                for i_ := I+1 to N-1 do
                begin
                    VC := C_Add(VC,C_Mul(Tmp[i_],X[i_]));
                end;
                Beta := C_Sub(X[I],VC);
            end
            else
            begin
                Beta := X[I];
            end;
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, VC);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := VC;
            Dec(I);
        end;
        Exit;
    end;
    if  not IsUpper and (Trans=0) then
    begin
        
        //
        // L*x = b
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_MulR(A[I,I],SA);
            end;
            if I>0 then
            begin
                for i_ := 0 to I-1 do
                begin
                    Tmp[i_] := C_MulR(A[I,i_],SA);
                end;
                VC := C_Complex(0.0);
                for i_ := 0 to I-1 do
                begin
                    VC := C_Add(VC,C_Mul(Tmp[i_],X[i_]));
                end;
                Beta := C_Sub(X[I],VC);
            end
            else
            begin
                Beta := X[I];
            end;
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, VC);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := VC;
            Inc(I);
        end;
        Exit;
    end;
    if IsUpper and (Trans=1) then
    begin
        
        //
        // U^T*x = b
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_MulR(A[I,I],SA);
            end;
            Beta := X[I];
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, VC);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := VC;
            
            //
            // update the rest of right part
            //
            if I<N-1 then
            begin
                for i_ := I+1 to N-1 do
                begin
                    Tmp[i_] := C_MulR(A[I,i_],SA);
                end;
                for i_ := I+1 to N-1 do
                begin
                    X[i_] := C_Sub(X[i_], C_Mul(VC, Tmp[i_]));
                end;
            end;
            Inc(I);
        end;
        Exit;
    end;
    if  not IsUpper and (Trans=1) then
    begin
        
        //
        // L^T*x = b
        //
        I:=N-1;
        while I>=0 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_MulR(A[I,I],SA);
            end;
            Beta := X[I];
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, VC);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := VC;
            
            //
            // update the rest of right part
            //
            if I>0 then
            begin
                for i_ := 0 to I-1 do
                begin
                    Tmp[i_] := C_MulR(A[I,i_],SA);
                end;
                for i_ := 0 to I-1 do
                begin
                    X[i_] := C_Sub(X[i_], C_Mul(VC, Tmp[i_]));
                end;
            end;
            Dec(I);
        end;
        Exit;
    end;
    if IsUpper and (Trans=2) then
    begin
        
        //
        // U^H*x = b
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_MulR(Conj(A[I,I]),SA);
            end;
            Beta := X[I];
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, VC);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := VC;
            
            //
            // update the rest of right part
            //
            if I<N-1 then
            begin
                for i_ := I+1 to N-1 do
                begin
                    Tmp[i_] := C_MulR(Conj(A[I,i_]),SA);
                end;
                for i_ := I+1 to N-1 do
                begin
                    X[i_] := C_Sub(X[i_], C_Mul(VC, Tmp[i_]));
                end;
            end;
            Inc(I);
        end;
        Exit;
    end;
    if  not IsUpper and (Trans=2) then
    begin
        
        //
        // L^T*x = b
        //
        I:=N-1;
        while I>=0 do
        begin
            
            //
            // Task is reduced to alpha*x[i] = beta
            //
            if IsUnit then
            begin
                Alpha := C_Complex(SA);
            end
            else
            begin
                Alpha := C_MulR(Conj(A[I,I]),SA);
            end;
            Beta := X[I];
            
            //
            // solve alpha*x[i] = beta
            //
            Result := CBasicSolveAndUpdate(Alpha, Beta, LnMax, NrmB, MaxGrowth, NrmX, VC);
            if  not Result then
            begin
                Exit;
            end;
            X[I] := VC;
            
            //
            // update the rest of right part
            //
            if I>0 then
            begin
                for i_ := 0 to I-1 do
                begin
                    Tmp[i_] := C_MulR(Conj(A[I,i_]),SA);
                end;
                for i_ := 0 to I-1 do
                begin
                    X[i_] := C_Sub(X[i_], C_Mul(VC, Tmp[i_]));
                end;
            end;
            Dec(I);
        end;
        Exit;
    end;
    Result := False;
end;


(*************************************************************************
complex basic solver-updater for reduced linear system

    alpha*x[i] = beta

solves this equation and updates it in overlfow-safe manner (keeping track
of relative growth of solution).

Parameters:
    Alpha   -   alpha
    Beta    -   beta
    LnMax   -   precomputed Ln(MaxRealNumber)
    BNorm   -   inf-norm of b (right part of original system)
    MaxGrowth-  maximum growth of norm(x) relative to norm(b)
    XNorm   -   inf-norm of other components of X (which are already processed)
                it is updated by CBasicSolveAndUpdate.
    X       -   solution

  -- ALGLIB routine --
     26.01.2009
     Bochkanov Sergey
*************************************************************************)
function CBasicSolveAndUpdate(Alpha : Complex;
     Beta : Complex;
     LnMax : Double;
     BNorm : Double;
     MaxGrowth : Double;
     var XNorm : Double;
     var X : Complex):Boolean;
var
    V : Double;
begin
    Result := False;
    if C_EqualR(Alpha,0) then
    begin
        Exit;
    end;
    if C_NotEqualR(Beta,0) then
    begin
        
        //
        // alpha*x[i]=beta
        //
        V := Ln(AbsComplex(Beta))-Ln(AbsComplex(Alpha));
        if AP_FP_Greater(V,LnMax) then
        begin
            Exit;
        end;
        X := C_Div(Beta,Alpha);
    end
    else
    begin
        
        //
        // alpha*x[i]=0
        //
        X := C_Complex(0);
    end;
    
    //
    // update NrmX, test growth limit
    //
    XNorm := Max(XNorm, AbsComplex(X));
    if AP_FP_Greater(XNorm,MaxGrowth*BNorm) then
    begin
        Exit;
    end;
    Result := True;
end;


end.