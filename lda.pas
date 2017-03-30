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
unit lda;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, evd, hqrnd, matgen, trfac, trlinsolve, safesolve, rcond, matinv;

procedure FisherLDA(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var W : TReal1DArray);
procedure FisherLDAN(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var W : TReal2DArray);

implementation

(*************************************************************************
Multiclass Fisher LDA

Subroutine finds coefficients of linear combination which optimally separates
training set on classes.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   linear combination coefficients, array[0..NVars-1]

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure FisherLDA(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var W : TReal1DArray);
var
    W2 : TReal2DArray;
    i_ : AlglibInteger;
begin
    FisherLDAN(XY, NPoints, NVars, NClasses, Info, W2);
    if Info>0 then
    begin
        SetLength(W, NVars-1+1);
        for i_ := 0 to NVars-1 do
        begin
            W[i_] := W2[i_,0];
        end;
    end;
end;


(*************************************************************************
N-dimensional multiclass Fisher LDA

Subroutine finds coefficients of linear combinations which optimally separates
training set on classes. It returns N-dimensional basis whose vector are sorted
by quality of training set separation (in descending order).

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   basis, array[0..NVars-1,0..NVars-1]
                    columns of matrix stores basis vectors, sorted by
                    quality of training set separation (in descending order)

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************)
procedure FisherLDAN(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     NClasses : AlglibInteger;
     var Info : AlglibInteger;
     var W : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    V : Double;
    C : TInteger1DArray;
    Mu : TReal1DArray;
    MuC : TReal2DArray;
    NC : TInteger1DArray;
    SW : TReal2DArray;
    ST : TReal2DArray;
    Z : TReal2DArray;
    Z2 : TReal2DArray;
    TM : TReal2DArray;
    SBRoot : TReal2DArray;
    A : TReal2DArray;
    XYProj : TReal2DArray;
    WProj : TReal2DArray;
    TF : TReal1DArray;
    D : TReal1DArray;
    D2 : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
begin
    
    //
    // Test data
    //
    if (NPoints<0) or (NVars<1) or (NClasses<2) then
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
    // Special case: NPoints<=1
    // Degenerate task.
    //
    if NPoints<=1 then
    begin
        Info := 2;
        SetLength(W, NVars-1+1, NVars-1+1);
        I:=0;
        while I<=NVars-1 do
        begin
            J:=0;
            while J<=NVars-1 do
            begin
                if I=J then
                begin
                    W[I,J] := 1;
                end
                else
                begin
                    W[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Prepare temporaries
    //
    SetLength(TF, NVars-1+1);
    SetLength(WORK, Max(NVars, NPoints)+1);
    
    //
    // Convert class labels from reals to integers (just for convenience)
    //
    SetLength(C, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        C[I] := Round(XY[I,NVars]);
        Inc(I);
    end;
    
    //
    // Calculate class sizes and means
    //
    SetLength(Mu, NVars-1+1);
    SetLength(MuC, NClasses-1+1, NVars-1+1);
    SetLength(NC, NClasses-1+1);
    J:=0;
    while J<=NVars-1 do
    begin
        Mu[J] := 0;
        Inc(J);
    end;
    I:=0;
    while I<=NClasses-1 do
    begin
        NC[I] := 0;
        J:=0;
        while J<=NVars-1 do
        begin
            MuC[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVAdd(@Mu[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        APVAdd(@MuC[C[I]][0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        NC[C[I]] := NC[C[I]]+1;
        Inc(I);
    end;
    I:=0;
    while I<=NClasses-1 do
    begin
        V := AP_Double(1)/NC[I];
        APVMul(@MuC[I][0], 0, NVars-1, V);
        Inc(I);
    end;
    V := AP_Double(1)/NPoints;
    APVMul(@Mu[0], 0, NVars-1, V);
    
    //
    // Create ST matrix
    //
    SetLength(ST, NVars-1+1, NVars-1+1);
    I:=0;
    while I<=NVars-1 do
    begin
        J:=0;
        while J<=NVars-1 do
        begin
            ST[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    K:=0;
    while K<=NPoints-1 do
    begin
        APVMove(@TF[0], 0, NVars-1, @XY[K][0], 0, NVars-1);
        APVSub(@TF[0], 0, NVars-1, @Mu[0], 0, NVars-1);
        I:=0;
        while I<=NVars-1 do
        begin
            V := TF[I];
            APVAdd(@ST[I][0], 0, NVars-1, @TF[0], 0, NVars-1, V);
            Inc(I);
        end;
        Inc(K);
    end;
    
    //
    // Create SW matrix
    //
    SetLength(SW, NVars-1+1, NVars-1+1);
    I:=0;
    while I<=NVars-1 do
    begin
        J:=0;
        while J<=NVars-1 do
        begin
            SW[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    K:=0;
    while K<=NPoints-1 do
    begin
        APVMove(@TF[0], 0, NVars-1, @XY[K][0], 0, NVars-1);
        APVSub(@TF[0], 0, NVars-1, @MuC[C[K]][0], 0, NVars-1);
        I:=0;
        while I<=NVars-1 do
        begin
            V := TF[I];
            APVAdd(@SW[I][0], 0, NVars-1, @TF[0], 0, NVars-1, V);
            Inc(I);
        end;
        Inc(K);
    end;
    
    //
    // Maximize ratio J=(w'*ST*w)/(w'*SW*w).
    //
    // First, make transition from w to v such that w'*ST*w becomes v'*v:
    //    v  = root(ST)*w = R*w
    //    R  = root(D)*Z'
    //    w  = (root(ST)^-1)*v = RI*v
    //    RI = Z*inv(root(D))
    //    J  = (v'*v)/(v'*(RI'*SW*RI)*v)
    //    ST = Z*D*Z'
    //
    //    so we have
    //
    //    J = (v'*v) / (v'*(inv(root(D))*Z'*SW*Z*inv(root(D)))*v)  =
    //      = (v'*v) / (v'*A*v)
    //
    if  not SMatrixEVD(ST, NVars, 1, True, D, Z) then
    begin
        Info := -4;
        Exit;
    end;
    SetLength(W, NVars-1+1, NVars-1+1);
    if AP_FP_Less_Eq(D[NVars-1],0) or AP_FP_Less_Eq(D[0],1000*MachineEpsilon*D[NVars-1]) then
    begin
        
        //
        // Special case: D[NVars-1]<=0
        // Degenerate task (all variables takes the same value).
        //
        if AP_FP_Less_Eq(D[NVars-1],0) then
        begin
            Info := 2;
            I:=0;
            while I<=NVars-1 do
            begin
                J:=0;
                while J<=NVars-1 do
                begin
                    if I=J then
                    begin
                        W[I,J] := 1;
                    end
                    else
                    begin
                        W[I,J] := 0;
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
            Exit;
        end;
        
        //
        // Special case: degenerate ST matrix, multicollinearity found.
        // Since we know ST eigenvalues/vectors we can translate task to
        // non-degenerate form.
        //
        // Let WG is orthogonal basis of the non zero variance subspace
        // of the ST and let WZ is orthogonal basis of the zero variance
        // subspace.
        //
        // Projection on WG allows us to use LDA on reduced M-dimensional
        // subspace, N-M vectors of WZ allows us to update reduced LDA
        // factors to full N-dimensional subspace.
        //
        M := 0;
        K:=0;
        while K<=NVars-1 do
        begin
            if AP_FP_Less_Eq(D[K],1000*MachineEpsilon*D[NVars-1]) then
            begin
                M := K+1;
            end;
            Inc(K);
        end;
        Assert(M<>0, 'FisherLDAN: internal error #1');
        SetLength(XYProj, NPoints-1+1, NVars-M+1);
        MatrixMatrixMultiply(XY, 0, NPoints-1, 0, NVars-1, False, Z, 0, NVars-1, M, NVars-1, False, Double(1.0), XYProj, 0, NPoints-1, 0, NVars-M-1, Double(0.0), WORK);
        I:=0;
        while I<=NPoints-1 do
        begin
            XYProj[I,NVars-M] := XY[I,NVars];
            Inc(I);
        end;
        FisherLDAN(XYProj, NPoints, NVars-M, NClasses, Info, WProj);
        if Info<0 then
        begin
            Exit;
        end;
        MatrixMatrixMultiply(Z, 0, NVars-1, M, NVars-1, False, WProj, 0, NVars-M-1, 0, NVars-M-1, False, Double(1.0), W, 0, NVars-1, 0, NVars-M-1, Double(0.0), WORK);
        K:=NVars-M;
        while K<=NVars-1 do
        begin
            for i_ := 0 to NVars-1 do
            begin
                W[i_,K] := Z[i_,K-(NVars-M)];
            end;
            Inc(K);
        end;
        Info := 2;
    end
    else
    begin
        
        //
        // General case: no multicollinearity
        //
        SetLength(TM, NVars-1+1, NVars-1+1);
        SetLength(A, NVars-1+1, NVars-1+1);
        MatrixMatrixMultiply(SW, 0, NVars-1, 0, NVars-1, False, Z, 0, NVars-1, 0, NVars-1, False, Double(1.0), TM, 0, NVars-1, 0, NVars-1, Double(0.0), WORK);
        MatrixMatrixMultiply(Z, 0, NVars-1, 0, NVars-1, True, TM, 0, NVars-1, 0, NVars-1, False, Double(1.0), A, 0, NVars-1, 0, NVars-1, Double(0.0), WORK);
        I:=0;
        while I<=NVars-1 do
        begin
            J:=0;
            while J<=NVars-1 do
            begin
                A[I,J] := A[I,J]/Sqrt(D[I]*D[J]);
                Inc(J);
            end;
            Inc(I);
        end;
        if  not SMatrixEVD(A, NVars, 1, True, D2, Z2) then
        begin
            Info := -4;
            Exit;
        end;
        K:=0;
        while K<=NVars-1 do
        begin
            I:=0;
            while I<=NVars-1 do
            begin
                TF[I] := Z2[I,K]/Sqrt(D[I]);
                Inc(I);
            end;
            I:=0;
            while I<=NVars-1 do
            begin
                V := APVDotProduct(@Z[I][0], 0, NVars-1, @TF[0], 0, NVars-1);
                W[I,K] := V;
                Inc(I);
            end;
            Inc(K);
        end;
    end;
    
    //
    // Post-processing:
    // * normalization
    // * converting to non-negative form, if possible
    //
    K:=0;
    while K<=NVars-1 do
    begin
        V := 0.0;
        for i_ := 0 to NVars-1 do
        begin
            V := V + W[i_,K]*W[i_,K];
        end;
        V := 1/Sqrt(V);
        for i_ := 0 to NVars-1 do
        begin
            W[i_,K] := V*W[i_,K];
        end;
        V := 0;
        I:=0;
        while I<=NVars-1 do
        begin
            V := V+W[I,K];
            Inc(I);
        end;
        if AP_FP_Less(V,0) then
        begin
            for i_ := 0 to NVars-1 do
            begin
                W[i_,K] := -1*W[i_,K];
            end;
        end;
        Inc(K);
    end;
end;


end.