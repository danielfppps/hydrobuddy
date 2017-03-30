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
unit pca;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd, descriptivestatistics;

procedure PCABuildBasis(const X : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var S2 : TReal1DArray;
     var V : TReal2DArray);

implementation

(*************************************************************************
Principal components analysis

Subroutine  builds  orthogonal  basis  where  first  axis  corresponds  to
direction with maximum variance, second axis maximizes variance in subspace
orthogonal to first axis and so on.

It should be noted that, unlike LDA, PCA does not use class labels.

INPUT PARAMETERS:
    X           -   dataset, array[0..NPoints-1,0..NVars-1].
                    matrix contains ONLY INDEPENDENT VARIABLES.
    NPoints     -   dataset size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1

ÂÛÕÎÄÍÛÅ ÏÀĞÀÌÅÒĞÛ:
    Info        -   return code:
                    * -4, if SVD subroutine haven't converged
                    * -1, if wrong parameters has been passed (NPoints<0,
                          NVars<1)
                    *  1, if task is solved
    S2          -   array[0..NVars-1]. variance values corresponding
                    to basis vectors.
    V           -   array[0..NVars-1,0..NVars-1]
                    matrix, whose columns store basis vectors.

  -- ALGLIB --
     Copyright 25.08.2008 by Bochkanov Sergey
*************************************************************************)
procedure PCABuildBasis(const X : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Info : AlglibInteger;
     var S2 : TReal1DArray;
     var V : TReal2DArray);
var
    A : TReal2DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    M : TReal1DArray;
    T : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    Mean : Double;
    Variance : Double;
    Skewness : Double;
    Kurtosis : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Check input data
    //
    if (NPoints<0) or (NVars<1) then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    
    //
    // Special case: NPoints=0
    //
    if NPoints=0 then
    begin
        SetLength(S2, NVars-1+1);
        SetLength(V, NVars-1+1, NVars-1+1);
        I:=0;
        while I<=NVars-1 do
        begin
            S2[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=NVars-1 do
        begin
            J:=0;
            while J<=NVars-1 do
            begin
                if I=J then
                begin
                    V[I,J] := 1;
                end
                else
                begin
                    V[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Calculate means
    //
    SetLength(M, NVars-1+1);
    SetLength(T, NPoints-1+1);
    J:=0;
    while J<=NVars-1 do
    begin
        for i_ := 0 to NPoints-1 do
        begin
            T[i_] := X[i_,J];
        end;
        CalculateMoments(T, NPoints, Mean, Variance, Skewness, Kurtosis);
        M[J] := Mean;
        Inc(J);
    end;
    
    //
    // Center, apply SVD, prepare output
    //
    SetLength(A, Max(NPoints, NVars)-1+1, NVars-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@A[I][0], 0, NVars-1, @X[I][0], 0, NVars-1);
        APVSub(@A[I][0], 0, NVars-1, @M[0], 0, NVars-1);
        Inc(I);
    end;
    I:=NPoints;
    while I<=NVars-1 do
    begin
        J:=0;
        while J<=NVars-1 do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    if  not RMatrixSVD(A, Max(NPoints, NVars), NVars, 0, 1, 2, S2, U, VT) then
    begin
        Info := -4;
        Exit;
    end;
    if NPoints<>1 then
    begin
        I:=0;
        while I<=NVars-1 do
        begin
            S2[I] := AP_Sqr(S2[I])/(NPoints-1);
            Inc(I);
        end;
    end;
    SetLength(V, NVars-1+1, NVars-1+1);
    CopyAndTranspose(VT, 0, NVars-1, 0, NVars-1, V, 0, NVars-1, 0, NVars-1);
end;


end.