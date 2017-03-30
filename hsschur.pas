{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

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
unit hsschur;
interface
uses Math, Sysutils, Ap, blas, reflections, rotations;

function UpperHessenbergSchurDecomposition(var H : TReal2DArray;
     N : AlglibInteger;
     var S : TReal2DArray):Boolean;
procedure InternalSchurDecomposition(var H : TReal2DArray;
     N : AlglibInteger;
     TNeeded : AlglibInteger;
     ZNeeded : AlglibInteger;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     var Z : TReal2DArray;
     var INFO : AlglibInteger);

implementation

procedure InternalAuxSchur(WANTT : Boolean;
     WANTZ : Boolean;
     N : AlglibInteger;
     ILO : AlglibInteger;
     IHI : AlglibInteger;
     var H : TReal2DArray;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     ILOZ : AlglibInteger;
     IHIZ : AlglibInteger;
     var Z : TReal2DArray;
     var WORK : TReal1DArray;
     var WORKV3 : TReal1DArray;
     var WORKC1 : TReal1DArray;
     var WORKS1 : TReal1DArray;
     var INFO : AlglibInteger);forward;
procedure Aux2X2Schur(var A : Double;
     var B : Double;
     var C : Double;
     var D : Double;
     var RT1R : Double;
     var RT1I : Double;
     var RT2R : Double;
     var RT2I : Double;
     var CS : Double;
     var SN : Double);forward;
function ExtSchurSign(a : Double; b : Double):Double;forward;
function ExtSchurSignToOne(b : Double):AlglibInteger;forward;


(*************************************************************************
Subroutine performing  the  Schur  decomposition  of  a  matrix  in  upper
Hessenberg form using the QR algorithm with multiple shifts.

The  source matrix  H  is  represented as  S'*H*S = T, where H - matrix in
upper Hessenberg form,  S - orthogonal matrix (Schur vectors),   T - upper
quasi-triangular matrix (with blocks of sizes  1x1  and  2x2  on  the main
diagonal).

Input parameters:
    H   -   matrix to be decomposed.
            Array whose indexes range within [1..N, 1..N].
    N   -   size of H, N>=0.


Output parameters:
    H   –   contains the matrix T.
            Array whose indexes range within [1..N, 1..N].
            All elements below the blocks on the main diagonal are equal
            to 0.
    S   -   contains Schur vectors.
            Array whose indexes range within [1..N, 1..N].

Note 1:
    The block structure of matrix T could be easily recognized: since  all
    the elements  below  the blocks are zeros, the elements a[i+1,i] which
    are equal to 0 show the block border.

Note 2:
    the algorithm  performance  depends  on  the  value  of  the  internal
    parameter NS of InternalSchurDecomposition  subroutine  which  defines
    the number of shifts in the QR algorithm (analog of  the  block  width
    in block matrix algorithms in linear algebra). If you require  maximum
    performance  on  your  machine,  it  is  recommended  to  adjust  this
    parameter manually.

Result:
    True, if the algorithm has converged and the parameters H and S contain
        the result.
    False, if the algorithm has not converged.

Algorithm implemented on the basis of subroutine DHSEQR (LAPACK 3.0 library).
*************************************************************************)
function UpperHessenbergSchurDecomposition(var H : TReal2DArray;
     N : AlglibInteger;
     var S : TReal2DArray):Boolean;
var
    WI : TReal1DArray;
    WR : TReal1DArray;
    INFO : AlglibInteger;
begin
    InternalSchurDecomposition(H, N, 1, 2, WR, WI, S, INFO);
    Result := INFO=0;
end;


procedure InternalSchurDecomposition(var H : TReal2DArray;
     N : AlglibInteger;
     TNeeded : AlglibInteger;
     ZNeeded : AlglibInteger;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     var Z : TReal2DArray;
     var INFO : AlglibInteger);
var
    WORK : TReal1DArray;
    I : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IERR : AlglibInteger;
    II : AlglibInteger;
    ITEMP : AlglibInteger;
    ITN : AlglibInteger;
    ITS : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    L : AlglibInteger;
    MAXB : AlglibInteger;
    NR : AlglibInteger;
    NS : AlglibInteger;
    NV : AlglibInteger;
    ABSW : Double;
    OVFL : Double;
    SMLNUM : Double;
    TAU : Double;
    TEMP : Double;
    TST1 : Double;
    ULP : Double;
    UNFL : Double;
    S : TReal2DArray;
    V : TReal1DArray;
    VV : TReal1DArray;
    WORKC1 : TReal1DArray;
    WORKS1 : TReal1DArray;
    WORKV3 : TReal1DArray;
    TmpWR : TReal1DArray;
    TmpWI : TReal1DArray;
    INITZ : Boolean;
    WANTT : Boolean;
    WANTZ : Boolean;
    CNST : Double;
    FailFlag : Boolean;
    P1 : AlglibInteger;
    P2 : AlglibInteger;
    VT : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Set the order of the multi-shift QR algorithm to be used.
    // If you want to tune algorithm, change this values
    //
    NS := 12;
    MAXB := 50;
    
    //
    // Now 2 < NS <= MAXB < NH.
    //
    MAXB := Max(3, MAXB);
    NS := Min(MAXB, NS);
    
    //
    // Initialize
    //
    CNST := Double(1.5);
    SetLength(WORK, Max(N, 1)+1);
    SetLength(S, NS+1, NS+1);
    SetLength(V, NS+1+1);
    SetLength(VV, NS+1+1);
    SetLength(WR, Max(N, 1)+1);
    SetLength(WI, Max(N, 1)+1);
    SetLength(WORKC1, 1+1);
    SetLength(WORKS1, 1+1);
    SetLength(WORKV3, 3+1);
    SetLength(TmpWR, Max(N, 1)+1);
    SetLength(TmpWI, Max(N, 1)+1);
    Assert(N>=0, 'InternalSchurDecomposition: incorrect N!');
    Assert((TNeeded=0) or (TNeeded=1), 'InternalSchurDecomposition: incorrect TNeeded!');
    Assert((ZNeeded=0) or (ZNeeded=1) or (ZNeeded=2), 'InternalSchurDecomposition: incorrect ZNeeded!');
    WANTT := TNeeded=1;
    INITZ := ZNeeded=2;
    WANTZ := ZNeeded<>0;
    INFO := 0;
    
    //
    // Initialize Z, if necessary
    //
    if INITZ then
    begin
        SetLength(Z, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                if I=J then
                begin
                    Z[I,J] := 1;
                end
                else
                begin
                    Z[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    if N=1 then
    begin
        WR[1] := H[1,1];
        WI[1] := 0;
        Exit;
    end;
    
    //
    // Set rows and columns 1 to N to zero below the first
    // subdiagonal.
    //
    J:=1;
    while J<=N-2 do
    begin
        I:=J+2;
        while I<=N do
        begin
            H[I,J] := 0;
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // Test if N is sufficiently small
    //
    if (NS<=2) or (NS>N) or (MAXB>=N) then
    begin
        
        //
        // Use the standard double-shift algorithm
        //
        InternalAuxSchur(WANTT, WANTZ, N, 1, N, H, WR, WI, 1, N, Z, WORK, WORKV3, WORKC1, WORKS1, INFO);
        
        //
        // fill entries under diagonal blocks of T with zeros
        //
        if WANTT then
        begin
            J := 1;
            while J<=N do
            begin
                if AP_FP_Eq(WI[J],0) then
                begin
                    I:=J+1;
                    while I<=N do
                    begin
                        H[I,J] := 0;
                        Inc(I);
                    end;
                    J := J+1;
                end
                else
                begin
                    I:=J+2;
                    while I<=N do
                    begin
                        H[I,J] := 0;
                        H[I,J+1] := 0;
                        Inc(I);
                    end;
                    J := J+2;
                end;
            end;
        end;
        Exit;
    end;
    UNFL := MinRealNumber;
    OVFL := 1/UNFL;
    ULP := 2*MachineEpsilon;
    SMLNUM := UNFL*(N/ULP);
    
    //
    // I1 and I2 are the indices of the first row and last column of H
    // to which transformations must be applied. If eigenvalues only are
    // being computed, I1 and I2 are set inside the main loop.
    //
    I1 := 1;
    I2 := N;
    
    //
    // ITN is the total number of multiple-shift QR iterations allowed.
    //
    ITN := 30*N;
    
    //
    // The main loop begins here. I is the loop index and decreases from
    // IHI to ILO in steps of at most MAXB. Each iteration of the loop
    // works with the active submatrix in rows and columns L to I.
    // Eigenvalues I+1 to IHI have already converged. Either L = ILO or
    // H(L,L-1) is negligible so that the matrix splits.
    //
    I := N;
    while True do
    begin
        L := 1;
        if I<1 then
        begin
            
            //
            // fill entries under diagonal blocks of T with zeros
            //
            if WANTT then
            begin
                J := 1;
                while J<=N do
                begin
                    if AP_FP_Eq(WI[J],0) then
                    begin
                        I:=J+1;
                        while I<=N do
                        begin
                            H[I,J] := 0;
                            Inc(I);
                        end;
                        J := J+1;
                    end
                    else
                    begin
                        I:=J+2;
                        while I<=N do
                        begin
                            H[I,J] := 0;
                            H[I,J+1] := 0;
                            Inc(I);
                        end;
                        J := J+2;
                    end;
                end;
            end;
            
            //
            // Exit
            //
            Exit;
        end;
        
        //
        // Perform multiple-shift QR iterations on rows and columns ILO to I
        // until a submatrix of order at most MAXB splits off at the bottom
        // because a subdiagonal element has become negligible.
        //
        FailFlag := True;
        ITS:=0;
        while ITS<=ITN do
        begin
            
            //
            // Look for a single small subdiagonal element.
            //
            K:=I;
            while K>=L+1 do
            begin
                TST1 := ABSReal(H[K-1,K-1])+ABSReal(H[K,K]);
                if AP_FP_Eq(TST1,0) then
                begin
                    TST1 := UpperHessenberg1Norm(H, L, I, L, I, WORK);
                end;
                if AP_FP_Less_Eq(ABSReal(H[K,K-1]),Max(ULP*TST1, SMLNUM)) then
                begin
                    Break;
                end;
                Dec(K);
            end;
            L := K;
            if L>1 then
            begin
                
                //
                // H(L,L-1) is negligible.
                //
                H[L,L-1] := 0;
            end;
            
            //
            // Exit from loop if a submatrix of order <= MAXB has split off.
            //
            if L>=I-MAXB+1 then
            begin
                FailFlag := False;
                Break;
            end;
            
            //
            // Now the active submatrix is in rows and columns L to I. If
            // eigenvalues only are being computed, only the active submatrix
            // need be transformed.
            //
            if (ITS=20) or (ITS=30) then
            begin
                
                //
                // Exceptional shifts.
                //
                II:=I-NS+1;
                while II<=I do
                begin
                    WR[II] := CNST*(ABSReal(H[II,II-1])+ABSReal(H[II,II]));
                    WI[II] := 0;
                    Inc(II);
                end;
            end
            else
            begin
                
                //
                // Use eigenvalues of trailing submatrix of order NS as shifts.
                //
                CopyMatrix(H, I-NS+1, I, I-NS+1, I, S, 1, NS, 1, NS);
                InternalAuxSchur(False, False, NS, 1, NS, S, TmpWR, TmpWI, 1, NS, Z, WORK, WORKV3, WORKC1, WORKS1, IERR);
                P1:=1;
                while P1<=NS do
                begin
                    WR[I-NS+P1] := TmpWR[P1];
                    WI[I-NS+P1] := TmpWI[P1];
                    Inc(P1);
                end;
                if IERR>0 then
                begin
                    
                    //
                    // If DLAHQR failed to compute all NS eigenvalues, use the
                    // unconverged diagonal elements as the remaining shifts.
                    //
                    II:=1;
                    while II<=IERR do
                    begin
                        WR[I-NS+II] := S[II,II];
                        WI[I-NS+II] := 0;
                        Inc(II);
                    end;
                end;
            end;
            
            //
            // Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
            // where G is the Hessenberg submatrix H(L:I,L:I) and w is
            // the vector of shifts (stored in WR and WI). The result is
            // stored in the local array V.
            //
            V[1] := 1;
            II:=2;
            while II<=NS+1 do
            begin
                V[II] := 0;
                Inc(II);
            end;
            NV := 1;
            J:=I-NS+1;
            while J<=I do
            begin
                if AP_FP_Greater_Eq(WI[J],0) then
                begin
                    if AP_FP_Eq(WI[J],0) then
                    begin
                        
                        //
                        // real shift
                        //
                        P1 := NV+1;
                        APVMove(@VV[0], 1, P1, @V[0], 1, P1);
                        MatrixVectorMultiply(H, L, L+NV, L, L+NV-1, False, VV, 1, NV, Double(1.0), V, 1, NV+1, -WR[J]);
                        NV := NV+1;
                    end
                    else
                    begin
                        if AP_FP_Greater(WI[J],0) then
                        begin
                            
                            //
                            // complex conjugate pair of shifts
                            //
                            P1 := NV+1;
                            APVMove(@VV[0], 1, P1, @V[0], 1, P1);
                            MatrixVectorMultiply(H, L, L+NV, L, L+NV-1, False, V, 1, NV, Double(1.0), VV, 1, NV+1, -2*WR[J]);
                            ITEMP := VectorIdxAbsMax(VV, 1, NV+1);
                            TEMP := 1/Max(ABSReal(VV[ITEMP]), SMLNUM);
                            P1 := NV+1;
                            APVMul(@VV[0], 1, P1, TEMP);
                            ABSW := Pythag2(WR[J], WI[J]);
                            TEMP := TEMP*ABSW*ABSW;
                            MatrixVectorMultiply(H, L, L+NV+1, L, L+NV, False, VV, 1, NV+1, Double(1.0), V, 1, NV+2, TEMP);
                            NV := NV+2;
                        end;
                    end;
                    
                    //
                    // Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
                    // reset it to the unit vector.
                    //
                    ITEMP := VectorIdxAbsMax(V, 1, NV);
                    TEMP := ABSReal(V[ITEMP]);
                    if AP_FP_Eq(TEMP,0) then
                    begin
                        V[1] := 1;
                        II:=2;
                        while II<=NV do
                        begin
                            V[II] := 0;
                            Inc(II);
                        end;
                    end
                    else
                    begin
                        TEMP := Max(TEMP, SMLNUM);
                        VT := 1/TEMP;
                        APVMul(@V[0], 1, NV, VT);
                    end;
                end;
                Inc(J);
            end;
            
            //
            // Multiple-shift QR step
            //
            K:=L;
            while K<=I-1 do
            begin
                
                //
                // The first iteration of this loop determines a reflection G
                // from the vector V and applies it from left and right to H,
                // thus creating a nonzero bulge below the subdiagonal.
                //
                // Each subsequent iteration determines a reflection G to
                // restore the Hessenberg form in the (K-1)th column, and thus
                // chases the bulge one step toward the bottom of the active
                // submatrix. NR is the order of G.
                //
                NR := Min(NS+1, I-K+1);
                if K>L then
                begin
                    P1 := K-1;
                    P2 := K+NR-1;
                    i1_ := (K) - (1);
                    for i_ := 1 to NR do
                    begin
                        V[i_] := H[i_+i1_,P1];
                    end;
                end;
                GenerateReflection(V, NR, Tau);
                if K>L then
                begin
                    H[K,K-1] := V[1];
                    II:=K+1;
                    while II<=I do
                    begin
                        H[II,K-1] := 0;
                        Inc(II);
                    end;
                end;
                V[1] := 1;
                
                //
                // Apply G from the left to transform the rows of the matrix in
                // columns K to I2.
                //
                ApplyReflectionFromTheLeft(H, Tau, V, K, K+NR-1, K, I2, WORK);
                
                //
                // Apply G from the right to transform the columns of the
                // matrix in rows I1 to min(K+NR,I).
                //
                ApplyReflectionFromTheRight(H, Tau, V, I1, Min(K+NR, I), K, K+NR-1, WORK);
                if WANTZ then
                begin
                    
                    //
                    // Accumulate transformations in the matrix Z
                    //
                    ApplyReflectionFromTheRight(Z, Tau, V, 1, N, K, K+NR-1, WORK);
                end;
                Inc(K);
            end;
            Inc(ITS);
        end;
        
        //
        // Failure to converge in remaining number of iterations
        //
        if FailFlag then
        begin
            INFO := I;
            Exit;
        end;
        
        //
        // A submatrix of order <= MAXB in rows and columns L to I has split
        // off. Use the double-shift QR algorithm to handle it.
        //
        InternalAuxSchur(WANTT, WANTZ, N, L, I, H, WR, WI, 1, N, Z, WORK, WORKV3, WORKC1, WORKS1, INFO);
        if INFO>0 then
        begin
            Exit;
        end;
        
        //
        // Decrement number of remaining iterations, and return to start of
        // the main loop with a new value of I.
        //
        ITN := ITN-ITS;
        I := L-1;
    end;
end;


procedure InternalAuxSchur(WANTT : Boolean;
     WANTZ : Boolean;
     N : AlglibInteger;
     ILO : AlglibInteger;
     IHI : AlglibInteger;
     var H : TReal2DArray;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     ILOZ : AlglibInteger;
     IHIZ : AlglibInteger;
     var Z : TReal2DArray;
     var WORK : TReal1DArray;
     var WORKV3 : TReal1DArray;
     var WORKC1 : TReal1DArray;
     var WORKS1 : TReal1DArray;
     var INFO : AlglibInteger);
var
    I : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    ITN : AlglibInteger;
    ITS : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    L : AlglibInteger;
    M : AlglibInteger;
    NH : AlglibInteger;
    NR : AlglibInteger;
    NZ : AlglibInteger;
    AVE : Double;
    CS : Double;
    DISC : Double;
    H00 : Double;
    H10 : Double;
    H11 : Double;
    H12 : Double;
    H21 : Double;
    H22 : Double;
    H33 : Double;
    H33S : Double;
    H43H34 : Double;
    H44 : Double;
    H44S : Double;
    OVFL : Double;
    S : Double;
    SMLNUM : Double;
    SN : Double;
    SUM : Double;
    T1 : Double;
    T2 : Double;
    T3 : Double;
    TST1 : Double;
    UNFL : Double;
    V1 : Double;
    V2 : Double;
    V3 : Double;
    FailFlag : Boolean;
    DAT1 : Double;
    DAT2 : Double;
    P1 : AlglibInteger;
    HIM1IM1 : Double;
    HIM1I : Double;
    HIIM1 : Double;
    HII : Double;
    WRIM1 : Double;
    WRI : Double;
    WIIM1 : Double;
    WII : Double;
    Ulp : Double;
begin
    INFO := 0;
    DAT1 := Double(0.75);
    DAT2 := -Double(0.4375);
    Ulp := MachineEpsilon;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    if ILO=IHI then
    begin
        WR[ILO] := H[ILO,ILO];
        WI[ILO] := 0;
        Exit;
    end;
    NH := IHI-ILO+1;
    NZ := IHIZ-ILOZ+1;
    
    //
    // Set machine-dependent constants for the stopping criterion.
    // If norm(H) <= sqrt(OVFL), overflow should not occur.
    //
    UNFL := MinRealNumber;
    OVFL := 1/UNFL;
    SMLNUM := UNFL*(NH/ULP);
    
    //
    // I1 and I2 are the indices of the first row and last column of H
    // to which transformations must be applied. If eigenvalues only are
    // being computed, I1 and I2 are set inside the main loop.
    //
    I1 := 1;
    I2 := N;
    
    //
    // ITN is the total number of QR iterations allowed.
    //
    ITN := 30*NH;
    
    //
    // The main loop begins here. I is the loop index and decreases from
    // IHI to ILO in steps of 1 or 2. Each iteration of the loop works
    // with the active submatrix in rows and columns L to I.
    // Eigenvalues I+1 to IHI have already converged. Either L = ILO or
    // H(L,L-1) is negligible so that the matrix splits.
    //
    I := IHI;
    while True do
    begin
        L := ILO;
        if I<ILO then
        begin
            Exit;
        end;
        
        //
        // Perform QR iterations on rows and columns ILO to I until a
        // submatrix of order 1 or 2 splits off at the bottom because a
        // subdiagonal element has become negligible.
        //
        FailFlag := True;
        ITS:=0;
        while ITS<=ITN do
        begin
            
            //
            // Look for a single small subdiagonal element.
            //
            K:=I;
            while K>=L+1 do
            begin
                TST1 := ABSReal(H[K-1,K-1])+ABSReal(H[K,K]);
                if AP_FP_Eq(TST1,0) then
                begin
                    TST1 := UpperHessenberg1Norm(H, L, I, L, I, WORK);
                end;
                if AP_FP_Less_Eq(ABSReal(H[K,K-1]),Max(ULP*TST1, SMLNUM)) then
                begin
                    Break;
                end;
                Dec(K);
            end;
            L := K;
            if L>ILO then
            begin
                
                //
                // H(L,L-1) is negligible
                //
                H[L,L-1] := 0;
            end;
            
            //
            // Exit from loop if a submatrix of order 1 or 2 has split off.
            //
            if L>=I-1 then
            begin
                FailFlag := False;
                Break;
            end;
            
            //
            // Now the active submatrix is in rows and columns L to I. If
            // eigenvalues only are being computed, only the active submatrix
            // need be transformed.
            //
            if (ITS=10) or (ITS=20) then
            begin
                
                //
                // Exceptional shift.
                //
                S := ABSReal(H[I,I-1])+ABSReal(H[I-1,I-2]);
                H44 := DAT1*S+H[I,I];
                H33 := H44;
                H43H34 := DAT2*S*S;
            end
            else
            begin
                
                //
                // Prepare to use Francis' double shift
                // (i.e. 2nd degree generalized Rayleigh quotient)
                //
                H44 := H[I,I];
                H33 := H[I-1,I-1];
                H43H34 := H[I,I-1]*H[I-1,I];
                S := H[I-1,I-2]*H[I-1,I-2];
                DISC := (H33-H44)*Double(0.5);
                DISC := DISC*DISC+H43H34;
                if AP_FP_Greater(DISC,0) then
                begin
                    
                    //
                    // Real roots: use Wilkinson's shift twice
                    //
                    DISC := SQRT(DISC);
                    AVE := Double(0.5)*(H33+H44);
                    if AP_FP_Greater(ABSReal(H33)-ABSReal(H44),0) then
                    begin
                        H33 := H33*H44-H43H34;
                        H44 := H33/(ExtSchurSign(DISC, AVE)+AVE);
                    end
                    else
                    begin
                        H44 := ExtSchurSign(DISC, AVE)+AVE;
                    end;
                    H33 := H44;
                    H43H34 := 0;
                end;
            end;
            
            //
            // Look for two consecutive small subdiagonal elements.
            //
            M:=I-2;
            while M>=L do
            begin
                
                //
                // Determine the effect of starting the double-shift QR
                // iteration at row M, and see if this would make H(M,M-1)
                // negligible.
                //
                H11 := H[M,M];
                H22 := H[M+1,M+1];
                H21 := H[M+1,M];
                H12 := H[M,M+1];
                H44S := H44-H11;
                H33S := H33-H11;
                V1 := (H33S*H44S-H43H34)/H21+H12;
                V2 := H22-H11-H33S-H44S;
                V3 := H[M+2,M+1];
                S := AbsReal(V1)+AbsReal(V2)+AbsReal(V3);
                V1 := V1/S;
                V2 := V2/S;
                V3 := V3/S;
                WORKV3[1] := V1;
                WORKV3[2] := V2;
                WORKV3[3] := V3;
                if M=L then
                begin
                    Break;
                end;
                H00 := H[M-1,M-1];
                H10 := H[M,M-1];
                TST1 := AbsReal(V1)*(AbsReal(H00)+AbsReal(H11)+AbsReal(H22));
                if AP_FP_Less_Eq(AbsReal(H10)*(ABSReal(V2)+ABSReal(V3)),ULP*TST1) then
                begin
                    Break;
                end;
                Dec(M);
            end;
            
            //
            // Double-shift QR step
            //
            K:=M;
            while K<=I-1 do
            begin
                
                //
                // The first iteration of this loop determines a reflection G
                // from the vector V and applies it from left and right to H,
                // thus creating a nonzero bulge below the subdiagonal.
                //
                // Each subsequent iteration determines a reflection G to
                // restore the Hessenberg form in the (K-1)th column, and thus
                // chases the bulge one step toward the bottom of the active
                // submatrix. NR is the order of G.
                //
                NR := Min(3, I-K+1);
                if K>M then
                begin
                    P1:=1;
                    while P1<=NR do
                    begin
                        WORKV3[P1] := H[K+P1-1,K-1];
                        Inc(P1);
                    end;
                end;
                GenerateReflection(WORKV3, NR, T1);
                if K>M then
                begin
                    H[K,K-1] := WORKV3[1];
                    H[K+1,K-1] := 0;
                    if K<I-1 then
                    begin
                        H[K+2,K-1] := 0;
                    end;
                end
                else
                begin
                    if M>L then
                    begin
                        H[K,K-1] := -H[K,K-1];
                    end;
                end;
                V2 := WORKV3[2];
                T2 := T1*V2;
                if NR=3 then
                begin
                    V3 := WORKV3[3];
                    T3 := T1*V3;
                    
                    //
                    // Apply G from the left to transform the rows of the matrix
                    // in columns K to I2.
                    //
                    J:=K;
                    while J<=I2 do
                    begin
                        SUM := H[K,J]+V2*H[K+1,J]+V3*H[K+2,J];
                        H[K,J] := H[K,J]-SUM*T1;
                        H[K+1,J] := H[K+1,J]-SUM*T2;
                        H[K+2,J] := H[K+2,J]-SUM*T3;
                        Inc(J);
                    end;
                    
                    //
                    // Apply G from the right to transform the columns of the
                    // matrix in rows I1 to min(K+3,I).
                    //
                    J:=I1;
                    while J<=Min(K+3, I) do
                    begin
                        SUM := H[J,K]+V2*H[J,K+1]+V3*H[J,K+2];
                        H[J,K] := H[J,K]-SUM*T1;
                        H[J,K+1] := H[J,K+1]-SUM*T2;
                        H[J,K+2] := H[J,K+2]-SUM*T3;
                        Inc(J);
                    end;
                    if WANTZ then
                    begin
                        
                        //
                        // Accumulate transformations in the matrix Z
                        //
                        J:=ILOZ;
                        while J<=IHIZ do
                        begin
                            SUM := Z[J,K]+V2*Z[J,K+1]+V3*Z[J,K+2];
                            Z[J,K] := Z[J,K]-SUM*T1;
                            Z[J,K+1] := Z[J,K+1]-SUM*T2;
                            Z[J,K+2] := Z[J,K+2]-SUM*T3;
                            Inc(J);
                        end;
                    end;
                end
                else
                begin
                    if NR=2 then
                    begin
                        
                        //
                        // Apply G from the left to transform the rows of the matrix
                        // in columns K to I2.
                        //
                        J:=K;
                        while J<=I2 do
                        begin
                            SUM := H[K,J]+V2*H[K+1,J];
                            H[K,J] := H[K,J]-SUM*T1;
                            H[K+1,J] := H[K+1,J]-SUM*T2;
                            Inc(J);
                        end;
                        
                        //
                        // Apply G from the right to transform the columns of the
                        // matrix in rows I1 to min(K+3,I).
                        //
                        J:=I1;
                        while J<=I do
                        begin
                            SUM := H[J,K]+V2*H[J,K+1];
                            H[J,K] := H[J,K]-SUM*T1;
                            H[J,K+1] := H[J,K+1]-SUM*T2;
                            Inc(J);
                        end;
                        if WANTZ then
                        begin
                            
                            //
                            // Accumulate transformations in the matrix Z
                            //
                            J:=ILOZ;
                            while J<=IHIZ do
                            begin
                                SUM := Z[J,K]+V2*Z[J,K+1];
                                Z[J,K] := Z[J,K]-SUM*T1;
                                Z[J,K+1] := Z[J,K+1]-SUM*T2;
                                Inc(J);
                            end;
                        end;
                    end;
                end;
                Inc(K);
            end;
            Inc(ITS);
        end;
        if FailFlag then
        begin
            
            //
            // Failure to converge in remaining number of iterations
            //
            INFO := I;
            Exit;
        end;
        if L=I then
        begin
            
            //
            // H(I,I-1) is negligible: one eigenvalue has converged.
            //
            WR[I] := H[I,I];
            WI[I] := 0;
        end
        else
        begin
            if L=I-1 then
            begin
                
                //
                // H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
                //
                //        Transform the 2-by-2 submatrix to standard Schur form,
                //        and compute and store the eigenvalues.
                //
                HIM1IM1 := H[I-1,I-1];
                HIM1I := H[I-1,I];
                HIIM1 := H[I,I-1];
                HII := H[I,I];
                Aux2X2Schur(HIM1IM1, HIM1I, HIIM1, HII, WRIM1, WIIM1, WRI, WII, CS, SN);
                WR[I-1] := WRIM1;
                WI[I-1] := WIIM1;
                WR[I] := WRI;
                WI[I] := WII;
                H[I-1,I-1] := HIM1IM1;
                H[I-1,I] := HIM1I;
                H[I,I-1] := HIIM1;
                H[I,I] := HII;
                if WANTT then
                begin
                    
                    //
                    // Apply the transformation to the rest of H.
                    //
                    if I2>I then
                    begin
                        WORKC1[1] := CS;
                        WORKS1[1] := SN;
                        ApplyRotationsFromTheLeft(True, I-1, I, I+1, I2, WORKC1, WORKS1, H, WORK);
                    end;
                    WORKC1[1] := CS;
                    WORKS1[1] := SN;
                    ApplyRotationsFromTheRight(True, I1, I-2, I-1, I, WORKC1, WORKS1, H, WORK);
                end;
                if WANTZ then
                begin
                    
                    //
                    // Apply the transformation to Z.
                    //
                    WORKC1[1] := CS;
                    WORKS1[1] := SN;
                    ApplyRotationsFromTheRight(True, ILOZ, ILOZ+NZ-1, I-1, I, WORKC1, WORKS1, Z, WORK);
                end;
            end;
        end;
        
        //
        // Decrement number of remaining iterations, and return to start of
        // the main loop with new value of I.
        //
        ITN := ITN-ITS;
        I := L-1;
    end;
end;


procedure Aux2X2Schur(var A : Double;
     var B : Double;
     var C : Double;
     var D : Double;
     var RT1R : Double;
     var RT1I : Double;
     var RT2R : Double;
     var RT2I : Double;
     var CS : Double;
     var SN : Double);
var
    MULTPL : Double;
    AA : Double;
    BB : Double;
    BCMAX : Double;
    BCMIS : Double;
    CC : Double;
    CS1 : Double;
    DD : Double;
    EPS : Double;
    P : Double;
    SAB : Double;
    SAC : Double;
    SCL : Double;
    SIGMA : Double;
    SN1 : Double;
    TAU : Double;
    TEMP : Double;
    Z : Double;
begin
    MULTPL := Double(4.0);
    EPS := MachineEpsilon;
    if AP_FP_Eq(C,0) then
    begin
        CS := 1;
        SN := 0;
    end
    else
    begin
        if AP_FP_Eq(B,0) then
        begin
            
            //
            // Swap rows and columns
            //
            CS := 0;
            SN := 1;
            TEMP := D;
            D := A;
            A := TEMP;
            B := -C;
            C := 0;
        end
        else
        begin
            if AP_FP_Eq(A-D,0) and (ExtSchurSignToOne(B)<>ExtSchurSignToOne(C)) then
            begin
                CS := 1;
                SN := 0;
            end
            else
            begin
                TEMP := A-D;
                P := Double(0.5)*TEMP;
                BCMAX := Max(ABSReal(B), ABSReal(C));
                BCMIS := Min(ABSReal(B), ABSReal(C))*ExtSchurSignToOne(B)*ExtSchurSignToOne(C);
                SCL := Max(ABSReal(P), BCMAX);
                Z := P/SCL*P+BCMAX/SCL*BCMIS;
                
                //
                // If Z is of the order of the machine accuracy, postpone the
                // decision on the nature of eigenvalues
                //
                if AP_FP_Greater_Eq(Z,MULTPL*EPS) then
                begin
                    
                    //
                    // Real eigenvalues. Compute A and D.
                    //
                    Z := P+ExtSchurSign(SQRT(SCL)*SQRT(Z), P);
                    A := D+Z;
                    D := D-BCMAX/Z*BCMIS;
                    
                    //
                    // Compute B and the rotation matrix
                    //
                    TAU := Pythag2(C, Z);
                    CS := Z/TAU;
                    SN := C/TAU;
                    B := B-C;
                    C := 0;
                end
                else
                begin
                    
                    //
                    // Complex eigenvalues, or real (almost) equal eigenvalues.
                    // Make diagonal elements equal.
                    //
                    SIGMA := B+C;
                    TAU := Pythag2(SIGMA, TEMP);
                    CS := SQRT(Double(0.5)*(1+ABSReal(SIGMA)/TAU));
                    SN := -P/(TAU*CS)*ExtSchurSign(1, SIGMA);
                    
                    //
                    // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
                    //         [ CC  DD ]   [ C  D ] [ SN  CS ]
                    //
                    AA := A*CS+B*SN;
                    BB := -A*SN+B*CS;
                    CC := C*CS+D*SN;
                    DD := -C*SN+D*CS;
                    
                    //
                    // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
                    //         [ C  D ]   [-SN  CS ] [ CC  DD ]
                    //
                    A := AA*CS+CC*SN;
                    B := BB*CS+DD*SN;
                    C := -AA*SN+CC*CS;
                    D := -BB*SN+DD*CS;
                    TEMP := Double(0.5)*(A+D);
                    A := TEMP;
                    D := TEMP;
                    if AP_FP_Neq(C,0) then
                    begin
                        if AP_FP_Neq(B,0) then
                        begin
                            if ExtSchurSignToOne(B)=ExtSchurSignToOne(C) then
                            begin
                                
                                //
                                // Real eigenvalues: reduce to upper triangular form
                                //
                                SAB := SQRT(ABSReal(B));
                                SAC := SQRT(ABSReal(C));
                                P := ExtSchurSign(SAB*SAC, C);
                                TAU := 1/SQRT(ABSReal(B+C));
                                A := TEMP+P;
                                D := TEMP-P;
                                B := B-C;
                                C := 0;
                                CS1 := SAB*TAU;
                                SN1 := SAC*TAU;
                                TEMP := CS*CS1-SN*SN1;
                                SN := CS*SN1+SN*CS1;
                                CS := TEMP;
                            end;
                        end
                        else
                        begin
                            B := -C;
                            C := 0;
                            TEMP := CS;
                            CS := -SN;
                            SN := TEMP;
                        end;
                    end;
                end;
            end;
        end;
    end;
    
    //
    // Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
    //
    RT1R := A;
    RT2R := D;
    if AP_FP_Eq(C,0) then
    begin
        RT1I := 0;
        RT2I := 0;
    end
    else
    begin
        RT1I := SQRT(ABSReal(B))*SQRT(ABSReal(C));
        RT2I := -RT1I;
    end;
end;


function ExtSchurSign(a : Double; b : Double):Double;
begin
    if AP_FP_Greater_Eq(b,0) then
    begin
        Result := AbsReal(a);
    end
    else
    begin
        Result := -AbsReal(a);
    end;
end;


function ExtSchurSignToOne(b : Double):AlglibInteger;
begin
    if AP_FP_Greater_Eq(b,0) then
    begin
        Result := 1;
    end
    else
    begin
        Result := -1;
    end;
end;


end.