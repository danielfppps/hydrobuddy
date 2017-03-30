{$MODESWITCH RESULT+}
{$GOTO ON}
unit xblas;
interface
uses Math, Sysutils, Ap;

procedure XDot(const A : TReal1DArray;
     const B : TReal1DArray;
     N : AlglibInteger;
     var Temp : TReal1DArray;
     var R : Double;
     var RErr : Double);
procedure XCDot(const A : TComplex1DArray;
     const B : TComplex1DArray;
     N : AlglibInteger;
     var Temp : TReal1DArray;
     var R : Complex;
     var RErr : Double);

implementation

procedure XSum(var W : TReal1DArray;
     MX : Double;
     N : AlglibInteger;
     var R : Double;
     var RErr : Double);forward;
function XFastPow(R : Double; N : AlglibInteger):Double;forward;


(*************************************************************************
More precise dot-product. Absolute error of  subroutine  result  is  about
1 ulp of max(MX,V), where:
    MX = max( |a[i]*b[i]| )
    V  = |(a,b)|

INPUT PARAMETERS
    A       -   array[0..N-1], vector 1
    B       -   array[0..N-1], vector 2
    N       -   vectors length, N<2^29.
    Temp    -   array[0..N-1], pre-allocated temporary storage

OUTPUT PARAMETERS
    R       -   (A,B)
    RErr    -   estimate of error. This estimate accounts for both  errors
                during  calculation  of  (A,B)  and  errors  introduced by
                rounding of A and B to fit in double (about 1 ulp).

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure XDot(const A : TReal1DArray;
     const B : TReal1DArray;
     N : AlglibInteger;
     var Temp : TReal1DArray;
     var R : Double;
     var RErr : Double);
var
    I : AlglibInteger;
    MX : Double;
    V : Double;
begin
    
    //
    // special cases:
    // * N=0
    //
    if N=0 then
    begin
        R := 0;
        RErr := 0;
        Exit;
    end;
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := A[I]*B[I];
        Temp[I] := V;
        MX := Max(MX, AbsReal(V));
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        R := 0;
        RErr := 0;
        Exit;
    end;
    XSum(Temp, MX, N, R, RErr);
end;


(*************************************************************************
More precise complex dot-product. Absolute error of  subroutine  result is
about 1 ulp of max(MX,V), where:
    MX = max( |a[i]*b[i]| )
    V  = |(a,b)|

INPUT PARAMETERS
    A       -   array[0..N-1], vector 1
    B       -   array[0..N-1], vector 2
    N       -   vectors length, N<2^29.
    Temp    -   array[0..2*N-1], pre-allocated temporary storage

OUTPUT PARAMETERS
    R       -   (A,B)
    RErr    -   estimate of error. This estimate accounts for both  errors
                during  calculation  of  (A,B)  and  errors  introduced by
                rounding of A and B to fit in double (about 1 ulp).

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure XCDot(const A : TComplex1DArray;
     const B : TComplex1DArray;
     N : AlglibInteger;
     var Temp : TReal1DArray;
     var R : Complex;
     var RErr : Double);
var
    I : AlglibInteger;
    MX : Double;
    V : Double;
    RErrX : Double;
    RErrY : Double;
begin
    
    //
    // special cases:
    // * N=0
    //
    if N=0 then
    begin
        R := C_Complex(0);
        RErr := 0;
        Exit;
    end;
    
    //
    // calculate real part
    //
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := A[I].X*B[I].X;
        Temp[2*I+0] := V;
        MX := Max(MX, AbsReal(V));
        V := -A[I].Y*B[I].Y;
        Temp[2*I+1] := V;
        MX := Max(MX, AbsReal(V));
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        R.X := 0;
        RErrX := 0;
    end
    else
    begin
        XSum(Temp, MX, 2*N, R.X, RErrX);
    end;
    
    //
    // calculate imaginary part
    //
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := A[I].X*B[I].Y;
        Temp[2*I+0] := V;
        MX := Max(MX, AbsReal(V));
        V := A[I].Y*B[I].X;
        Temp[2*I+1] := V;
        MX := Max(MX, AbsReal(V));
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        R.Y := 0;
        RErrY := 0;
    end
    else
    begin
        XSum(Temp, MX, 2*N, R.Y, RErrY);
    end;
    
    //
    // total error
    //
    if AP_FP_Eq(RErrX,0) and AP_FP_Eq(RErrY,0) then
    begin
        RErr := 0;
    end
    else
    begin
        RErr := Max(RErrX, RErrY)*Sqrt(1+AP_Sqr(Min(RErrX, RErrY)/Max(RErrX, RErrY)));
    end;
end;


(*************************************************************************
Internal subroutine for extra-precise calculation of SUM(w[i]).

INPUT PARAMETERS:
    W   -   array[0..N-1], values to be added
            W is modified during calculations.
    MX  -   max(W[i])
    N   -   array size
    
OUTPUT PARAMETERS:
    R   -   SUM(w[i])
    RErr-   error estimate for R

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure XSum(var W : TReal1DArray;
     MX : Double;
     N : AlglibInteger;
     var R : Double;
     var RErr : Double);
var
    I : AlglibInteger;
    K : AlglibInteger;
    KS : AlglibInteger;
    V : Double;
    S : Double;
    LN2 : Double;
    Chunk : Double;
    InvChunk : Double;
    AllZeros : Boolean;
begin
    
    //
    // special cases:
    // * N=0
    // * N is too large to use integer arithmetics
    //
    if N=0 then
    begin
        R := 0;
        RErr := 0;
        Exit;
    end;
    if AP_FP_Eq(MX,0) then
    begin
        R := 0;
        RErr := 0;
        Exit;
    end;
    Assert(N<536870912, 'XDot: N is too large!');
    
    //
    // Prepare
    //
    LN2 := Ln(2);
    RErr := MX*MachineEpsilon;
    
    //
    // 1. find S such that 0.5<=S*MX<1
    // 2. multiply W by S, so task is normalized in some sense
    // 3. S:=1/S so we can obtain original vector multiplying by S
    //
    K := Round(Ln(MX)/LN2);
    S := XFastPow(2, -K);
    while AP_FP_Greater_Eq(S*MX,1) do
    begin
        S := Double(0.5)*S;
    end;
    while AP_FP_Less(S*MX,Double(0.5)) do
    begin
        S := 2*S;
    end;
    APVMul(@W[0], 0, N-1, S);
    S := 1/S;
    
    //
    // find Chunk=2^M such that N*Chunk<2^29
    //
    // we have chosen upper limit (2^29) with enough space left
    // to tolerate possible problems with rounding and N's close
    // to the limit, so we don't want to be very strict here.
    //
    K := Trunc(Ln(AP_Double(536870912)/N)/LN2);
    Chunk := XFastPow(2, K);
    if AP_FP_Less(Chunk,2) then
    begin
        Chunk := 2;
    end;
    InvChunk := 1/Chunk;
    
    //
    // calculate result
    //
    R := 0;
    APVMul(@W[0], 0, N-1, Chunk);
    while True do
    begin
        S := S*InvChunk;
        AllZeros := True;
        KS := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := W[I];
            K := Trunc(V);
            if AP_FP_Neq(V,K) then
            begin
                AllZeros := False;
            end;
            W[I] := Chunk*(V-K);
            KS := KS+K;
            Inc(I);
        end;
        R := R+S*KS;
        V := AbsReal(R);
        if AllZeros or AP_FP_Eq(S*N+MX,MX) then
        begin
            Break;
        end;
    end;
    
    //
    // correct error
    //
    RErr := Max(RErr, AbsReal(R)*MachineEpsilon);
end;


(*************************************************************************
Fast Pow

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
function XFastPow(R : Double; N : AlglibInteger):Double;
begin
    if N>0 then
    begin
        if N mod 2=0 then
        begin
            Result := AP_Sqr(XFastPow(R, N div 2));
        end
        else
        begin
            Result := R*XFastPow(R, N-1);
        end;
        Exit;
    end;
    if N=0 then
    begin
        Result := 1;
    end;
    if N<0 then
    begin
        Result := XFastPow(1/R, -N);
    end;
end;


end.