{$MODESWITCH RESULT+}
{$GOTO ON}
unit nearunityunit;
interface
uses Math, Sysutils, Ap;

function Log1P(x : Double):Double;
function ExpM1(x : Double):Double;
function CosM1(x : Double):Double;

implementation

function Log1P(x : Double):Double;
var
    z : Double;
    LP : Double;
    LQ : Double;
begin
    z := Double(1.0)+x;
    if AP_FP_Less(z,Double(0.70710678118654752440)) or AP_FP_Greater(z,Double(1.41421356237309504880)) then
    begin
        Result := Ln(z);
        Exit;
    end;
    z := x*x;
    LP := Double(4.5270000862445199635215E-5);
    LP := LP*x+Double(4.9854102823193375972212E-1);
    LP := LP*x+Double(6.5787325942061044846969E0);
    LP := LP*x+Double(2.9911919328553073277375E1);
    LP := LP*x+Double(6.0949667980987787057556E1);
    LP := LP*x+Double(5.7112963590585538103336E1);
    LP := LP*x+Double(2.0039553499201281259648E1);
    LQ := Double(1.0000000000000000000000E0);
    LQ := LQ*x+Double(1.5062909083469192043167E1);
    LQ := LQ*x+Double(8.3047565967967209469434E1);
    LQ := LQ*x+Double(2.2176239823732856465394E2);
    LQ := LQ*x+Double(3.0909872225312059774938E2);
    LQ := LQ*x+Double(2.1642788614495947685003E2);
    LQ := LQ*x+Double(6.0118660497603843919306E1);
    z := -Double(0.5)*z+x*(z*LP/LQ);
    Result := x+z;
end;


function ExpM1(x : Double):Double;
var
    r : Double;
    xx : Double;
    EP : Double;
    EQ : Double;
begin
    if AP_FP_Less(x,-Double(0.5)) or AP_FP_Greater(x,Double(0.5)) then
    begin
        Result := exp(x)-Double(1.0);
        Exit;
    end;
    xx := x*x;
    EP := Double(1.2617719307481059087798E-4);
    EP := EP*xx+Double(3.0299440770744196129956E-2);
    EP := EP*xx+Double(9.9999999999999999991025E-1);
    EQ := Double(3.0019850513866445504159E-6);
    EQ := EQ*xx+Double(2.5244834034968410419224E-3);
    EQ := EQ*xx+Double(2.2726554820815502876593E-1);
    EQ := EQ*xx+Double(2.0000000000000000000897E0);
    r := x*EP;
    r := r/(EQ-r);
    Result := r+r;
end;


function CosM1(x : Double):Double;
var
    xx : Double;
    C : Double;
begin
    if AP_FP_Less(x,-Double(0.25)*Pi) or AP_FP_Greater(x,Double(0.25)*Pi) then
    begin
        Result := Cos(x)-1;
        Exit;
    end;
    xx := x*x;
    C := Double(4.7377507964246204691685E-14);
    C := C*xx-Double(1.1470284843425359765671E-11);
    C := C*xx+Double(2.0876754287081521758361E-9);
    C := C*xx-Double(2.7557319214999787979814E-7);
    C := C*xx+Double(2.4801587301570552304991E-5);
    C := C*xx-Double(1.3888888888888872993737E-3);
    C := C*xx+Double(4.1666666666666666609054E-2);
    Result := -Double(0.5)*xx+xx*xx*C;
end;


end.