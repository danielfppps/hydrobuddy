unit hb_commercialnutrient;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  ComCtrls, StdCtrls, Menus, ExtCtrls, Buttons, hb_comparison;

type

  { TForm5 }

  TForm5 = class(TForm)
    Button1: TBitBtn;
    Button2: TButton;
    Button3: TButton;
    ComboBox1: TComboBox;
    ComboBox2: TComboBox;
    ComboBox3: TComboBox;
    Edit1: TEdit;
    Edit10: TEdit;
    Edit11: TEdit;
    Edit12: TEdit;
    Edit13: TEdit;
    Edit14: TEdit;
    Edit18: TEdit;
    Edit15: TEdit;
    Edit16: TEdit;
    Edit17: TEdit;
    Edit19: TEdit;
    Edit2: TEdit;
    Edit3: TEdit;
    Edit4: TEdit;
    Edit5: TEdit;
    Edit6: TEdit;
    Edit7: TEdit;
    Edit8: TEdit;
    Edit9: TEdit;
    Label1: TLabel;
    Label10: TLabel;
    Label11: TLabel;
    Label12: TLabel;
    Label13: TLabel;
    Label14: TLabel;
    Label29: TLabel;
    Label15: TLabel;
    Label16: TLabel;
    Label19: TLabel;
    Label2: TLabel;
    Label20: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label5: TLabel;
    Label6: TLabel;
    Label7: TLabel;
    Label8: TLabel;
    Label9: TLabel;
    Panel1: TPanel ;
    CheckBox1: TCheckBox;
    Panel2: TPanel;
    RadioButton1: TRadioButton;
    RadioButton2: TRadioButton;
    RadioButton3: TRadioButton;
    RadioButton4: TRadioButton;
    RadioButton5: TRadioButton;
    RadioButton6: TRadioButton;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
    procedure Button3Click(Sender: TObject);
    procedure RadioButton1Change(Sender: TObject);
    procedure RadioButton2Change(Sender: TObject);
    procedure RadioButton3Change(Sender: TObject);
    procedure RadioButton5Change(Sender: TObject);
    procedure RadioButton6Change(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form5: TForm5; 

implementation

uses HB_Main ;

{ TForm5 }


procedure TForm5.Button1Click(Sender: TObject);
var
i : integer ;
j : integer ;
varnames : array of string ;
result : array of double ;
test : double ;
Volume : double ;
currentValP: integer;
currentValK: Integer;
currentValSi: Integer;
begin

currentValP := ComboBox1.ItemIndex;
currentValK := ComboBox2.ItemIndex;
currentValSi := ComboBox3.ItemIndex;

SetLength (varnames, 16) ;
SetLength (result, 16) ;

if  RadioButton3.Checked then
Volume := 1      ;

if RadioButton4.Checked then
Volume := 3.78541178 ;

for i := 1 to 16 do
 begin

 // load all element names (this time we don't need to discriminate as
 // we simply calculate for everyone
 varnames[ i - 1 ] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption ;

 end;


    for i := 1 to 16 do

    begin

    test := StrtoFloat(((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)) ;

    if RadioButton2.Checked and HB_Main.Form1.RadioButton8.Checked then
    result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit18.Text) * 1000 )/Volume, 2) ;

     if RadioButton2.Checked and HB_Main.Form1.RadioButton9.Checked then
    result[i - 1] := Form1.round2(((1/0.0352739619)*0.01 * test * StrtoFloat(Edit18.Text) * 1000 )/Volume, 2) ;

    if (RadioButton1.Checked) and (RadioButton5.Checked) then
    result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit18.Text) * StrtoFloat(Edit17.Text) * 1000 )/Volume, 2) ;

    if (RadioButton1.Checked) and (RadioButton6.Checked)then
    result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit18.Text) * 1000 )/Volume, 2) ;

    end;

    HB_Main.Form1.cleanresults ;
    HB_Main.Form1.RadioButton10.Checked := true ;


    // finally copy values to edit boxes



           if CheckBox1.Checked = false then

           begin

          HB_Main.Form1.Edit1.Text := FloattoStr(result[0]) ;

          if currentValP = 0 then
          HB_Main.Form1.Edit3.Text := FloattoStr(result[2]) ;

          if currentValP = 1 then
          HB_Main.Form1.Edit3.Text := FloattoStr(Form1.round2(result[2]*0.4364, 2))  ;


          if currentValK = 0 then
          HB_Main.Form1.Edit2.Text := FloattoStr(result[1]) ;

          if currentValK = 1 then
          HB_Main.Form1.Edit2.Text := FloattoStr(Form1.round2(result[1]*0.8301, 2)) ;


          if currentValSi = 0 then
          HB_Main.Form1.Edit13.Text := FloattoStr(result[13]) ;

          if currentValSi = 1 then
          HB_Main.Form1.Edit13.Text := FloattoStr(Form1.round2(result[13]*0.4684, 2)) ;


          HB_Main.Form1.Edit4.Text := FloattoStr(result[3]) ;
          HB_Main.Form1.Edit5.Text := FloattoStr(result[4]) ;
          HB_Main.Form1.Edit6.Text := FloattoStr(result[5]) ;
          HB_Main.Form1.Edit7.Text := FloattoStr(result[6]) ;
          HB_Main.Form1.Edit8.Text := FloattoStr(result[7]) ;
          HB_Main.Form1.Edit9.Text := FloattoStr(result[8]) ;
          HB_Main.Form1.Edit10.Text := FloattoStr(result[9]) ;
          HB_Main.Form1.Edit11.Text := FloattoStr(result[10]) ;
          HB_Main.Form1.Edit12.Text := FloattoStr(result[11]) ;
          HB_Main.Form1.Edit13.Text := FloattoStr(result[12]) ;
          HB_Main.Form1.Edit15.Text := FloattoStr(result[14]) ;
          HB_Main.Form1.Edit16.Text := FloattoStr(result[15]) ;

          end ;




          if CheckBox1.Checked then

          begin

          HB_Main.Form1.Edit1.Text := FloattoStr(result[0]+ StrtoFloat(HB_Main.Form1.Edit1.Text)) ;

          if currentValP = 0 then
          HB_Main.Form1.Edit3.Text := FloattoStr(result[2]+ StrtoFloat(HB_Main.Form1.Edit3.Text)) ;

          if currentValP = 1 then
          HB_Main.Form1.Edit3.Text := FloattoStr(Form1.round2(result[2]*0.4364, 2))  ;


          if currentValK = 0 then
          HB_Main.Form1.Edit2.Text := FloattoStr(result[1]+ StrtoFloat(HB_Main.Form1.Edit2.Text)) ;  ;

          if currentValK = 1 then
          HB_Main.Form1.Edit2.Text := FloattoStr(Form1.round2(result[1]*0.8301+ StrtoFloat(HB_Main.Form1.Edit2.Text), 2)) ;


          if currentValSi = 0 then
          HB_Main.Form1.Edit13.Text := FloattoStr(result[13]+ StrtoFloat(HB_Main.Form1.Edit13.Text)) ;

          if currentValSi = 1 then
          HB_Main.Form1.Edit13.Text := FloattoStr(Form1.round2(result[13]*0.4684+ StrtoFloat(HB_Main.Form1.Edit13.Text), 2)) ;



          HB_Main.Form1.Edit4.Text := FloattoStr(result[3]+ StrtoFloat(HB_Main.Form1.Edit4.Text)) ;
          HB_Main.Form1.Edit5.Text := FloattoStr(result[4]+ StrtoFloat(HB_Main.Form1.Edit5.Text)) ;
          HB_Main.Form1.Edit6.Text := FloattoStr(result[5]+ StrtoFloat(HB_Main.Form1.Edit6.Text)) ;
          HB_Main.Form1.Edit7.Text := FloattoStr(result[6]+ StrtoFloat(HB_Main.Form1.Edit7.Text)) ;
          HB_Main.Form1.Edit8.Text := FloattoStr(result[7]+ StrtoFloat(HB_Main.Form1.Edit8.Text)) ;
          HB_Main.Form1.Edit9.Text := FloattoStr(result[8]+ StrtoFloat(HB_Main.Form1.Edit9.Text)) ;
          HB_Main.Form1.Edit10.Text := FloattoStr(result[9]+ StrtoFloat(HB_Main.Form1.Edit10.Text)) ;
          HB_Main.Form1.Edit11.Text := FloattoStr(result[10]+ StrtoFloat(HB_Main.Form1.Edit11.Text));
          HB_Main.Form1.Edit12.Text := FloattoStr(result[11]+ StrtoFloat(HB_Main.Form1.Edit12.Text));
          HB_Main.Form1.Edit13.Text := FloattoStr(result[12]+ StrtoFloat(HB_Main.Form1.Edit13.Text));
          HB_Main.Form1.Edit14.Text := FloattoStr(result[13]+ StrtoFloat(HB_Main.Form1.Edit14.Text));
          HB_Main.Form1.Edit15.Text := FloattoStr(result[14]+ StrtoFloat(HB_Main.Form1.Edit15.Text));
          HB_Main.Form1.Edit16.Text :=FloattoStr(result[15]+ StrtoFloat(HB_Main.Form1.Edit16.Text));

          end ;


          Form5.Visible := false ;


          end ;

procedure TForm5.Button2Click(Sender: TObject);
var
i : integer ;
j : integer ;
varnames : array of string ;
result : array of double ;
test : double ;
Volume : double ;
currentValP: integer;
currentValK: Integer;
currentValSi: Integer;
colCount: integer;
begin

currentValP := ComboBox1.ItemIndex;
currentValK := ComboBox2.ItemIndex;
currentValSi := ComboBox3.ItemIndex;

SetLength (varnames, 16) ;
SetLength (result, 16) ;

if  RadioButton3.Checked then
Volume := 1      ;

if RadioButton4.Checked then
Volume := 3.78541178 ;

for i := 1 to 16 do
 begin

 // load all element names (this time we don't need to discriminate as
 // we simply calculate for everyone
 varnames[ i - 1 ] := (FindComponent('Label' + IntToStr(i)) as TLabel).Caption ;

 end;


    for i := 1 to 16 do

    begin

    test := StrtoFloat(((FindComponent('Edit' + IntToStr(i)) as TEdit).Text)) ;

    if RadioButton2.Checked and HB_Main.Form1.RadioButton8.Checked then
    result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit18.Text) * 1000 )/Volume, 2) ;

     if RadioButton2.Checked and HB_Main.Form1.RadioButton9.Checked then
    result[i - 1] := Form1.round2(((1/0.0352739619)*0.01 * test * StrtoFloat(Edit18.Text) * 1000 )/Volume, 2) ;

    if RadioButton1.Checked then
    result[i - 1] := Form1.round2((0.01 * test * StrtoFloat(Edit18.Text) * StrtoFloat(Edit17.Text) * 1000 )/Volume, 2) ;

    end;


    // finally copy values for comparison

    hb_comparison.Form15.StringGrid1.ColCount:= hb_comparison.Form15.StringGrid1.ColCount + 1 ;

    colCount := hb_comparison.Form15.StringGrid1.ColCount ;

    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 0] := Edit19.Text  ;
    hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 1] := FloattoStr(result[0]) ;

          if currentValP = 0 then
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 2] := FloattoStr(result[2]) ;

          if currentValP = 1 then
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 2] := FloattoStr(Form1.round2(result[2]*0.4364, 2))  ;


          if currentValK = 0 then
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 3] := FloattoStr(result[1]) ;

          if currentValK = 1 then
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 3] := FloattoStr(Form1.round2(result[1]*0.8301, 2)) ;


          if currentValSi = 0 then
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 14] := FloattoStr(result[13]) ;

          if currentValSi = 1 then
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 14] := FloattoStr(Form1.round2(result[13]*0.4684, 2)) ;


          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 4] := FloattoStr(result[3]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 5] := FloattoStr(result[4]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 6] := FloattoStr(result[5]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 7] := FloattoStr(result[6]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 8] := FloattoStr(result[7]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 9] := FloattoStr(result[8]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 10] := FloattoStr(result[9]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 11] := FloattoStr(result[10]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 12] := FloattoStr(result[11]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 13] := FloattoStr(result[12]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 15] := FloattoStr(result[14]) ;
          hb_comparison.Form15.StringGrid1.Cells[ColCount-1, 16] := FloattoStr(result[15]) ;;

          ShowMessage('Product final ppm values added to comparison chart') ;

end;

procedure TForm5.Button3Click(Sender: TObject);
begin
  hb_comparison.Form15.Visible := true ;
end;


procedure TForm5.RadioButton1Change(Sender: TObject);
begin

Label20.Caption :=  'Volume of addition (mL) ' ;
Edit17.Enabled := true ;

if RadioButton1.Checked then
  Panel2.Enabled := true ;

if RadioButton6.Checked then
  Edit17.Enabled := false ;

end;

procedure TForm5.RadioButton2Change(Sender: TObject);
begin

if HB_Main.Form1.RadioButton8.Checked then
Label20.Caption :=  'Mass of addition (g)' ;

if HB_Main.Form1.RadioButton9.Checked then
Label20.Caption :=  'Mass of addition (oz)' ;

if RadioButton2.Checked then
Panel2.Enabled := false ;

Edit17.Enabled := false ;

end;

procedure TForm5.RadioButton3Change(Sender: TObject);
begin

end;

procedure TForm5.RadioButton5Change(Sender: TObject);
begin

  if RadioButton5.Checked then
  Edit17.Enabled := true ;

end;

procedure TForm5.RadioButton6Change(Sender: TObject);
begin

if RadioButton6.Checked then
Edit17.Enabled := false ;

end;

initialization
  {$I hb_commercialnutrient.lrs}

end.

