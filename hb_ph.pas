unit hb_ph;

{$mode objfpc}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls;

type

  { TForm13 }

  TForm13 = class(TForm)
    Button1: TButton;
    Edit17: TEdit;
    Edit18: TEdit;
    Edit19: TEdit;
    Label17: TLabel;
    Label18: TLabel;
    Label19: TLabel;
    Label20: TLabel;
    RadioButton1: TRadioButton;
    RadioButton2: TRadioButton;
    procedure Button1Click(Sender: TObject);
    procedure RadioButton1Change(Sender: TObject);
    procedure RadioButton2Change(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form13: TForm13; 

implementation

{ TForm13 }

procedure TForm13.Button1Click(Sender: TObject);
begin
  Form13.Visible := false ;
end;

procedure TForm13.RadioButton1Change(Sender: TObject);
  var
GH : double ;
KH : double ;
begin

GH := StrToFloat(Edit18.Text) ;
KH := StrToFloat(Edit19.Text) ;

if RadioButton1.Checked then

   begin

   Edit18.Text := FloatToStr(GH/17.86) ;
   Edit19.Text := FloatToStr(KH/17.86) ;

   end ;
end;

procedure TForm13.RadioButton2Change(Sender: TObject);
var
GH : double ;
KH : double ;
begin

GH := StrToFloat(Edit18.Text) ;
KH := StrToFloat(Edit19.Text) ;

if RadioButton2.Checked then

   begin

   Edit18.Text := FloatToStr(GH*17.86) ;
   Edit19.Text := FloatToStr(KH*17.86) ;

   end ;

end;

initialization
  {$I hb_ph.lrs}

end.

