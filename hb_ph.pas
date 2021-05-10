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
    ComboBox1: TComboBox;
    Edit18: TEdit;
    Label1: TLabel;
    Label18: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form13: TForm13; 

implementation

uses hb_waterquality;

{ TForm13 }

procedure TForm13.Button1Click(Sender: TObject);
var
  total_alkalinity: double;
  ppm_contribution: double;
begin

  //see here https://ag.umass.edu/greenhouse-floriculture/fact-sheets/adjusting-alkalinity-with-acids

  total_alkalinity := StrToFloat(Edit18.Text);

  if ComboBox1.ItemIndex = 0 then
  begin
        ppm_contribution := total_alkalinity*0.7*(0.033814/1)*(25.6);
        Form6.Edit3.Text := FloatToStr(StrToFloat(Form6.Edit3.Text)+ppm_contribution);
  end;

  if ComboBox1.ItemIndex = 1 then
  begin
       ppm_contribution := total_alkalinity*0.23*(0.033814/1)*(43.6);
       Form6.Edit6.Text := FloatToStr(StrToFloat(Form6.Edit6.Text)+ppm_contribution);
  end;

  if ComboBox1.ItemIndex = 2 then
  begin
       ppm_contribution := total_alkalinity*0.56*(0.033814/1)*(14.6);
       Form6.Edit1.Text := FloatToStr(StrToFloat(Form6.Edit1.Text)+ppm_contribution);
  end;


  Form13.Visible := false ;
end;

initialization
  {$I hb_ph.lrs}

end.

