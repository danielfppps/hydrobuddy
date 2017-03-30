unit hb_freedom;

{$mode objfpc}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls;

type

  { TForm12 }

  TForm12 = class(TForm)
    Button1: TButton;
    ComboBox1: TComboBox;
    Label1: TLabel;
    Label2: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form12: TForm12; 

implementation

{ TForm12 }

procedure TForm12.Button1Click(Sender: TObject);
begin
  Form12.Visible := false ;
end;

initialization
  {$I hb_freedom.lrs}

end.

