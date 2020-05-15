unit hb_insprecision;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls;

type

  { TForm7 }

  TForm7 = class(TForm)
    Button1: TButton;
    Edit1: TEdit;
    Edit2: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    procedure Button1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form7: TForm7; 

implementation

{ TForm7 }

procedure TForm7.Button1Click(Sender: TObject);
begin

Form7.Visible := false ;

end;

procedure TForm7.FormCreate(Sender: TObject);
begin

end;

initialization
  {$I hb_insprecision.lrs}

end.

