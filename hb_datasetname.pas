unit hb_datasetname;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls;

type

  { TForm10 }

  TForm10 = class(TForm)
    Button1: TButton;
    Edit1: TEdit;
    Label1: TLabel;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form10: TForm10; 

implementation

uses HB_Main ;

{ TForm10 }

procedure TForm10.Button1Click(Sender: TObject);

begin

datasetname := Edit1.Text ;

ShowMessage('Data Set name set to ' + Edit1.Text + '. Database creation will finish when you add a data point');

Form10.Visible := False ;

end;

initialization
  {$I hb_datasetname.lrs}

end.

