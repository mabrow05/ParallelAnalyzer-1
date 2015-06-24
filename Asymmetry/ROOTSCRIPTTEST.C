{

#include <string>;
  string dbAddress = "localhost";
  string dbname = "UCNADB_full_quad1";
  string port = "3306";
  string dbAddressFull = "mysql://"+dbAddress+":"+port+"/"+dbname;
  string dbUser = "ucn";
  string dbPass = "T1td4ctc";
  TSQLServer *db = new TMySQLServer(dbAddressFull.c_str(),dbUser.c_str(),dbPass.c_str());

  TSQLResult *res = db->Query("SELECT * from posmap_set where posmap_set_id=283");
  Int_t fieldCount = res->GetFieldCount(); 
  
  TSQLRow *row;
  while (row = res->Next())
    {
      for (int i=0; i<fieldCount; i++)
	{
	  cout << row->GetField(i) << " ";
	}
      cout << endl;
    }
}
 
