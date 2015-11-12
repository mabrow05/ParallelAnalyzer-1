#include <map>
using namespace std;


unsigned int getSrcRunPeriod(int runNumber) {
  unsigned int calibrationPeriod=0;
  if (runNumber <= 17297) {
    calibrationPeriod=1;
  }
  else if (runNumber <= 17439) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 17734) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 17955) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 18386) {
    calibrationPeriod=5;
  }
  else if (runNumber <= 18683) {
    calibrationPeriod=6;
  }
  else if (runNumber <= 18994) {
    calibrationPeriod=7;
  }
  else if (runNumber <= 19239) {
    calibrationPeriod=8;
  }
  else if (runNumber <= 19544) {
    calibrationPeriod=9;
  }
  else if (runNumber <= 20000) {
    calibrationPeriod=11;
  }
  return calibrationPeriod;
}

unsigned int getXeRunPeriod(int runNumber) {
  unsigned int calibrationPeriod=0;
  if (runNumber <= 18080) {
    calibrationPeriod=2;
  }
  else if (runNumber <= 18389) {
    calibrationPeriod=3;
  }
  else if (runNumber <= 18711) {
    calibrationPeriod=4;
  }
  else if (runNumber <= 19238) {
    calibrationPeriod=5;
  }
  //else if (runNumber <= 19872) {
  //calibrationPeriod=6;
  //}
  else if (runNumber <= 20000) {
    calibrationPeriod=7;
  }
  else {
    cout << "Bad run number\n";
    exit(0);}

  return calibrationPeriod;
}

unsigned int getGainReferenceRun(int runNumber) {
  map <unsigned int,unsigned int> refRun;
  refRun[1]=17238;
  refRun[2]=17370;
  refRun[3]=17521;
  refRun[4]=17925;
  refRun[5]=18361;
  refRun[6]=18621;
  refRun[7]=18749;
  refRun[8]=19232;
  refRun[9]=19359;
  refRun[10]=19511;
  refRun[11]=19857;
  refRun[12]=19899;
  unsigned int calPeriod = getSrcRunPeriod(runNumber);
  map <unsigned int, unsigned int>::iterator it = refRun.find(calPeriod);
  return it->second;
}
