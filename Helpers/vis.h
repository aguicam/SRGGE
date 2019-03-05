#ifndef VIS_H
#define VIS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include <cstdlib>
#include <algorithm>

class Vis
{
public:
    Vis();
    void LoadIN();
    void LoadWALL();
    void LoadT();
    void LoadB();
    void LoadR();
    void LoadL();
    void CalculateVisibility();
    bool checkWall(int c);
    bool checkIN(int c);
    void saveCluster(std::vector<int> cluster);
    void makeUnique();
    void writeToFile();
    void LoadDic();
    void CalculateBresenham();

    struct CellIN{
      int cellID;
      int i;
      int k;
      float x;
      float z;
      std::vector<int> visibility;
    };
    struct CellWALL{
      int cellID;
      int i;
      int k;
      float x;
      float z;
    };
    struct CellBorder{
      int cellID;
      int i;
      int k;
      float x;
      float z;
    };

    std::vector<CellIN> cellsIN;
    std::vector<CellWALL> cellsWALL;
    std::vector<CellBorder> cellsBoardT;
    std::vector<CellBorder> cellsBoardB;
    std::vector<CellBorder> cellsBoardR;
    std::vector<CellBorder> cellsBoardL;
    std::vector<int> dic;
    int PosCurrentCell;
    bool axisChanged =false;
    int lastAxis=0;  //x=0 , z=1;
    int lastK;
    int lastI;
  //  bool oneWall=false;
    bool lastIn=false;
};

#endif // VIS_H
