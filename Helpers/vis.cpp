#include "vis.h"

Vis::Vis()
{
    cellsIN.clear();
    cellsWALL.clear();
    cellsBoardB.clear();
    cellsBoardL.clear();
    cellsBoardT.clear();
    cellsBoardR.clear();
    LoadIN();
    LoadWALL();
    LoadB();
    LoadL();
    LoadR();
    LoadT();
    LoadDic();

}

void Vis::LoadIN(){
// We read the positions of the grid that are interior from a txt file
    std::ifstream infile("Files/IN.txt");
    std::string line;
    int L;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss>>L;
    std::cout<<"In Long"<<std::to_string(L)<<std::endl;
    for(int l=0;l<L;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;

        int k =fmod(cellR,61);
        int i = (cellR-k)/61;
        // We store the cell with its characteristics
        CellIN c;
        c.cellID=cellR;
        c.i=i;
        c.k=k;
        c.x=(i+fmod(i,2))/2*cos(M_PI*i);
        c.z=(k+fmod(k,2))/2*cos(M_PI*k);
        cellsIN.push_back(c);
    }

}

void Vis::LoadWALL(){
// We read the positions of the grid that are walls from a txt file
    std::ifstream infile("Files/WALL.txt");
    std::string line;
    int L;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss>>L;

    for(int l=0;l<L;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;

        int k =fmod(cellR,61);
        int i = (cellR-k)/61;
        // We store the cell with its characteristics
        CellWALL c;
        c.cellID=cellR;
        c.i=i;
        c.k=k;
        c.x=(i+fmod(i,2))/2*cos(M_PI*i);
        c.z=(k+fmod(k,2))/2*cos(M_PI*k);
        cellsWALL.push_back(c);
    }

}
// We read the positions of the grid that are borders from a txt file : T-Top, B-Bottom, R-Right, L-Left
void Vis::LoadT(){
    std::ifstream infile("Files/T.txt");
    std::string line;
    int L;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss>>L;

    for(int l=0;l<L;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;

        int k =fmod(cellR,61);
        int i = (cellR-k)/61;

        CellBorder c;
        c.cellID=cellR;
        c.i=i;
        c.k=k;
        c.x=(i+fmod(i,2))/2*cos(M_PI*i);
        c.z=(k+fmod(k,2))/2*cos(M_PI*k);
        cellsBoardT.push_back(c);
    }
}
void Vis::LoadB(){
    std::ifstream infile("Files/B.txt");
    std::string line;
    int L;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss>>L;

    for(int l=0;l<L;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;

        int k =fmod(cellR,61);
        int i = (cellR-k)/61;

        CellBorder c;
        c.cellID=cellR;
        c.i=i;
        c.k=k;
        c.x=(i+fmod(i,2))/2*cos(M_PI*i);
        c.z=(k+fmod(k,2))/2*cos(M_PI*k);
        cellsBoardB.push_back(c);
    }
}
void Vis::LoadR(){
    std::ifstream infile("Files/R.txt");
    std::string line;
    int L;
    std::getline(infile, line);
    std::istringstream iss(line);
    iss>>L;

    for(int l=0;l<L;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;

        int k =fmod(cellR,61);
        int i = (cellR-k)/61;

        CellBorder c;
        c.cellID=cellR;
        c.i=i;
        c.k=k;
        c.x=(i+fmod(i,2))/2*cos(M_PI*i);
        c.z=(k+fmod(k,2))/2*cos(M_PI*k);
        cellsBoardR.push_back(c);
    }
}

void Vis::LoadL(){
        std::ifstream infile("Files/L.txt");
        std::string line;
        int L;
        std::getline(infile, line);
        std::istringstream iss(line);
        iss>>L;

        for(int l=0;l<L;l++){
            std::getline(infile, line);
            std::istringstream iss(line);
            int cellR;
            iss>>cellR;

            int k =fmod(cellR,61);
            int i = (cellR-k)/61;

            CellBorder c;
            c.cellID=cellR;
            c.i=i;
            c.k=k;
            c.x=(i+fmod(i,2))/2*cos(M_PI*i);
            c.z=(k+fmod(k,2))/2*cos(M_PI*k);
            cellsBoardL.push_back(c);
        }
    }
//We need a dictionary to because the distributions of cellID is different for the map and the visibility grid.
void Vis::LoadDic(){
    std::ifstream infile("Files/DIC.txt");
    std::string line;
    for(int l=0;l<61*61;l++){
        std::getline(infile, line);
        std::istringstream iss(line);
        int cellR;
        iss>>cellR;
        dic.push_back(cellR);
    }
}

void Vis::CalculateVisibility(){

    int numRays=100000;
       //we'll throw 100000 rays
    for(int r=0;r<numRays;r++){
        std::cout<<"------------------------------------"<<std::endl;
      //  std::cout<<"RAYS: "<<std::to_string(r)<<std::endl;
        //We decide two borders randomly
        int b1 =rand()%4;
        int b2= b1;
        while(b2==b1) b2=rand()%4;

        //We decide which cell from each border the ray is going to go through
        int c1=rand()%61;
        int c2=rand()%61;

        // we make sure that the rays always go from left to right
        if(b2==3){
            b2=b1;
            b1=3;
        }
        if(b1==2){
            b1=b2;
            b2=2;
        }
        if(b1<2&&b2<2){
            if(c1>c2){
                int temp=b1;
                b1=b2;
                b2=temp;
                temp=c1;
                c1=c2;
                c2=temp;
            }

        }
        std::cout<<"RAYS: "<<std::to_string(r)<<" L1: "<<std::to_string(b1)<<" L2: "<<std::to_string(b2)<<" c1: "<<std::to_string(c1)<<" c2: "<<std::to_string(c2)<<std::endl;


        //We take the two corresponding cells
        CellBorder CB1;
        CellBorder CB2;
        switch(b1){
        case 0:
            CB1 = cellsBoardT[c1];
            break;
        case 1:
            CB1 = cellsBoardB[c1];
            break;
        case 2:
            CB1 = cellsBoardR[c1];
            break;
        case 3:
            CB1 = cellsBoardL[c1];
            break;
        }
        switch(b2){
        case 0:
            CB2 = cellsBoardT[c2];
            break;
        case 1:
            CB2 = cellsBoardB[c2];
            break;
        case 2:
            CB2 = cellsBoardR[c2];
            break;
        case 3:
            CB2 = cellsBoardL[c2];
            break;
        }

        //We make our current cell the starting cell
        int curr_cell;
        for(int i=0;i<(int)dic.size();i++){
            if(dic[i]==CB1.cellID)curr_cell=i;
        }
        //we take the current k and i
        int curr_k=fmod(curr_cell,61);
        int curr_i=(curr_cell-curr_k)/61;
        //We take the ending cell
        int end_cell;
        for(int i=0;i<(int)dic.size();i++){
            if(dic[i]==CB2.cellID)end_cell=i;
        }
        int end_k=fmod(end_cell,61);
        int end_i=(end_cell-end_k)/61;

        //We compute the linear expression of the ray with the two cells
        float pos[2]={0,0};

        pos[0]=(end_i-curr_i);
        pos[1]=(end_k-curr_k);

        //z=m*x+n
        float m=pos[1]/pos[0];  //slope
        int sign=1;
        if(m<0)sign=-1; // the ray going up or down on the z axis

        float n=curr_k-m*curr_i;  //bias

        pos[0]=(float)curr_i;
        pos[1]=(float)curr_k;

        std::vector<int> cluster;
        // We go through all the cells that the ray crosses.
        while(curr_i<61&&curr_k<61&&curr_i>-1&&curr_k>-1){
        // Calculate distance between the current position and the intersection with the next x and z axis (cell borders)
        float h_l= curr_i+0.5 -pos[0];
        float v_l= curr_k+sign*0.5-pos[1];

        if(h_l<v_l){ // if the distance in x is lower than in z
            pos[0]=curr_i+0.5; // we move x (we'll be in the border of the cell)
            pos[1]=m*(curr_i+0.5)+n; // we compute z
            curr_i+=1;// we increase the current i
        }
        else{// if the distance in z is lower than in x
            pos[0]=(curr_k+sign*0.5-n)/m; // we compute x
            pos[1]=curr_k+sign*0.5; // we move z (we'll be in the border of the cell)
            curr_k+=sign; // we increase the current i according to the direction of the slope
        }
        curr_cell=61*curr_i+curr_k; // we compute the current cell

        std::cout<<std::to_string(dic[curr_cell])<<"   current i "<<std::to_string(curr_i)<<"   current k "<<std::to_string(curr_k)<<std::endl;
            //we check if the current cell is an interior cell
            if(checkIN(dic[curr_cell])){
                lastIn=true;
                cluster.push_back(PosCurrentCell);// we add the cell to the cluster
            }
            //we check if the current cell is a wall
            if(checkWall(dic[curr_cell])){
                if(lastIn&&cluster.size()>0){
                    saveCluster(cluster);//if the last cell visited was interior, we save the cluster, and start a new one
                    cluster.clear();
                }
                lastIn=false;
            }
        }
    }
}
void Vis::CalculateBresenham(){

    int numRays=100000;
//we'll throw 100000 rays
    for(int r=0;r<numRays;r++){
        std::cout<<"------------------------------------"<<std::endl;
        //We decide two borders randomly
        int b1 =rand()%4;
        int b2= b1;
        while(b2==b1) b2=rand()%4;

         //We decide which cell from each border the ray is going to go through
        int c1=rand()%61;
        int c2=rand()%61;
           // we make sure that the rays always go from left to right
        if(b2==3){
            b2=b1;
            b1=3;
        }
        if(b1==2){
            b1=b2;
            b2=2;
        }
        if(b1<2&&b2<2){
            if(c1>c2){
                int temp=b1;
                b1=b2;
                b2=temp;
                temp=c1;
                c1=c2;
                c2=temp;
            }

        }
        std::cout<<"RAYS: "<<std::to_string(r)<<" L1: "<<std::to_string(b1)<<" L2: "<<std::to_string(b2)<<" c1: "<<std::to_string(c1)<<" c2: "<<std::to_string(c2)<<std::endl;
        //We take the two corresponding cells
        CellBorder CB1;
        CellBorder CB2;
        switch(b1){
        case 0:
            CB1 = cellsBoardT[c1];
            break;
        case 1:
            CB1 = cellsBoardB[c1];
            break;
        case 2:
            CB1 = cellsBoardR[c1];
            break;
        case 3:
            CB1 = cellsBoardL[c1];
            break;
        }
        switch(b2){
        case 0:
            CB2 = cellsBoardT[c2];
            break;
        case 1:
            CB2 = cellsBoardB[c2];
            break;
        case 2:
            CB2 = cellsBoardR[c2];
            break;
        case 3:
            CB2 = cellsBoardL[c2];
            break;
        }


        //We make our current cell the starting cell
        int curr_cell;
        for(int i=0;i<(int)dic.size();i++){
            if(dic[i]==CB1.cellID)curr_cell=i;
        }
         //we take the current k and i
        int curr_k=fmod(curr_cell,61);
        int curr_i=(curr_cell-curr_k)/61;
        //We take the ending cell
        int end_cell;
        for(int i=0;i<(int)dic.size();i++){
            if(dic[i]==CB2.cellID)end_cell=i;
        }
        int end_k=fmod(end_cell,61);
        int end_i=(end_cell-end_k)/61;

        //We compute the linear expression of the ray with the two cells
        float pos[2]={0,0};

        pos[0]=(end_i-curr_i);
        pos[1]=(end_k-curr_k);

        //z=m*x+n
        float m=pos[1]/pos[0];//slope
        int sign=1;
        if(m<0)sign=-1; // the ray going up or down on the z axis

        float n=curr_k-m*curr_i; //bias

        pos[0]=(float)curr_i;
        pos[1]=(float)curr_k;

        lastI=curr_i;
        lastK=curr_k;

        std::vector<int> cluster;
        // We go through all the cells that the ray crosses.
        while(curr_i<61&&curr_k<61&&curr_i>-1&&curr_k>-1){

        // Calculate distance between the current position and the intersection with the next x and z axis (cell borders)
        float h_l= curr_i+0.5 -pos[0];
        float v_l= curr_k+sign*0.5-pos[1];

        if(h_l<v_l){// if the distance in x is lower than in z
            pos[0]=curr_i+0.5;// we move x (we'll be in the border of the cell)
            pos[1]=m*(curr_i+0.5)+n;// we compute z
            if(lastAxis==1){axisChanged=true; //We check if there has been a change of direction/axis
            lastAxis=0;}
            else{axisChanged=false;}
            curr_i+=1;// we increase the current i

        }
        else{// if the distance in z is lower than in x
            pos[0]=(curr_k+sign*0.5-n)/m;// we compute x
            pos[1]=curr_k+sign*0.5;// we move z (we'll be in the border of the cell)
            if(lastAxis==0){//We check if there has been a change of direction/axis
                axisChanged=true;
            lastAxis=1;}
            else{axisChanged=false;}
            lastK=curr_k;
            curr_k+=sign;
        }
        curr_cell=61*curr_i+curr_k;// we compute the current cell

        std::cout<<std::to_string(dic[curr_cell])<<"   current i "<<std::to_string(curr_i)<<"   current k "<<std::to_string(curr_k)<<std::endl;
            //we check if the current cell is an interior cell
            if(checkIN(dic[curr_cell])){
                lastIn=true;
                cluster.push_back(PosCurrentCell);//we save that position to the current cluster
            }
            //If the direction/ axis has changed we compute the supercover
            if(axisChanged){
                //Depending on what axis was the previous axis, we check for one or another cell contained in the corner formed
                // and we check if it is an interior cell
                if(lastAxis==1&&(61*(curr_i-1)+curr_k)<61*61&&(61*(curr_i-1)+curr_k)>=0){
                    std::cout<<std::to_string(dic[61*(curr_i-1)+curr_k])<<"   SUPERCOVER  last 0"<<std::endl;
                    if(checkIN(dic[61*(curr_i-1)+curr_k])){
                        lastIn=true;
                        cluster.push_back(PosCurrentCell); //we save that position to the current cluster
                    }
                }
                if(lastAxis==0&&(61*(curr_i)+curr_k-sign)<61*61&&(61*(curr_i)+curr_k-sign)>=0){
                    std::cout<<std::to_string(dic[61*(curr_i)+curr_k-sign])<<"   SUPERCOVER last 1"<<std::endl;
                    if(checkIN(dic[61*(curr_i)+curr_k-sign])){
                        lastIn=true;
                        cluster.push_back(PosCurrentCell);//we save that position to the current cluster

                    }

                }

            }
            //we check if the current cell is a wall
            if(checkWall(dic[curr_cell])){
                if(lastIn&&cluster.size()>0){
                    saveCluster(cluster);
                    cluster.clear();
                }
                lastIn=false;
            }
        }
    }
}

bool Vis::checkWall(int c){
    // we run all the positions where there is a wall, and check if coincides or not
    for(int t=0;t<(int)cellsWALL.size();t++){
        if(c==cellsWALL[t].cellID){
           return true;
        }
    }
    return false;
}
bool Vis::checkIN(int c){
    // we run all the positions considered as intirior of the museu, and check if coincides or not
    for(int t=0;t<(int)cellsIN.size();t++){
        if(c==cellsIN[t].cellID){
            PosCurrentCell=t;// If coincides, we save that position to the current cluster
           return true;
        }
    }
    return false;
}

void Vis::saveCluster(std::vector<int> cluster){
    // For all the cells in a cluster we add all the cells to their visibility vector
    for(int i=0;i<(int)cluster.size();i++){
        for(int j=0;j<(int)cluster.size();j++){
                cellsIN[cluster[i]].visibility.push_back(cellsIN[cluster[j]].cellID);
        }
    }
}
void Vis::makeUnique(){
// We delete all the repeated values
for(int i =0;i<(int)cellsIN.size();i++){
    std::sort( cellsIN[i].visibility.begin(), cellsIN[i].visibility.end() );
    cellsIN[i].visibility.erase(std::unique( cellsIN[i].visibility.begin(), cellsIN[i].visibility.end() ), cellsIN[i].visibility.end() );

}

}

void Vis::writeToFile(){
    //We write into a file,the interior cells and their visibility vectors
    std::ofstream myfile;
    myfile.open ("Map/Visibility.txt");
    for(int i =0;i<(int)cellsIN.size();i++){
        myfile <<std::to_string(cellsIN[i].cellID);
        for(int j =0;j<(int)cellsIN[i].visibility.size();j++){
            myfile <<" ";
            myfile <<std::to_string(cellsIN[i].visibility[j]);
        }
        myfile << "\n";
    }
    myfile.close();


}
