#include<cstdlib>
#include<cassert>
#include<sys/types.h>
#include<time.h>

#include<fstream>
#include<map>
#include<vector>
#include<string>
#include<iostream>


using namespace std;

size_t points_count;
typedef pair<double, double> Point;
vector<Point> points1;
vector<Point> points2;

string smallerFileName;
string largerFileName;

void write(ofstream& file, const vector<Point>& points){
    file<<"NAME: ulysses22.tsp"<<endl;
    file<<"TYPE: TSP"<<endl;
    file<<"COMMENT: Shapley Experiment"<<endl;
    file<<"DIMENSION: "<<points.size()<<endl;
    file<<"EDGE_WEIGHT_TYPE: EUC_2D"<<endl;
    file<<"DISPLAY_DATA_TYPE: COORD_DISPLAY"<<endl;
    file<<"NODE_COORD_SECTION"<<endl;
    for(size_t i = 1; i <= points.size(); i++){
        file<<i<<" "<<points[i-1].first<<" "<<points[i-1].second<<endl;
    }
    file<<"EOF"<<endl;
}

void write_all(){
    
    ofstream file;
    file.open (smallerFileName.c_str());
    write(file, points1);
    file.close();
    file.open (largerFileName.c_str());
    write(file, points2);
    file.close();
}

int main(int argc, char** argv){
    assert(argc==4);
    points_count = atoi(argv[1]);
    smallerFileName = string(argv[2]);
    largerFileName = string(argv[3]);
    
    time_t t1; (void) time(&t1);
    srand48((long) t1);
    
    double x;
    double y;
    
    for (size_t i =0 ; i < points_count; i++){
        x = drand48() * 10000;
        y = drand48() * 10000;
        points1.push_back(Point(x, y));
        points2.push_back(Point(x, y));
    }
    x = drand48() * 10000;
    y = drand48() * 10000;
    points2.push_back(Point(x, y));

    write_all();

    return 0;
}
