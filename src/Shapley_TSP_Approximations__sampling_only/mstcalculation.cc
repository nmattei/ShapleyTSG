#include<cstdlib>
#include<cassert>
#include<sys/types.h>
#include<time.h>
#include<cmath>

#include<fstream>
#include<map>
#include<vector>
#include<string>
#include<iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/one_bit_color_map.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>


using namespace std;

string primaryFileName;
string tourFileName;


typedef pair<double, double> Point;
vector<Point> points;

double mst_tour_length = 0;
double shorter_mst_tour_length = 0;
double tsp_tour_length = 0;

struct edge_t
{
  size_t first;
  size_t second;
};


double calc_distance(int i, int j){
    
    double x1 = points[i].first;
    double y1 = points[i].second;
          
    double x2 = points[j].first;
    double y2 = points[j].second;
          
    double dist = sqrt((x1 - x2) * (x1 - x2)  + (y1 - y2) * (y1 - y2)  );

    return dist;
}


void make_circuit(
    vector<int>& answer,
    const map<int, set<int> >::const_iterator& tree_node,
    const map<int, set<int> >& forward_tree, 
    set<int>& included_in_tour)
{
    assert(included_in_tour.find(tree_node->first) == included_in_tour.end());

    answer.push_back(tree_node->first) ;
    included_in_tour.insert(tree_node->first);
    
    if (tree_node->second.size() != 0 ){ // Has  children
        const set<int>& children = tree_node->second;

        for (set<int>::const_iterator child = children.begin() ; 
             child != children.end(); 
             child ++){
            assert(forward_tree.find(*child) != forward_tree.end());
            if(forward_tree.find(*child) == forward_tree.end()) {
                cerr<<"Unrecoverable error!\n";
                exit(0);
            }

            if(included_in_tour.find(forward_tree.find(*child)->first) == included_in_tour.end()){
                map<int, set<int> >::const_iterator next_node = forward_tree.find(*child);
                make_circuit(answer, 
                             next_node,
                             forward_tree,
                             included_in_tour);
            }
        }
    } 
}

vector<int> make_circuit(const map<int, set<int> >& forward_tree)
{
  set<int> included_in_tour;
  vector<int> answer;

  map<int, set<int> >::const_iterator tree_node = forward_tree.begin();
  make_circuit(answer, tree_node, forward_tree, included_in_tour);
  return answer;
}


void process(){

    
  using namespace boost;

  typedef adjacency_list < vecS, vecS, undirectedS,
    no_property, property < edge_weight_t, double > > Graph;
  typedef graph_traits < Graph >::edge_descriptor Edge;
  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::pair<int, int> E;

  int number_of_edges = ceil( static_cast<double>(points.size() * points.size()) / 2.0);

  cerr<<"number of edges is :: "<<number_of_edges<<endl;

  E edges[number_of_edges];
  double ws[number_of_edges];
  
  size_t count = 0;
  for(size_t i = 0 ; i < points.size(); i++){
      for(size_t j = i+1 ; j < points.size(); j++){
          cerr<<count<<endl;;
          E tmp;
          tmp.first = i;
          tmp.second = j;
          edges[count] = tmp;
          
          double dist = calc_distance(i, j);
          
          cerr<<"Distance from "<<i<<" to "<<j<<" is :"<<dist<<endl;

          ws[count] = dist;

          count++;
      }
  }

  map<int, set<int> > forward_tree;

  Graph g(edges, edges + number_of_edges, ws, points.size());//, number_of_edges);
  property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
  std::vector < Edge > spanning_tree;
  kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
  double total_weight = 0;
  for (std::vector < Edge >::iterator ei = spanning_tree.begin();
       ei != spanning_tree.end(); ++ei) 
  {
      total_weight += weight[*ei];
      
      int ids = source(*ei, g);
      int idt = target(*ei, g);

      cerr<<"tree from "<<ids<<" to "<<idt<<endl;

      if(forward_tree.find(ids) == forward_tree.end()){forward_tree[ids]=set<int>();}
      if(forward_tree.find(idt) == forward_tree.end()){forward_tree[idt]=set<int>();}
      forward_tree[ids].insert(idt);
      forward_tree[idt].insert(ids);
  }
  std::cout << "MST Tour length is = " << total_weight * 2 << std::endl;
  
  mst_tour_length = total_weight * 2 ; 

  
  vector<int> circuit = make_circuit(forward_tree);//, included_in_tour, not_included_in_tour, pending);  
  total_weight = 0;
  total_weight += calc_distance(circuit[0], circuit.back());
  for(size_t i =0 ; i < circuit.size() - 1; i++){
      total_weight += calc_distance(circuit[i], circuit[i+1]);
  }
  shorter_mst_tour_length =  total_weight;
 
  cout<<"TMST Tour length is = "<<shorter_mst_tour_length<<endl;

  assert(mst_tour_length >= shorter_mst_tour_length);

}

void parse(ifstream& in)
{
    if(in.eof())return;
    string tmp;
    do{
        getline(in, tmp);
        if(!tmp.size())continue;
        string str;
        istringstream iss(tmp);
        iss>>str;
        if(str == "NAME:")continue;
        if(str == "TYPE:")continue;
        if(str == "COMMENT:")continue;
        if(str == "DIMENSION:")continue;
        if(str == "EDGE_WEIGHT_TYPE:")continue;
        if(str == "DISPLAY_DATA_TYPE:")continue;
        if(str == "NODE_COORD_SECTION")continue;
        if(str == "EOF")break;

        double x;
        double y;
        iss>>x;
        iss>>y;
        points.push_back(Point(x, y));
        
    }while(!in.eof());
}

void parse()
{    
    ifstream file;
    file.open (primaryFileName.c_str());
    parse(file);
    file.close();
}

void parse_tour(ifstream& in)
{
    cerr<<"Parsing tour"<<endl;
    
    if(in.eof())return;
    string tmp;
    getline(in, tmp);

    if(in.eof())return;
    
    size_t first_vertex_id = 1000000;
    double total_distance=0;
    size_t old_vertex_id = 1000000;
    size_t vertex_id;
    do{
        getline(in, tmp);
        if(!tmp.size())continue;
        cerr<<tmp<<endl;;

        istringstream iss(tmp);
        
        do {
            iss>>vertex_id;
            if(old_vertex_id == 1000000){
                old_vertex_id = vertex_id;

                first_vertex_id = vertex_id; 
                continue;
            } else {
                cerr<<vertex_id<<" to "<<vertex_id<<endl;

                Point p1 = points[vertex_id];
                Point p2 = points[old_vertex_id];

                
                double x1 = p1.first;
                double y1 = p1.second;
                
                double x2 = p2.first;
                double y2 = p2.second;
                
                double distance = sqrt((x1 - x2) * (x1 - x2)  + (y1 - y2) * (y1 - y2)  );
                total_distance+=distance;
            }

            old_vertex_id = vertex_id;
        }while(!iss.eof());
        
        
    }while(!in.eof());

    Point p1 = points[vertex_id];
    Point p2 = points[first_vertex_id];

                
    double x1 = p1.first;
    double y1 = p1.second;
                
    double x2 = p2.first;
    double y2 = p2.second;
                
    double distance = sqrt((x1 - x2) * (x1 - x2)  + (y1 - y2) * (y1 - y2)  );
    total_distance+=distance;
    
    tsp_tour_length = total_distance; 
    
    cout<<"TSP Tour length is = "<<total_distance<<endl;
}

void parse_tour()
{    
    ifstream file;
    file.open (tourFileName.c_str());
    parse_tour(file);
    file.close();
}

int main(int argc, char** argv){
    assert(argc>=2);
    primaryFileName = string(argv[1]);
    parse();
    
    // calculate_mst();
    // report_mst_length();
    
    process();

    if(argc==3){
        cerr<<"Processing Euclidean tour.\n";
        tourFileName = string(argv[2]);
        parse_tour();
        assert(mst_tour_length / 2 < tsp_tour_length);
        assert(mst_tour_length < 2 * tsp_tour_length);
        assert(tsp_tour_length <= mst_tour_length);
        assert(shorter_mst_tour_length >= tsp_tour_length);
        assert(shorter_mst_tour_length <= mst_tour_length);

    }

    return 0;
}
