
vector<mpf_class>  moat_margin_values()
{
    set<int> all_points;
    for(uint i = 0 ; i < points.size(); i++){
        all_points.insert(i);
    }

    WRITING_TO_CONCORDE__FUNCTION(all_points);

    
    ostringstream command;
    

    command<<"matlab -r \"cd('/home/cgretton/Work2015/vrg/src/MMP/') ; [fval, w] = maxmoatpacking('"<<INSTALLED_PREFIX
        //January 9th 2014//command<<"matlab -r \"cd('/home/cgretton/Work2013/Software/Casey/MMP/') ; [fval, w] = maxmoatpacking('"<<INSTALLED_PREFIX
           <<concorde_file_name<<"'); exit; \" -nosplash -nojvm -nodisplay > matlab.data"; 

    cerr<<command.str()<<endl;
    //exit(-1);

    //cerr<<command.str().c_str()<<flush;
    //exit(0);
    int res  = system(command.str().c_str());

    assert(res == 0);
    ifstream file;
    file.open("matlab.data");//{assert(0);cerr<<"Failing;\n";exit(0);}
    if(file.fail()){
        assert(0);cerr<<"Failing;\n";exit(0);
    } else if (file.eof()){
        assert(0);cerr<<"Failing;\n";exit(0);
    }

    vector<mpf_class> answer(all_points.size()); for (int i =0 ; i < all_points.size(); i++){answer[i] = 0.0;}

    while(!file.eof()){
        string tmp;
        getline(file, tmp);
        
        istringstream iss(tmp);
        string first_word;
        iss>>first_word;
        //cerr<<first_word<<endl;

        
        if("Allocation" == first_word){
            double multiplicative_factor = 1.0;
            //getline(file, tmp);
            for (int i = 0 ; i < all_points.size(); i++){
                
                getline(file, tmp);
                //cerr<<tmp<<endl;
                istringstream iss(tmp);
                double number;
                string possible_factor;
                iss>>number;
                iss>>possible_factor;

                if ("*" == possible_factor){
                    i--;
                    multiplicative_factor = number;
                    getline(file, tmp);
                    continue;
                }

                
                if(i==0 && number != 0.0) {
                    cerr<<"Depot was given a value : "<<number<<endl;
                    exit(0);
                }

                answer[i] = number * multiplicative_factor;
            }
            //exit(0);
            break;
        }
        
    }
    file.close();

    return answer;

}

map<set<int>, double> TSP__calculate_value_of__cache;
map<set<int>, double> HELD_KARP__calculate_value_of__cache;
// double TSP__calculate_value_of_matrix(const set<int>& local_points){
//     if(local_points.size() == 1){
//         return 0.0;
//     }
//     if(local_points.size() == 2){
//         set<int>::const_iterator p1 = local_points.begin();
//         set<int>::const_iterator p2 = local_points.begin();
//         p2++;
//         return calc_distance_matrix(*p1, *p2, true) + calc_distance_matrix(*p2, *p1, true);
//     }
//     map<set<int>, double>::const_iterator cached_value = TSP__calculate_value_of__cache.find(local_points);
//     if( cached_value != TSP__calculate_value_of__cache.end()){
//         //cerr<<"Using cache.\n";
//         //cerr<< cached_value->second<<" -- ";
//         return  cached_value->second;
//     }

//     WRITING_TO_CONCORDE__FUNCTION__matrix(local_points);

//     ostringstream command;
//     command<<INSTALLED_PREFIX
//            <<"/concorde/TSP/concorde -o tour.output "
//            <<concorde_file_name
//            <<" 2> /dev/null | grep \"Optimal Solution:\" | awk '{print $3}' > concorde.data";
    
//     int res  = system(command.str().c_str());
//     assert(res == 0);
    
//     //cerr<<"Done exec.\n";
//     ifstream file;
//     file.open("concorde.data");//{assert(0);cerr<<"Failing;\n";exit(0);}
//     if(file.fail()){
//         assert(0);cerr<<"Failing;\n";exit(0);
//     } else if (file.eof()){
//         assert(0);cerr<<"Failing;\n";exit(0);
//     }
    
//     double concorde_value;
//   {
//       //cerr<<"Done exec.\n";
//       ifstream file;
//       file.open("concorde.data");//{assert(0);cerr<<"Failing;\n";exit(0);}
//       if(file.fail()){
//           assert(0);cerr<<"Failing;\n";exit(0);
//       } else if (file.eof()){
//           assert(0);cerr<<"Failing;\n";exit(0);
//       }
//       string tmp;
//       getline(file, tmp);
//       file.close();
//       istringstream iss(tmp);
//       iss>>concorde_value;
//   }


//     double tour_length;
//   {
      
//     ifstream file;
//     file.open("tour.output");
//     if(file.fail()){
//         assert(0);cerr<<"Failing;\n";exit(0);
//     } else if (file.eof()){
//         assert(0);cerr<<"Failing;\n";exit(0);
//     }


//     pair<double, vector<double> > parsed_tour = parse_tour(file, calc_distance_matrix__sym);//calc_distance);
//     file.close();
//     tour_length = parsed_tour.first;
//   }
    
//     TSP__calculate_value_of__cache[local_points] = value;

//     return value;
// }


double HELD_KARP__calculate_value_of(const set<int>& local_points){

    //assert(0 != local_points.size());

    if(local_points.size() == 1){
        return 0;
    }
    
    if(local_points.size() == 2){
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return COST_OF_EDGE__FUNCTION(*first, *second) * 2;
    }

    if(local_points.size() == 3){
        set<int>::const_iterator third = local_points.begin()++;third++;
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return COST_OF_EDGE__FUNCTION(*first, *second) + COST_OF_EDGE__FUNCTION(*second, *third) +  COST_OF_EDGE__FUNCTION(*third, *first) ;
    }


    map<set<int>, double>::const_iterator cached_value = HELD_KARP__calculate_value_of__cache.find(local_points);
    if( cached_value != HELD_KARP__calculate_value_of__cache.end()){
        //cerr<<"Using cache.\n";
        //cerr<< cached_value->second<<" -- ";
        return  cached_value->second;
    }


    // ostringstream debugging;
    // debugging<<"cp "<<concorde_file_name<<" backup.concorde";
    // int   res2  = system(debugging.str().c_str());
    // assert(res2 == 0);

    WRITING_TO_CONCORDE__FUNCTION(local_points);
    ostringstream command;
    command<<INSTALLED_PREFIX
           <<"/concorde/HELDKARP/heldkarp -o tour.output "
           <<concorde_file_name
           <<" 2> /dev/null | grep \"Optimal Solution:\" | awk '{print $3}' > concorde.data";
    
    cerr<<command.str()<<endl;

    int res  = system(command.str().c_str());
    assert(res == 0);
    
    
    double concorde_value;
  {
      //cerr<<"Done exec.\n";
      ifstream file;
      file.open("concorde.data");//{assert(0);cerr<<"Failing;\n";exit(0);}
      if(file.fail()){
          assert(0);cerr<<"Failing;\n";exit(0);
      } else if (file.eof()){
          assert(0);cerr<<"Failing;\n";exit(0);
      }
      string tmp;
      getline(file, tmp);
      file.close();
      istringstream iss(tmp);
      iss>>concorde_value;
  }

    double tour_length;
  {
      
    ifstream file;
    file.open("tour.output");
    if(file.fail()){
        assert(0);cerr<<"Failing;\n";exit(0);
    } else if (file.eof()){
        assert(0);cerr<<"Failing;\n";exit(0);
    }

    pair<double, vector<double> > parsed_tour = parse_tour(file, COST_OF_EDGE__FUNCTION);
    file.close();
    tour_length = parsed_tour.first;
  }
    
  HELD_KARP__calculate_value_of__cache[local_points] = tour_length;//value;

  return tour_length;//value;
}

double TSP__calculate_value_of__indigo(const set<int>& local_points);

double TSP__calculate_value_of__concorde(const set<int>& local_points);


double(*TSP__calculate_value_of)(const set<int>& local_points);
// double *TSP__calculate_value_of(const set<int>& local_points){
//     //return TSP__calculate_value_of__concorde(local_points);
//     return TSP__calculate_value_of__indigo(local_points);
// }


double TSP__calculate_value_of__indigo(const set<int>& local_points){

    //assert(0 != local_points.size());

    if(local_points.size() == 1){
        return 0;
    }
    
    if(local_points.size() == 2){
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();

        return COST_OF_EDGE__FUNCTION(*first, *second) + COST_OF_EDGE__FUNCTION(*second, *first);
    }

    if(local_points.size() == 3){
        set<int>::const_iterator third = local_points.begin()++;third++;
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return std::min(COST_OF_EDGE__FUNCTION(*first, *second) + COST_OF_EDGE__FUNCTION(*second, *third) +  COST_OF_EDGE__FUNCTION(*third, *first)
                        , COST_OF_EDGE__FUNCTION(*first, *third) + COST_OF_EDGE__FUNCTION(*third, *second) +  COST_OF_EDGE__FUNCTION(*second, *first));
    }

    map<set<int>, double>::const_iterator cached_value = TSP__calculate_value_of__cache.find(local_points);
    if( cached_value != TSP__calculate_value_of__cache.end()){
        //cerr<<"Using cache.\n";
        //cerr<< cached_value->second<<" -- ";
        return  cached_value->second;
    }


    // ostringstream debugging;
    // debugging<<"cp "<<concorde_file_name<<" backup.concorde";
    // int   res2  = system(debugging.str().c_str());
    // assert(res2 == 0);

    WRITING_TO_CONCORDE__FUNCTION(local_points);
    ostringstream command;
    command<<"/home/cgretton/bin/indigo -i lns-3-1000-5 -r tour.output "
           <<concorde_file_name
           <<" 2> /dev/null";
    
    cerr<<command.str()<<endl;

    int res  = system(command.str().c_str());
    assert(res == 0);

    double tour_length;
  {
      
    ifstream file;
    file.open("tour.output");
    if(file.fail()){
        assert(0);cerr<<"Failing;\n";exit(0);
    } else if (file.eof()){
        assert(0);cerr<<"Failing;\n";exit(0);
    }

    pair<double, vector<double> > parsed_tour = parse_tour__indigo(file, COST_OF_EDGE__FUNCTION);
    file.close();
    tour_length = parsed_tour.first;
  }
  
  // cerr<<"I calculate the tour length to be :: "<<tour_length<<endl
  //     <<"And concorde reports the your length to be :: "<<concorde_value<<endl;

    
  TSP__calculate_value_of__cache[local_points] = tour_length;//value;

  return tour_length;//value;
}


double TSP__calculate_value_of__concorde(const set<int>& local_points){

    //assert(0 != local_points.size());

    if(local_points.size() == 1){
        return 0;
    }
    
    if(local_points.size() == 2){
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return COST_OF_EDGE__FUNCTION(*first, *second) * 2;
    }

    if(local_points.size() == 3){
        set<int>::const_iterator third = local_points.begin()++;third++;
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return COST_OF_EDGE__FUNCTION(*first, *second) + COST_OF_EDGE__FUNCTION(*second, *third) +  COST_OF_EDGE__FUNCTION(*third, *first) ;
    }

    map<set<int>, double>::const_iterator cached_value = TSP__calculate_value_of__cache.find(local_points);
    if( cached_value != TSP__calculate_value_of__cache.end()){
        //cerr<<"Using cache.\n";
        //cerr<< cached_value->second<<" -- ";
        return  cached_value->second;
    }


    // ostringstream debugging;
    // debugging<<"cp "<<concorde_file_name<<" backup.concorde";
    // int   res2  = system(debugging.str().c_str());
    // assert(res2 == 0);

    WRITING_TO_CONCORDE__FUNCTION(local_points);
    ostringstream command;
    command<<INSTALLED_PREFIX
           <<"/concorde/TSP/concorde -o tour.output "
           <<concorde_file_name
           <<" 2> /dev/null | grep \"Optimal Solution:\" | awk '{print $3}' > concorde.data";
    
    //cerr<<command.str()<<endl;

    int res  = system(command.str().c_str());
    assert(res == 0);
    
    
    double concorde_value;
  {
      //cerr<<"Done exec.\n";
      ifstream file;
      file.open("concorde.data");//{assert(0);cerr<<"Failing;\n";exit(0);}
      if(file.fail()){
          assert(0);cerr<<"Failing;\n";exit(0);
      } else if (file.eof()){
          assert(0);cerr<<"Failing;\n";exit(0);
      }
      string tmp;
      getline(file, tmp);
      file.close();
      istringstream iss(tmp);
      iss>>concorde_value;
  }

    double tour_length;
  {
      
    ifstream file;
    file.open("tour.output");
    if(file.fail()){
        assert(0);cerr<<"Failing;\n";exit(0);
    } else if (file.eof()){
        assert(0);cerr<<"Failing;\n";exit(0);
    }

    pair<double, vector<double> > parsed_tour = parse_tour(file, COST_OF_EDGE__FUNCTION);
    file.close();
    tour_length = parsed_tour.first;
  }
  
  // cerr<<"I calculate the tour length to be :: "<<tour_length<<endl
  //     <<"And concorde reports the your length to be :: "<<concorde_value<<endl;

    
  TSP__calculate_value_of__cache[local_points] = tour_length;//value;

  return tour_length;//value;
}


// /* Recursive function. Given an MST described in
//  * \argument{forward_tree}, generate a circuit which traverses that
//  * tree without visiting a node twice. On first call, tree is rooted
//  * at \argument{tree_node}. In subsequent calls, \argument{tree_node}
//  * is the current node of the tour being generated, and
//  * \argument{included_in_tour} contains all nodes visited in the tour
//  * so far. The list of nodes in the order visited by the generated
//  * tour is given in \argument{answer} (initially empty). */
// void make_circuit(
//     vector<int>& answer,
//     const map<int, set<int> >::const_iterator& tree_node,
//     const map<int, set<int> >& forward_tree, 
//     set<int>& included_in_tour)
// {
//     assert(included_in_tour.find(tree_node->first) == included_in_tour.end());

//     answer.push_back(tree_node->first) ;
//     included_in_tour.insert(tree_node->first);
    
//     if (tree_node->second.size() != 0 ){ // Has children
//         const set<int>& _children = tree_node->second;
//         vector<int> children(_children.begin(), _children.end());

//         while(children.size()){
//             int index = random() % children.size();
//             int child = children[index];
            
//             children[index] = children.back();
//             children.resize(children.size() - 1);
            
            
//             assert(forward_tree.find(child) != forward_tree.end());
//             if(forward_tree.find(child) == forward_tree.end()) {
//                 cerr<<"Unrecoverable error!\n";
//                 exit(0);
//             }

//             if(included_in_tour.find(forward_tree.find(child)->first) == included_in_tour.end()){
//                 map<int, set<int> >::const_iterator next_node = forward_tree.find(child);
//                 make_circuit(answer, 
//                              next_node,
//                              forward_tree,
//                              included_in_tour);
//             }
//         }
//     } 
// }

/* (see overloaded \function{make_circuit}) */
// vector<int> make_circuit(const map<int, set<int> >& forward_tree)
// {
//   set<int> included_in_tour;
//   vector<int> answer;

//   assert(forward_tree.size());
//   int index = random() % forward_tree.size();
  
//   map<int, set<int> >::const_iterator tree_node = forward_tree.begin();
//   for(int i = 0; i < index; i++){
//       tree_node++;
//   }

  
//   make_circuit(answer, tree_node, forward_tree, included_in_tour);
//   return answer;
// }


// map<set<int>, double> TMST__calculate_value_of__cache;

// /* Calculates a series of \global{MST_CIRCUIT_SAMPLES} tours from an
//  * MST over \argument{local_points}, and returns the resulting tour
//  * with the shortest distance.*/
// double TMST__calculate_value_of_matrix(const set<int>& local_points)
// {    
//     if(local_points.size() == 2){
//         set<int>::const_iterator second = local_points.begin()++;
//         set<int>::const_iterator first =  local_points.begin();
        
//         return COST_OF_EDGE__FUNCTION(*first, *second) * 2;
//     }


//     map<set<int>, double>::const_iterator cached_value = TMST__calculate_value_of__cache.find(local_points);
//     if( cached_value != TMST__calculate_value_of__cache.end()){
//         //cerr<<"Using cache.\n";
//         //cerr<< cached_value->second<<" -- ";
//         return  cached_value->second;
//     }


//     vector<int> local_game_points_index;
    
//     for(set<int>::const_iterator p = local_points.begin()
//             ; p != local_points.end()
//             ; p++){
//         local_game_points_index.push_back(*p);
//     }
    
//     using namespace boost;

//     typedef adjacency_list < vecS, vecS, undirectedS,
//                              no_property, property < edge_weight_t, double > > Graph;
//     typedef graph_traits < Graph >::edge_descriptor Edge;
//     typedef graph_traits < Graph >::vertex_descriptor Vertex;
//     typedef std::pair<int, int> E;

//     int number_of_edges = ceil( static_cast<double>(local_game_points_index.size() * local_game_points_index.size()) / 2.0);

  
//     //cerr<<"number of edges is :: "<<number_of_edges<<endl;

//     E edges[number_of_edges];
//     double ws[number_of_edges];
  
//     assert(local_game_points_index.size());

//     size_t count = 0;
//     for(size_t i = 0 ; i < local_game_points_index.size(); i++){
//         for(size_t j = i+1 ; j < local_game_points_index.size(); j++){
//             //cerr<<count<<endl;;
//             E tmp;
//             tmp.first = i;
//             tmp.second = j;
//             edges[count] = tmp;
          
//             double dist = calc_distance_matrix(local_game_points_index[i], local_game_points_index[j]);

//             assert( dist > 0.0);
//             ws[count] = dist;

//             count++;
//         }
//     }
    
//     map<int, set<int> > forward_tree;

//     Graph g(edges, edges + number_of_edges, ws, points.size());//, number_of_edges);
//     property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
//     std::vector < Edge > spanning_tree;
//     kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
//     double total_weight = 0;
//     assert(spanning_tree.size());
//     for (std::vector < Edge >::iterator ei = spanning_tree.begin();
//          ei != spanning_tree.end(); ++ei) 
//     {
//         total_weight += weight[*ei];
      
//         int ids = source(*ei, g);
//         int idt = target(*ei, g);

//         //cerr<<"tree from "<<ids<<" to "<<idt<<endl;

//         if(forward_tree.find(ids) == forward_tree.end()){forward_tree[ids]=set<int>();}
//         if(forward_tree.find(idt) == forward_tree.end()){forward_tree[idt]=set<int>();}
//         forward_tree[ids].insert(idt);
//         forward_tree[idt].insert(ids);
//     }

    
//     //sample for shortest circuit.
//     //int MST_CIRCUIT_SAMPLES = 500;
//     double current_best_weight = 1000000000;
//     for(int i = 0; i < MST_CIRCUIT_SAMPLES; i++){
//         vector<int> circuit = make_circuit(forward_tree); 
//         total_weight = 0;
//         total_weight += calc_distance_matrix(local_game_points_index[circuit[0]], local_game_points_index[circuit.back()], true);
//         for(size_t i = 0 ; i < circuit.size() - 1; i++){
//             total_weight += calc_distance_matrix(local_game_points_index[circuit[i]], local_game_points_index[circuit[i+1]], true);
//         }    
//         if(total_weight < current_best_weight){
//             current_best_weight = total_weight;
//         }
        
//     }

//     TMST__calculate_value_of__cache[local_points] = current_best_weight;
//     return current_best_weight;
// }

// /* Calculates a series of \global{MST_CIRCUIT_SAMPLES} tours from an
//  * MST over \argument{local_points}, and returns the resulting tour
//  * with the shortest distance.*/
// double TMST__calculate_value_of(const set<int>& local_points)
// {    
//     map<set<int>, double>::const_iterator cached_value = TMST__calculate_value_of__cache.find(local_points);
//     if( cached_value != TMST__calculate_value_of__cache.end()){
//         //cerr<<"Using cache.\n";
//         //cerr<< cached_value->second<<" -- ";
//         return  cached_value->second;
//     }

//     vector<Point> local_game_points;
//     vector<int> local_game_points_index;
    
//     for(set<int>::const_iterator p = local_points.begin()
//             ; p != local_points.end()
//             ; p++){
//         local_game_points.push_back(points[*p]);
//         local_game_points_index.push_back(*p);
//     }
    
//     using namespace boost;

//     typedef adjacency_list < vecS, vecS, undirectedS,
//                              no_property, property < edge_weight_t, double > > Graph;
//     typedef graph_traits < Graph >::edge_descriptor Edge;
//     typedef graph_traits < Graph >::vertex_descriptor Vertex;
//     typedef std::pair<int, int> E;

//     int number_of_edges = ceil( static_cast<double>(local_game_points.size() * local_game_points.size()) / 2.0);

  
//     //cerr<<"number of edges is :: "<<number_of_edges<<endl;

//     E edges[number_of_edges];
//     double ws[number_of_edges];
  
//     size_t count = 0;
//     for(size_t i = 0 ; i < local_game_points.size(); i++){
//         for(size_t j = i+1 ; j < local_game_points.size(); j++){
//             //cerr<<count<<endl;;
//             E tmp;
//             tmp.first = i;
//             tmp.second = j;
//             edges[count] = tmp;
          
//             double dist = calc_distance(local_game_points_index[i], local_game_points_index[j]);
          
//             //cerr<<"Distance from "<<i<<" to "<<j<<" is :"<<dist<<endl;

//             ws[count] = dist;

//             count++;
//         }
//     }
    
//     map<int, set<int> > forward_tree;

//     Graph g(edges, edges + number_of_edges, ws, points.size());//, number_of_edges);
//     property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
//     std::vector < Edge > spanning_tree;
//     kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
//     double total_weight = 0;
//     for (std::vector < Edge >::iterator ei = spanning_tree.begin();
//          ei != spanning_tree.end(); ++ei) 
//     {
//         total_weight += weight[*ei];
      
//         int ids = source(*ei, g);
//         int idt = target(*ei, g);

//         //cerr<<"tree from "<<ids<<" to "<<idt<<endl;

//         if(forward_tree.find(ids) == forward_tree.end()){forward_tree[ids]=set<int>();}
//         if(forward_tree.find(idt) == forward_tree.end()){forward_tree[idt]=set<int>();}
//         forward_tree[ids].insert(idt);
//         forward_tree[idt].insert(ids);
//     }

    
//     //sample for shortest circuit.
//     //int MST_CIRCUIT_SAMPLES = 500;
//     double current_best_weight = 1000000000;
//     for(int i = 0; i < MST_CIRCUIT_SAMPLES; i++){
//         vector<int> circuit = make_circuit(forward_tree); 
//         total_weight = 0;
//         total_weight += calc_distance(local_game_points_index[circuit[0]], local_game_points_index[circuit.back()]);
//         for(size_t i = 0 ; i < circuit.size() - 1; i++){
//             total_weight += calc_distance(local_game_points_index[circuit[i]], local_game_points_index[circuit[i+1]]);
//         }    
//         if(total_weight < current_best_weight){
//             current_best_weight = total_weight;
//         }
        
//     }

//     TMST__calculate_value_of__cache[local_points] = current_best_weight;
//     return current_best_weight;
// }


#include "asp.hh"


/* We have a multigraph data structure \type{Colored_Edges} to
 * generate Eulerian tour from MSTs according to Christofides' TSP
 * heuristic. Here, there can be multiple edges between two integer
 * index vertices. In order to distinguish between edges, we given
 * them an integer color.  */
class Colored_Edge{
public:
    int source;
    int destination;
    int color;

    bool operator==( const Colored_Edge& in) const{
        return source == in.source && destination == in.destination && color == in.color;
    }

    bool operator<( const Colored_Edge& in) const{
        if(source < in.source){
            return true;
        } else if (source == in.source){
            if(destination < in.destination){
                return true;
            } else if (destination == in.destination){
                if(color < in.color){
                    return true;
                }
            }
        }

        return false;
    }
};


typedef set<Colored_Edge> Colored_Edges; /* Multi-graph with colored edges.*/
typedef map<int,  Colored_Edges> Forward_Graph;/* For each vertex in a
                                                * multigraph, keeps
                                                * track of which
                                                * colored edges the
                                                * vertex participates
                                                * in.*/


/* Recursive function. Given a multigraph described in
 * \argument{forward_graph}, generate an Eulerian tour of that
 * graph. On first call, circuit begins at \argument{tree_node}. In
 * subsequent calls, \argument{tree_node} is the current node of the
 * tour being generated, and \argument{included_in_tour} contains all
 * nodes visited in the tour so far. Here,
 * \argument{edge_included_in_tour} is the list of edges included in
 * the tour so far. Algorithm is completed when all edges are
 * included. The list of nodes in the order visited by the generated
 * tour is given in \argument{answer} (initially empty). */
void make_circuit_Euler(
    vector<int>& answer,
    const Forward_Graph::const_iterator& graph_node,
    const Forward_Graph& forward_graph, 
    set<Colored_Edge>& edge_included_in_tour, 
    set<int>& node_included_in_tour)
{
    //cerr<<"Pushing back "<<graph_node->first<<endl;
    answer.push_back(graph_node->first) ;
    node_included_in_tour.insert(graph_node->first);

    /*Pick the first edge that is not already included in the Eulerian circuit.*/
    const set<Colored_Edge>& edges = graph_node->second;
    set<Colored_Edge>::const_iterator edge = edges.begin();
    for(
        ; edge != edges.end()
            ; edge++){

        //cerr<<"Trying edge :: "<<edge->source<<" to "<<edge->destination<<endl;

        if(edge_included_in_tour.find(*edge) == edge_included_in_tour.end()){
             
            edge_included_in_tour.insert(*edge);
    
            int selected_node = 0;
            if(edge->source == graph_node->first){
                selected_node = edge->destination;
            } else{
                selected_node = edge->source;
            }
            assert(forward_graph.find(selected_node) != forward_graph.end());

            const Forward_Graph::const_iterator next_node = forward_graph.find(selected_node);


            make_circuit_Euler(answer, 
                               next_node,
                               forward_graph,
                               edge_included_in_tour,
                               node_included_in_tour);
        }

         
    }
} 

/* (see overloaded \function{make_circuit_Euler}). */
vector<int> make_circuit_Euler(const Forward_Graph& forward_graph)
{
  set<Colored_Edge> edge_included_in_tour;
  set<int> node_included_in_tour;

  vector<int> answer;

  assert(forward_graph.size());
  Forward_Graph::const_iterator graph_node = forward_graph.begin();
  
  make_circuit_Euler(answer, graph_node, forward_graph, edge_included_in_tour, node_included_in_tour);
  
  set<int> seen;
  vector<int> shortcut_answer;
  for(size_t i = 0 ; i < answer.size(); i++){
      if (seen.find(answer[i]) == seen.end()){
          shortcut_answer.push_back(answer[i]);
          seen.insert(answer[i]);
      }
  }   

  assert(shortcut_answer.size() >= forward_graph.size());
  return shortcut_answer;
}


map<set<int>, double> Christofides__calculate_value_of__cache;


/* Calculates the value of the Christofides heuristic for the TSP over
 * the points \argument{local_points}.*/
double Christofides__calculate_value_of(const set<int>& local_points)
{    
    if(local_points.size() == 1){
        return 0;
    }

    if(local_points.size() == 2){
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return calc_distance(*first, *second) * 2;
    }


    map<set<int>, double>::const_iterator cached_value = Christofides__calculate_value_of__cache.find(local_points);
    if( cached_value != Christofides__calculate_value_of__cache.end()){
        return  cached_value->second;
    }

    // Insert caching mechanism here.

    vector<Point> local_game_points;
    vector<int> local_game_points_index;
    
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        local_game_points.push_back(points[*p]);
        local_game_points_index.push_back(*p);
    }
    
    using namespace boost;

    typedef adjacency_list < vecS, vecS, undirectedS,
                             no_property, property < edge_weight_t, double > > Graph;
    typedef graph_traits < Graph >::edge_descriptor Edge;
    typedef graph_traits < Graph >::vertex_descriptor Vertex;
    typedef std::pair<int, int> E;

    int number_of_edges = ceil( static_cast<double>(local_game_points.size() * local_game_points.size()) / 2.0);

  
    //cerr<<"number of edges is :: "<<number_of_edges<<endl;

    E edges[number_of_edges];
    double ws[number_of_edges];
  
    size_t count = 0;
    for(size_t i = 0 ; i < local_game_points.size(); i++){
        for(size_t j = i+1 ; j < local_game_points.size(); j++){
            //cerr<<count<<endl;;
            E tmp;
            tmp.first = i;
            tmp.second = j;
            edges[count] = tmp;
          
            double dist = COST_OF_EDGE__FUNCTION(local_game_points_index[i], local_game_points_index[j]);//calc_distance(local_game_points_index[i], local_game_points_index[j]);
          
            //cerr<<"Distance from "<<i<<" to "<<j<<" is :"<<dist<<endl;

            ws[count] = dist;

            count++;
        }
    }
    
    map<int, set<int> > forward_tree;

    Forward_Graph forward_graph;

    Graph g(edges, edges + number_of_edges, ws, points.size());//, number_of_edges);
    property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
    std::vector < Edge > spanning_tree;
    kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
         ei != spanning_tree.end(); ++ei) 
    {
      
        int ids = source(*ei, g);
        int idt = target(*ei, g);

        Colored_Edge colored_edge;
        colored_edge.source = ids;
        colored_edge.destination = idt;
        colored_edge.color = 0;


        //cerr<<"tree from "<<ids<<" to "<<idt<<endl;

        if(forward_tree.find(ids) == forward_tree.end()){forward_tree[ids]=set<int>();}
        if(forward_tree.find(idt) == forward_tree.end()){forward_tree[idt]=set<int>();}
        forward_tree[ids].insert(idt);
        forward_tree[idt].insert(ids);

        
        if(forward_graph.find(ids) == forward_graph.end()){forward_graph[ids]=Colored_Edges();}
        if(forward_graph.find(idt) == forward_graph.end()){forward_graph[idt]=Colored_Edges();}
        forward_graph[ids].insert(colored_edge);
        forward_graph[idt].insert(colored_edge);
    }

    set<int> odd_degree;
    vector<int> index_to_odd_degree;
    for( map<int, set<int> >:: const_iterator vertex = forward_tree.begin()
             ; vertex != forward_tree.end()
             ; vertex++){
        if(vertex->second.size() % 2){// Collect only odd degree vertices. 
            odd_degree.insert(vertex->first);
            index_to_odd_degree.push_back(vertex->first);
        }
    }
    
    /* The number of elements of odd degree must be even.*/
    assert(!(odd_degree.size() %2));

    int dim = odd_degree.size();
    double **Array = new double*[dim];
    for(int i =0 ; i < odd_degree.size(); i++){
        Array[i] =  new double[dim];
    }    

    /* Form a complete graph over odd degree vertices, and then
     * compute the minimum weight matching of that.*/
    int _i = 0;
    int _j = 0;
    for(set<int>::const_iterator i =  odd_degree.begin()
            ; i != odd_degree.end()
            ; i++, _i++){
        set<int>::const_iterator j = i;
        _j = _i;
        for(; j != odd_degree.end(); j++, _j++){

            assert(*i < local_game_points_index.size());
            assert(*j < local_game_points_index.size());
            
            if(i == j){
                Array[_i][_j] = 1000000;/*Big number.*/
                continue;
            }

            Array[_i][_j] = COST_OF_EDGE__FUNCTION(local_game_points_index[*i], local_game_points_index[*j]);
            Array[_j][_i] = COST_OF_EDGE__FUNCTION(local_game_points_index[*j], local_game_points_index[*i]);
        }
    }

    long *col_mate =  new long[dim];//(long*)malloc(sizeof(long) * odd_degree.size());;
    long *row_mate =  new long[dim];//(long*)malloc(sizeof(long) * odd_degree.size());;
    
    asp(odd_degree.size(), Array, col_mate, row_mate);

    
    for(int i =0 ; i < odd_degree.size() ; i++){
        assert(col_mate[i] < index_to_odd_degree.size());
        assert(i < index_to_odd_degree.size());


        int ids = index_to_odd_degree[i];
        int idt = index_to_odd_degree[col_mate[i]];

        
        Colored_Edge colored_edge;
        colored_edge.source = ids;
        colored_edge.destination = idt;
        colored_edge.color = 1;
        
        //cerr<<"ADDING  from "<<ids<<" to "<<idt<<endl;
        
        if(forward_graph.find(ids) == forward_graph.end()){forward_graph[ids]=Colored_Edges();}
        if(forward_graph.find(idt) == forward_graph.end()){forward_graph[idt]=Colored_Edges();}
        forward_graph[ids].insert(colored_edge);
        forward_graph[idt].insert(colored_edge);
    }
    

    vector<int> circuit = make_circuit_Euler(forward_graph);

    double total_weight = 0;
    total_weight += COST_OF_EDGE__FUNCTION(local_game_points_index[circuit[0]], local_game_points_index[circuit.back()]);
    //cerr<<"Euler circuit from  (number of points is "<<local_points.size()<<") :: \n";
    for(size_t i = 0 ; i < circuit.size() - 1; i++){
        //cerr<<local_game_points_index[circuit[i]]<<" to "<<local_game_points_index[circuit[i+1]]<<endl;
        total_weight += COST_OF_EDGE__FUNCTION(local_game_points_index[circuit[i]], local_game_points_index[circuit[i+1]]);
    }

    Christofides__calculate_value_of__cache[local_points] = total_weight;

    
    for(int i =0 ; i < odd_degree.size(); i++){
        delete[] Array[i];
    }    
    delete[] Array;
    delete[] col_mate;
    delete[] row_mate;
    
    return total_weight;

}

// /* Calculates the value of the Christofides heuristic for the TSP over
//  * the points \argument{local_points}.*/
// double Christofides__calculate_value_of(const set<int>& local_points)
// {   
//     if(local_points.size() == 1){
//         return 0;
//     }

//     if(local_points.size() == 2){
//         set<int>::const_iterator second = local_points.begin()++;
//         set<int>::const_iterator first =  local_points.begin();
        
//         return COST_OF_EDGE__FUNCTION(*first, *second) * 2;
//     }

//     if(local_points.size() == 3){
//         set<int>::const_iterator third = local_points.begin()++;third++;
//         set<int>::const_iterator second = local_points.begin()++;
//         set<int>::const_iterator first =  local_points.begin();
        
//         return COST_OF_EDGE__FUNCTION(*first, *second) + COST_OF_EDGE__FUNCTION(*second, *third) +  COST_OF_EDGE__FUNCTION(*third, *first) ;
//     }

//     map<set<int>, double>::const_iterator cached_value = Christofides__calculate_value_of__cache.find(local_points);
//     if( cached_value != Christofides__calculate_value_of__cache.end()){
//         return  cached_value->second;
//     }

//     // Insert caching mechanism here.

//     vector<Point> local_game_points;
//     vector<int> local_game_points_index;
    
//     for(set<int>::const_iterator p = local_points.begin()
//             ; p != local_points.end()
//             ; p++){
//         local_game_points.push_back(points[*p]);
//         local_game_points_index.push_back(*p);
//     }
    
//     using namespace boost;

//     typedef adjacency_list < vecS, vecS, undirectedS,
//                              no_property, property < edge_weight_t, double > > Graph;
//     typedef graph_traits < Graph >::edge_descriptor Edge;
//     typedef graph_traits < Graph >::vertex_descriptor Vertex;
//     typedef std::pair<int, int> E;

//     int number_of_edges = ceil( static_cast<double>(local_game_points.size() * local_game_points.size()) / 2.0);

  
//     //cerr<<"number of edges is :: "<<number_of_edges<<endl;

//     E edges[number_of_edges];
//     double ws[number_of_edges];
  
//     size_t count = 0;
//     for(size_t i = 0 ; i < local_game_points.size(); i++){
//         for(size_t j = i+1 ; j < local_game_points.size(); j++){
//             //cerr<<count<<endl;;
//             E tmp;
//             tmp.first = i;
//             tmp.second = j;
//             edges[count] = tmp;
          
//             double dist = COST_OF_EDGE__FUNCTION(local_game_points_index[i], local_game_points_index[j]);
          
//             //cerr<<"Distance from "<<i<<" to "<<j<<" is :"<<dist<<endl;

//             ws[count] = dist;

//             count++;
//         }
//     }
    
//     map<int, set<int> > forward_tree;

//     Forward_Graph forward_graph;

//     Graph g(edges, edges + number_of_edges, ws, points.size());//, number_of_edges);
//     property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
//     std::vector < Edge > spanning_tree;
//     kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

//     for (std::vector < Edge >::iterator ei = spanning_tree.begin();
//          ei != spanning_tree.end(); ++ei) 
//     {
      
//         int ids = source(*ei, g);
//         int idt = target(*ei, g);

//         Colored_Edge colored_edge;
//         colored_edge.source = ids;
//         colored_edge.destination = idt;
//         colored_edge.color = 0;


//         //cerr<<"tree from "<<ids<<" to "<<idt<<endl;

//         if(forward_tree.find(ids) == forward_tree.end()){forward_tree[ids]=set<int>();}
//         if(forward_tree.find(idt) == forward_tree.end()){forward_tree[idt]=set<int>();}
//         forward_tree[ids].insert(idt);
//         forward_tree[idt].insert(ids);

        
//         if(forward_graph.find(ids) == forward_graph.end()){forward_graph[ids]=Colored_Edges();}
//         if(forward_graph.find(idt) == forward_graph.end()){forward_graph[idt]=Colored_Edges();}
//         forward_graph[ids].insert(colored_edge);
//         forward_graph[idt].insert(colored_edge);
//     }

//     set<int> odd_degree;
//     vector<int> index_to_odd_degree;
//     for( map<int, set<int> >:: const_iterator vertex = forward_tree.begin()
//              ; vertex != forward_tree.end()
//              ; vertex++){
//         if(vertex->second.size() % 2){// Collect only odd degree vertices. 
//             odd_degree.insert(vertex->first);
//             index_to_odd_degree.push_back(vertex->first);
//         }
//     }
    
//     /* The number of elements of odd degree must be even.*/
//     assert(!(odd_degree.size() %2));

//     int dim = odd_degree.size();
//     double **Array = new double*[dim];
//     for(int i =0 ; i < odd_degree.size(); i++){
//         Array[i] =  new double[dim];
//     }    

//     /* Form a complete graph over odd degree vertices, and then
//      * compute the minimum weight matching of that.*/
//     int _i = 0;
//     int _j = 0;
//     for(set<int>::const_iterator i =  odd_degree.begin()
//             ; i != odd_degree.end()
//             ; i++, _i++){
//         set<int>::const_iterator j = i;
//         _j = _i;
//         for(; j != odd_degree.end(); j++, _j++){

//             assert(*i < local_game_points_index.size());
//             assert(*j < local_game_points_index.size());
            
//             if(i == j){
//                 Array[_i][_j] = 1000000;/*Big number.*/
//                 continue;
//             }

//             Array[_i][_j] = COST_OF_EDGE__FUNCTION(local_game_points_index[*i], local_game_points_index[*j]);
//             Array[_j][_i] = COST_OF_EDGE__FUNCTION(local_game_points_index[*j], local_game_points_index[*i]);
//         }
//     }

//     long *col_mate =  new long[dim];//(long*)malloc(sizeof(long) * odd_degree.size());;
//     long *row_mate =  new long[dim];//(long*)malloc(sizeof(long) * odd_degree.size());;
    
//     asp(odd_degree.size(), Array, col_mate, row_mate);

    
//     for(int i =0 ; i < odd_degree.size() ; i++){
//         assert(col_mate[i] < index_to_odd_degree.size());
//         assert(i < index_to_odd_degree.size());


//         int ids = index_to_odd_degree[i];
//         int idt = index_to_odd_degree[col_mate[i]];

        
//         Colored_Edge colored_edge;
//         colored_edge.source = ids;
//         colored_edge.destination = idt;
//         colored_edge.color = 1;
        
//         //cerr<<"ADDING  from "<<ids<<" to "<<idt<<endl;
        
//         if(forward_graph.find(ids) == forward_graph.end()){forward_graph[ids]=Colored_Edges();}
//         if(forward_graph.find(idt) == forward_graph.end()){forward_graph[idt]=Colored_Edges();}
//         forward_graph[ids].insert(colored_edge);
//         forward_graph[idt].insert(colored_edge);
//     }
    

//     vector<int> circuit = make_circuit_Euler(forward_graph);

//     double total_weight = 0;
//     total_weight += COST_OF_EDGE__FUNCTION(local_game_points_index[circuit[0]], local_game_points_index[circuit.back()]);
//     //cerr<<"Euler circuit from  (number of points is "<<local_points.size()<<") :: \n";
//     for(size_t i = 0 ; i < circuit.size() - 1; i++){
//         //cerr<<local_game_points_index[circuit[i]]<<" to "<<local_game_points_index[circuit[i+1]]<<endl;
//         total_weight += COST_OF_EDGE__FUNCTION(local_game_points_index[circuit[i]], local_game_points_index[circuit[i+1]]);
//     }

//     Christofides__calculate_value_of__cache[local_points] = total_weight;

    
//     for(int i =0 ; i < odd_degree.size(); i++){
//         delete[] Array[i];
//     }    
//     delete[] Array;
//     delete[] col_mate;
//     delete[] row_mate;
    
//     return total_weight;

// }



/* Calculates the value of the your which traverses the MST over the
 * points \argument{local_points}. That is 2 times the weight of the
 * MST.*/
double MST__calculate_value_of(const set<int>& local_points)
{
    if(local_points.size() == 1){
        return 0;
    }

    if(local_points.size() == 2){
        set<int>::const_iterator second = local_points.begin()++;
        set<int>::const_iterator first =  local_points.begin();
        
        return COST_OF_EDGE__FUNCTION(*first, *second);
    }

    vector<Point> local_game_points;
    vector<int> local_game_points_index;
    
    for(set<int>::const_iterator p = local_points.begin()
            ; p != local_points.end()
            ; p++){
        local_game_points.push_back(points[*p]);
        local_game_points_index.push_back(*p);
    }
   


    using namespace boost;

    typedef adjacency_list < vecS, vecS, undirectedS,
                             no_property, property < edge_weight_t, double > > Graph;
    typedef graph_traits < Graph >::edge_descriptor Edge;
    typedef graph_traits < Graph >::vertex_descriptor Vertex;
    typedef std::pair<int, int> E;

    int number_of_edges = ceil( static_cast<double>(local_game_points.size() * local_game_points.size()) / 2.0);

  
    //cerr<<"number of edges is :: "<<number_of_edges<<endl;

    E edges[number_of_edges];
    double ws[number_of_edges];
  
    size_t count = 0;
    for(size_t i = 0 ; i < local_game_points.size(); i++){
        for(size_t j = i+1 ; j < local_game_points.size(); j++){
            //cerr<<count<<endl;;
            E tmp;
            tmp.first = i;
            tmp.second = j;
            edges[count] = tmp;
          
            double dist = COST_OF_EDGE__FUNCTION(local_game_points_index[i], local_game_points_index[j]);
          
            //cerr<<"Distance from "<<i<<" to "<<j<<" is :"<<dist<<endl;

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

        //cerr<<"tree from "<<ids<<" to "<<idt<<endl;

        if(forward_tree.find(ids) == forward_tree.end()){forward_tree[ids]=set<int>();}
        if(forward_tree.find(idt) == forward_tree.end()){forward_tree[idt]=set<int>();}
        forward_tree[ids].insert(idt);
        forward_tree[idt].insert(ids);
    }


    return total_weight * 2;
  
}




vector<mpf_class> distance_from_depot_values()
{
    assert(are_we_solving_a_problem_with_a_single_depot);
    
    vector<mpf_class> answer(points.size());
    answer[0] = 0.0;
    for(uint i = 1 ; i < points.size(); i++){
        answer[i] = COST_OF_EDGE__FUNCTION(i, 0);
    }

    return answer;
}


vector<mpf_class> Classical_margin_values()
{
    assert(are_we_solving_a_problem_with_a_single_depot);

    set<int> all_points;
    for(uint i = 0 ; i < points.size(); i++){
        all_points.insert(i);
    }
    
    double grand_tour_value = TSP__calculate_value_of(all_points);

    vector<mpf_class> answer(points.size());
    answer[0] = 0.0;
    for(uint i = 1 ; i < points.size(); i++){
        
        all_points.erase(i);
        double local_tour_value = TSP__calculate_value_of(all_points);
        all_points.insert(i);

        //assert(grand_tour_value >= local_tour_value);

        answer[i] = grand_tour_value - local_tour_value;
    }

    return answer;
}

vector<mpf_class> Phil_Kilby_margin_values__indigo(){
    
    set<int> all_points;
    for(uint i = 0 ; i < points.size(); i++){
        all_points.insert(i);
    }

    WRITING_TO_CONCORDE__FUNCTION(all_points);
    ostringstream command;

    command<<"/home/cgretton/bin/indigo -i lns-3-10000-5 -r tour.output "
           <<concorde_file_name
           <<" 2> /dev/null";
    
    exit(0);

    cerr<<command.str()<<endl;

    int res  = system(command.str().c_str());
    
    assert(res == 0);
    
    ifstream file;
    file.open("tour.output");
    if(file.fail()){
        assert(0);cerr<<"Failing;\n";exit(0);
    } else if (file.eof()){
        assert(0);cerr<<"Failing;\n";exit(0);
    }


    pair<double, vector<double> > parsed_tour = parse_tour__indigo(file, COST_OF_EDGE__FUNCTION);
    
    file.close();
    
    double real_total =  parsed_tour.first;
    double actual_total = std::accumulate(parsed_tour.second.begin(), parsed_tour.second.end(), 0);
    mpf_class ratio = actual_total / real_total; 
    const vector<double>& margins = parsed_tour.second; 
    vector<mpf_class> answer(margins.begin(), margins.end());
    for(size_t i = 0; i < answer.size(); i++){
        answer[i] = answer[i] * ratio;
    }
    

    return answer;
}

vector<mpf_class> Phil_Kilby_margin_values(){
    
    set<int> all_points;
    for(uint i = 0 ; i < points.size(); i++){
        all_points.insert(i);
    }

    WRITING_TO_CONCORDE__FUNCTION(all_points);
    ostringstream command;
    command<<INSTALLED_PREFIX
           <<"/concorde/TSP/concorde  -o tour.output  "
           <<concorde_file_name
           <<" 2> /dev/null";


    int res  = system(command.str().c_str());
    
    assert(res == 0);
    
    ifstream file;
    file.open("tour.output");
    if(file.fail()){
        assert(0);cerr<<"Failing;\n";exit(0);
    } else if (file.eof()){
        assert(0);cerr<<"Failing;\n";exit(0);
    }


    pair<double, vector<double> > parsed_tour = parse_tour(file, COST_OF_EDGE__FUNCTION);
    
    file.close();
    
    double real_total =  parsed_tour.first;
    double actual_total = std::accumulate(parsed_tour.second.begin(), parsed_tour.second.end(), 0);
    mpf_class ratio = actual_total / real_total; 
    const vector<double>& margins = parsed_tour.second; 
    vector<mpf_class> answer(margins.begin(), margins.end());
    for(size_t i = 0; i < answer.size(); i++){
        answer[i] = answer[i] * ratio;
    }
    

    return answer;
}

// vector<mpf_class> Phil_Kilby_margin_values_matrix(){
    
//     set<int> all_points;
//     for(uint i = 0 ; i < points.size(); i++){
//         all_points.insert(i);
//     }

//     WRITING_TO_CONCORDE__FUNCTION(all_points);
//     ostringstream command;
//     command<<INSTALLED_PREFIX
//            <<"/concorde/TSP/concorde  -o tour.output  "
//            <<concorde_file_name
//            <<" 2> /dev/null";

//     int res  = system(command.str().c_str());
    
//     assert(res == 0);
    
//     ifstream file;
//     file.open("tour.output");
//     if(file.fail()){
//         assert(0);cerr<<"Failing;\n";exit(0);
//     } else if (file.eof()){
//         assert(0);cerr<<"Failing;\n";exit(0);
//     }


//     pair<double, vector<double> > parsed_tour = parse_tour(file, calc_distance_matrix__asym);
    
//     file.close();
    
//     double real_total =  parsed_tour.first;
//     double actual_total = std::accumulate(parsed_tour.second.begin(), parsed_tour.second.end(), 0);
//     mpf_class ratio = actual_total / real_total; 
//     const vector<double>& margins = parsed_tour.second; 
//     vector<mpf_class> answer(margins.begin(), margins.end());
//     for(size_t i = 0; i < answer.size(); i++){
//         answer[i] = answer[i] * ratio;
//     }
    

//     return answer;
// }
