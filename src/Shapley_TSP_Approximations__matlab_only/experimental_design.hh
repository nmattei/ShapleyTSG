
class ComparisonValues {
public:
    double tsp_value;
    double mst_value;
    double customer_margin_value;
    double moat_margin_value;
    double tour_margin_value;
    double depot_distance_value;
};

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
double _get_cpu_time(){ 
	struct rusage usage;
  	getrusage(RUSAGE_SELF, &usage);
  	return (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) * 
	(1e-6) + (usage.ru_utime.tv_sec + usage.ru_stime.tv_sec); 
}

// vector<ComparisonValues > calculate_all_shapley_values_matrix()
// {
    
//     cerr<<"MST Values are : \n";
//     vector<mpf_class> MST_tour_values 
//         = ApproShapley(matrix.size(), num_samples_approShapley,  TMST__calculate_value_of_matrix);// MST__calculate_value_of


//     vector<mpf_class> Phil_Kilby_tour_values = 
//         Phil_Kilby_margin_values_matrix();

//     cerr<<"TSP Values are : \n";
//     TSP__calculate_value_of__cache.clear();
//     TMST__calculate_value_of__cache.clear();
//     vector<mpf_class> TSP_tour_values
//         = ApproShapley(matrix.size(), num_samples_approShapley,  TSP__calculate_value_of_matrix);

//     assert(TSP_tour_values.size() == MST_tour_values.size());

//     vector<ComparisonValues > answer;
//     for(size_t i =0 ; i != TSP_tour_values.size(); i++){
//         ComparisonValues comparisonValues;
//         comparisonValues.tsp_value = 
//             TSP_tour_values[i].get_d();
//         comparisonValues.mst_value = 
//             MST_tour_values[i].get_d();
//         comparisonValues.customer_margin_value = 
//             Phil_Kilby_tour_values[i].get_d();
//         comparisonValues.moat_margin_value = 
//             moat_tour_values[i].get_d();
//         comparisonValues.tour_margin_value = 
//             classical_margin_values[i].get_d();
//         comparisonValues.depot_distance_value = 
//             Distance_from_Depot_tour_values[i].get_d();
//         answer.push_back(comparisonValues);//pair<double, double>(mst_val, tsp_val));
//     }

//     return answer;
// }

vector<ComparisonValues > calculate_all_shapley_values()
{
    cerr<<"Calculating all allocations (sample count set at "<<num_samples_approShapley<<").\n";

    Christofides__calculate_value_of__cache.clear();
    TSP__calculate_value_of__cache.clear();
    HELD_KARP__calculate_value_of__cache.clear();
    //Used to be doing TMST here.. //TMST__calculate_value_of__cache.clear();

    TSP__calculate_value_of = TSP__calculate_value_of__concorde;
    
    // TSP__calculate_value_of__cache = map<set<int>, double>();
    // HELD_KARP__calculate_value_of__cache = map<set<int>, double>();
    // TSP__calculate_value_of = TSP__calculate_value_of__indigo;
    // COST_OF_EDGE__FUNCTION = calc_distance_matrix__asym;
    // WRITING_TO_CONCORDE__FUNCTION = write__for_indigo__matrix;

    double start = _get_cpu_time();
    /*Classical margin values, this is the margin between optimal tours of size N and N-1.*/
    cerr<<"Classical tour margin values being calculated : \n";
    vector<mpf_class> classical_margin_values = 
        Classical_margin_values();

    double classical_margin__time =  _get_cpu_time() - start;

    
    // TSP__calculate_value_of__cache = map<set<int>, double>();
    // HELD_KARP__calculate_value_of__cache = map<set<int>, double>();
    // TSP__calculate_value_of = TSP__calculate_value_of__concorde;
    // COST_OF_EDGE__FUNCTION = calc_distance_matrix__sym;
    // WRITING_TO_CONCORDE__FUNCTION = write__for_concorde__matrix;
    

    start = _get_cpu_time();

    /*Difference from depot proxy values.*/
    cerr<<"Distance-from-depot heuristic values being calculated : \n";
    vector<mpf_class> Distance_from_Depot_tour_values  
        = distance_from_depot_values();

    double distance_from_depot__time =  _get_cpu_time() - start;
    start = _get_cpu_time();
    
    cerr<<"Christofides heuristic values being calculated : \n";
    vector<mpf_class> MST_tour_values 
      //= ApproShapley(points.size(), num_samples_approShapley, Christofides__calculate_value_of);
      = Distance_from_Depot_tour_values;
//        = ApproShapley(points.size(), num_samples_approShapley, HELD_KARP__calculate_value_of);

    double Christofides__time =  _get_cpu_time() - start;
    start = _get_cpu_time();


    
    
    // TSP__calculate_value_of__cache = map<set<int>, double>();
    // HELD_KARP__calculate_value_of__cache = map<set<int>, double>();
    // TSP__calculate_value_of = TSP__calculate_value_of__indigo;
    // COST_OF_EDGE__FUNCTION = calc_distance_matrix__asym;
    // WRITING_TO_CONCORDE__FUNCTION = write__for_indigo__matrix;


    cerr<<"Phil Kilby heuristic values being calculated : \n";
    vector<mpf_class> Phil_Kilby_tour_values 
      = Phil_Kilby_margin_values();
    //vector<mpf_class> Phil_Kilby_tour_values = Phil_Kilby_margin_values__indigo();

    double  Phil_Kilby_tour__time =  _get_cpu_time() - start;
    start = _get_cpu_time();


    
    // TSP__calculate_value_of__cache = map<set<int>, double>();
    // HELD_KARP__calculate_value_of__cache = map<set<int>, double>();
    // TSP__calculate_value_of = TSP__calculate_value_of__concorde;
    // COST_OF_EDGE__FUNCTION = calc_distance_matrix__sym;
    // WRITING_TO_CONCORDE__FUNCTION = write__for_concorde__matrix;
    
    cerr<<"Moat-Packing heuristic values being calculated : \n";
    vector<mpf_class>  moat_tour_values  
        //= MST_tour_values ;
        = moat_margin_values();


    double  Moat__time =  _get_cpu_time() - start;
    start = _get_cpu_time();

    cerr<<"TSP Values are : \n";
    tspmode = true;
    vector<mpf_class> TSP_tour_values 
      = Distance_from_Depot_tour_values;
//        = MST_tour_values ;
//        = ApproShapley(points.size(), num_samples_approShapley,  TSP__calculate_value_of);
    tspmode = false;

    double  TSP_Shapley__time =  _get_cpu_time() - start;
    

    cout<<"Times are :: \n "
        <<TSP_Shapley__time<<endl
        <<Christofides__time<<endl
        <<Phil_Kilby_tour__time<<endl
        <<Moat__time<<endl
        <<classical_margin__time<<endl
        <<distance_from_depot__time<<endl<<endl;

    assert(TSP_tour_values.size() == MST_tour_values.size());

    assert(TSP_tour_values.size());

    vector<ComparisonValues > answer;
    for(size_t i =0 ; i != TSP_tour_values.size(); i++){
        ComparisonValues comparisonValues;
        comparisonValues.tsp_value = 
            TSP_tour_values[i].get_d();
        comparisonValues.mst_value = 
            MST_tour_values[i].get_d();
        comparisonValues.customer_margin_value = 
            Phil_Kilby_tour_values[i].get_d();
        comparisonValues.moat_margin_value = 
            moat_tour_values[i].get_d();
        comparisonValues.tour_margin_value = 
            classical_margin_values[i].get_d();
        comparisonValues.depot_distance_value = 
            Distance_from_Depot_tour_values[i].get_d();
        answer.push_back(comparisonValues);//pair<double, double>(mst_val, tsp_val));
    }

    return answer;
}


vector<vector<ComparisonValues > > shapley_samples;
vector<vector<int > > removal_prescription_values;
vector<vector<int > > group_removal_prescription_values;

void one_point_experiment__Euclidean()
{
    vector<ComparisonValues > shapley_values 
        = calculate_all_shapley_values();
            
    shapley_samples.push_back(vector<ComparisonValues >());
    vector<ComparisonValues >& existing_values 
        = shapley_samples.back();//.push_back(shapley_value);

    for(size_t k =0; k < shapley_values.size(); k++)
    {
        existing_values.push_back(shapley_values[k]);
    }
}
void one_point_experiment__matrix()
{
    vector<ComparisonValues > shapley_values 
        = calculate_all_shapley_values();
            
    shapley_samples.push_back(vector<ComparisonValues >());
    vector<ComparisonValues >& existing_values 
        = shapley_samples.back();//.push_back(shapley_value);

    for(size_t k =0; k < shapley_values.size(); k++)
    {
        existing_values.push_back(shapley_values[k]);
    }
}

// void main_loop_matrix(){
    
//     for(size_t i = lower_range ; i <= upper_range; i++){
//         points_count = i;
//         shapley_samples.push_back(vector<ComparisonValues >());
//         for(size_t j = 0 ; j < num_samples; j++){
//             generate_model_matrix();

//             vector<ComparisonValues > shapley_values 
//                 = calculate_all_shapley_values_matrix();
            
//             vector<ComparisonValues >& existing_values 
//                 = shapley_samples.back();//.push_back(shapley_value);

//             for(size_t k =0; k < shapley_values.size(); k++)
//             {
//                 existing_values.push_back(shapley_values[k]);
//             }
//         }
//     }
// }

#define vector_to_set_expansion(data,X,Y) {     \
        if(data.find(X) == data.end()){         \
            data[X] = set<int>();               \
        }                                       \
        data[X].insert(Y);                      \
    }                                           \
 

double scenario_value(const set<int>& game, int index, double(*calculate_value_of)(const set<int> &)){
    set<int> _game = game;
    _game.erase(index);
    if(!_game.size()){
        return 0;
    }
    if(_game.size() == 1){
        return 0;
    }
    return calculate_value_of(_game);
}

/* For the transport \argument{scenario} where players have cost
 * allocations from \argument{prescriptions}, how many times (integer
 * return value) is the cost-allocation prescriptive, in terms of
 * telling us what customer to remove to get the most profitable
 * scenario with one less player.*/
int count_correct_prescriptions(map<double, set<int> >& prescriptions, 
                                set<int>& scenario, 
                                double(*calculate_value_of)(const set<int> &))
{
    if(!scenario.size()) return 0;

    int answer = 0 ;

    double best_prescription_value = 10000000000;
    int best_prescription_index = 0;

    /* This loop finds the most profitable customer to delete. */
    for(set<int>::const_iterator p = scenario.begin()
            ; p != scenario.end()
            ; p++){
        double value = scenario_value(scenario, *p, calculate_value_of);
        if(value < best_prescription_value){
            best_prescription_value = value;
            best_prescription_index = *p;
        }
    }

    assert(prescriptions.size());
    assert((prescriptions.end()--)->first >= prescriptions.begin()->first);

    map<double, set<int> >::iterator p = --(prescriptions.end());
    set<int> prescribed_points = p->second;
    double prescribed_index = p->first;
    
    /* This condition checks if the cost allocation procedure is prescribing the optimal customer for deletion.*/
    if(prescribed_points.find(best_prescription_index) != prescribed_points.end()){
        prescribed_points.erase(best_prescription_index);// This prescription was correct, and has now been made.
        if(!prescribed_points.size()){// If there are no more prescribed points with the $p->first$ value...
            prescriptions.erase(p);
        } else {
            prescriptions[prescribed_index] = prescribed_points;
        }

        scenario.erase(best_prescription_index);
        
        answer += 1;
        
        int add_to_answer = count_correct_prescriptions(prescriptions, scenario, calculate_value_of);

        answer += add_to_answer;
    }
    
    return answer;
}

/* Second argument is a function that calculates a solution for the
 * input TSP (should be optimal). First argument is a set of
 * prescriptions, for each unique allocation-value to a player, we
 * have the set of players associated with that value.*/
int count_correct_prescriptions(const map<double, set<int> >& prescriptions, 
                                double(*calculate_value_of)(const set<int> &))
{
    
    set<int> initial_scenario;
    for(size_t i = 0 ; i < points.size(); i++){
        initial_scenario.insert(i);
    }
    
    map<double, set<int> > _prescriptions = prescriptions;
    return count_correct_prescriptions(_prescriptions, initial_scenario, calculate_value_of);
}


/* This function does three things. First, it uses the cost
 * allocations as a policy, evaluates the consequence of deleting
 * those customers, and then tests if those deleted customers were
 * actually optimal to delete. */
void prescribe_test_store(const vector<ComparisonValues >& shapley_values)
{
    
    map<double, set<int> > mst;
    map<double, set<int> > tsp;
    map<double, set<int> > phil;
    map<double, set<int> > moat;
    map<double, set<int> > tour_margin;
    map<double, set<int> > depot_dist;
    map<int, double> _mst;
    map<int, double> _tsp;
    map<int, double> _phil;
    map<int, double> _moat;
    map<int, double> _tour_margin;
    map<int, double> _depot_dist;

    size_t index = 0; 
    while(index < shapley_values.size()){
        vector_to_set_expansion(mst, shapley_values[index].mst_value, index);
        _mst[index] = shapley_values[index].mst_value;
        vector_to_set_expansion(tsp, shapley_values[index].tsp_value, index);
        _tsp[index] = shapley_values[index].tsp_value;
        vector_to_set_expansion(phil, shapley_values[index].customer_margin_value, index);
        _phil[index] = shapley_values[index].customer_margin_value;
        vector_to_set_expansion(moat, shapley_values[index].moat_margin_value, index);
        _moat[index] = shapley_values[index].moat_margin_value;
        vector_to_set_expansion(tour_margin, shapley_values[index].tour_margin_value, index);
        _tour_margin[index] = shapley_values[index].tour_margin_value;
        vector_to_set_expansion(depot_dist, shapley_values[index].depot_distance_value, index);
        _depot_dist[index] = shapley_values[index].depot_distance_value;

        index++;
    }

    cerr<<"Counting prescriptions TSP.\n";
    int tsp_prescriptions = count_correct_prescriptions(tsp, TSP__calculate_value_of);

    cerr<<"Counting prescriptions MST.\n";
    int mst_prescriptions = count_correct_prescriptions(mst, TSP__calculate_value_of);
    
    cerr<<"Counting prescriptions PHIL.\n";
    int phil_prescriptions = count_correct_prescriptions(phil, TSP__calculate_value_of);

    cerr<<"Counting prescriptions MOAT.\n";
    int moat_prescriptions = count_correct_prescriptions(moat, TSP__calculate_value_of);

    cerr<<"Counting prescriptions TOUR_MARGIN.\n";
    int tour_margin_prescriptions = count_correct_prescriptions(tour_margin, TSP__calculate_value_of);

    cerr<<"Counting prescriptions DEPOT_DIST.\n";
    int depot_dist_prescriptions = count_correct_prescriptions(depot_dist, TSP__calculate_value_of);

    vector<int> values;
    values.push_back(tsp_prescriptions);
    values.push_back(mst_prescriptions);
    values.push_back(phil_prescriptions);
    values.push_back(moat_prescriptions);
    values.push_back(tour_margin_prescriptions);
    values.push_back(depot_dist_prescriptions);
    
    removal_prescription_values.push_back(values);
}

//first is sum of cost allocations to members of \argument{evaluation_group} 
//second is cost of the transport scenario without elements from \argument{evaluation_group} 
pair<double, double>
group_prescribe_test_store(map<int, double>& assignments,
                           set<int> all_players,
                           set<int> evaluation_group)
{
    pair<double, double> answer;
    
    double value_of_group = 0;
    for(set<int>::iterator p =  evaluation_group.begin()
            ; p != evaluation_group.end()
            ; p++){
        all_players.erase(*p);

        value_of_group += assignments[*p]; 
    }

    double scenario_value = TSP__calculate_value_of(all_players);
    

    return pair<double, double>(value_of_group, scenario_value);
}

/* Group the players into 3 groups. Then test if the cost
 * allocation is prescriptive, in terms of telling you which group to
 * delete to get a transport scenario with the lowest cost. */
void group_prescribe_test_store(const vector<ComparisonValues >& shapley_values)
{
    /*This experiments fixes the number of test groups to 3. */
    vector<set<int> > groups(3);
    size_t num_groups = groups.size();
    for(size_t i =1 ; i< points.size(); i++){//do not count the depot.
        size_t index = random()%num_groups;
        groups[index].insert(i);
    }
    vector<double > groups_values(3);

    map<double, set<int> > mst;
    map<double, set<int> > tsp;
    map<double, set<int> > phil;
    map<double, set<int> > moat;
    map<double, set<int> > tour_margin;
    map<double, set<int> > depot_dist;
    map<int, double> _mst;
    map<int, double> _tsp;
    map<int, double> _phil;
    map<int, double> _moat;
    map<int, double> _tour_margin;
    map<int, double> _depot_dist;

    set<int> all_players;
    size_t index = 0; 
    while(index < shapley_values.size()){
        all_players.insert(index);

        vector_to_set_expansion(mst, shapley_values[index].mst_value, index);
        _mst[index] = shapley_values[index].mst_value;
        vector_to_set_expansion(tsp, shapley_values[index].tsp_value, index);
        _tsp[index] = shapley_values[index].tsp_value;
        vector_to_set_expansion(phil, shapley_values[index].customer_margin_value, index);
        _phil[index] = shapley_values[index].customer_margin_value;
        vector_to_set_expansion(moat, shapley_values[index].moat_margin_value, index);
        _moat[index] = shapley_values[index].moat_margin_value;
        vector_to_set_expansion(tour_margin, shapley_values[index].tour_margin_value, index);
        _tour_margin[index] = shapley_values[index].tour_margin_value;
        vector_to_set_expansion(depot_dist, shapley_values[index].depot_distance_value, index);
        _depot_dist[index] = shapley_values[index].depot_distance_value;
        index++;
    }

//first is sum of allocations
//second is score of game
    pair<double, double> tsp_pair_G1 = group_prescribe_test_store(_tsp, all_players, groups[0]);
    pair<double, double> mst_pair_G1 = group_prescribe_test_store(_mst, all_players, groups[0]);
    pair<double, double> phil_pair_G1  = group_prescribe_test_store(_phil, all_players, groups[0]);
    pair<double, double> moat_pair_G1  = group_prescribe_test_store(_moat, all_players, groups[0]);
    pair<double, double> tour_margin_pair_G1  = group_prescribe_test_store(_tour_margin, all_players, groups[0]);
    pair<double, double> depot_dist_pair_G1  = group_prescribe_test_store(_depot_dist, all_players, groups[0]);
    
    pair<double, double> tsp_pair_G2 = group_prescribe_test_store(_tsp, all_players, groups[1]);
    pair<double, double> mst_pair_G2 = group_prescribe_test_store(_mst, all_players, groups[1]);
    pair<double, double> phil_pair_G2  = group_prescribe_test_store(_phil, all_players, groups[1]);
    pair<double, double> moat_pair_G2  = group_prescribe_test_store(_moat, all_players, groups[1]);
    pair<double, double> tour_margin_pair_G2  = group_prescribe_test_store(_tour_margin, all_players, groups[1]);
    pair<double, double> depot_dist_pair_G2  = group_prescribe_test_store(_depot_dist, all_players, groups[1]);
    
    pair<double, double> tsp_pair_G3 = group_prescribe_test_store(_tsp, all_players, groups[2]);
    pair<double, double> mst_pair_G3 = group_prescribe_test_store(_mst, all_players, groups[2]);
    pair<double, double> phil_pair_G3  = group_prescribe_test_store(_phil, all_players, groups[2]);
    pair<double, double> moat_pair_G3  = group_prescribe_test_store(_moat, all_players, groups[2]);
    pair<double, double> tour_margin_pair_G3  = group_prescribe_test_store(_tour_margin, all_players, groups[2]);
    pair<double, double> depot_dist_pair_G3  = group_prescribe_test_store(_depot_dist, all_players, groups[2]);

    map<double, double> TSP_SCENARIO;
    map<double, double> MST_SCENARIO;
    map<double, double> PHIL_SCENARIO;
    map<double, double> MOAT_SCENARIO;
    map<double, double> TOUR_MARGIN_SCENARIO;
    map<double, double> DEPOT_DIST_SCENARIO;

    map<double, double> _TSP_SCENARIO;
    map<double, double> _MST_SCENARIO;
    map<double, double> _PHIL_SCENARIO;
    map<double, double> _MOAT_SCENARIO;
    map<double, double> _TOUR_MARGIN_SCENARIO;
    map<double, double> _DEPOT_DIST_SCENARIO;
    
    TSP_SCENARIO[tsp_pair_G1.first] = tsp_pair_G1.second;
    TSP_SCENARIO[tsp_pair_G2.first] = tsp_pair_G2.second;
    TSP_SCENARIO[tsp_pair_G3.first] = tsp_pair_G3.second;

    MST_SCENARIO[mst_pair_G1.first] = mst_pair_G1.second;
    MST_SCENARIO[mst_pair_G2.first] = mst_pair_G2.second;
    MST_SCENARIO[mst_pair_G3.first] = mst_pair_G3.second;

    PHIL_SCENARIO[phil_pair_G1.first] = phil_pair_G1.second;
    PHIL_SCENARIO[phil_pair_G2.first] = phil_pair_G2.second;
    PHIL_SCENARIO[phil_pair_G3.first] = phil_pair_G3.second;

    MOAT_SCENARIO[moat_pair_G1.first] = moat_pair_G1.second;
    MOAT_SCENARIO[moat_pair_G2.first] = moat_pair_G2.second;
    MOAT_SCENARIO[moat_pair_G3.first] = moat_pair_G3.second;

    TOUR_MARGIN_SCENARIO[tour_margin_pair_G1.first] = tour_margin_pair_G1.second;
    TOUR_MARGIN_SCENARIO[tour_margin_pair_G2.first] = tour_margin_pair_G2.second;
    TOUR_MARGIN_SCENARIO[tour_margin_pair_G3.first] = tour_margin_pair_G3.second;

    DEPOT_DIST_SCENARIO[depot_dist_pair_G1.first] = depot_dist_pair_G1.second;
    DEPOT_DIST_SCENARIO[depot_dist_pair_G2.first] = depot_dist_pair_G2.second;
    DEPOT_DIST_SCENARIO[depot_dist_pair_G3.first] = depot_dist_pair_G3.second;

    _TSP_SCENARIO[tsp_pair_G1.second] = tsp_pair_G1.first;
    _TSP_SCENARIO[tsp_pair_G2.second] = tsp_pair_G2.first;
    _TSP_SCENARIO[tsp_pair_G3.second] = tsp_pair_G3.first;

    _MST_SCENARIO[mst_pair_G1.second] = mst_pair_G1.first;
    _MST_SCENARIO[mst_pair_G2.second] = mst_pair_G2.first;
    _MST_SCENARIO[mst_pair_G3.second] = mst_pair_G3.first;

    _PHIL_SCENARIO[phil_pair_G1.second] = phil_pair_G1.first;
    _PHIL_SCENARIO[phil_pair_G2.second] = phil_pair_G2.first;
    _PHIL_SCENARIO[phil_pair_G3.second] = phil_pair_G3.first;

    _MOAT_SCENARIO[moat_pair_G1.second] = moat_pair_G1.first;
    _MOAT_SCENARIO[moat_pair_G2.second] = moat_pair_G2.first;
    _MOAT_SCENARIO[moat_pair_G3.second] = moat_pair_G3.first;

    _TOUR_MARGIN_SCENARIO[tour_margin_pair_G1.second] = tour_margin_pair_G1.first;
    _TOUR_MARGIN_SCENARIO[tour_margin_pair_G2.second] = tour_margin_pair_G2.first;
    _TOUR_MARGIN_SCENARIO[tour_margin_pair_G3.second] = tour_margin_pair_G3.first;

    _DEPOT_DIST_SCENARIO[depot_dist_pair_G1.second] = depot_dist_pair_G1.first;
    _DEPOT_DIST_SCENARIO[depot_dist_pair_G2.second] = depot_dist_pair_G2.first;
    _DEPOT_DIST_SCENARIO[depot_dist_pair_G3.second] = depot_dist_pair_G3.first;


    bool tsp_is_correct =false;
    //Lowest game value to highest allocation of cost.
    if(_TSP_SCENARIO.begin()->second == (--(TSP_SCENARIO.end()))->first){
        tsp_is_correct = true;
    } 
    bool mst_is_correct =false;
    //Lowest game value to highest allocation of cost.
    if(_MST_SCENARIO.begin()->second == (--(MST_SCENARIO.end()))->first){
        mst_is_correct = true;
    } 
    bool phil_is_correct =false;
    //Lowest game value to highest allocation of cost.
    if(_PHIL_SCENARIO.begin()->second == (--(PHIL_SCENARIO.end()))->first){
        phil_is_correct = true;
    } 
    bool moat_is_correct =false;
    //Lowest game value to highest allocation of cost.
    if(_MOAT_SCENARIO.begin()->second == (--(MOAT_SCENARIO.end()))->first){
        moat_is_correct = true;
    } 
    bool tour_margin_is_correct =false;
    //Lowest game value to highest allocation of cost.
    if(_TOUR_MARGIN_SCENARIO.begin()->second == (--(TOUR_MARGIN_SCENARIO.end()))->first){
        tour_margin_is_correct = true;
    } 
    bool depot_dist_is_correct =false;
    //Lowest game value to highest allocation of cost.
    if(_DEPOT_DIST_SCENARIO.begin()->second == (--(DEPOT_DIST_SCENARIO.end()))->first){
        depot_dist_is_correct = true;
    } 

    vector<int > results;
    results.push_back(tsp_is_correct);
    results.push_back(mst_is_correct);
    results.push_back(phil_is_correct);
    results.push_back(moat_is_correct);
    results.push_back(tour_margin_is_correct);
    results.push_back(depot_dist_is_correct);

    group_removal_prescription_values.push_back(results);
}


void main_loop(){
    
    for(size_t i = lower_range ; i <= upper_range; i++){
        points_count = i;
        shapley_samples.push_back(vector<ComparisonValues >());
        for(size_t j = 0 ; j < num_samples; j++){
            generate_model();
            ofstream file;
            ostringstream oss;
            oss<<"euctsp_size="<<i<<"_sample="<<j<<".tsp";
            file.open (oss.str().c_str());
            set<int> all_points;
            for(uint i = 0 ; i < points.size(); i++){
                all_points.insert(i);
            }
            write__for_concorde(file, all_points);//local_points);
            file.close();
            // exit(0);
            vector<ComparisonValues > shapley_values 
                = calculate_all_shapley_values();
            assert(shapley_values.size());
            vector<ComparisonValues >& existing_values 
                = shapley_samples.back();//.push_back(shapley_value);



            for(size_t k =0; k < shapley_values.size(); k++)
            {
                existing_values.push_back(shapley_values[k]);
            }
            
            // /* Group customers by type, prescribe removal of type, test
            //  * profitability of removal of type, report that
            //  * statistic. */
            cerr<<"Checking individual prescriptions\n";
            prescribe_test_store(shapley_values);
            cerr<<"Checking group prescriptions\n";
            group_prescribe_test_store(shapley_values);
        }
    }
}

void write_raw_samples(ofstream& file){
    size_t size = lower_range;
    for(vector<vector<ComparisonValues > >::iterator p = shapley_samples.begin()
            ; p != shapley_samples.end()
            ; p++){

        vector<ComparisonValues >& existing_values  = *p;
        
        for (size_t i = 0 ; i < existing_values.size(); i++){
            file<<size<<" "<<existing_values[i].tsp_value
                <<" "<<existing_values[i].mst_value
                <<" "<<existing_values[i].customer_margin_value
                <<" "<<existing_values[i].moat_margin_value
                <<" "<<existing_values[i].tour_margin_value
                <<" "<<existing_values[i].depot_distance_value<<endl;
        }
        size++;
    }
}

void write_raw_samples(){
    
    ostringstream oss;
    oss<<"raw_samples_"<<lower_range<<"_"<<upper_range<<".data";

    ofstream file;
    file.open (oss.str().c_str());
    write_raw_samples(file);
    file.close();

}

void write_ratio_data(ofstream& file){
    
    size_t size = lower_range;
    for(vector<vector<ComparisonValues > >::iterator p = shapley_samples.begin()
            ; p != shapley_samples.end()
            ; p++){

        vector<ComparisonValues >& existing_values  = *p;
        for (size_t i = 0 ; i < existing_values.size(); i++){
            file<<size<<" "<<existing_values[i].mst_value       / existing_values[i].tsp_value
                <<" "<<existing_values[i].customer_margin_value / existing_values[i].tsp_value
                <<" "<<existing_values[i].moat_margin_value / existing_values[i].tsp_value
                <<" "<<existing_values[i].tour_margin_value / existing_values[i].tsp_value
                <<" "<<existing_values[i].depot_distance_value / existing_values[i].tsp_value<<endl;
        }
        size++;
    }
}

void write_ratio_data(){
    
    ostringstream oss;
    oss<<"ratio_data_"<<lower_range<<"_"<<upper_range<<".data";

    ofstream file;
    file.open (oss.str().c_str());
    write_ratio_data(file);
    file.close();
}

// void write_ranking_data(ofstream& file){
    
//     size_t size = lower_range;
//     for(vector<vector<ComparisonValues > >::iterator p = shapley_samples.begin()
//             ; p != shapley_samples.end()
//             ; p++){

//         vector<ComparisonValues >& existing_values  = *p;
//         size_t index = 0; 
        
//         while(index < existing_values.size()){
//             map<double, int> mst;
//             map<double, int> tsp;
//             map<double, int> phil;
//             map<double, int> moat;
//             map<int, double> _mst;
//             map<int, double> _tsp;
//             map<int, double> _phil;
//             map<int, double> _moat;
//             for(size_t i =0; i < size; i++,index++){
//                 mst[existing_values[index].mst_value] = i;
//                 tsp[existing_values[index].tsp_value] = i;
//                 phil[existing_values[index].customer_margin_value] = i;
//                 moat[existing_values[index].moat_margin_value] = i;
//                 _mst[i] = existing_values[index].mst_value;
//                 _tsp[i] = existing_values[index].tsp_value;
//                 _phil[i] = existing_values[index].customer_margin_value;
//                 _moat[i] = existing_values[index].moat_margin_value;
//             }

//             int sum_of_mst_rank_errors = 0;
//             int sum_of_phil_rank_errors = 0;
//             int sum_of_moat_rank_errors = 0;
//             int ranking = 0;
//             for(map<double, int>::iterator q = tsp.begin(); q != tsp.end(); q++, ranking++){
                
//                 double mst_value = _mst[q->second];
//                 map<double, int>::iterator iter = mst.find(mst_value);
//                 int mst_error = 0; 
//                 for(; iter != mst.begin(); iter--){
//                     mst_error++;
//                 }
//                 sum_of_mst_rank_errors+=mst_error;

//                 double phil_value = _phil[q->second];
//                 iter = phil.find(phil_value);
//                 int phil_error = 0; 
//                 for(; iter != phil.begin(); iter--){
//                     phil_error++;
//                 }
//                 sum_of_phil_rank_errors+=phil_error;

                
//                 double moat_value = _moat[q->second];
//                 iter = moat.find(moat_value);

//                 int moat_error = 0; 
//                 for(; iter != moat.begin(); iter--){
//                     moat_error++;
//                 }
//                 sum_of_moat_rank_errors+=moat_error;
//             }
//             file<<size<<" "<<sum_of_mst_rank_errors / tsp.size()<<" "<<sum_of_phil_rank_errors / tsp.size()
//                 <<" "<<sum_of_moat_rank_errors / tsp.size()<<endl;
//         }
        
//         size++;
//     }
// }

// void write_ranking_data(){
//     ostringstream oss;
//     oss<<"ranking_data_"<<lower_range<<"_"<<upper_range<<".data";

//     ofstream file;
//     file.open (oss.str().c_str());
//     write_ranking_data(file);
//     file.close();
// }



void write_apportionment(ofstream& file){
    
    size_t size = lower_range;
    for(vector<vector<ComparisonValues > >::iterator p = shapley_samples.begin()
            ; p != shapley_samples.end()
            ; p++){

        vector<ComparisonValues >& existing_values  = *p;
        size_t index = 0; 
        
        while(index < existing_values.size()){
            map<double, int> mst;
            map<double, int> tsp;
            map<double, int> phil;
            map<double, int> moat;
            map<double, int> tour_margin;
            map<double, int> depot_distance;

            map<int, double> _mst;
            map<int, double> _tsp;
            map<int, double> _phil;
            map<int, double> _moat;
            map<int, double> _tour_margin;
            map<int, double> _depot_distance;

            double tsp_total = 0;
            double mst_total = 0;
            double phil_total = 0;
            double moat_total = 0;
            double tour_margin_total = 0;
            double depot_distance_total = 0;

            for(size_t i =0; i < size; i++,index++){
                mst[existing_values[index].mst_value] = i;
                tsp[existing_values[index].tsp_value] = i;
                phil[existing_values[index].customer_margin_value] = i;
                moat[existing_values[index].moat_margin_value] = i;
                tour_margin[existing_values[index].tour_margin_value] = i;
                depot_distance[existing_values[index].depot_distance_value] = i;

                _mst[i] = existing_values[index].mst_value;
                _tsp[i] = existing_values[index].tsp_value;
                _moat[i] = existing_values[index].moat_margin_value;
                _tour_margin[i] = existing_values[index].tour_margin_value;
                _depot_distance[i] = existing_values[index].depot_distance_value;

                tsp_total+=existing_values[index].tsp_value;
                mst_total+=existing_values[index].mst_value;
                phil_total+=existing_values[index].customer_margin_value;
                moat_total+=existing_values[index].moat_margin_value;
                tour_margin_total+=existing_values[index].tour_margin_value;
                depot_distance_total+=existing_values[index].depot_distance_value;
            }

            
            for(size_t i =0; i < size; i++){
                _mst[i] /= mst_total;
                _tsp[i] /= tsp_total;
                _phil[i] /= phil_total;
                _moat[i] /= moat_total;
                _tour_margin[i] /= tour_margin_total;
                _depot_distance[i] /= depot_distance_total;

                /* Deltas in order are : Christofides, phil-margin, moat-packing, classical-margin, depot-distance.*/
                file<<size<<" "<<_mst[i] - _tsp[i]<<" "<<_phil[i] - _tsp[i]<<" "<<_moat[i] - _tsp[i]<<" "<<_tour_margin[i] - _tsp[i]<<" "<<_depot_distance[i] - _tsp[i]<<endl;
            }
        }
        
        size++;
    }
}



void write_apportionment(){
    ostringstream oss;
    oss<<"apportionment_data_"<<lower_range<<"_"<<upper_range<<".data";

    ofstream file;
    file.open (oss.str().c_str());
    write_apportionment(file);
    file.close();
}



void write_prescription_data(ofstream& file)
{
    size_t size = lower_range;
    size_t index =0 ; 
    for(vector<vector<ComparisonValues > >::iterator p = shapley_samples.begin()
            ; p != shapley_samples.end()
            ; p++){
        for(size_t j = 0 ; j < num_samples; j++){
            file<<size<<" ";
            for (int k =0  ; k < removal_prescription_values[index].size(); k++){
                file<<removal_prescription_values[index][k]<<" ";
            }
            file<<endl<<endl;
            index++;
        }
        
        size++;
    }
    

}

void write_prescription_data()
{
    ostringstream oss;
    oss<<"prescription_data_"<<lower_range<<"_"<<upper_range<<".data";

    ofstream file;
    file.open (oss.str().c_str());
    write_prescription_data(file);
    file.close();
}

void write_group_prescription_data(ofstream& file)
{
    size_t size = lower_range;
    size_t index =0 ; 
    for(vector<vector<ComparisonValues > >::iterator p = shapley_samples.begin()
            ; p != shapley_samples.end()
            ; p++){
        for(size_t j = 0 ; j < num_samples; j++){
            file<<size<<" ";
            for (int k =0  ; k < group_removal_prescription_values[index].size(); k++){
                file<<group_removal_prescription_values[index][k]<<" ";
            }
            file<<endl<<endl;
            index++;
        }
        size++;
    }
}

void write_group_prescription_data()
{
    ostringstream oss;
    oss<<"group_prescription_data_"<<lower_range<<"_"<<upper_range<<".data";

    ofstream file;
    file.open (oss.str().c_str());
    write_group_prescription_data(file);
    file.close();
}
