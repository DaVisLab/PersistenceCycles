#define FG_PERSISTENTHOMOLOGY_TEMPLATE_H

#include "FG_PersistentHomology.h"

using namespace std;
using namespace ttk;

template <class dataType> void ttk::FG_PersistentHomology::computeIndexing(double& maxF, double& minF){


    dataType* field = (dataType*)inputData_;

    maxF = 0;
    minF = double(field[0]);
    // building the indexing

    int vertices = triangulation_->getNumberOfVertices();

    vector<pair<dataType,SimplexId> > thepairs(vertices);
    for(int i=0; i<vertices; i++){
        maxF = std::max(maxF, (double)field[i]);
        minF = std::min(minF, (double)field[i]);
      thepairs[i] = pair<dataType,SimplexId>(field[i],i);
    }

    sort(thepairs.begin(), thepairs.end());


    indexing_ = new vector<SimplexId>(vertices);
    for(SimplexId i=0; i<vertices; i++){
      (*indexing_)[thepairs[i].second]=i;
    }

    computeGradient(indexing_);
}

template <class dataType> void ttk::FG_PersistentHomology::computeBoundayMatrix(){

    BoundaryMatrix<dataType>* bdmatrix = new BoundaryMatrix<dataType>();
    map<Simplex, int> simplex_to_index;


    int global_index=0;
    for(int d=0; d <= dimensionality_; d++){
        for(auto i : criticalSimplices[d]){
            Simplex simplex = Simplex(d,i);
            simplex_to_index[simplex] = global_index++; //add the simplex to the set
        }
    }



    for(int d=0; d <= dimensionality_; d++){
        vector<SimplexId> criticals(criticalSimplices[d].begin(),criticalSimplices[d].end());
        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for shared(bdmatrix,criticals)
        #endif
        for(int i=0; i<criticals.size(); i++){
            Simplex simplex = Simplex(d,criticals[i]);
            vector<Simplex> boundary;
            if(simplex.first == 0){
                #ifdef TTK_ENABLE_OPENMP
                #pragma omp critical
                {
                #endif
                bdmatrix->addValue(simplex_to_index[simplex],-1,getIndex(simplex));
                #ifdef TTK_ENABLE_OPENMP
                }
                #endif
            }
            else{
                getBoundarySimplices(simplex,boundary);

                if(boundary.size() == 0){
                    #ifdef TTK_ENABLE_OPENMP
                    #pragma omp critical
                    {
                    #endif
                    bdmatrix->addValue(simplex_to_index[simplex],-1,getIndex(simplex));
                    #ifdef TTK_ENABLE_OPENMP
                    }
                    #endif

                }
                else{
                    #ifdef TTK_ENABLE_OPENMP
                    #pragma omp critical
                    {
                    #endif
                    for(auto s : boundary){
                        bdmatrix->addValue(simplex_to_index[simplex],simplex_to_index[s],getIndex(simplex));
                    }
                    #ifdef TTK_ENABLE_OPENMP
                    }
                    #endif
                }
            }
        }
    }


    index_to_simplex = vector<Simplex>(global_index);
    for(auto p : simplex_to_index)
        index_to_simplex[p.second] = p.first;

    Timer t;


    bdmatrix->sort();
    bdmatrix->reduce();
    bdmatrix->getMatrix(matrix);
    bdmatrix->getPairs(allpairs);
    bdmatrix->getHomology(homology);

}

template <class dataType>
void FG_PersistentHomology::findPersistenceInterval(double& minPers, double& maxPers){

    dataType minP = 0;
    dataType maxP = 0;

    int count_real_pairs = 0;
    for(auto ppair :allpairs){

        Simplex simpl1 = index_to_simplex[ppair.first];
        Simplex simpl2 = index_to_simplex[ppair.second];

        dataType f1 = getFiltration<dataType>(simpl1);
        dataType f2 = getFiltration<dataType>(simpl2);

        dataType filtr = f1-f2;
        if(filtr < 0) filtr = -filtr;

        if(filtr != 0)
            count_real_pairs++;

        if (minP > filtr)
            minP = filtr;

        if (maxP < filtr)
            maxP = filtr;

    }

    minPers = double(minP);
    maxPers = double(maxP);

    cout << allpairs.size() << " total number of pairs " << count_real_pairs << " with non-zero persistence" << endl;
    cout << " Minimum Persistence: " << minPers << endl;
    cout << " Maximum Persistence: " << maxPers << endl;
}


template <class dataType>
void FG_PersistentHomology::readPersistencePairs(vector<float>& points, vector<double>& fvalues, vector<char>& dims, double realMinPers, double realMaxPers){

  vector<float> coords;

  for(auto ppair :allpairs){

    Simplex simpl1 = index_to_simplex[ppair.first];
    Simplex simpl2 = index_to_simplex[ppair.second];

    dataType f1 = getFiltration<dataType>(simpl1);
    dataType f2 = getFiltration<dataType>(simpl2);

    double filtr = double(f1-f2);
    if(filtr < 0) filtr = -filtr;

    if(filtr >= realMinPers && filtr <= realMaxPers){

      computeBarycenter(simpl1, coords);
      points.push_back(coords[0]);
      points.push_back(coords[1]);
      points.push_back(coords[2]);
      dims.push_back(simpl1.first);

      computeBarycenter(simpl2, coords);
      points.push_back(coords[0]);
      points.push_back(coords[1]);
      points.push_back(coords[2]);
      dims.push_back(simpl2.first);

      fvalues.push_back(f1);
      fvalues.push_back(f2);
    }

  }

}


template <class dataType>
dataType FG_PersistentHomology::getFiltration(Simplex& simplex){

  dataType *inputData = (dataType *) inputData_;

  vector<SimplexId> vertices = vector<SimplexId>();
  simplexToVertices(simplex, vertices);

  //from the vertices select the field value of the simplex
  dataType fieldValue = inputData[vertices[0]];
  for(SimplexId v=1; v<vertices.size(); v++){
    fieldValue = fieldValue > inputData[vertices[v]] ? fieldValue : inputData[vertices[v]];
  }

  return fieldValue;
}


template <class dataType>
void FG_PersistentHomology::computerPersistenceCycles(Simplex simpl1, Simplex simpl2, map<Simplex,int>& simplexToIndex, vector<int>& generators, bool save_hole, vector<int>& holes){

  int dim = simpl2.first;
  int pers_val = getIndex(simpl1);

  set<SimplexId> growingHole;

  set<SimplexId>  visited_simpl;
  visited_simpl.insert(simpl1.second);

  set<SimplexId> growing_generator;

  queue<int> seeds;
  seeds.push(simpl2.second);

  while(!seeds.empty()){
    SimplexId triangleIdx = seeds.front();
    seeds.pop();

      if(save_hole)
          growingHole.insert(triangleIdx);

    for(int i=0; i<dim+1; i++){
      SimplexId face = extractBoundary(Simplex(dim,triangleIdx),dim-1,i);
      if(growing_generator.find(face) == growing_generator.end()){
        growing_generator.insert(face);
      }
      else{
        growing_generator.erase(face);
      }
    }

    for(int i=0; i<dim+1; i++){
      SimplexId face = extractBoundary(Simplex(dim,triangleIdx),dim-1,i);

      Simplex old = Simplex(dim-1,face);

      if(pers_val > getIndex(old))
        continue;

      Simplex pair;
      bool isPaired = getPair(old, pair);

      if(isPaired && pair.first == dim && pair.second != triangleIdx){

        if(pers_val < getIndex(pair)) {
            seeds.push(pair.second);
        }
      }
      else if(!isPaired){

          if(visited_simpl.find(old.second) == visited_simpl.end() && pers_val < getIndex(old)){
              int index = simplexToIndex[old];
              if(allpairs.find(index) != allpairs.end()){
                  seeds.push(index_to_simplex[allpairs[index]].second);
              }

              visited_simpl.insert(old.second);
          }
      }
    }
  }

  generators.clear();
  generators.insert(generators.end(), growing_generator.begin(), growing_generator.end());

  if(save_hole){
      holes.clear();
      holes.insert(holes.end(), growingHole.begin(), growingHole.end());
  }

}


// template <class dataType>
// void FG_PersistentHomology::computeGeneratorForman(Simplex simpl1, Simplex simpl2, vector<int>& column, vector<int>& generators){


//     for(int cSaddle : column){
//        Simplex simpl = index_to_simplex[cSaddle];

//         queue<int> seeds;
//         seeds.push(simpl.second);
//         generators.push_back(simpl.second);

//         while(!seeds.empty()){
//             SimplexId simplId = seeds.front();
//             seeds.pop();

//             for(int i=0; i<simpl.first+1; i++) {
//                 SimplexId face = extractBoundary(Simplex(simpl.first, simplId), simpl.first - 1, i);
//                 Simplex pair;
//                 bool isPaired = getPair(Simplex(simpl.first-1,face), pair);

//                 if(isPaired && pair.first == simpl.first && pair.second != simplId){
//                     seeds.push(pair.second);
//                     generators.push_back(pair.second);
//                 }
//             }
//         }
//     }

// }

template <class dataType>
void FG_PersistentHomology::readCycle(vector<float>& coordinates_gen,
                                      vector<vector<int> >& simplices_gen,
                                      list<int>& indices_gen,
                                      vector<int>& vertices_hole,
                                      vector<vector<int> >& simplices_hole,
                                      list<int>& indices_hole,
                                      vector<double>& filtration, double minPers, double maxPers, bool formanCycles, bool compute_hole, int dim){

    list<vector<int> > list_generators;
    list<vector<int> > list_holes;

    map<Simplex,int> simplexToIndex;
    if(!formanCycles){
        for(int i=0; i<index_to_simplex.size(); i++){
            simplexToIndex[index_to_simplex[i]]=i;
        }
    }

    vector<SimplexId> origin_homology = vector<SimplexId>();
    for(auto ppair :allpairs) {
        Simplex simpl1 = index_to_simplex[ppair.first];
        if(simpl1.first != dim)
            continue;
        origin_homology.push_back(ppair.first);
    }

    int tot_simplices=0;
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for shared(filtration, list_generators, tot_simplices)
    #endif
    for(int i=0; i<origin_homology.size(); i++){

        Simplex simpl1 = index_to_simplex[origin_homology[i]];
        Simplex simpl2 = index_to_simplex[allpairs[origin_homology[i]]];

        if(simpl1.first != dim)
            continue;

        dataType f1 = getFiltration<dataType>(simpl1);
        dataType f2 = getFiltration<dataType>(simpl2);

        double filtr = double(f2-f1);
        if(filtr < 0) filtr = -filtr;

        if( filtr >= minPers && filtr <= maxPers){
            vector<int> generator;
            vector<int> hole;
            computerPersistenceCycles<dataType>(simpl1, simpl2, simplexToIndex, generator, compute_hole, hole);


            #ifdef TTK_ENABLE_OPENMP
            #pragma omp critical
            {
            #endif
                list_generators.push_back(generator);
                filtration.push_back(filtr);
                if(compute_hole)
                    list_holes.push_back(hole);
                tot_simplices += generator.size();
            #ifdef TTK_ENABLE_OPENMP
            }
            #endif
        }
    }


    //prepare vertices from the list of generators
    map<int,int> unique_vertices;
    int count = 0;
    vector<float> coords;

    for(auto gen : list_generators){

        for(auto simpl : gen){
            for(int i=0; i<dim+1; i++){
                SimplexId v = extractBoundary(Simplex(dim,simpl),0,i);
                if(unique_vertices.find(v) == unique_vertices.end()){
                    unique_vertices[v] = count++;
                    Simplex vertex = Simplex(0,v);
                    computeBarycenter(vertex, coords);

                    for(int k=0; k<coords.size(); k++)
                        coordinates_gen.push_back(coords[k]);

                }
            }
        }
    }

    for(auto gen : list_generators){

        for(auto simpl : gen){
            vector<int> vertices(dim+1);
            for(int i=0; i<dim+1; i++){
                vertices[i] = unique_vertices[extractBoundary(Simplex(dim,simpl),0,i)];
            }
            simplices_gen.push_back(vertices);
        }
        indices_gen.push_back(gen.size());

    }

    //prepare vertices from the list of holes
    if(compute_hole){
        unique_vertices.clear();
        count = 0;

        for (auto gen : list_holes) {
            for (auto simpl : gen) {
                for (int i = 0; i < dim + 2; i++) {
                    SimplexId v = extractBoundary(Simplex(dim + 1, simpl), 0, i);
                    if (unique_vertices.find(v) == unique_vertices.end()) {
                        unique_vertices[v] = count++;
                        vertices_hole.push_back(v);

                    }
                }
            }
        }

        for (auto gen : list_holes) {
            for (auto simpl : gen) {
                vector<int> vertices(dim + 2);
                for (int i = 0; i < dim + 2; i++) {
                    vertices[i] = unique_vertices[extractBoundary(Simplex(dim + 1, simpl), 0, i)];
                }
                simplices_hole.push_back(vertices);
            }
            indices_hole.push_back(gen.size());

        }
    }

}



