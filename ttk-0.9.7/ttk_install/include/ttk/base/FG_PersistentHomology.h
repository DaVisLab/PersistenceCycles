/// \ingroup base
/// \class ttk::FG_PersistentHomology
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK %fG_PersistentHomology processing package.
///
/// %FG_PersistentHomology is a TTK processing package that takes a scalar field on the input
/// and produces a scalar field on the output.
///
/// \sa ttk::Triangulation
/// \sa ttkFG_PersistentHomology.cpp %for a usage example.

#pragma once

// base code includes
#include                  <Triangulation.h>
#include                  <Wrapper.h>
#include                  <FormanGradient.h>
#include                  <boundarymatrix.h>

#include <queue>

namespace ttk{

  class FG_PersistentHomology : public FormanGradient{

    public:

      FG_PersistentHomology();
      ~FG_PersistentHomology();

      template <class dataType>
      void computeIndexing(double&,double&);

      template <class dataType>
      void computeBoundayMatrix();

      template <class dataType>
      void findPersistenceInterval(double& minPers, double& maxPers);

      void readHomology(vector<float>&, vector<char>&);

      template <class dataType>
      void readPersistencePairs(vector<float>&, vector<double>&, vector<char>&, double realMinPers, double realMaxPers);

      template <class dataType>
      void computerPersistenceCycles(Simplex simpl1, Simplex simpl2, map<Simplex,int>& simplexToIndex, vector<int>& generators, bool save_hole, vector<int>& holes);

      template <class dataType>
      void readCycle(vector<float>& coordinates_gen,
                     vector<vector<int> >& simplices_gen,
                     list<int>& indices_gen,
                     vector<int>& vertices_hole,
                     vector<vector<int> >& simplices_hole,
                     list<int>& indices_hole,
                     vector<double>& filtration, double minPers, double maxPers, bool formanCycles, bool compute_hole, int dim);

      void readGradientVectors(vector<float>& points, vector<float>& vects);

      inline void setInputDataPointer(const void *data){
        inputData_ = (const void*)data;
      }

      template <class dataType>
      void printPersistenceDiagram(double realMinPers, double realMaxPers);

    protected:
      void getBoundarySimplices(Simplex,vector<Simplex>&);
      void collectSeeds(int dim, SimplexId face, vector<SimplexId>& cofaces);

      template <class dataType>
      dataType getFiltration(Simplex& simplex);


    protected:
      vector<vector<int> >        matrix;
      map<int,int>                allpairs;
      vector<int>                 homology;

      vector<Simplex>             index_to_simplex; //mapping each critical simplex to the corresponding index in the boundary matrix

      const void*                       inputData_;
  };
}

#include                  <FG_PersistentHomology_template.h>
