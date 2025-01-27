#include                  <ttkFG_PersistentHomology.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkFG_PersistentHomology)

int ttkFG_PersistentHomology::doIt(vector<vtkDataSet *> &inputs, vector<vtkDataSet *> &outputs){

    MemoryUsage m;
    vtkDataArray *inputScalarField = NULL;
    vtkDataSet *input = inputs[0];

    if(ScalarField.length()){
        inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
    }
    else{
        inputScalarField = input->GetPointData()->GetArray(0);
    }

    if(!inputScalarField)
        return -2;

    Triangulation *triangulation = ttkTriangulation::getTriangulation(input);

    if(fG_PersistentHomology_ == NULL){
        fG_PersistentHomology_ = new FG_PersistentHomology();

        if(!triangulation)
          return -1;

        triangulation->setWrapper(this);
        fG_PersistentHomology_->setupTriangulation(triangulation);
        fG_PersistentHomology_->setWrapper(this);

        fG_PersistentHomology_->setInputDataPointer(inputScalarField->GetVoidPointer(0));

        Timer time;

        switch(inputScalarField->GetDataType()){
          ttkTemplateMacro(fG_PersistentHomology_->computeIndexing<VTK_TT>(maxFunctValue, minFunctionValue ))
        }

        cout << "Gradient computed in " << time.getElapsedTime() << " seconds" << endl;
        cout << "Memory used for discrete gradient: " << m.getValue_in_MB(false) << " MB" <<endl;
        //generally the template for this function should agree with the input scalar field. But here we are working on the indexing
        //so we should use integers
        time.reStart();
        fG_PersistentHomology_->computeBoundayMatrix<int>();
        cout << "Persistent homology computed in " << time.getElapsedTime() << " seconds" << endl;
        cout << "Memory used for boundary matrix: " << m.getValue_in_MB(false) << " MB"<< endl;

      switch(inputScalarField->GetDataType()){
          ttkTemplateMacro(fG_PersistentHomology_->findPersistenceInterval<VTK_TT>(realMinPers, realMaxPers));
      }
  }
  else{
    triangulation->setWrapper(this);
    fG_PersistentHomology_->setupTriangulation(triangulation);
    fG_PersistentHomology_->setWrapper(this);

    fG_PersistentHomology_->setInputDataPointer(inputScalarField->GetVoidPointer(0));
  }

    outPersistencePairs(inputScalarField,outputs[0],outputs[1]);
    outHomology(outputs[2]);

    Timer time;
    if(cycles1){
        time.reStart();
        if(!holes1)
            out1Cycles(inputScalarField, triangulation, outputs[3], NULL);
        else
            out1Cycles(inputScalarField, triangulation, outputs[3], outputs[4]);
        cout << "1-cycles computed in " << time.getElapsedTime() << " seconds" << endl;
    }
    if(cycles2){
        time.reStart();
        if(!holes2)
            out2Cycles(inputScalarField, triangulation, outputs[5], NULL);
        else
            out2Cycles(inputScalarField, triangulation, outputs[5], outputs[6]);
        cout << "2-cycles computed in " << time.getElapsedTime() << " seconds" << endl;
    }

    // if(Pairs){
    //     switch(inputScalarField->GetDataType()){
    //         ttkTemplateMacro(fG_PersistentHomology_->printPersistenceDiagram<VTK_TT>(realMaxPers*MinPers, realMaxPers*MaxPers));
    //     }
    // }

  return 0;
}

void ttkFG_PersistentHomology::outPersistencePairs(vtkDataArray *inputScalarField, vtkDataSet* output, vtkDataSet* outputDiag){

    vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(output);
    vtkUnstructuredGrid* outputPersDiagrams=vtkUnstructuredGrid::SafeDownCast(outputDiag);

    vector<float> criticalPoints = vector<float>();
    vector<char> criticalPointsCellDimension = vector<char>();
    vector<double> criticalPointsFunctionValue = vector<double>();

    switch(inputScalarField->GetDataType()){
        ttkTemplateMacro(fG_PersistentHomology_->readPersistencePairs<VTK_TT>(criticalPoints, criticalPointsFunctionValue, criticalPointsCellDimension, realMaxPers*MinPers, realMaxPers*MaxPers));
    }

    vtkSmartPointer<vtkPoints> pointsDomain=vtkSmartPointer<vtkPoints>::New();
    //prepare array of critical points dimension
    vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
    cellDimensions->SetNumberOfComponents(1);
    cellDimensions->SetName("CellDimension");

    vtkSmartPointer<vtkDoubleArray> step=vtkSmartPointer<vtkDoubleArray>::New();
    step->SetNumberOfComponents(1);
    step->SetName("Filtration");

    //Embed the persistence diagram into the domain
    for(int i=0; i<criticalPointsCellDimension.size(); i++){
        //add new point to list
        pointsDomain->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);
        //add new point dimension (corresponding to the original cell originating it)
        cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);
        step->InsertNextTuple1(criticalPointsFunctionValue[i]);

    }
    outputCriticalPoints->SetPoints(pointsDomain);
    vtkPointData* pointData=outputCriticalPoints->GetPointData();
    pointData->AddArray(cellDimensions);
    pointData->AddArray(step);


    vtkSmartPointer<vtkDoubleArray> filtration=vtkSmartPointer<vtkDoubleArray>::New();
    filtration->SetNumberOfComponents(1);
    filtration->SetName("Persistence");

    vtkSmartPointer<vtkDoubleArray> pairtype=vtkSmartPointer<vtkDoubleArray>::New();
    pairtype->SetNumberOfComponents(1);
    pairtype->SetName("Type");

    outputCriticalPoints->Allocate(criticalPointsCellDimension.size()/2.0);
    for(int i=0; i<criticalPointsCellDimension.size(); i=i+2){
        vtkIdType line[2];
        line[0]=i;
        line[1]=i+1;
        outputCriticalPoints->InsertNextCell(VTK_LINE, 2, line);
        filtration->InsertNextTuple1(criticalPointsFunctionValue[i+1] - criticalPointsFunctionValue[i]);
        pairtype->InsertNextTuple1(criticalPointsCellDimension[i]);
    }

    vtkCellData* cellData = outputCriticalPoints->GetCellData();
    cellData->AddArray(filtration);
    cellData->AddArray(pairtype);

    //Draw classic persistence diagram
    vtkSmartPointer<vtkPoints> pointsPD=vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCharArray> pairType=vtkSmartPointer<vtkCharArray>::New();
    pairType->SetNumberOfComponents(1);
    pairType->SetName("Pair Type");

    vtkSmartPointer<vtkDoubleArray> filtrationPD=vtkSmartPointer<vtkDoubleArray>::New();
    filtrationPD->SetNumberOfComponents(1);
    filtrationPD->SetName("Filtration");

    for(int i=0; i<criticalPointsFunctionValue.size(); i+=2){
        //add new point to list
        pointsPD->InsertNextPoint(criticalPointsFunctionValue[i],criticalPointsFunctionValue[i+1],0);
        pairType->InsertNextTuple1(criticalPointsCellDimension[i]);
        filtrationPD->InsertNextTuple1(criticalPointsFunctionValue[i+1]-criticalPointsFunctionValue[i]);
    }

    pointsPD->InsertNextPoint((float)minFunctionValue,(float)minFunctionValue,0);
    pairType->InsertNextTuple1(-1);
    filtrationPD->InsertNextTuple1(maxFunctValue+1);

    pointsPD->InsertNextPoint((float)maxFunctValue,(float)maxFunctValue,0);
    pairType->InsertNextTuple1(-1);
    filtrationPD->InsertNextTuple1(maxFunctValue+1);

    vtkIdType line[2];
    line[0]=criticalPointsFunctionValue.size()/2;
    line[1]=criticalPointsFunctionValue.size()/2+1;
    outputPersDiagrams->Allocate(1);
    outputPersDiagrams->InsertNextCell(VTK_LINE, 2, line);


    outputPersDiagrams->SetPoints(pointsPD);
    outputPersDiagrams->GetPointData()->AddArray(pairType);
    outputPersDiagrams->GetPointData()->AddArray(filtrationPD);

}

void ttkFG_PersistentHomology::outHomology(vtkDataSet* output){

  vtkUnstructuredGrid* outputCriticalPoints=vtkUnstructuredGrid::SafeDownCast(output);

  vector<float> criticalPoints = vector<float>();
  vector<char> criticalPointsCellDimension = vector<char>();

  fG_PersistentHomology_->readHomology(criticalPoints,criticalPointsCellDimension);

  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
  //prepare array of critical points dimension
  vtkSmartPointer<vtkCharArray> cellDimensions=vtkSmartPointer<vtkCharArray>::New();
  cellDimensions->SetNumberOfComponents(1);
  cellDimensions->SetName("CellDimension");

  for(int i=0; i<criticalPointsCellDimension.size(); i++){
    //add new point to list
    points->InsertNextPoint(criticalPoints[3*i],criticalPoints[3*i+1],criticalPoints[3*i+2]);

    //add new point dimension (corresponding to the original cell originating it)
    cellDimensions->InsertNextTuple1(criticalPointsCellDimension[i]);

  }
  outputCriticalPoints->SetPoints(points);

  vtkPointData* pointData=outputCriticalPoints->GetPointData();
  pointData->AddArray(cellDimensions);

}

void ttkFG_PersistentHomology::out1Cycles(vtkDataArray *inputScalarField, Triangulation* triangulation, vtkDataSet* output_gen, vtkDataSet* output_hole){

  vtkUnstructuredGrid* outputCycles = vtkUnstructuredGrid::SafeDownCast(output_gen);

  vector<float> coordinates_gen;
  vector<vector<int> > simplices_gen;
  list<int> indices_gen;

  vector<int> vertices_hole;
  vector<vector<int> > simplices_hole;
  list<int> indices_hole;

    vector<double> filtration;


  switch(inputScalarField->GetDataType()){
      ttkTemplateMacro(fG_PersistentHomology_->readCycle<VTK_TT>(coordinates_gen,simplices_gen,indices_gen,vertices_hole,simplices_hole,indices_hole, filtration, realMaxPers*MinPers, realMaxPers*MaxPers, formanCycles, holes1, 1));
    }
  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

  for(int i=0; i<coordinates_gen.size(); i+=3){
    //add new point to list
    points->InsertNextPoint(coordinates_gen[i],coordinates_gen[i+1],coordinates_gen[i+2]);
  }
  outputCycles->SetPoints(points);

  outputCycles->Allocate(simplices_gen.size());
  for(int i=0; i<simplices_gen.size(); i++){
    vtkIdType line[2];
    line[0]=simplices_gen[i][0];
    line[1]=simplices_gen[i][1];
    outputCycles->InsertNextCell(VTK_LINE, 2, line);
  }

  vtkSmartPointer<vtkIntArray> cycleId=vtkSmartPointer<vtkIntArray>::New();
  cycleId->SetNumberOfComponents(1);
  cycleId->SetName("CycleId");

    vtkSmartPointer<vtkDoubleArray> filtr=vtkSmartPointer<vtkDoubleArray>::New();
    filtr->SetNumberOfComponents(1);
    filtr->SetName("Filtration");

  int count=0;
  int f=0;
  for(auto val : indices_gen){
    for(int i=0; i<val; i++){
      cycleId->InsertNextTuple1(count);
      filtr->InsertNextTuple1(filtration[f]);
    }
    f++;
    count++;
  }

  vtkCellData* cellData=outputCycles->GetCellData();
  cellData->AddArray(cycleId);
  cellData->AddArray(filtr);


    if(holes1){

        vtkUnstructuredGrid* outputHoles = vtkUnstructuredGrid::SafeDownCast(output_hole);
        vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkDoubleArray> funct=vtkSmartPointer<vtkDoubleArray>::New();
        funct->SetNumberOfComponents(1);
        funct->SetName("FunctValues");

        for(int i=0; i<vertices_hole.size(); i++){
            //add new point to list
            float x,y,z;
            triangulation->getVertexPoint(vertices_hole[i],x,y,z);
            points->InsertNextPoint(x,y,z);
            funct->InsertNextTuple1(inputScalarField->GetTuple1(vertices_hole[i]));


        }
        outputHoles->SetPoints(points);

        vtkPointData* pointData=outputHoles->GetPointData();
        pointData->AddArray(funct);

        outputHoles->Allocate(simplices_hole.size());
        for(int i=0; i<simplices_hole.size(); i++){
            vtkIdType triangle[3];
            triangle[0]=simplices_hole[i][0];
            triangle[1]=simplices_hole[i][1];
            triangle[2]=simplices_hole[i][2];
            outputHoles->InsertNextCell(VTK_TRIANGLE, 3, triangle);
        }

        vtkSmartPointer<vtkIntArray> cycleId=vtkSmartPointer<vtkIntArray>::New();
        cycleId->SetNumberOfComponents(1);
        cycleId->SetName("CycleId");

        vtkSmartPointer<vtkDoubleArray> filtr=vtkSmartPointer<vtkDoubleArray>::New();
        filtr->SetNumberOfComponents(1);
        filtr->SetName("Filtration");

        int count=0;
        int f=0;
        for(auto val : indices_hole){
            for(int i=0; i<val; i++){
                cycleId->InsertNextTuple1(count);
                filtr->InsertNextTuple1(filtration[f]);
            }
            f++;
            count++;
        }

        vtkCellData* cellData=outputHoles->GetCellData();
        cellData->AddArray(cycleId);
        cellData->AddArray(filtr);
    }

}

void ttkFG_PersistentHomology::out2Cycles(vtkDataArray *inputScalarField, Triangulation* triangulation, vtkDataSet* output_gen, vtkDataSet* output_hole){

   vtkUnstructuredGrid* outputCycles = vtkUnstructuredGrid::SafeDownCast(output_gen);

    vector<float> coordinates_gen;
    vector<vector<int> > simplices_gen;
    list<int> indices_gen;

    vector<int> vertices_hole;
    vector<vector<int> > simplices_hole;
    list<int> indices_hole;

    vector<double> f_values;

    vector<double> filtration;


    switch(inputScalarField->GetDataType()){
        ttkTemplateMacro(fG_PersistentHomology_->readCycle<VTK_TT>(coordinates_gen,simplices_gen,indices_gen,vertices_hole,simplices_hole,indices_hole, filtration, realMaxPers*MinPers, realMaxPers*MaxPers, formanCycles, holes2, 2));
    }

   vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

   for(int i=0; i<coordinates_gen.size(); i+=3){
     //add new point to list
     points->InsertNextPoint(coordinates_gen[i],coordinates_gen[i+1],coordinates_gen[i+2]);
   }
   outputCycles->SetPoints(points);

   outputCycles->Allocate(simplices_gen.size());
   for(int i=0; i<simplices_gen.size(); i++){
     vtkIdType triangle[3];
     triangle[0]=simplices_gen[i][0];
     triangle[1]=simplices_gen[i][1];
     triangle[2]=simplices_gen[i][2];
     outputCycles->InsertNextCell(VTK_TRIANGLE, 3, triangle);
   }

    vtkSmartPointer<vtkIntArray> cycleId=vtkSmartPointer<vtkIntArray>::New();
    cycleId->SetNumberOfComponents(1);
    cycleId->SetName("CycleId");

    vtkSmartPointer<vtkDoubleArray> filtr=vtkSmartPointer<vtkDoubleArray>::New();
    filtr->SetNumberOfComponents(1);
    filtr->SetName("Filtration");

    int count=0;
    int f=0;
    for(auto val : indices_gen){
        for(int i=0; i<val; i++){
            cycleId->InsertNextTuple1(count);
            filtr->InsertNextTuple1(filtration[f]);
        }
        f++;
        count++;
    }

    vtkCellData* cellData=outputCycles->GetCellData();
    cellData->AddArray(cycleId);
    cellData->AddArray(filtr);


    if(holes2){

        vtkUnstructuredGrid* outputHoles = vtkUnstructuredGrid::SafeDownCast(output_hole);
        vtkSmartPointer<vtkDoubleArray> funct=vtkSmartPointer<vtkDoubleArray>::New();
        funct->SetNumberOfComponents(1);
        funct->SetName("FunctValues");
        vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

        for(int i=0; i<vertices_hole.size(); i++){
            //add new point to list
            float x,y,z;
            triangulation->getVertexPoint(vertices_hole[i],x,y,z);
            points->InsertNextPoint(x,y,z);
            funct->InsertNextTuple1(inputScalarField->GetTuple1(vertices_hole[i]));
        }
        outputHoles->SetPoints(points);

        vtkPointData* pointData=outputHoles->GetPointData();
        pointData->AddArray(funct);

        outputHoles->Allocate(simplices_hole.size());
        for(int i=0; i<simplices_hole.size(); i++){
            vtkIdType tetrahedron[4];
            tetrahedron[0]=simplices_hole[i][0];
            tetrahedron[1]=simplices_hole[i][1];
            tetrahedron[2]=simplices_hole[i][2];
            tetrahedron[3]=simplices_hole[i][3];
            outputHoles->InsertNextCell(VTK_TETRA, 4, tetrahedron);
        }

        vtkSmartPointer<vtkIntArray> cycleId=vtkSmartPointer<vtkIntArray>::New();
        cycleId->SetNumberOfComponents(1);
        cycleId->SetName("CycleId");

        vtkSmartPointer<vtkDoubleArray> filtr=vtkSmartPointer<vtkDoubleArray>::New();
        filtr->SetNumberOfComponents(1);
        filtr->SetName("Filtration");

        int count=0;
        int f=0;
        for(auto val : indices_hole){
            for(int i=0; i<val; i++){
                cycleId->InsertNextTuple1(count);
                filtr->InsertNextTuple1(filtration[f]);
            }
            f++;
            count++;
        }

        vtkCellData* cellData=outputHoles->GetCellData();
        cellData->AddArray(cycleId);
        cellData->AddArray(filtr);
    }

}

void ttkFG_PersistentHomology::outGradient(vtkDataSet* output){

    vtkUnstructuredGrid* outputGradient = vtkUnstructuredGrid::SafeDownCast(output);

    vector<float> pts = vector<float>();
    vector<float> vects = vector<float>();

    fG_PersistentHomology_->readGradientVectors(pts, vects);


    vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkDoubleArray> vectglyph=vtkSmartPointer<vtkDoubleArray>::New();
    vectglyph->SetNumberOfComponents(3);
    vectglyph->SetName("Gradient");

    for(int i=0; i<pts.size(); i+=3){
        //add new point to list
        points->InsertNextPoint(pts[i],pts[i+1],pts[i+2]);
        vectglyph->InsertNextTuple3(vects[i],vects[i+1],vects[i+2]);
    }
    outputGradient->SetPoints(points);
    vtkPointData* pointData=outputGradient->GetPointData();
    pointData->AddArray(vectglyph);

}
