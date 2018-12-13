#include "ascent.hpp"
#include "conduit.hpp"
#include "conduit_blueprint.hpp"

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetReader.h>
#include <vtkDataSet.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>

#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <cstring>
#include <string.h>

using namespace conduit;
using namespace ascent;
using namespace std;

void GetVTKData(string input_path, string filename, 
                float *x, float *y, float *z, int *dims,
                vtkSmartPointer<vtkDataSetReader> rdr, vtkSmartPointer<vtkDataSet> ds,
                vtkSmartPointer<vtkAbstractArray> aa_x, vtkSmartPointer<vtkFloatArray> da_x,
                vtkSmartPointer<vtkAbstractArray> aa_y, vtkSmartPointer<vtkFloatArray> da_y,
                vtkSmartPointer<vtkAbstractArray> aa_z, vtkSmartPointer<vtkFloatArray> da_z,
                vtkSmartPointer<vtkAbstractArray> aa_d, vtkSmartPointer<vtkFloatArray> da_d,
                int cycle, double *time)
{
  string xprefix = "xmom";  
  string yprefix = "ymom";  
  string zprefix = "zmom";
  string dprefix = "density";  

  stringstream ss;
  ss << input_path << filename << setfill('0') << setw(5) << cycle << ".vtk";
  cout << ss.str() << endl;
  rdr->SetFileName(ss.str().c_str());
  rdr->Update();
  ds = rdr->GetOutput();

  vtkSmartPointer<vtkDoubleArray> t = vtkDoubleArray::SafeDownCast(ds->GetFieldData()->GetAbstractArray("TIME"));
  *time = t->GetValue(0);

  aa_x = ds->GetPointData()->GetAbstractArray(xprefix.c_str());
  da_x = vtkFloatArray::SafeDownCast(aa_x);
  aa_y = ds->GetPointData()->GetAbstractArray(yprefix.c_str());
  da_y = vtkFloatArray::SafeDownCast(aa_y);
  aa_z = ds->GetPointData()->GetAbstractArray(zprefix.c_str());
  da_z = vtkFloatArray::SafeDownCast(aa_z);
  aa_d = ds->GetPointData()->GetAbstractArray(dprefix.c_str());
  da_d = vtkFloatArray::SafeDownCast(aa_d);

  for(int i = 0 ; i < dims[2] ; i++)
  {
    for(int j = 0 ; j < dims[1] ; j++)
    {
      for(int k = 0; k < dims[0]; k++)
      {
        int index = i*dims[1]*dims[0] + j*dims[0] + k;
        x[index] = (da_x->GetValue(index))/(da_d->GetValue(index));
        y[index] = (da_y->GetValue(index))/(da_d->GetValue(index));
        z[index] = (da_z->GetValue(index))/(da_d->GetValue(index));
      }
    }
  }
}

void GetVTKTime(string input_path, string filename, 
                vtkSmartPointer<vtkDataSetReader> rdr, vtkSmartPointer<vtkDataSet> ds,
                int cycle, double *time)
{
  stringstream ss;
  ss << input_path << filename << setfill('0') << setw(5) << cycle << ".vtk";
  rdr->SetFileName(ss.str().c_str());
  rdr->Update();
  ds = rdr->GetOutput();
  vtkSmartPointer<vtkDoubleArray> t = vtkDoubleArray::SafeDownCast(ds->GetFieldData()->GetAbstractArray("TIME"));
  *time = t->GetValue(0);
}


int main(int argc, char *argv[])
{
  if(argc != 18)
    cout << "Too few parameters" << endl;  

  cout << "Starting Nyx Simulation" << endl;

  int dims[3];
  dims[0] = atoi(argv[1]);
  dims[1] = atoi(argv[2]);
  dims[2] = atoi(argv[3]);
  
  cout << "DIMENSIONS : " << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;

  int start_cycle = atoi(argv[4]);
  int end_cycle = atoi(argv[5]);

  int x_res = atoi(argv[6]);
  int y_res = atoi(argv[7]);
  int z_res = atoi(argv[8]);

  int write_frequency = atoi(argv[9]);

  string input_path = argv[10];
  string filename = argv[11];

  double xmin = atof(argv[12]);
  double ymin = atof(argv[13]);
  double zmin = atof(argv[14]);
  double xmax = atof(argv[15]);
  double ymax = atof(argv[16]);
  double zmax = atof(argv[17]);

  double spacing_x = (xmax-xmin)/(dims[0]-1);
  double spacing_y = (ymax-ymin)/(dims[1]-1);
  double spacing_z = (zmax-zmin)/(dims[2]-1);

  int npts = dims[0]*dims[1]*dims[2];

  float *x_vec = (float*)malloc(sizeof(float)*npts);
  float *y_vec = (float*)malloc(sizeof(float)*npts);
  float *z_vec = (float*)malloc(sizeof(float)*npts);

  Ascent ascent;

  Node ascent_opts;
  ascent_opts["runtime/type"] = "ascent";
  ascent.open(ascent_opts);

  double time_current, time_next, time_step_diff;
    
  vtkSmartPointer<vtkDataSetReader> rdr = vtkSmartPointer<vtkDataSetReader>::New();
  vtkSmartPointer<vtkDataSet> ds;
  vtkSmartPointer<vtkAbstractArray> aa_x;
  vtkSmartPointer<vtkAbstractArray> aa_y;
  vtkSmartPointer<vtkAbstractArray> aa_z;
  vtkSmartPointer<vtkAbstractArray> aa_d;
  vtkSmartPointer<vtkFloatArray> da_x;
  vtkSmartPointer<vtkFloatArray> da_y;
  vtkSmartPointer<vtkFloatArray> da_z;
  vtkSmartPointer<vtkFloatArray> da_d;

  for(int cycle = start_cycle; cycle < end_cycle; cycle++)
  {
    GetVTKData(input_path, filename, x_vec, y_vec, z_vec, dims, 
             rdr, ds, aa_x, da_x, aa_y, da_y, aa_z, da_z, aa_d, da_d,
             cycle, &time_current);
    
    GetVTKTime(input_path, filename, rdr, ds, cycle+1, &time_next); 
    time_step_diff = time_next - time_current;

    cout << "Setting step size to " << time_step_diff << endl;
    
    Conduit:Node mesh_data;
    // TODO Populate mesh_data;
    
    mesh_data["state/time"] = time_current;
    mesh_data["state/cycle"] = cycle;
    mesh_data["state/domain_id"] = 0;
    mesh_data["state/info"] = "Nyx";

    mesh_data["coordsets/coords/type"] = "uniform";
    mesh_data["coordsets/coords/dims/i"] = dims[0];
    mesh_data["coordsets/coords/dims/j"] = dims[1];
    mesh_data["coordsets/coords/dims/k"] = dims[2];

    mesh_data["coordsets/coords/origin/x"] = 0; 
    mesh_data["coordsets/coords/origin/y"] = 0; 
    mesh_data["coordsets/coords/origin/z"] = 0; 

    mesh_data["coordsets/coords/spacing/dx"] = spacing_x;
    mesh_data["coordsets/coords/spacing/dy"] = spacing_y;
    mesh_data["coordsets/coords/spacing/dz"] = spacing_z;
  
    mesh_data["topologies/mesh/type"] = "uniform";
    mesh_data["topologies/mesh/coordset"] = "coords";

    mesh_data["fields/vel/association"] = "vertex";
    mesh_data["fields/vel/type"] = "vector";
    mesh_data["fields/vel/topology"] = "mesh";
    mesh_data["fields/vel/values/u"].set_external(x_vec, npts);
    mesh_data["fields/vel/values/v"].set_external(y_vec, npts);
    mesh_data["fields/vel/values/w"].set_external(z_vec, npts);
    
    conduit::Node pipelines;
    // pipeline 1
    pipelines["pl1/f1/type"] = "lagrangian";
    // filter knobs
    conduit::Node &lagrangian_params = pipelines["pl1/f1/params"];
    lagrangian_params["field"] = "vel";
    lagrangian_params["step_size"] = time_step_diff;
    lagrangian_params["write_frequency"] = write_frequency;
    lagrangian_params["cust_res"] = 1;
    lagrangian_params["x_res"] = x_res;
    lagrangian_params["y_res"] = y_res;
    lagrangian_params["z_res"] = z_res;
    
    conduit::Node actions;
    // add the pipeline
    conduit::Node &add_pipelines = actions.append();
    add_pipelines["action"] = "add_pipelines";
    add_pipelines["pipelines"] = pipelines;

    // execute
    conduit::Node &execute  = actions.append();
    execute["action"] = "execute";
    // reset
    conduit::Node &reset  = actions.append();
    reset["action"] = "reset";

    ascent.publish(mesh_data);
    ascent.execute(actions);
  }

  ascent.close();

}
