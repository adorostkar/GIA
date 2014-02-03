#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/numerics/data_out_faces.h>

#include "parameters.h"

#include <fstream>
#include <sstream>
#include <string>
#include <typeinfo>

#ifndef SURFACEDATAOUT_H
#define SURFACEDATAOUT_H

template<int dim>
class SurfaceDataOut : public DataOutFaces<dim> {
public:
    SurfaceDataOut(){
        par = parameters::getInstance();
    }

    virtual typename DataOutFaces<dim>::FaceDescriptor
    first_face (){
        typename DataOutFaces<dim>::active_cell_iterator
                cell = this->dofs->begin_active();
        for (; cell != this->dofs->end(); ++cell){
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
                if (cell->face(f)->at_boundary() && (
                        cell->face(f)->boundary_indicator() == par->b_ice ||
                            cell->face(f)->boundary_indicator() == par->b_up))
                {
                    return typename DataOutFaces<dim>::FaceDescriptor(cell, f);
                }
            }
        }
        return typename DataOutFaces<dim>::FaceDescriptor();
    }

    virtual typename DataOutFaces<dim>::FaceDescriptor
    next_face (const typename DataOutFaces<dim>::FaceDescriptor &old_face){
        typename DataOutFaces<dim>::FaceDescriptor face = old_face;

        for (unsigned int f=face.second+1; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (face.first->face(f)->at_boundary() &&
                    (face.first->face(f)->boundary_indicator() == par->b_ice ||
                        face.first->face(f)->boundary_indicator() == par->b_up))
            {
                face.second = f;
                return face;
            };
        typename DataOutFaces<dim>::active_cell_iterator active_cell = face.first;

        // increase face pointer by one
        ++active_cell;

        // while there are active cells
        while (active_cell != this->dofs->end()) {
            // check all the faces of this
            // active cell
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (active_cell->face(f)->at_boundary() &&
                        (active_cell->face(f)->boundary_indicator() == par->b_ice ||
                            active_cell->face(f)->boundary_indicator() == par->b_up)){
                    face.first  = active_cell;
                    face.second = f;
                    return face;
                };
            // the present cell had no
            // faces on the boundary, so
            // check next cell
            ++active_cell;
        };

        // we fell off the edge, so return
        // with invalid pointer
        face.first  = this->dofs->end();
        face.second = 0;
        return face;
    }

private:
    parameters* par;
};

#endif // SURFACEDATAOUT_H
