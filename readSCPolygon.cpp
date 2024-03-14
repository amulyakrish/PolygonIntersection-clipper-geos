#include "readSCPolygon.h"

void gpc_read_polygon(FILE *fp, gpc_polygon *p)
{
    int v;

    //fscanf(fp, "%d", &(p->num_contours));
  
    MALLOC(p->contour, 1
         * sizeof(gpc_vertex_list), "contour creation", gpc_vertex_list);
  
    fscanf(fp, "%d", &(p->contour[0].num_vertices));

    MALLOC(p->contour[0].vertex, p->contour[0].num_vertices
           * sizeof(gpc_vertex), "vertex creation", gpc_vertex);
           
    for (v= 0; v < p->contour[0].num_vertices; v++)
    {
      fscanf(fp, "%lf,%lf", &(p->contour[0].vertex[v].x),
                            &(p->contour[0].vertex[v].y));
     // printf(" %f %f, ", p->contour[0].vertex[v].x,
     //                    p->contour[0].vertex[v].y );
    }
}
