/**
 *
 */
#ifndef MESH_H
#define MESH_H

#include <vector>
#include <list>

#define MESH_CLASS_ID 125486

template <typename T>
struct vec3d
{
  T x,y,z;

  explicit vec3d() : x(0), y(0), z(0) {}
  vec3d(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  template <typename S>
  vec3d(const vec3d<S>& s) : x(T(s.x)), y(T(s.y)), z(T(s.z)) {}
  
  vec3d<T> operator-(const vec3d<T>& p) const
  {
	 return vec3d<T>(x - p.x, y - p.y, z - p.z);
  }
};

template <typename T>
inline T dot_product(const vec3d<T>& v1, const vec3d<T>& v2) 
{
	return (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z);
}

class mesh_t
{
public:
	
   //const int mesh_class_id = MESH_CLASS_ID;
   
   struct vertex_data {
      vertex_data() : n_tris(0) {}
      vec3d<double> p;
      size_t n_tris; // no. of triangles this vertex belongs to
      std::list<int> tri_indices; // triangles this vertex belongs to
   };

   struct triangle_data {
      int p0, p1, p2; // indices of the 3 vertices
      double E,F,G;
   };
   
   std::vector<vertex_data> vertices;
   std::vector<triangle_data> triangles;
   
   void mesh_t::put_vertices(const std::vector< vec3d<double> >& points)
	{
	   const size_t n_points = points.size();
	   vertices.resize(n_points);

	   // Put the points in the appropriate structure for the kd-tree
	   // and in the vertices collection for future indexed reference

	   for (size_t k=0; k<n_points; ++k)
	   {
		  // indexed vertices
		  vertex_data v;
		  v.p = vec3d<double>(points[k]);
		  v.n_tris = 0;
		  vertices[k] = v;
	   }
	}
	
	void mesh_t::add_triangle(int p0, int p1, int p2)
	{
	   triangle_data td;

	   td.p0 = p0;
	   td.p1 = p1;
	   td.p2 = p2;

	   triangles.push_back(td);
	   const int tri_idx = static_cast<int>(triangles.size() - 1);

	   vertices.at(p0).tri_indices.push_back(tri_idx);
	   vertices.at(p0).n_tris++;
	   vertices.at(p1).tri_indices.push_back(tri_idx);
	   vertices.at(p1).n_tris++;
	   vertices.at(p2).tri_indices.push_back(tri_idx);
	   vertices.at(p2).n_tris++;
	}
	
	/**
	 * Get the indices of the triangles neighboring to a given triangle,
	 * i.e., the triangles that share an edge with it (at most 3 if the mesh is manifold).
	 *
	 * WARNING: - This is currently returning repeating indices (FIXME!), please keep this in mind!
	 *          - This code assumes the mesh is manifold.
	 *
	 * @param tri
	 * @param neighs
	 */
	 void mesh_t::get_tri_imm_neighbors(int tri, std::vector<int>& neighs) const
	 {
		const triangle_data& t = triangles.at(tri);
		const vertex_data& p0 = vertices.at(t.p0), p1 = vertices.at(t.p1), p2 = vertices.at(t.p2);
		
		neighs.reserve(3);
		
		// p0,p1 and p0,p2
		for (std::list<int>::const_iterator it=p0.tri_indices.begin(); it!=p0.tri_indices.end(); ++it)
		{
			if (*it != tri)
			{
				const triangle_data& neigh_t = triangles.at(*it);
				if (neigh_t.p0 == t.p1 || neigh_t.p1 == t.p1 || neigh_t.p2 == t.p1)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
				if (neigh_t.p0 == t.p2 || neigh_t.p1 == t.p2 || neigh_t.p2 == t.p2)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
			}
		}
		
		// p1,p0 and p1,p2
		for (std::list<int>::const_iterator it=p1.tri_indices.begin(); it!=p1.tri_indices.end(); ++it)
		{
			if (*it != tri)
			{
				const triangle_data& neigh_t = triangles.at(*it);
				if (neigh_t.p0 == t.p0 || neigh_t.p1 == t.p0 || neigh_t.p2 == t.p0)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
				if (neigh_t.p0 == t.p2 || neigh_t.p1 == t.p2 || neigh_t.p2 == t.p2)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
			}
		}
		
		// p2,p0 and p2,p1
		for (std::list<int>::const_iterator it=p2.tri_indices.begin(); it!=p2.tri_indices.end(); ++it)
		{
			if (*it != tri)
			{
				const triangle_data& neigh_t = triangles.at(*it);
				if (neigh_t.p0 == t.p0 || neigh_t.p1 == t.p0 || neigh_t.p2 == t.p0)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
				if (neigh_t.p0 == t.p1 || neigh_t.p1 == t.p1 || neigh_t.p2 == t.p1)
					neighs.push_back(*it);
				//if (neighs.size()==3) 
				//	break;
			}
		}
	 }
};
   
#endif
