#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <vector>
#include <cmath>
#include <utility>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;
typedef boost::tuple<Point_2, int> Point_and_int;
typedef CGAL::Search_traits_2<Kernel> Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
                                    CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
                                    Traits_base> Traits;
typedef CGAL::Orthogonal_k_neighbor_search<Traits> K_neighbor_search;
typedef K_neighbor_search::Tree Tree;
typedef K_neighbor_search::Distance Distance;
int main()
{
  const unsigned int N = 2;
  std::vector<Point_2> points;
  points.push_back(Point_2(0,0));
  points.push_back(Point_2(1,1));
  points.push_back(Point_2(2,2));
  points.push_back(Point_2(3,3));
  points.push_back(Point_2(4,4));
  points.push_back(Point_2(1.2,0.8));

  points[0] = Point_2(0.1, 0);
  std::cout << points[0].x() << std::endl;

  std::vector<int> idx;
  idx.push_back(0);
  idx.push_back(1);
  idx.push_back(2);
  idx.push_back(3);
  idx.push_back(4);
  idx.push_back(5);

  // Insert number_of_data_points in the tree
  Tree tree(boost::make_zip_iterator(boost::make_tuple( points.begin(),idx.begin())),
            boost::make_zip_iterator(boost::make_tuple( points.end(),idx.end())));

  // Tree tree(points.begin(), points.end());
  // Initialize the search structure, and search all N points
  Point_2 query(0.9,0.9);
  K_neighbor_search search(tree, query, N);
   // report the N nearest neighbors and their distance
  // This should sort all N points by increasing distance from origin
  for(K_neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    std::cout << "i=" << boost::get<1>(it->first) << ", pos="<< boost::get<0>(it->first) << " dis="<< std::sqrt(it->second) << std::endl;
  return 0;
}