// [[Rcpp::depends(BH)]]

// Enable C++11 via this plugin to suppress 'long long' errors
// [[Rcpp::plugins("cpp11")]]

// Some of this code based on http://gallery.rcpp.org/articles/Rtree-examples/

#include <vector>
#include <Rcpp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
using namespace Rcpp;

// Mnemonics
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, 2, bg::cs::cartesian> point_t;
typedef bg::model::box<point_t> box;
typedef std::pair<point_t, unsigned int> value_t;

class RTreeCpp2 {
public:

  // Constructor, creates R-Tree on points in mat
  // mat must have 2 columns with x and y coordinates of points
  RTreeCpp2(NumericMatrix mat) {

    coords = clone(mat);
    int size = mat.nrow();
    NumericVector x = mat(_,0);
    NumericVector y = mat(_,1);
    std::vector<value_t> point_value_pairs;

    // Create vector of point-index pairs
    for(int i = 0; i < size; ++i) {
      point_t point_t(x[i], y[i]);
      point_value_pairs.push_back(std::make_pair(point_t, static_cast<unsigned int>(i)));
    }

    // Create rtree using packing algorithm
    bgi::rtree<value_t, bgi::quadratic<16> > rtree_new(point_value_pairs);
    rtree_ = rtree_new;
  }

  // Within-box lookup
  // box_vec is [lrx, lry, urx, ury] coordinates
  std::vector<int> intersects(NumericVector box_vec) {
    box query_box(point_t(box_vec[0], box_vec[1]), point_t(box_vec[2], box_vec[3]));
    std::vector<value_t> result_n;
    rtree_.query(bgi::intersects(query_box), std::back_inserter(result_n));
    std::vector<int> indexes;
    std::vector<value_t>::iterator itr;
    for (itr = result_n.begin(); itr != result_n.end(); ++itr) {
      value_t value = *itr;
      indexes.push_back(value.second);
    }
    return indexes;
  }

  // Get euclidean distances between point and points referenced by their indexes
  std::vector<double> get_distances(NumericVector point_vec, std::vector<int> indexes) {
    point_t start = point_t(point_vec[0], point_vec[1]);
    std::vector<double> distances;
    for (int i = 0; i < indexes.size(); i++) {
      point_t end = point_t(coords(indexes[i],0), coords(indexes[i],1));
      double dist = bg::distance(start, end);
      distances.push_back(dist);
    }
    return(distances);
  }

  // Get indices of points within distance of point
  // Note: Returns R indices (starting at 1, not 0!)
  std::vector<int> within_distance(NumericVector point_vec, double distance) {
    std::vector<int> indist_indexes;
    NumericVector box_vec = NumericVector::create(point_vec[0]-distance, point_vec[1]-distance,
                                                  point_vec[0]+distance, point_vec[1]+distance);

    std::vector<int> inbox_indexes = intersects(box_vec);
    if (inbox_indexes.size()==0) {
      return(indist_indexes);
    }

    std::vector<double> inbox_distances = get_distances(point_vec, inbox_indexes);
    for (int i=0; i<inbox_distances.size(); i++) {
      if (inbox_distances[i] <= distance) {
        indist_indexes.push_back(inbox_indexes[i] + 1);
      }
    }
    return(indist_indexes);
  }

  // Multi point version of within_distance
  // Note: Returns R indices (starting at 1, not 0!)
  List within_distance_list(NumericMatrix point_mat, double distance) {
    List indist_indexes_ls(point_mat.nrow());
    for (int i=0; i < point_mat.nrow(); i++){
      NumericVector point_vec = point_mat(i,_);
      std::vector<int> indist_indexes = within_distance(point_vec, distance);
      indist_indexes_ls(i) = wrap(indist_indexes);
    }
    return(indist_indexes_ls);
  }

  // Get indices of points within boxes defined yb ditance of a point.
  // Note: Returns R indices (starting at 1, not 0!)
  std::vector<int> within_box(NumericVector point_vec, double dx, double dy) {
    NumericVector box_vec = NumericVector::create(point_vec[0]-dx, point_vec[1]-dy,
                                                  point_vec[0]+dx, point_vec[1]+dy);
    std::vector<int> inbox_indexes = intersects(box_vec);
    return(inbox_indexes);
  }

  // Note: Returns R indices (starting at 1, not 0!)
  List within_box_list(NumericMatrix point_mat, double dx, double dy) {
    List inbox_indexes_ls(point_mat.nrow());
    for (int i=0; i < point_mat.nrow(); i++){
      NumericVector point_vec = point_mat(i,_);
      std::vector<int> indist_indexes = within_box(point_vec, dx, dy);
      inbox_indexes_ls(i) = wrap(indist_indexes);
    }
    return(inbox_indexes_ls);
  }


  // KNN
  // Note: Returns R indices (starting at 1, not 0!)
  std::vector<int> knn(NumericVector point, unsigned int n) {
    std::vector<value_t> result_n;
    rtree_.query(bgi::nearest(point_t(point[0], point[1]), n), std::back_inserter(result_n));
    std::vector<int> indexes;
    std::vector<value_t>::iterator itr;
    for ( itr = result_n.begin(); itr != result_n.end(); ++itr ) {
      value_t value = *itr;
      indexes.push_back(value.second + 1);
    }
    return indexes;
  }

  // Multi point version of KNN
  // Note: Returns R indices (starting at 1, not 0!)
  List knn_list(NumericMatrix point_mat, unsigned int n) {
    List knn_indexes_ls(point_mat.nrow());
    for (int i=0; i < point_mat.nrow(); i++){
      NumericVector point_vec = point_mat(i,_);
      std::vector<int> knn_indexes = knn(point_vec, n);
      knn_indexes_ls(i) = wrap(knn_indexes);
    }
    return(knn_indexes_ls);
  }

private:
  bgi::rtree<value_t, bgi::quadratic<16> > rtree_;
  NumericMatrix coords;
};

RCPP_MODULE(rtreecpp2) {
  class_<RTreeCpp2>( "RTreeCpp2" )

  .constructor<NumericMatrix>()

  .method( "intersects", &RTreeCpp2::intersects)
  .method( "within_distance", &RTreeCpp2::within_distance)
  .method( "within_distance_list", &RTreeCpp2::within_distance_list)
  .method( "within_box", &RTreeCpp2::within_box)
  .method( "within_box_list", &RTreeCpp2::within_box_list)
  .method( "knn", &RTreeCpp2::knn)
  .method( "knn_list", &RTreeCpp2::knn_list)
  ;
}
