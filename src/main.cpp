// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/cameras/Camera_Pinhole.hpp"
#include "openMVG/cameras/Camera_Pinhole_Radial.hpp"
#include "openMVG/features/feature.hpp"
#include "openMVG/features/sift/SIFT_Anatomy_Image_Describer.hpp"
#include "openMVG/features/svg_features.hpp"
#include "openMVG/geometry/pose3.hpp"
#include "openMVG/image/image_io.hpp"
#include "openMVG/image/image_concat.hpp"
#include "openMVG/matching/indMatchDecoratorXY.hpp"
#include "openMVG/matching/regions_matcher.hpp"
#include "openMVG/matching/svg_matches.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "openMVG/sfm/sfm_data_BA.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"
#include "openMVG/sfm/sfm_data_io.hpp"

#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <iostream>
#include <string>
#include <utility>

using namespace openMVG;
using namespace openMVG::image;
using namespace openMVG::matching;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::sfm;
using namespace std;

/// Read intrinsic K matrix from a file (ASCII)
/// F 0 ppx
/// 0 F ppy
/// 0 0 1
bool readIntrinsic(const std::string & fileName, Mat3 & K);

/// Show :
///  how computing an essential with know internal calibration matrix K
///  how refine the camera motion, focal and structure with Bundle Adjustment
///   way 1: independent cameras [R|t|f] and structure
///   way 2: independent cameras motion [R|t], shared focal [f] and structure
int main() {

  const std::string sInputDir = "./imageData/SceauxCastle/";

  Image<RGBColor> image;


  const string jpg_filenameL = sInputDir + "100_7101.jpg";
  const string jpg_filenameR = sInputDir + "100_7102.jpg";

  Image<unsigned char> imageL, imageR;
  ReadImage(jpg_filenameL.c_str(), &imageL);
  ReadImage(jpg_filenameR.c_str(), &imageR);

  //--
  // Detect regions thanks to an image_describer
  //--
  using namespace openMVG::features;
  std::unique_ptr<Image_describer> image_describer(new SIFT_Anatomy_Image_describer);
  std::map<IndexT, std::unique_ptr<features::Regions>> regions_perImage;

  image_describer->Describe(imageL, regions_perImage[0]);
  image_describer->Describe(imageR, regions_perImage[1]);



//常用特征的 Regions
/*
SIFT Keypoint
using SIFT_Regions = Scalar_Regions<SIOPointFeature, unsigned char, 128>;

AKAZE Keypoint (with a float descriptor)
using AKAZE_Float_Regions = Scalar_Regions<SIOPointFeature, float, 64>;

 AKAZE Keypoint (with a LIOP descriptor)
using AKAZE_Liop_Regions = Scalar_Regions<SIOPointFeature, unsigned char, 144>;

 AKAZE Keypoint (with a binary descriptor saved in an uchar array)
using AKAZE_Binary_Regions = Binary_Regions<SIOPointFeature, 64>;

Regions 定义, 包括所有的features 和descriptors
Describe an image a set of regions (position, ...) + attributes
Each region is described by a set of attributes (descriptor)

SIOPointFeature： Add scale and orientation description to basis PointFeature.

 */


  const SIFT_Regions* regionsL = dynamic_cast<SIFT_Regions*>(regions_perImage.at(0).get());
  const SIFT_Regions* regionsR = dynamic_cast<SIFT_Regions*>(regions_perImage.at(1).get());

  const PointFeatures
    featsL = regions_perImage.at(0)->GetRegionsPositions(),
    featsR = regions_perImage.at(1)->GetRegionsPositions();

  // Show both images side by side
  {
    Image<unsigned char> concat;
    ConcatH(imageL, imageR, concat);
    string out_filename = "01_concat.jpg";
    WriteImage(out_filename.c_str(), concat);
  }

  //- Draw features on the two image (side by side)
  {
    Features2SVG
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regionsL->Features(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regionsR->Features(),
      "02_features.svg",
      false
    );
  }

  std::vector<IndMatch> vec_PutativeMatches;
  //-- Perform matching -> find Nearest neighbor, filtered with Distance ratio
  {
    // Find corresponding points
    matching::DistanceRatioMatch(
      0.8, matching::BRUTE_FORCE_L2,
      *regions_perImage.at(0).get(),
      *regions_perImage.at(1).get(),
      vec_PutativeMatches);

    //对重复的特征点进行筛选
    IndMatchDecorator<float> matchDeduplicator( vec_PutativeMatches, featsL, featsR);
    matchDeduplicator.getDeduplicated(vec_PutativeMatches);

    std::cout
      << regions_perImage.at(0)->RegionCount() << " #Features on image A" << std::endl
      << regions_perImage.at(1)->RegionCount() << " #Features on image B" << std::endl
      << vec_PutativeMatches.size() << " #matches with Distance Ratio filter" << std::endl;

    // Draw correspondences after Nearest Neighbor ratio filter
    const bool bVertical = false;
    Matches2SVG
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regionsL->GetRegionsPositions(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regionsR->GetRegionsPositions(),
      vec_PutativeMatches,
      "03_Matches.svg",
      bVertical
    );
  }


  // Essential geometry filtering of putative matches
  {
    Mat3 K;
    //read K from file
    if (!readIntrinsic(stlplus::create_filespec(sInputDir,"K","txt"), K))
    {
      std::cerr << "Cannot read intrinsic parameters." << std::endl;
      return EXIT_FAILURE;
    }

    const Pinhole_Intrinsic
      camL(imageL.Width(), imageL.Height(), K(0,0), K(0,2), K(1,2)),
      camR(imageR.Width(), imageR.Height(), K(0,0), K(0,2), K(1,2));

    //A. prepare the corresponding putatives points
    Mat xL(2, vec_PutativeMatches.size());
    Mat xR(2, vec_PutativeMatches.size());  //xL[i]与xR[i]相对应
    for (size_t k = 0; k < vec_PutativeMatches.size(); ++k)  {
      const PointFeature & imaL = featsL[vec_PutativeMatches[k].i_];
      const PointFeature & imaR = featsR[vec_PutativeMatches[k].j_];
      xL.col(k) = imaL.coords().cast<double>();
      xR.col(k) = imaR.coords().cast<double>();
    }

    //B. Compute the relative pose thanks to a essential matrix estimation
    const std::pair<size_t, size_t> size_imaL(imageL.Width(), imageL.Height());
    const std::pair<size_t, size_t> size_imaR(imageR.Width(), imageR.Height());
    RelativePose_Info relativePose_info;
    if (!robustRelativePose(&camL, &camR, xL, xR, relativePose_info, size_imaL, size_imaR, 256))
        //最后一个参数是最大迭代次数
    {
      std::cerr << " /!\\ Robust relative pose estimation failure."
        << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "\nFound an Essential matrix:\n"
      << "\tprecision: " << relativePose_info.found_residual_precision << " pixels\n"
      << "\t#inliers: " << relativePose_info.vec_inliers.size() << "\n"
      << "\t#matches: " << vec_PutativeMatches.size()
      << std::endl;

    // Show Essential validated point
    const bool bVertical = false;
    InlierMatches2SVG    //用本征矩阵的inliers 属性剔除错误的匹配
    (
      jpg_filenameL,
      {imageL.Width(), imageL.Height()},
      regionsL->GetRegionsPositions(),
      jpg_filenameR,
      {imageR.Width(), imageR.Height()},
      regionsR->GetRegionsPositions(),
      vec_PutativeMatches,
      relativePose_info.vec_inliers,    //存放的是inlier在 vec_PutativeMatches[] 中的index
      "04_ACRansacEssential.svg",
      bVertical
    );

    std::cout << std::endl
      << "-- Rotation|Translation matrices: --" << "\n"
      << relativePose_info.relativePose.rotation() << "\n\n"
      << relativePose_info.relativePose.translation() << "\n" << std::endl;

    //C. Triangulate and check valid points
    // invalid points that do not respect cheirality are discarded (removed
    //  from the list of inliers).

    std::cout << "Which BA do you want ?\n"
      << "\t 1: Refine [X],[f,ppx,ppy,R|t] (individual cameras)\n"
      << "\t 2: Refine [X],[R|t], shared [f, ppx, ppy]\n"
      << "\t 3: Refine [X],[R|t], shared brown K3 distortion model [f,ppx,ppy,k1,k2,k3]\n" << std::endl;
    int iBAType = -1;
    std::cin >> iBAType;
    const bool bSharedIntrinsic = (iBAType == 2 || iBAType == 3) ? true : false;


    //搭建SFM进行全局的BA

    // Setup a SfM scene with two view corresponding the pictures
    SfM_Data tiny_scene;
    tiny_scene.views[0].reset(new View("", 0, bSharedIntrinsic ? 0 : 1,    0,   imageL.Width(), imageL.Height()));
    tiny_scene.views[1].reset(new View("", 1, bSharedIntrinsic ? 0 : 1,    1,   imageR.Width(), imageR.Height()));
    /* View(
    const std::string & sImgPath = "",
    IndexT view_id = UndefinedIndexT,
    IndexT intrinsic_id = UndefinedIndexT,
    IndexT pose_id = UndefinedIndexT,
    IndexT width = UndefinedIndexT, IndexT height ） */


    // Setup intrinsics camera data
    switch (iBAType)
    {
      case 1: // Each view use it's own pinhole camera intrinsic
        tiny_scene.intrinsics[0].reset(new Pinhole_Intrinsic(imageL.Width(), imageL.Height(), K(0, 0), K(0, 2), K(1, 2)));
        tiny_scene.intrinsics[1].reset(new Pinhole_Intrinsic(imageR.Width(), imageR.Height(), K(0, 0), K(0, 2), K(1, 2)));
        break;
      case 2: // Shared pinhole camera intrinsic
        tiny_scene.intrinsics[0].reset(new Pinhole_Intrinsic(imageL.Width(), imageL.Height(), K(0, 0), K(0, 2), K(1, 2)));
        break;
      case 3: // Shared pinhole camera intrinsic with radial K3 distortion
        tiny_scene.intrinsics[0].reset(new Pinhole_Intrinsic_Radial_K3(imageL.Width(), imageL.Height(),
                                                                                                             K(0, 0), K(0, 2), K(1, 2)));
        break;
      default:
        std::cerr << "Invalid input number" << std::endl;
        return EXIT_FAILURE;
    }

    // Setup poses camera data
    const Pose3 pose0 = tiny_scene.poses[tiny_scene.views[0]->id_pose] = Pose3(Mat3::Identity(), Vec3::Zero());
    const Pose3 pose1 = tiny_scene.poses[tiny_scene.views[1]->id_pose] = relativePose_info.relativePose;


    // Init structure by inlier triangulation     分别得到两帧的投影矩阵
    const Mat34 P1 = tiny_scene.intrinsics[tiny_scene.views[0]->id_intrinsic]->get_projective_equivalent(pose0);
    const Mat34 P2 = tiny_scene.intrinsics[tiny_scene.views[1]->id_intrinsic]->get_projective_equivalent(pose1);
    Landmarks & landmarks = tiny_scene.structure;


    for (size_t i = 0; i < relativePose_info.vec_inliers.size(); ++i)  {
      const SIOPointFeature & LL = regionsL->Features()[vec_PutativeMatches[relativePose_info.vec_inliers[i]].i_];
      const SIOPointFeature & RR = regionsR->Features()[vec_PutativeMatches[relativePose_info.vec_inliers[i]].j_];

      // Point triangulation   对每一对正确的匹配，进行三角化
      Vec3 X;
      TriangulateDLT(P1, LL.coords().cast<double>().homogeneous(),
                                P2, RR.coords().cast<double>().homogeneous(), &X);


      // Reject point that is behind the camera
      if (pose0.depth(X) < 0 && pose1.depth(X) < 0)   //不应该取两个深度值都大于0吗？？
        continue;

     // Landmarks structure;     //std::map<IndexT,  Landmark>    indexed by their TrackId  点云的id
      /* struct Landmark
      {
        Vec3 X;      a 3D point, with its 2d observations
        Observations obs;  //std::map<IndexT, Observation>   /// Observations are indexed by their View_id
           struct Observation
          {
             Vec2 x;
             IndexT id_feat;
           }
      }*/

      // Add a new landmark (3D point with it's 2d observations)    添加到点云中
      landmarks[i].obs[tiny_scene.views[0]->id_view] =
              Observation(LL.coords().cast<double>(),   vec_PutativeMatches[relativePose_info.vec_inliers[i]].i_);

      landmarks[i].obs[tiny_scene.views[1]->id_view] =
              Observation(RR.coords().cast<double>(),  vec_PutativeMatches[relativePose_info.vec_inliers[i]].j_);

      landmarks[i].X = X;
    }
    Save(tiny_scene, "EssentialGeometry_start.ply", ESfM_Data(ALL));



    //D. Perform Bundle Adjustment of the scene

    Bundle_Adjustment_Ceres  bundle_adjustment_obj;
    bundle_adjustment_obj.Adjust(tiny_scene,
      Optimize_Options(
        Intrinsic_Parameter_Type::ADJUST_ALL,
        Extrinsic_Parameter_Type::ADJUST_ALL,
        Structure_Parameter_Type::ADJUST_ALL));

    Save(tiny_scene, "EssentialGeometry_refined.ply", ESfM_Data(ALL));
  }

  return EXIT_SUCCESS;
}



bool readIntrinsic(const std::string & fileName, Mat3 & K)
{
  // Load the K matrix
  ifstream in;
  in.open( fileName.c_str(), ifstream::in);
  if (in.is_open())  {
    for (int j=0; j < 3; ++j)
      for (int i=0; i < 3; ++i)
        in >> K(j,i);
  }
  else  {
    std::cerr << std::endl
      << "Invalid input K.txt file" << std::endl;
    return false;
  }
  return true;
}
