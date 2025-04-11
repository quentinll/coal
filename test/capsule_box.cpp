/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2014-2016, CNRS-LAAS and AIST
 *  Copyright (c) 2025, INRIA
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of CNRS-LAAS and AIST nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/** \authors Florent Lamiraux*/

#define BOOST_TEST_MODULE COAL_GEOMETRIC_SHAPES
#include <boost/test/included/unit_test.hpp>

#define CHECK_CLOSE_TO_0(x, eps) BOOST_CHECK_CLOSE((x + 1.0), (1.0), (eps))

#include <cmath>
#include "coal/distance.h"
#include "coal/math/transform.h"
#include "coal/collision.h"
#include "coal/collision_object.h"
#include "coal/shape/geometric_shapes.h"

#include "utility.h"

using coal::generateRandomTransform;
using coal::Scalar;
using coal::Transform3s;

template <int VecSize>
Eigen::Matrix<Scalar, VecSize, 1> generateRandomVector(Scalar min, Scalar max) {
  typedef Eigen::Matrix<Scalar, VecSize, 1> VecType;
  // Generate a random vector in the [min, max] range
  VecType v = VecType::Random() * (max - min) * 0.5 +
              VecType::Ones() * (max + min) * 0.5;
  return v;
}

Scalar generateRandomNumber(Scalar min, Scalar max) {
  Scalar r = static_cast<Scalar>(rand()) / static_cast<Scalar>(RAND_MAX);
  r = 2 * r - 1;
  const Scalar half(0.5);
  return r * (max - min) * half + (max + min) * half;
}

BOOST_AUTO_TEST_CASE(distance_capsule_box) {
  using coal::CollisionGeometryPtr_t;
  // Capsule of radius 2 and of height 4
  CollisionGeometryPtr_t capsuleGeometry(new coal::Capsule(2., 4.));
  // Box of size 1 by 2 by 4
  CollisionGeometryPtr_t boxGeometry(new coal::Box(1., 2., 4.));

  // Enable computation of nearest points
  coal::DistanceRequest distanceRequest(true, 0, 0);
  coal::DistanceResult distanceResult;

  // Test case 1: Capsule to the right of box
  coal::Transform3s tf1(coal::Vec3s(3., 0, 0));
  coal::Transform3s tf2;
  coal::CollisionObject capsule(capsuleGeometry, tf1);
  coal::CollisionObject box(boxGeometry, tf2);

  // test distance
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  // Nearest point on capsule
  coal::Vec3s o1(distanceResult.nearest_points[0]);
  // Nearest point on box
  coal::Vec3s o2(distanceResult.nearest_points[1]);
  BOOST_CHECK_CLOSE(distanceResult.min_distance, 0.5, 1e-1);
  BOOST_CHECK_CLOSE(o1[0], 1.0, 1e-1);
  CHECK_CLOSE_TO_0(o1[1], 1e-1);
  BOOST_CHECK_CLOSE(o2[0], 0.5, 1e-1);
  CHECK_CLOSE_TO_0(o2[1], 1e-1);

  // Test case 2: Capsule above box
  tf1 = coal::Transform3s(coal::Vec3s(0., 0., 8.));
  capsule.setTransform(tf1);

  // test distance
  distanceResult.clear();
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  o1 = distanceResult.nearest_points[0];
  o2 = distanceResult.nearest_points[1];

  BOOST_CHECK_CLOSE(distanceResult.min_distance, 2.0, 1e-1);
  CHECK_CLOSE_TO_0(o1[0], 1e-1);
  CHECK_CLOSE_TO_0(o1[1], 1e-1);
  BOOST_CHECK_CLOSE(o1[2], 4.0, 1e-1);

  CHECK_CLOSE_TO_0(o2[0], 1e-1);
  CHECK_CLOSE_TO_0(o2[1], 1e-1);
  BOOST_CHECK_CLOSE(o2[2], 2.0, 1e-1);

  // Test case 3: Rotated capsule behind box
  tf1.setTranslation(coal::Vec3s(-10., 0., 0.));
  tf1.setQuatRotation(
      coal::makeQuat(sqrt(Scalar(2)) / 2, 0, sqrt(Scalar(2)) / 2, 0));
  capsule.setTransform(tf1);

  // test distance
  distanceResult.clear();
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  o1 = distanceResult.nearest_points[0];
  o2 = distanceResult.nearest_points[1];

  BOOST_CHECK_CLOSE(distanceResult.min_distance, 5.5, 1e-1);
  BOOST_CHECK_CLOSE(o1[0], -6, 1e-2);
  CHECK_CLOSE_TO_0(o1[1], 1e-1);
  CHECK_CLOSE_TO_0(o1[2], 1e-1);
  BOOST_CHECK_CLOSE(o2[0], -0.5, 1e-2);
  CHECK_CLOSE_TO_0(o2[1], 1e-1);
  CHECK_CLOSE_TO_0(o2[2], 1e-1);

  // Test case 4: Rotated capsule behind box with offset
  tf1 = coal::Transform3s(
      coal::makeQuat(sqrt(Scalar(2)) / 2, 0, sqrt(Scalar(2)) / 2, 0),
      coal::Vec3s(Scalar(-10.), Scalar(0.8), Scalar(1.5)));
  capsule.setTransform(tf1);

  // test distance
  distanceResult.clear();
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  o1 = distanceResult.nearest_points[0];
  o2 = distanceResult.nearest_points[1];

  BOOST_CHECK_CLOSE(distanceResult.min_distance, 5.5, 1e-2);
  BOOST_CHECK_CLOSE(o1[0], -6, 1e-2);
  BOOST_CHECK_CLOSE(o1[1], 0.8, 1e-1);
  BOOST_CHECK_CLOSE(o1[2], 1.5, 1e-2);
  BOOST_CHECK_CLOSE(o2[0], -0.5, 1e-2);
  BOOST_CHECK_CLOSE(o2[1], 0.8, 1e-1);
  BOOST_CHECK_CLOSE(o2[2], 1.5, 1e-2);
}

BOOST_AUTO_TEST_CASE(distance_segment_box) {
  using coal::CollisionGeometryPtr_t;
  using coal::CollisionObject;
  using coal::DistanceRequest;
  using coal::DistanceResult;
  using coal::Vec3s;
  CollisionGeometryPtr_t segmentGeometry(new coal::Capsule(0.1, 1));
  CollisionGeometryPtr_t boxGeometry(new coal::Box(1., 1., 1.));
  coal::Transform3s tf1, tf2;
  // Test case 1: Segment to the right of box is colliding
  tf1.setTranslation(coal::Vec3s(0.9, 0., 0.));
  tf1.setQuatRotation(
      coal::makeQuat(sqrt(Scalar(2)) / 2, 0, sqrt(Scalar(2)) / 2, 0));
  tf2.setIdentity();
  CollisionObject segment(segmentGeometry, tf1);
  CollisionObject box(boxGeometry, tf2);
  coal::DistanceRequest distanceRequest(true, 0, 0);
  coal::DistanceResult distanceResult;
  coal::distance(&segment, &box, distanceRequest, distanceResult);
  BOOST_CHECK_CLOSE(distanceResult.min_distance, -0.2, 1e-1);
  BOOST_CHECK_CLOSE(distanceResult.nearest_points[0][0], 0.3, 1e-1);
  BOOST_CHECK_CLOSE(distanceResult.nearest_points[1][0], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(distanceResult.normal[0], -1, 1e-1);
  CHECK_CLOSE_TO_0(distanceResult.normal[1], 1e-1);
  CHECK_CLOSE_TO_0(distanceResult.normal[2], 1e-1);
  // Test case 2: Segment to the right of box is not colliding but capsule is
  // colliding
  tf1.setTranslation(coal::Vec3s(1.01, 0., 0.));
  segment.setTransform(tf1);
  distanceResult.clear();
  coal::distance(&segment, &box, distanceRequest, distanceResult);
  BOOST_CHECK_CLOSE(distanceResult.min_distance, -0.09, 1e-1);
  BOOST_CHECK_CLOSE(distanceResult.nearest_points[0][0], 0.41, 1e-1);
  BOOST_CHECK_CLOSE(distanceResult.nearest_points[1][0], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(distanceResult.normal[0], -1, 1e-1);
  CHECK_CLOSE_TO_0(distanceResult.normal[1], 1e-1);
  CHECK_CLOSE_TO_0(distanceResult.normal[2], 1e-1);
}

BOOST_AUTO_TEST_CASE(contact_points_sphere_box) {
  using coal::CollisionGeometryPtr_t;
  // We use a capsule of radius 1 and height 0 to represent a sphere
  CollisionGeometryPtr_t capsuleGeometry(new coal::Capsule(1., 0.));
  // Box of size 1 by 1 by 1
  CollisionGeometryPtr_t boxGeometry(new coal::Box(1., 1., 1.));

  // Enable computation of nearest points
  coal::DistanceRequest distanceRequest(true, 0, 0);
  coal::DistanceResult distanceResult;

  // Test case 1: Capsule collides to the right of box
  coal::Transform3s tf1(coal::Vec3s(1.4, 0, 0));
  coal::Transform3s tf2;
  coal::CollisionObject capsule(capsuleGeometry, tf1);
  coal::CollisionObject box(boxGeometry, tf2);

  // test distance
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  // Nearest point on capsule
  coal::Vec3s o1(distanceResult.nearest_points[0]);
  // Nearest point on box
  coal::Vec3s o2(distanceResult.nearest_points[1]);
  coal::Vec3s normal = distanceResult.normal;
  BOOST_CHECK_CLOSE(distanceResult.min_distance, -0.1, 1e-1);
  BOOST_CHECK_CLOSE(o1[0], 0.4, 1e-1);
  CHECK_CLOSE_TO_0(o1[1], 1e-1);
  BOOST_CHECK_CLOSE(o2[0], 0.5, 1e-1);
  CHECK_CLOSE_TO_0(o2[1], 1e-1);
  BOOST_CHECK_CLOSE(normal[0], -1, 1e-1);
  CHECK_CLOSE_TO_0(normal[1], 1e-1);
  CHECK_CLOSE_TO_0(normal[2], 1e-1);

  // Test case 2: Capsule above and collide with box
  tf1 = coal::Transform3s(coal::Vec3s(0., 0., 1.4));
  capsule.setTransform(tf1);

  // test distance
  distanceResult.clear();
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  o1 = distanceResult.nearest_points[0];
  o2 = distanceResult.nearest_points[1];
  normal = distanceResult.normal;

  BOOST_CHECK_CLOSE(distanceResult.min_distance, -0.1, 1e-1);
  CHECK_CLOSE_TO_0(o1[0], 1e-1);
  CHECK_CLOSE_TO_0(o1[1], 1e-1);
  BOOST_CHECK_CLOSE(o1[2], 0.4, 1e-1);
  CHECK_CLOSE_TO_0(o2[0], 1e-1);
  CHECK_CLOSE_TO_0(o2[1], 1e-1);
  BOOST_CHECK_CLOSE(o2[2], 0.5, 1e-1);
  CHECK_CLOSE_TO_0(normal[0], 1e-1);
  CHECK_CLOSE_TO_0(normal[1], 1e-1);
  BOOST_CHECK_CLOSE(normal[2], -1, 1e-1);

  // Test case 3: Capsule collides with box on an edge
  tf1 = coal::Transform3s(
      coal::Vec3s(0.5 + std::sqrt(0.405), 0.5 + std::sqrt(0.405), 0.));
  capsule.setTransform(tf1);

  // test distance
  distanceResult.clear();
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  o1 = distanceResult.nearest_points[0];
  o2 = distanceResult.nearest_points[1];
  normal = distanceResult.normal;
  BOOST_CHECK_CLOSE(distanceResult.min_distance, -0.1, 1e-1);
  BOOST_CHECK_CLOSE(o1[0], 0.5 + std::sqrt(0.405) - std::sqrt(0.5), 1e-1);
  BOOST_CHECK_CLOSE(o1[1], 0.5 + std::sqrt(0.405) - std::sqrt(0.5), 1e-1);
  BOOST_CHECK_CLOSE(o2[0], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(o2[1], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(normal[0], -1 / std::sqrt(2), 1e-1);
  BOOST_CHECK_CLOSE(normal[1], -1 / std::sqrt(2), 1e-1);
  CHECK_CLOSE_TO_0(normal[2], 1e-1);

  // Test case 4: Capsule collides with box on a corner
  tf1 = coal::Transform3s(coal::Vec3s(
      0.5 + std::sqrt(0.27), 0.5 + std::sqrt(0.27), 0.5 + std::sqrt(0.27)));
  capsule.setTransform(tf1);

  // test distance
  distanceResult.clear();
  coal::distance(&capsule, &box, distanceRequest, distanceResult);
  o1 = distanceResult.nearest_points[0];
  o2 = distanceResult.nearest_points[1];
  normal = distanceResult.normal;

  BOOST_CHECK_CLOSE(distanceResult.min_distance, -0.1, 1e-1);
  BOOST_CHECK_CLOSE(o1[0], 0.5 + std::sqrt(0.27) - std::sqrt(1. / 3.), 1e-1);
  BOOST_CHECK_CLOSE(o1[1], 0.5 + std::sqrt(0.27) - std::sqrt(1. / 3.), 1e-1);
  BOOST_CHECK_CLOSE(o1[2], 0.5 + std::sqrt(0.27) - std::sqrt(1. / 3.), 1e-1);
  BOOST_CHECK_CLOSE(o2[0], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(o2[1], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(o2[2], 0.5, 1e-1);
  BOOST_CHECK_CLOSE(normal[0], -1 / std::sqrt(3.), 1e-1);
  BOOST_CHECK_CLOSE(normal[1], -1 / std::sqrt(3.), 1e-1);
  BOOST_CHECK_CLOSE(normal[2], -1 / std::sqrt(3.), 1e-1);
}