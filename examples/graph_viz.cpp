#include <geometry_msgs/Point.h>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

#include <srrg_solver/solver_core/factor_graph.h>
#include <srrg_solver/solver_core/solver.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_point3.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_se3.h>
#include <srrg_solver/variables_and_factors/types_3d/variable_se3_ad.h>

#include <srrg_solver/solver_core/instances.h>
#include <srrg_system_utils/parse_command_line.h>

#include "variables_and_factors/imu_preintegration_factor.h"
#include "variables_and_factors/instances.h"

void initTypes() {
  srrg2_solver::variables_and_factors_imu_registerTypes();
}

using BAPoseVariableType  = srrg2_solver::VariableSE3QuaternionRight;
using EBAPoseVariableType = srrg2_solver::VariableSE3QuaternionRightAD;
using PointVariableType   = srrg2_solver::VariablePoint3;

const char* banner[] = {"graph visualizer", 0};

std::map<std::string, srrg2_solver::FactorGraphPtr> graphs;

int main(int argc, char** argv) {
  initTypes();

  ros::init(argc, argv, "graph_viz");
  ros::NodeHandle n;
  ros::Rate r(1);
  ros::Publisher marker_pub =
    n.advertise<visualization_msgs::MarkerArray>("visualization_marker_array", 1);

  std::set<srrg2_solver::VariableBase::Id> gt_vars;
  std::set<srrg2_solver::FactorBase::Id> gt_factors;

  srrg2_core::ParseCommandLine cmd_line(argv, banner);
  cmd_line.parse();

  if (cmd_line.lastParsedArgs().empty()) {
    throw std::runtime_error("no input file");
  }

  for (const auto& input_file : cmd_line.lastParsedArgs()) {
    srrg2_solver::FactorGraphPtr graph = srrg2_solver::FactorGraph::read(input_file);
    if (!graph)
      continue;

    graphs[input_file] = graph;
  }

  if (graphs.empty())
    return -1;

  uint32_t poses_shape  = visualization_msgs::Marker::CUBE_LIST;
  uint32_t points_shape = visualization_msgs::Marker::SPHERE_LIST;

  srrg2_solver::VariableSE3ExpMapRight* pose_var;
  srrg2_solver::Isometry3f pose(srrg2_solver::Isometry3f::Identity());
  srrg2_solver::Vector3f point(srrg2_solver::Vector3f::Zero());
  int num_graph = 0;

  while (marker_pub.getNumSubscribers() < 1) {
    ROS_WARN_ONCE("Please create a subscriber to the marker");
    sleep(1);
  }

  visualization_msgs::MarkerArray graphs_list;

  std::cerr << graphs.size() << " graphs read" << std::endl;

  for (const auto& it_graph : graphs) {
    visualization_msgs::Marker poses_list;
    poses_list.header.frame_id = "map";
    poses_list.header.stamp    = ros::Time::now();
    poses_list.type            = poses_shape;
    poses_list.ns              = it_graph.first;
    poses_list.id              = num_graph;

    poses_list.action = visualization_msgs::Marker::MODIFY;

    poses_list.pose.orientation.w = 1;
    poses_list.scale.x            = 0.5;
    poses_list.scale.y            = 0.5;
    poses_list.scale.z            = 0.5;
    poses_list.color.a            = 1.0;

    visualization_msgs::Marker points_list;
    points_list.header.frame_id = "map";
    points_list.header.stamp    = ros::Time::now();
    points_list.type            = points_shape;
    points_list.ns              = it_graph.first + "_points";
    points_list.id              = num_graph + 1;

    points_list.action = visualization_msgs::Marker::MODIFY;

    points_list.pose.orientation.w = 1;
    points_list.scale.x            = 0.1;
    points_list.scale.y            = 0.1;
    points_list.scale.z            = 0.1;
    points_list.color.a            = 1.0;

    switch (num_graph / 2) {
      case 0:
        poses_list.color.r  = 1;
        points_list.color.r = 0.5;
        break;

      case 1:
        poses_list.color.g  = 1;
        points_list.color.g = 0.5;
        break;

      case 2:
        poses_list.color.b  = 1;
        points_list.color.b = 0.5;
        break;

      case 3:
        points_list.color.r = 1;
        points_list.color.g = 1;
        points_list.color.r = 0.5;
        points_list.color.g = 0.5;
        break;

      case 4:
        points_list.color.g = 1;
        points_list.color.b = 1;
        points_list.color.g = 0.5;
        points_list.color.b = 0.5;
        break;

      case 5:
        points_list.color.r = 1;
        points_list.color.b = 1;
        points_list.color.r = 0.5;
        points_list.color.b = 0.5;
        break;

      default:
        break;
    }

    const auto& graph = it_graph.second;
    for (const auto& it_var : graph->variables()) {
      pose_var = dynamic_cast<srrg2_solver::VariableSE3ExpMapRight*>(it_var.second);
      if (pose_var) {
        pose = pose_var->estimate();

        geometry_msgs::Point p;

        p.x = pose.translation().x();
        p.y = pose.translation().y();
        p.z = pose.translation().z();
        poses_list.points.push_back(p);
        continue;
      }
    }
    if (poses_list.points.size())
      graphs_list.markers.push_back(poses_list);
    if (points_list.points.size())
      graphs_list.markers.push_back(points_list);
    num_graph += 2;
  }

  while (ros::ok()) {
    marker_pub.publish(graphs_list);
  }

  return 0;
}