#include <thread>
#include <unistd.h>
// ia system utils
#include <srrg_system_utils/parse_command_line.h>
// ia solver stuff
#include "srrg_solver/solver_core/solver.h"
#include "srrg_solver/variables_and_factors/types_3d/all_types.h"
#include "variables_and_factors/imu_preintegration_factor.h"
#include "variables_and_factors/variable_se3_expmap_right.h"

#include <srrg_solver/solver_core/instances.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>

// ia viewer stuff
#include <srrg_qgl_viewport/viewer_core_shared_qgl.h>
#include <srrg_solver/solver_core/factor_graph.h>

#include "variables_and_factors/instances.h"

const std::string exe_name = "test_graph_viewer|";
#define LOG std::cerr << exe_name

using namespace srrg2_core;
using namespace srrg2_solver;
using namespace srrg2_qgl_viewport;

void initTypes() {
  variables_and_factors_3d_registerTypes();
  variables_and_factors_imu_registerTypes();
}

const char* banner[] = {"showing a standard 3d PoseGraph.", "usage: <exe> <filename>", 0};

std::map<std::string, FactorGraphPtr> graphs;

bool gl_list_generated = false;
void viewGraph(ViewerCanvasPtr canvas_) {
  std::cerr << "canvas_ok" << std::endl;
  while (ViewerCoreSharedQGL::isRunning()) {
    if (!canvas_->_setup()) {
      usleep(10000);
      continue;
    }
    if (!gl_list_generated) {
      for (auto& g_it : graphs) {
        std::string name = g_it.first;
        auto graph       = g_it.second;
        canvas_->createList(name);
        canvas_->beginList(name);
        for (auto v : graph->variables()) {
          v.second->_drawImpl(canvas_);
        }

        for (auto f : graph->factors()) {
          f.second->_drawImpl(canvas_);
        }
        canvas_->endList();
      }
      gl_list_generated = true;
    } else {
      for (auto& g_it : graphs) {
        std::string name = g_it.first;
        canvas_->callList(name);
      }
    }
    canvas_->flush();
    sleep(1);
  }
}

// to read gtsam's shittings
void parseCSVFile(const std::string& csv_file, FactorGraphPtr graph) {
  std::ifstream file(csv_file);

  std::string line;

  VariableSE3ExpMapRightAD* prev_pose = new VariableSE3ExpMapRightAD();

  prev_pose->setGraphId(0);
  prev_pose->setStatus(VariableBase::Status::Fixed);
  graph->addVariable(VariableBasePtr(prev_pose));

  // dummy variables
  VariableVector3AD* bias_acc_from  = new VariableVector3AD();
  VariableVector3AD* bias_gyro_from = new VariableVector3AD();
  VariableVector3AD* bias_acc_to    = new VariableVector3AD();
  VariableVector3AD* bias_gyro_to   = new VariableVector3AD();
  VariableVector3AD* vel_from       = new VariableVector3AD();
  VariableVector3AD* vel_to         = new VariableVector3AD();

  bias_acc_from->setGraphId(1);
  bias_gyro_from->setGraphId(2);
  bias_acc_to->setGraphId(3);
  bias_gyro_to->setGraphId(4);
  vel_from->setGraphId(5);
  vel_to->setGraphId(6);
  graph->addVariable(VariableBasePtr(bias_acc_from));
  graph->addVariable(VariableBasePtr(bias_gyro_from));
  graph->addVariable(VariableBasePtr(bias_acc_to));
  graph->addVariable(VariableBasePtr(bias_gyro_to));
  graph->addVariable(VariableBasePtr(vel_from));
  graph->addVariable(VariableBasePtr(vel_to));

  int graph_id = 7;
  int i        = 0;
  while (std::getline(file, line)) {
    std::cout << graph_id << "\n";
    // Skip lines starting with '#'
    if (line[0] == '#') {
      continue;
    }
    // Tokenize the line using a stringstream
    std::istringstream ss(line);
    std::vector<std::string> tokens;
    std::string token;

    while (std::getline(ss, token, ',')) {
      tokens.push_back(token);
    }

    // Check if there are at least 8 columns
    if (tokens.size() < 8) {
      std::cerr << "Error: CSV file does not have enough columns." << std::endl;
      return;
    }

    // Extract the relevant columns
    float x  = std::stod(tokens[1]);
    float y  = std::stod(tokens[2]);
    float z  = std::stod(tokens[3]);
    float qx = std::stod(tokens[4]);
    float qy = std::stod(tokens[5]);
    float qz = std::stod(tokens[6]);
    float qw = std::stod(tokens[7]);

    Eigen::Isometry3f pose = Eigen::Isometry3f::Identity();

    pose.translation() = Vector3f(x, y, z);
    pose.linear()      = Eigen::Quaternionf(qw, qx, qy, qz).normalized().toRotationMatrix();
    if (i == 0) {
      prev_pose->setEstimate(pose);
      ++i;
      continue;
    }
    VariableSE3ExpMapRightAD* new_pose = new VariableSE3ExpMapRightAD();
    new_pose->setEstimate(pose);
    new_pose->setGraphId(++graph_id);
    new_pose->setStatus(VariableBase::Status::Fixed);
    graph->addVariable(VariableBasePtr(new_pose));

    // put a dummy factor between two variables for better visualization
    ImuPreintegrationFactorAD* factor = new ImuPreintegrationFactorAD();
    factor->setVariableId(0, prev_pose->graphId());
    factor->setVariableId(1, vel_from->graphId());
    factor->setVariableId(2, new_pose->graphId());
    factor->setVariableId(3, vel_to->graphId());
    factor->setVariableId(4, bias_acc_from->graphId());
    factor->setVariableId(5, bias_acc_to->graphId());
    factor->setVariableId(6, bias_gyro_from->graphId());
    factor->setVariableId(7, bias_gyro_to->graphId());

    graph->addFactor(FactorBasePtr(factor));

    prev_pose = new_pose;
    ++i;
  }
  file.close();

  graph->setSerializationLevel(-1);
  graph->write("/workspace/src/test_imu/prova.boss");
  std::cerr << "graph written in "
            << "/workspace/src/test_imu/prova.boss" << std::endl;
  return;
}

int main(int argc, char** argv) {
  initTypes();

  ParseCommandLine cmd_line(argv);
  cmd_line.parse();

  if (cmd_line.lastParsedArgs().empty()) {
    throw std::runtime_error("you forgot the input file, exiting");
  }

  for (const auto& input_file : cmd_line.lastParsedArgs()) {
    LOG << "loading graph from file [" << FG_YELLOW(input_file) << "]" << std::endl;

    // Check the file extension
    size_t dotPosition = input_file.find_last_of(".");
    if (dotPosition != std::string::npos) {
      std::string extension = input_file.substr(dotPosition + 1);

      if (extension == "csv") {
        // If it's a CSV file, call the parseCSVFile function
        FactorGraphPtr graph(new FactorGraph);
        parseCSVFile(input_file, graph);
        if (graph) {
          graphs[input_file] = graph;
        }
      } else if (extension == "boss") {
        FactorGraphPtr graph = FactorGraph::read(input_file);
        if (!graph) {
          continue;
        }
        graphs[input_file] = graph;
      }
    }
  }

  if (graphs.empty()) {
    return -1;
  }

  QApplication qapp(argc, argv);
  ViewerCoreSharedQGL viewer_core(argc, argv, &qapp);

  std::thread graph_t(viewGraph, viewer_core.getCanvas("viewer_core_shared_canvas"));
  viewer_core.startViewerServer();

  graph_t.join();
  return 0;
}
