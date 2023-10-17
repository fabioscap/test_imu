#include <thread>
#include <unistd.h>
// ia system utils
#include <srrg_system_utils/parse_command_line.h>
// ia solver stuff
#include <srrg_solver/solver_core/instances.h>
#include <srrg_solver/solver_core/internals/linear_solvers/instances.h>
#include <srrg_solver/solver_core/internals/linear_solvers/sparse_block_linear_solver_cholesky_csparse.h>
#include <srrg_solver/variables_and_factors/types_2d/instances.h>
#include <srrg_solver/variables_and_factors/types_3d/instances.h>
#include <srrg_solver/variables_and_factors/types_projective/instances.h>
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
  variables_and_factors_2d_registerTypes();
  variables_and_factors_3d_registerTypes();
  variables_and_factors_projective_registerTypes();
  solver_registerTypes();
  linear_solver_registerTypes();
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

int main(int argc, char** argv) {
  initTypes();

  ParseCommandLine cmd_line(argv);
  cmd_line.parse();

  if (cmd_line.lastParsedArgs().empty()) {
    throw std::runtime_error("you forgot the input file, exiting");
  }

  for (const auto& input_file : cmd_line.lastParsedArgs()) {
    LOG << "loading graph from file [" << FG_YELLOW(input_file) << "]" << std::endl;
    FactorGraphPtr graph = FactorGraph::read(input_file);
    if (!graph) {
      continue;
    }
    graphs[input_file] = graph;
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
